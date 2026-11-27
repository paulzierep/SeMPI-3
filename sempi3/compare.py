"""Comparing a predicted backbone with a known structure.

A single similarity number hides what a backbone predictor gets right. Three
complementary measures are reported:

``tanimoto``
    Morgan fingerprint similarity, the standard whole-molecule measure.

``mcs_coverage_pred``
    Fraction of the *prediction* found in the known structure. High values mean
    "everything predicted is really there", which is what a backbone predictor
    should achieve even when tailoring enzymes add more atoms later.

``mcs_coverage_true``
    Fraction of the *known structure* explained by the prediction. Low values
    point at missing tailoring or a truncated backbone.
"""

from __future__ import annotations

from dataclasses import (
    asdict,
    dataclass,
)
from typing import (
    Iterable,
    List,
    Optional,
)

from rdkit import (
    DataStructs,
    RDLogger,
)
from rdkit.Chem import (
    rdFMCS,
    rdFingerprintGenerator,
)
from rdkit import Chem

RDLogger.DisableLog("rdApp.*")

#: Morgan fingerprint generator, radius 2 (ECFP4-like), 2048 bits.
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

#: Seconds allowed for one maximum-common-substructure search.
MCS_TIMEOUT = 10


@dataclass
class Similarity:
    """Scores for one predicted/known pair."""

    tanimoto: float = 0.0
    mcs_atoms: int = 0
    mcs_coverage_pred: float = 0.0
    mcs_coverage_true: float = 0.0
    #: Same, but with bond orders required to match. A single bond standing
    #: where the product has a double one breaks the match, so this is the
    #: measure that notices an oxidation level.
    mcs_strict_atoms: int = 0
    mcs_strict_coverage_pred: float = 0.0
    #: Largest fraction of the prediction that is literally a substructure of
    #: the product, element and bond order included.
    substructure_match: float = 0.0
    pred_atoms: int = 0
    true_atoms: int = 0
    compound: Optional[str] = None
    release: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


#: A pyranose or furanose ring, the usual glycosylation product.
SUGAR_RING = Chem.MolFromSmarts("[OX2;R]1[CX4;R][CX4;R][CX4;R][CX4;R][CX4;R]1")
SUGAR_RING_5 = Chem.MolFromSmarts("[OX2;R]1[CX4;R][CX4;R][CX4;R][CX4;R]1")

#: Halogens added by tailoring enzymes.
HALOGENS = (9, 17, 35, 53)


def strip_tailoring(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Remove sugars and halogens, leaving something closer to the backbone.

    SeMPI-3 predicts the chain that comes off the assembly line, while the
    reference is the finished natural product. Comparing the two directly
    charges the prediction for every decoration added afterwards, which hides
    whether the backbone itself got better. Stripping the most common
    decorations off the *reference* makes the comparison answer the backbone
    question instead.

    This is deliberately crude: sugar rings and halogen atoms only. It is a
    second view on the same prediction, never a replacement for the
    whole-molecule score.
    """
    if mol is None:
        return None
    editable = Chem.RWMol(mol)
    doomed = set()
    for pattern in (SUGAR_RING, SUGAR_RING_5):
        if pattern is None:
            continue
        for match in mol.GetSubstructMatches(pattern):
            ring = set(match)
            # Take the ring and whatever hangs off it that is not the aglycone:
            # a sugar is attached through exactly one bond.
            attachments = [
                atom.GetIdx()
                for index in ring
                for atom in mol.GetAtomWithIdx(index).GetNeighbors()
                if atom.GetIdx() not in ring
            ]
            linking = [a for a in attachments if mol.GetAtomWithIdx(a).GetDegree() > 1]
            if len(linking) <= 1:
                doomed |= ring
                doomed |= {a for a in attachments if mol.GetAtomWithIdx(a).GetDegree() == 1}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in HALOGENS:
            doomed.add(atom.GetIdx())
    for index in sorted(doomed, reverse=True):
        editable.RemoveAtom(index)
    try:
        stripped = editable.GetMol()
        Chem.SanitizeMol(stripped)
    except Exception:  # noqa: BLE001
        return mol
    fragments = Chem.GetMolFrags(stripped, asMols=True, sanitizeFrags=False)
    if not fragments:
        return mol
    # Keep the largest piece: the aglycone.
    return max(fragments, key=lambda f: f.GetNumHeavyAtoms())


def backbone_tanimoto(pred_smiles: str, true_smiles: str) -> Optional[float]:
    """Similarity against the reference with its decorations stripped off."""
    pred = parse(pred_smiles)
    true = parse(true_smiles)
    if pred is None or true is None:
        return None
    core = strip_tailoring(true)
    if core is None or core.GetNumHeavyAtoms() == 0:
        return None
    return tanimoto(pred, core)


def parse(smiles: Optional[str]) -> Optional[Chem.Mol]:
    """Parse SMILES, returning ``None`` instead of raising."""
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def tanimoto(pred: Chem.Mol, true: Chem.Mol) -> float:
    """Morgan fingerprint Tanimoto similarity, ignoring stereochemistry."""
    fp_a = _MORGAN.GetFingerprint(pred)
    fp_b = _MORGAN.GetFingerprint(true)
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def _mcs_atoms(pred, true, timeout: int, strict: bool) -> int:
    """Atoms in the largest common substructure of two molecules.

    With ``strict`` the bond orders have to agree, which makes the result
    sensitive to oxidation levels; without it a chain matches whether it is
    saturated or not.
    """
    try:
        result = rdFMCS.FindMCS(
            [pred, true],
            timeout=timeout,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact
            if strict
            else rdFMCS.BondCompare.CompareAny,
        )
        return result.numAtoms
    except Exception:  # noqa: BLE001 - RDKit raises broadly on odd inputs
        return 0


def compare(
    pred_smiles: str,
    true_smiles: str,
    timeout: int = MCS_TIMEOUT,
    with_mcs: bool = True,
    strict: bool = False,
) -> Optional[Similarity]:
    """Score one prediction against one known structure.

    Args:
        with_mcs: Run the maximum-common-substructure search. It dominates the
            runtime, so callers that only need the fingerprint score — such as
            :func:`best_match` while it is still choosing a pair — switch it off.
        strict: Also run the search with bond orders required to match, and
            test whether the prediction is a literal substructure.
    """
    pred = parse(pred_smiles)
    true = parse(true_smiles)
    if pred is None or true is None:
        return None

    result = Similarity(
        tanimoto=tanimoto(pred, true),
        pred_atoms=pred.GetNumHeavyAtoms(),
        true_atoms=true.GetNumHeavyAtoms(),
    )
    if not with_mcs:
        return result
    result.mcs_atoms = _mcs_atoms(pred, true, timeout, strict=False)
    if result.pred_atoms:
        result.mcs_coverage_pred = result.mcs_atoms / result.pred_atoms
    if result.true_atoms:
        result.mcs_coverage_true = result.mcs_atoms / result.true_atoms

    if strict:
        result.mcs_strict_atoms = _mcs_atoms(pred, true, timeout, strict=True)
        if result.pred_atoms:
            result.mcs_strict_coverage_pred = result.mcs_strict_atoms / result.pred_atoms
        result.substructure_match = 1.0 if true.HasSubstructMatch(pred) else 0.0
    return result


def best_match(
    predictions: dict,
    true_smiles_list: Iterable[str],
    compound_names: Optional[List[str]] = None,
    timeout: int = MCS_TIMEOUT,
    strict: bool = False,
) -> Optional[Similarity]:
    """Best score over all release variants and all known compounds.

    A cluster often makes several congeners, and the release mechanism is
    rarely annotated, so the kindest pairing is used. The chosen variant and
    compound are recorded on the result so this stays visible.
    """
    names = list(compound_names or [])
    best: Optional[Similarity] = None
    best_pair: Optional[tuple] = None
    # Choosing the pair needs only the cheap fingerprint score. The
    # maximum-common-substructure search runs once, on the winner.
    for index, true_smiles in enumerate(true_smiles_list):
        for release, pred_smiles in predictions.items():
            score = compare(pred_smiles, true_smiles, timeout=timeout, with_mcs=False)
            if score is None:
                continue
            score.release = release
            score.compound = names[index] if index < len(names) else None
            if best is None or score.tanimoto > best.tanimoto:
                best = score
                best_pair = (pred_smiles, true_smiles)

    if best is None or best_pair is None:
        return best
    full = compare(best_pair[0], best_pair[1], timeout=timeout, with_mcs=True, strict=strict)
    if full is not None:
        full.release = best.release
        full.compound = best.compound
        return full
    return best
