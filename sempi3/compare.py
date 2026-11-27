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
    pred_atoms: int = 0
    true_atoms: int = 0
    compound: Optional[str] = None
    release: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


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


def compare(
    pred_smiles: str,
    true_smiles: str,
    timeout: int = MCS_TIMEOUT,
    with_mcs: bool = True,
) -> Optional[Similarity]:
    """Score one prediction against one known structure.

    Args:
        with_mcs: Run the maximum-common-substructure search. It dominates the
            runtime, so callers that only need the fingerprint score — such as
            :func:`best_match` while it is still choosing a pair — switch it off.
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
    try:
        mcs = rdFMCS.FindMCS(
            [pred, true],
            timeout=timeout,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            bondCompare=rdFMCS.BondCompare.CompareAny,
        )
        result.mcs_atoms = mcs.numAtoms
    except Exception:  # noqa: BLE001 - RDKit raises broadly on odd inputs
        result.mcs_atoms = 0
    if result.pred_atoms:
        result.mcs_coverage_pred = result.mcs_atoms / result.pred_atoms
    if result.true_atoms:
        result.mcs_coverage_true = result.mcs_atoms / result.true_atoms
    return result


def best_match(
    predictions: dict,
    true_smiles_list: Iterable[str],
    compound_names: Optional[List[str]] = None,
    timeout: int = MCS_TIMEOUT,
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
    full = compare(best_pair[0], best_pair[1], timeout=timeout, with_mcs=True)
    if full is not None:
        full.release = best.release
        full.compound = best.compound
        return full
    return best
