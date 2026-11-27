"""Post-assembly tailoring.

Most of the difference between a predicted backbone and the finished natural
product is put there after the chain leaves the assembly line: sugars,
halogens, methyl groups and hydroxyls added by enzymes encoded elsewhere in
the cluster.

Which enzymes are present can be read from the cluster annotation. *Where*
they act cannot, so each reaction here is applied at a fixed, chemically
sensible site. That is a guess, and a wrong decoration costs similarity, so
every reaction is measured on its own before being switched on
(``scripts/experiment_tailoring.py``).
"""

from __future__ import annotations

import re
from dataclasses import (
    dataclass,
    field,
)
from typing import (
    Dict,
    List,
    Optional,
)

from rdkit import Chem

#: Enzyme families recognised in the cluster annotation.
#:
#: The patterns are matched against the ``product``, ``gene`` and
#: ``sec_met_domain`` text of every CDS in the record.
ENZYME_PATTERNS: Dict[str, re.Pattern] = {
    "glycosyltransferase": re.compile(r"glycosyl\s*transferase|glycosyltransferase|\bgtf\b", re.I),
    "halogenase": re.compile(r"halogenase|chlorinase|brominase", re.I),
    "o_methyltransferase": re.compile(r"\bo-methyltransferase|o_methyltransferase|omt\b", re.I),
    "hydroxylase": re.compile(r"hydroxylase|monooxygenase|cytochrome\s*p450|\bp450\b", re.I),
}

#: beta-D-glucose, attached through its anomeric carbon.
GLUCOSE = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"

#: A free hydroxyl that is not part of an acid: the usual acceptor.
FREE_HYDROXYL = Chem.MolFromSmarts("[OX2H][CX4]")

#: An aromatic carbon carrying a hydrogen: where halogenases act.
AROMATIC_CH = Chem.MolFromSmarts("[cH]")


@dataclass
class TailoringGenes:
    """Which tailoring enzymes a cluster encodes, and how many of each."""

    counts: Dict[str, int] = field(default_factory=dict)

    def has(self, name: str) -> bool:
        return self.counts.get(name, 0) > 0

    def to_dict(self) -> Dict[str, int]:
        return dict(self.counts)


def detect(record) -> TailoringGenes:
    """Count tailoring enzymes in a parsed GenBank record.

    The MIBiG files carry their original annotation, so a cluster whose genes
    are all "hypothetical protein" yields nothing. That is a limit of the
    input, not of the method.
    """
    counts: Dict[str, int] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        text = " ".join(
            value
            for key in ("product", "gene", "sec_met_domain", "note")
            for value in feature.qualifiers.get(key, [])
        )
        if not text:
            continue
        for name, pattern in ENZYME_PATTERNS.items():
            if pattern.search(text):
                counts[name] = counts.get(name, 0) + 1
    return TailoringGenes(counts=counts)


def _free_hydroxyls(mol: Chem.Mol) -> List[int]:
    """Oxygens of hydroxyl groups that are not part of a carboxylic acid."""
    acid = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_oxygens = {match[2] for match in mol.GetSubstructMatches(acid)}
    return [
        match[0]
        for match in mol.GetSubstructMatches(FREE_HYDROXYL)
        if match[0] not in acid_oxygens
    ]


def methylate(smiles: str) -> Optional[str]:
    """Add a methyl group to the last free hydroxyl."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    sites = _free_hydroxyls(mol)
    if not sites:
        return None
    editable = Chem.RWMol(mol)
    methyl = editable.AddAtom(Chem.Atom(6))
    editable.AddBond(sites[-1], methyl, Chem.BondType.SINGLE)
    editable.GetAtomWithIdx(sites[-1]).SetNumExplicitHs(0)
    return _finish(editable)


def glycosylate(smiles: str, sugar: str = GLUCOSE) -> Optional[str]:
    """Attach a sugar to the last free hydroxyl through a glycosidic bond."""
    mol = Chem.MolFromSmiles(smiles)
    sugar_mol = Chem.MolFromSmiles(sugar)
    if mol is None or sugar_mol is None:
        return None
    sites = _free_hydroxyls(mol)
    if not sites:
        return None
    anomeric = _anomeric_carbon(sugar_mol)
    if anomeric is None:
        return None
    anomeric_carbon, leaving_oxygen = anomeric

    offset = mol.GetNumAtoms()
    editable = Chem.RWMol(Chem.CombineMols(mol, sugar_mol))
    editable.AddBond(sites[-1], offset + anomeric_carbon, Chem.BondType.SINGLE)
    editable.GetAtomWithIdx(sites[-1]).SetNumExplicitHs(0)
    editable.RemoveAtom(offset + leaving_oxygen)
    return _finish(editable)


def _anomeric_carbon(sugar: Chem.Mol):
    """The ring carbon bearing both the ring oxygen and a free hydroxyl."""
    pattern = Chem.MolFromSmarts("[C;R]([OX2H])[OX2;R]")
    matches = sugar.GetSubstructMatches(pattern)
    if not matches:
        return None
    carbon, hydroxyl, _ring_oxygen = matches[0]
    return carbon, hydroxyl


def halogenate(smiles: str, halogen: str = "Cl") -> Optional[str]:
    """Put a halogen on the first aromatic carbon that carries a hydrogen."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    matches = mol.GetSubstructMatches(AROMATIC_CH)
    if not matches:
        return None
    editable = Chem.RWMol(mol)
    atom = editable.AddAtom(Chem.Atom(halogen))
    editable.AddBond(matches[0][0], atom, Chem.BondType.SINGLE)
    return _finish(editable)


def hydroxylate(smiles: str) -> Optional[str]:
    """Add a hydroxyl to an unsubstituted aliphatic carbon."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    pattern = Chem.MolFromSmarts("[CX4H2]")
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return None
    editable = Chem.RWMol(mol)
    oxygen = editable.AddAtom(Chem.Atom(8))
    editable.GetAtomWithIdx(oxygen).SetNumExplicitHs(1)
    editable.AddBond(matches[-1][0], oxygen, Chem.BondType.SINGLE)
    return _finish(editable)


#: Reaction name -> the function that applies it once.
REACTIONS = {
    "glycosyltransferase": glycosylate,
    "halogenase": halogenate,
    "o_methyltransferase": methylate,
    "hydroxylase": hydroxylate,
}


def _finish(editable: Chem.RWMol) -> Optional[str]:
    mol = editable.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001 - a failed decoration is simply dropped
        return None
    return Chem.MolToSmiles(mol)


def apply(smiles: str, genes: TailoringGenes, enabled: Optional[List[str]] = None, max_each: int = 1) -> str:
    """Decorate a backbone with the reactions its cluster encodes.

    Args:
        smiles: The assembled backbone.
        genes: Enzymes detected in the cluster.
        enabled: Which reactions to apply; all of them when ``None``.
        max_each: Upper bound on how often one reaction is applied, however
            many copies of the gene the cluster carries.
    """
    current = smiles
    for name, reaction in REACTIONS.items():
        if enabled is not None and name not in enabled:
            continue
        repeats = min(genes.counts.get(name, 0), max_each)
        for _ in range(repeats):
            result = reaction(current)
            if result is None:
                break
            current = result
    return current
