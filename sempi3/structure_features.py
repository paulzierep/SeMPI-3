"""Properties of a product, predicted from the cluster rather than copied.

Finding the most similar known cluster and reporting its molecule is not
prediction: on MIBiG it scores well only because near-identical clusters are
in the database, and it collapses as soon as they are excluded.

What a cluster vector can honestly support is weaker but real: *how the
product is likely to look*. Does it carry a halogen, a sugar, a macrocycle?
Roughly how large is it, how many nitrogens does it hold? Those are learned
here as separate targets, each with a probability.

The predictions are then used to choose between structures SeMPI-3 has itself
assembled (`rank.py`), so the chemistry still comes from the assembly line and
the model only decides which of its outputs fits the cluster best.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import (
    Callable,
    Dict,
    List,
    Optional,
)

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    rdMolDescriptors,
)

#: Smallest ring treated as a macrocycle.
MACROCYCLE_MIN_SIZE = 8

#: A pyranose or furanose ring: the usual sign of glycosylation.
SUGAR = Chem.MolFromSmarts("[OX2;R][CX4;R][CX4;R][CX4;R][CX4;R]")

#: Ester and amide bonds, counted as a read on how the chain is linked.
ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
AMIDE = Chem.MolFromSmarts("[CX3](=O)[NX3]")


def _macrocycle(mol: Chem.Mol) -> float:
    return float(any(len(ring) >= MACROCYCLE_MIN_SIZE for ring in mol.GetRingInfo().AtomRings()))


def _halogens(mol: Chem.Mol) -> float:
    return float(sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ("F", "Cl", "Br", "I")))


def _element(symbol: str) -> Callable[[Chem.Mol], float]:
    def count(mol: Chem.Mol) -> float:
        return float(sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == symbol))

    return count


def _matches(pattern: Chem.Mol) -> Callable[[Chem.Mol], float]:
    def count(mol: Chem.Mol) -> float:
        return float(len(mol.GetSubstructMatches(pattern)))

    return count


#: Targets the model learns. Binary ones are probabilities, the rest counts.
TARGETS: Dict[str, Callable[[Chem.Mol], float]] = {
    "has_macrocycle": _macrocycle,
    "has_halogen": lambda mol: float(_halogens(mol) > 0),
    "has_sugar": lambda mol: float(mol.HasSubstructMatch(SUGAR)),
    "has_aromatic_ring": lambda mol: float(rdMolDescriptors.CalcNumAromaticRings(mol) > 0),
    "n_heavy_atoms": lambda mol: float(mol.GetNumHeavyAtoms()),
    "n_nitrogen": _element("N"),
    "n_oxygen": _element("O"),
    "n_sulfur": _element("S"),
    "n_rings": lambda mol: float(rdMolDescriptors.CalcNumRings(mol)),
    "n_aromatic_rings": lambda mol: float(rdMolDescriptors.CalcNumAromaticRings(mol)),
    "n_amide": _matches(AMIDE),
    "n_ester": _matches(ESTER),
    "n_rotatable": lambda mol: float(Descriptors.NumRotatableBonds(mol)),
    "mol_weight": lambda mol: float(Descriptors.MolWt(mol)),
}

#: Targets that are yes/no rather than counts.
BINARY_TARGETS = ("has_macrocycle", "has_halogen", "has_sugar", "has_aromatic_ring")


@dataclass
class StructureProfile:
    """Measured or predicted properties of one product."""

    values: Dict[str, float]

    def get(self, name: str, default: float = 0.0) -> float:
        return self.values.get(name, default)

    def to_dict(self) -> Dict[str, float]:
        return dict(self.values)


def describe(smiles: str) -> Optional[StructureProfile]:
    """Measure every target on a real structure."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    return StructureProfile({name: function(mol) for name, function in TARGETS.items()})


def describe_many(smiles_list: List[str]) -> List[StructureProfile]:
    profiles = []
    for smiles in smiles_list:
        profile = describe(smiles)
        if profile is not None:
            profiles.append(profile)
    return profiles
