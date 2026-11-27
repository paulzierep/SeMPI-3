"""Tests for the chemistry primitives.

Products are checked against independently written reference SMILES and
against CIP codes, so that a wrong stereo descriptor cannot pass unnoticed.
"""

import pytest
from rdkit import Chem

from sempi3 import chem

L_ALANINE = "C[C@@H](C(=O)O)N"
GLYCINE = "NCC(=O)O"
ACETIC_ACID = "CC(=O)O"


def canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


def cip_codes(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return [a.GetPropsAsDict()["_CIPCode"] for a in mol.GetAtoms() if a.HasProp("_CIPCode")]


def test_find_backbone_alpha_amino_acid():
    mol = chem.parse_monomer(L_ALANINE)
    n_idx, c_idx, o_idx = chem.find_backbone(mol)
    assert mol.GetAtomWithIdx(n_idx).GetSymbol() == "N"
    assert mol.GetAtomWithIdx(c_idx).GetSymbol() == "C"
    assert mol.GetAtomWithIdx(o_idx).GetTotalNumHs() == 1


def test_find_backbone_ignores_side_chain_amide():
    # Asparagine's side-chain amide must not be taken for the backbone amine.
    mol = chem.parse_monomer("N[C@@H](CC(N)=O)C(=O)O")
    n_idx, _c_idx, _o_idx = chem.find_backbone(mol)
    nitrogen = mol.GetAtomWithIdx(n_idx)
    assert all(nb.GetSymbol() == "C" and not _is_amide_carbon(nb) for nb in nitrogen.GetNeighbors())


def _is_amide_carbon(atom) -> bool:
    return any(nb.GetSymbol() == "O" for nb in atom.GetNeighbors())


def test_find_backbone_acid_without_amine():
    mol = chem.parse_monomer("O=C(O)c1cccc(O)c1O")
    n_idx, c_idx, _o_idx = chem.find_backbone(mol)
    assert n_idx is None
    assert mol.GetAtomWithIdx(c_idx).GetSymbol() == "C"


def test_monomer_without_acid_is_rejected():
    with pytest.raises(chem.ChemError):
        chem.start_chain("CCO")


def test_dipeptide_keeps_both_configurations():
    chain = chem.start_chain(L_ALANINE)
    chain = chem.append_amino_acid(chain, L_ALANINE)
    product = chem.to_smiles(chem.terminate(chain, "hydrolysis"))
    assert product == canon("C[C@H](N)C(=O)N[C@@H](C)C(=O)O")
    assert cip_codes(product) == ["S", "S"]


def test_tripeptide_length():
    chain = chem.start_chain(GLYCINE)
    for _ in range(2):
        chain = chem.append_amino_acid(chain, GLYCINE)
    product = chem.to_smiles(chem.terminate(chain, "hydrolysis"))
    assert product == canon("NCC(=O)NCC(=O)NCC(=O)O")


def test_n_methylation():
    chain = chem.start_chain(GLYCINE)
    chain = chem.append_amino_acid(chain, GLYCINE, n_methylation=True)
    product = chem.to_smiles(chem.terminate(chain, "hydrolysis"))
    assert product == canon("NCC(=O)N(C)CC(=O)O")


@pytest.mark.parametrize(
    "state,expected",
    [
        ("keto", "CC(=O)CC(=O)O"),
        ("hydroxyl", "CC(O)CC(=O)O"),
        ("ene", "CC=CC(=O)O"),
        ("alkyl", "CCCC(=O)O"),
    ],
)
def test_pks_reduction_states(state, expected):
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, state)
    assert chem.to_smiles(chem.terminate(chain, "hydrolysis")) == canon(expected)


def test_pks_methylmalonyl_branch():
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, "keto", alpha_substituent="C")
    assert chem.to_smiles(chem.terminate(chain, "hydrolysis")) == canon("CC(=O)C(C)C(=O)O")


def test_two_pks_extensions():
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, "alkyl")
    chain = chem.append_pks_unit(chain, "keto")
    assert chem.to_smiles(chem.terminate(chain, "hydrolysis")) == canon("CCCC(=O)CC(=O)O")


def test_unknown_reduction_state_is_rejected():
    chain = chem.start_chain(ACETIC_ACID)
    with pytest.raises(chem.ChemError):
        chem.append_pks_unit(chain, "fully-reduced")


def test_macrolactam_closes_head_to_tail():
    chain = chem.start_chain(L_ALANINE)
    chain = chem.append_amino_acid(chain, L_ALANINE)
    product = chem.to_smiles(chem.terminate(chain, "macrolactam"))
    mol = Chem.MolFromSmiles(product)
    assert mol is not None
    assert mol.GetRingInfo().NumRings() == 1
    assert cip_codes(product) == ["S", "S"]


def test_macrolactone_closes_onto_hydroxyl():
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, "hydroxyl")
    chain = chem.append_pks_unit(chain, "alkyl")
    product = chem.to_smiles(chem.terminate(chain, "macrolactone"))
    mol = Chem.MolFromSmiles(product)
    assert mol is not None
    assert mol.GetRingInfo().NumRings() == 1
    # An ester, not a ketone: the ring oxygen is still there.
    assert product.count("O") >= 2


def test_reductive_release_gives_aldehyde():
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, "keto")
    product = chem.to_smiles(chem.terminate(chain, "reduction"))
    assert product == canon("CC(=O)CC=O")


def test_macrolactam_without_nitrogen_is_rejected():
    chain = chem.start_chain(ACETIC_ACID)
    chain = chem.append_pks_unit(chain, "keto")
    with pytest.raises(chem.ChemError):
        chem.terminate(chain, "macrolactam")
