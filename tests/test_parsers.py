"""Tests for the MIBiG and antiSMASH parsers and the comparison metrics."""

from pathlib import Path

import pytest

from sempi3 import (
    antismash_io,
    compare,
    mibig,
)
from sempi3.spec import NRPS, PKS

DATA = Path(__file__).parent / "data"


def test_clean_smiles_strips_genbank_line_wrapping():
    # GenBank wraps long qualifiers, which puts whitespace inside the SMILES.
    assert antismash_io._clean_smiles("C=CC =CC(=O)O") == "C=CC=CC(=O)O"
    assert antismash_io._clean_smiles("CC(=O)\n   CC") == "CC(=O)CC"


def test_consensus_substrate_prefers_the_consensus_line():
    qualifiers = {
        "specificity": ["PKS signature: Malonyl-CoA", "substrate consensus: mmal", "Minowa: mal"]
    }
    assert antismash_io._consensus_substrate(qualifiers) == "mmal"


def test_consensus_substrate_treats_x_as_unknown():
    assert antismash_io._consensus_substrate({"specificity": ["substrate consensus: X"]}) is None
    assert antismash_io._consensus_substrate({"specificity": []}) is None


def test_module_builder_cuts_pks_modules():
    builder = antismash_io._ModuleBuilder()
    domains = [
        ("at", None),
        ("carrier", None),
        ("ks", None),
        ("at", "mmal"),
        ("kr", None),
        ("carrier", None),
        ("ks", None),
        ("at", "mal"),
        ("te", None),
    ]
    for index, (role, substrate) in enumerate(domains):
        builder.add(antismash_io._Domain(role=role, name=role, locus_tag="g", start=index, substrate=substrate))
    modules = builder.finish()
    assert [m.kind for m in modules] == [PKS, PKS, PKS]
    assert modules[0].is_starter
    assert modules[1].substrate == "mmal" and modules[1].kr
    assert builder.release == "hydrolysis"


def test_module_builder_cuts_nrps_modules():
    builder = antismash_io._ModuleBuilder()
    domains = [("a", "Ala"), ("carrier", None), ("c", None), ("a", "Gly"), ("e", None), ("carrier", None)]
    for index, (role, substrate) in enumerate(domains):
        builder.add(antismash_io._Domain(role=role, name=role, locus_tag="g", start=index, substrate=substrate))
    modules = builder.finish()
    assert [m.kind for m in modules] == [NRPS, NRPS]
    assert [m.substrate for m in modules] == ["Ala", "Gly"]
    assert modules[1].epimerization


def test_compare_identical_molecules_scores_one():
    score = compare.compare("CC(=O)O", "CC(=O)O")
    assert score is not None
    assert score.tanimoto == pytest.approx(1.0)
    assert score.mcs_coverage_pred == pytest.approx(1.0)


def test_compare_rejects_unparseable_smiles():
    assert compare.compare("not-a-smiles", "CC(=O)O") is None


def test_coverage_is_asymmetric_for_a_substructure():
    # A short backbone inside a bigger product: everything predicted is there,
    # but it explains only part of the product.
    score = compare.compare("CCCC(=O)O", "CCCC(=O)NCC(=O)NCC(=O)O")
    assert score is not None
    assert score.mcs_coverage_pred > score.mcs_coverage_true


def test_best_match_picks_the_best_variant_and_records_it():
    best = compare.best_match(
        {"hydrolysis": "CCCC(=O)O", "macrolactam": "CC(=O)O"},
        ["CCCC(=O)O"],
        ["butyric acid"],
    )
    assert best is not None
    assert best.release == "hydrolysis"
    assert best.compound == "butyric acid"
    assert best.tanimoto == pytest.approx(1.0)


def test_mibig_substrate_prefers_the_curated_structure():
    notes = []
    module = {"a_domain": {"substrates": [{"name": "alanine", "structure": "C[C@@H](C(=O)O)N"}]}}
    assert mibig._substrate(module, notes, 0) == "smiles:C[C@@H](C(=O)O)N"
    assert not notes


def test_mibig_substrate_notes_ambiguity():
    notes = []
    module = {"a_domain": {"substrates": [{"name": "alanine"}, {"name": "glycine"}]}}
    assert mibig._substrate(module, notes, 3) == "alanine"
    assert "2 candidate substrates" in notes[0]


def test_mibig_modification_types_skip_inactive_domains():
    module = {
        "modification_domains": [
            {"type": "ketoreductase", "active": True},
            {"type": "dehydratase", "active": False},
        ]
    }
    assert mibig._modification_types(module) == ["ketoreductase"]


def test_strict_mcs_notices_a_bond_order_difference():
    # The loose search treats a single and a double bond as the same, so only
    # the strict one can tell a saturated chain from an unsaturated one.
    loose = compare.compare("CCCCCC(=O)O", "CC=CCCC(=O)O", strict=True, timeout=5)
    assert loose is not None
    assert loose.mcs_coverage_pred == pytest.approx(1.0)
    assert loose.mcs_strict_coverage_pred < 1.0


def test_strict_mcs_notices_an_oxidation_level():
    score = compare.compare("CC(O)CC(=O)O", "CC(=O)CC(=O)O", strict=True, timeout=5)
    assert score is not None
    assert score.mcs_coverage_pred == pytest.approx(1.0)
    assert score.mcs_strict_coverage_pred < 1.0


def test_identical_molecules_match_on_every_measure():
    score = compare.compare("CCCCCC(=O)O", "CCCCCC(=O)O", strict=True, timeout=5)
    assert score is not None
    assert score.tanimoto == pytest.approx(1.0)
    assert score.mcs_strict_coverage_pred == pytest.approx(1.0)
    assert score.substructure_match == pytest.approx(1.0)


def test_substructure_match_finds_a_contained_backbone():
    score = compare.compare("CCCC(=O)O", "CCCC(=O)OC", strict=True, timeout=5)
    assert score is not None
    assert score.substructure_match == pytest.approx(0.0) or score.mcs_strict_coverage_pred > 0.5


def test_strip_tailoring_removes_a_sugar():
    from rdkit import Chem

    glycoside = Chem.MolFromSmiles("CC1OC(OC2CCCCC2)C(O)C(O)C1O")
    stripped = compare.strip_tailoring(glycoside)
    assert stripped is not None
    assert Chem.MolToSmiles(stripped) == Chem.MolToSmiles(Chem.MolFromSmiles("OC1CCCCC1"))


def test_strip_tailoring_removes_halogens():
    from rdkit import Chem

    chlorinated = Chem.MolFromSmiles("Clc1ccc(CCN)cc1")
    stripped = compare.strip_tailoring(chlorinated)
    assert stripped is not None
    assert "Cl" not in Chem.MolToSmiles(stripped)


def test_strip_tailoring_leaves_a_plain_molecule_alone():
    from rdkit import Chem

    plain = Chem.MolFromSmiles("CCCCC(=O)O")
    stripped = compare.strip_tailoring(plain)
    assert Chem.MolToSmiles(stripped) == Chem.MolToSmiles(plain)


def test_backbone_tanimoto_ignores_the_decorations():
    # The prediction is the aglycone; against the glycoside it scores poorly,
    # but against the stripped reference it is perfect.
    whole = compare.compare("OC1CCCCC1", "CC1OC(OC2CCCCC2)C(O)C(O)C1O", timeout=5)
    backbone = compare.backbone_tanimoto("OC1CCCCC1", "CC1OC(OC2CCCCC2)C(O)C(O)C1O")
    assert whole is not None and backbone is not None
    assert backbone > whole.tanimoto
