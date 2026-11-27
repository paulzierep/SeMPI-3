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
