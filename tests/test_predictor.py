"""Tests for the end-to-end predictor used downstream of antiSMASH."""

from pathlib import Path

import pytest

from sempi3 import (
    predictor,
    stachelhaus as st,
)

DATA = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
EXAMPLE = DATA / "BGC0000290.gbk"

needs_data = pytest.mark.skipif(not EXAMPLE.exists(), reason="run scripts/fetch_data.sh first")


def test_signature_table_is_packaged():
    rows = predictor._reference()
    assert len(rows) > 1000
    signature, substrate, support = rows[0]
    assert len(signature) == len(st.EXTENDED_POSITIONS)
    assert substrate
    assert support >= 1


def test_predict_substrate_recognises_the_reference_pocket():
    # The packaged reference contains real domains; the phenylalanine pocket
    # of the GrsA reference has to come back as phenylalanine.
    substrate, confidence = predictor.predict_substrate(st.reference_sequence())
    assert substrate == "phenylalanine"
    assert confidence == pytest.approx(1.0, abs=0.2)


def test_predict_substrate_abstains_on_nonsense():
    substrate, _confidence = predictor.predict_substrate("M" * 400)
    assert substrate is None


def test_predict_substrate_abstains_on_a_short_sequence():
    substrate, confidence = predictor.predict_substrate("MKV")
    assert substrate is None
    assert confidence == 0.0


@needs_data
def test_predict_on_a_real_cluster():
    result = predictor.predict(EXAMPLE)
    assert result.bgc_id == "BGC0000290"
    assert result.smiles
    assert result.modules
    assert all(m["kind"] in ("nrps", "pks") for m in result.modules)
    # antiSMASH's own prediction is carried along for comparison.
    assert result.antismash_smiles


@needs_data
def test_tailoring_can_be_switched_off():
    decorated = predictor.predict(EXAMPLE, apply_tailoring=True)
    bare = predictor.predict(EXAMPLE, apply_tailoring=False)
    assert decorated.smiles != bare.smiles
    assert "Cl" in decorated.smiles
    assert "Cl" not in bare.smiles


@needs_data
def test_substrate_filling_can_be_switched_off():
    filled = predictor.predict(EXAMPLE, fill_substrates=True)
    unfilled = predictor.predict(EXAMPLE, fill_substrates=False)
    assert filled.filled_substrates >= unfilled.filled_substrates
