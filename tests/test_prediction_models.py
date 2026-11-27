"""Tests for the binding-pocket model, the datasets and the tailoring steps."""

import pytest
from rdkit import Chem

from sempi3 import (
    datasets,
    stachelhaus as st,
    substrates,
    tailoring,
)


def test_reference_reproduces_the_published_code():
    # The numbering is only right if the reference itself reads DAWTIAAIC.
    assert st.self_check()


def test_extract_code_on_the_reference():
    assert st.extract_code(st.reference_sequence()) == st.REFERENCE_CODE


def test_extract_code_handles_a_fragment():
    # A domain is usually a slice of a larger protein; the alignment must
    # still land on the right positions.
    fragment = st.reference_sequence()[180:420]
    assert st.extract_code(fragment) == st.REFERENCE_CODE


def test_extract_code_rejects_a_short_sequence():
    assert st.extract_code("MKV") == st.GAP * len(st.POCKET_POSITIONS)


def test_signature_contains_the_code():
    signature = st.extract_signature(st.reference_sequence())
    assert len(signature) == len(st.EXTENDED_POSITIONS)
    assert st.code_from_signature(signature) == st.REFERENCE_CODE


def test_every_pocket_position_is_in_the_extended_set():
    assert set(st.POCKET_POSITIONS) <= set(st.EXTENDED_POSITIONS)


def test_code_from_signature_rejects_the_wrong_length():
    assert st.code_from_signature("TOOSHORT") == st.GAP * len(st.POCKET_POSITIONS)


def test_predictor_returns_an_exact_match_first():
    predictor = st.StachelhausPredictor().fit(["DAWTIAAIC", "DVWHFSLVD"], ["phenylalanine", "serine"])
    match = predictor.predict_one("DAWTIAAIC")
    assert match.substrate == "phenylalanine"
    assert match.identity == 1.0


def test_predictor_abstains_when_nothing_is_close():
    predictor = st.StachelhausPredictor(min_identity=0.9).fit(["DAWTIAAIC"], ["phenylalanine"])
    assert predictor.predict_one("GGGGGGGGG").substrate is None


def test_predictor_takes_the_nearest_code():
    predictor = st.StachelhausPredictor(min_identity=0.5).fit(
        ["DAWTIAAIC", "GGGGGGGGG"], ["phenylalanine", "glycine"]
    )
    # One substitution away from the phenylalanine code.
    assert predictor.predict_one("DAWTIAAIV").substrate == "phenylalanine"


@pytest.mark.parametrize(
    "name,expected",
    [
        ("Leu", "leucine"),
        ("leucine", "leucine"),
        ("D-leucine", "leucine"),
        ("Dhpg", "3,5-dihydroxyphenylglycine"),
        ("X", None),
        ("", None),
        (None, None),
    ],
)
def test_normalise_substrate_names(name, expected):
    assert substrates.normalise(name) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC1CCCCCCCCCC(=O)O1", "macrolactone"),
        ("O=C1NCCCCCCCCC1", "macrolactam"),
        ("CCCCCC(=O)O", "linear"),
        # A five-membered lactam such as proline is not a macrocycle.
        ("O=C(O)[C@@H]1CCCN1", "linear"),
        ("CC(=O)OC", "linear"),
        ("not a molecule", "unknown"),
    ],
)
def test_release_label(smiles, expected):
    assert datasets.release_label(smiles) == expected


def _valid(smiles):
    return smiles is not None and Chem.MolFromSmiles(smiles) is not None


BACKBONE = "CC(O)CC(=O)NCc1ccccc1"


def test_methylation_adds_one_carbon():
    product = tailoring.methylate(BACKBONE)
    assert _valid(product)
    before = Chem.MolFromSmiles(BACKBONE).GetNumHeavyAtoms()
    assert Chem.MolFromSmiles(product).GetNumHeavyAtoms() == before + 1


def test_glycosylation_attaches_a_sugar():
    product = tailoring.glycosylate(BACKBONE)
    assert _valid(product)
    # Glucose contributes eleven heavy atoms once its anomeric OH has left.
    before = Chem.MolFromSmiles(BACKBONE).GetNumHeavyAtoms()
    assert Chem.MolFromSmiles(product).GetNumHeavyAtoms() == before + 11


def test_halogenation_puts_chlorine_on_a_ring():
    product = tailoring.halogenate(BACKBONE)
    assert _valid(product)
    assert "Cl" in product


def test_hydroxylation_adds_an_oxygen():
    product = tailoring.hydroxylate(BACKBONE)
    assert _valid(product)
    before = Chem.MolFromSmiles(BACKBONE).GetNumHeavyAtoms()
    assert Chem.MolFromSmiles(product).GetNumHeavyAtoms() == before + 1


def test_reactions_leave_a_molecule_without_a_site_alone():
    # No hydroxyl, no aromatic ring: nothing to decorate.
    assert tailoring.methylate("CCC(=O)N") is None
    assert tailoring.halogenate("CCC(=O)N") is None


def test_apply_only_runs_the_enabled_reactions():
    genes = tailoring.TailoringGenes(counts={"halogenase": 1, "o_methyltransferase": 1})
    only_halogen = tailoring.apply(BACKBONE, genes, enabled=["halogenase"])
    assert "Cl" in only_halogen
    assert only_halogen.count("OC") == BACKBONE.count("OC")


def test_apply_respects_the_repeat_cap():
    genes = tailoring.TailoringGenes(counts={"halogenase": 5})
    once = tailoring.apply(BACKBONE, genes, enabled=["halogenase"], max_each=1)
    twice = tailoring.apply(BACKBONE, genes, enabled=["halogenase"], max_each=2)
    assert once.count("Cl") == 1
    assert twice.count("Cl") == 2


def test_apply_without_genes_changes_nothing():
    assert tailoring.apply(BACKBONE, tailoring.TailoringGenes()) == BACKBONE


class _Feature:
    def __init__(self, feature_type, qualifiers):
        self.type = feature_type
        self.qualifiers = qualifiers


class _Record:
    def __init__(self, features):
        self.features = features


def test_detect_counts_enzyme_families():
    record = _Record(
        [
            _Feature("CDS", {"product": ["putative glycosyltransferase"]}),
            _Feature("CDS", {"product": ["tryptophan halogenase"]}),
            _Feature("CDS", {"product": ["hypothetical protein"]}),
            _Feature("gene", {"product": ["glycosyltransferase"]}),  # not a CDS
        ]
    )
    counts = tailoring.detect(record).to_dict()
    assert counts == {"glycosyltransferase": 1, "halogenase": 1}
