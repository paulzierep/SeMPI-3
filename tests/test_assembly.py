"""Tests for the assembly rules and the monomer library."""

import pytest
from rdkit import Chem

from sempi3 import (
    assembly,
    monomers,
)
from sempi3.spec import (
    ClusterSpec,
    Module,
    NRPS,
    PKS,
)

L_ALANINE = "C[C@@H](C(=O)O)N"


def canon(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def nrps_module(substrate, **kwargs):
    return Module(kind=NRPS, substrate=substrate, **kwargs)


def test_every_monomer_parses():
    for name, smiles in monomers.AMINO_ACIDS.items():
        assert Chem.MolFromSmiles(smiles) is not None, name


def test_every_alias_points_at_a_structure():
    for alias, target in monomers.ALIASES.items():
        assert target in monomers.AMINO_ACIDS, alias


def test_resolve_is_case_insensitive_and_handles_isomer_prefix():
    assert monomers.resolve("Ala") == monomers.AMINO_ACIDS["alanine"]
    assert monomers.resolve("D-alanine") == monomers.AMINO_ACIDS["alanine"]
    # An isomer prefix in front of a short name has to work too, or the
    # residue is silently dropped when the chain is built.
    assert monomers.resolve("D-Ala") == monomers.AMINO_ACIDS["alanine"]
    assert monomers.resolve("D-Leu") == monomers.AMINO_ACIDS["leucine"]
    assert monomers.resolve("not-a-monomer") is None


def test_invert_alpha_carbon_flips_configuration():
    inverted = monomers.invert_alpha_carbon(L_ALANINE)
    assert canon(inverted) != canon(L_ALANINE)
    assert canon(monomers.invert_alpha_carbon(inverted)) == canon(L_ALANINE)


def test_invert_alpha_carbon_leaves_glycine_alone():
    glycine = monomers.AMINO_ACIDS["glycine"]
    assert monomers.invert_alpha_carbon(glycine) == glycine


def test_dipeptide_from_spec():
    spec = ClusterSpec(
        bgc_id="TEST",
        modules=[nrps_module("alanine", order=0), nrps_module("glycine", order=1)],
    )
    prediction = assembly.assemble(spec)
    assert not prediction.failed
    assert prediction.n_modules == 2
    assert prediction.variants["hydrolysis"] == canon("C[C@H](N)C(=O)NCC(=O)O")


def test_substrate_given_as_smiles():
    spec = ClusterSpec(bgc_id="TEST", modules=[nrps_module("smiles:" + L_ALANINE)])
    prediction = assembly.assemble(spec)
    assert prediction.variants["hydrolysis"] == canon(L_ALANINE)


def test_epimerisation_inverts_the_residue():
    plain = assembly.assemble(ClusterSpec(bgc_id="T", modules=[nrps_module("alanine")]))
    epimerised = assembly.assemble(
        ClusterSpec(bgc_id="T", modules=[nrps_module("alanine", epimerization=True)])
    )
    assert plain.variants["hydrolysis"] != epimerised.variants["hydrolysis"]


def test_unknown_substrate_fails_by_default():
    spec = ClusterSpec(bgc_id="T", modules=[nrps_module(None)])
    prediction = assembly.assemble(spec)
    assert prediction.failed


def test_unknown_substrate_can_be_skipped():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[nrps_module("alanine", order=0), nrps_module(None, order=1)],
    )
    prediction = assembly.assemble(spec, skip_unknown_substrates=True)
    assert not prediction.failed
    assert prediction.n_modules == 1
    assert any("no substrate" in note for note in prediction.notes)


def test_pks_starter_module_is_not_an_extension():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[
            Module(kind=PKS, substrate="mal", is_starter=True, order=0),
            Module(kind=PKS, substrate="mal", kr=True, order=1),
        ],
    )
    prediction = assembly.assemble(spec)
    # Acetyl starter plus one reduced extension: 3-hydroxybutyric acid.
    assert prediction.variants["hydrolysis"] == canon("CC(O)CC(=O)O")


def test_first_pks_module_assumes_an_acetyl_starter():
    spec = ClusterSpec(bgc_id="T", modules=[Module(kind=PKS, substrate="mal", order=0)])
    prediction = assembly.assemble(spec)
    assert prediction.variants["hydrolysis"] == canon("CC(=O)CC(=O)O")
    assert any("acetyl starter" in note for note in prediction.notes)


def test_assume_acetyl_starter_can_be_switched_off():
    spec = ClusterSpec(bgc_id="T", modules=[Module(kind=PKS, substrate="mal", order=0)])
    prediction = assembly.assemble(spec, assume_acetyl_starter=False)
    assert prediction.failed


def test_methylmalonyl_adds_a_branch():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[
            Module(kind=PKS, substrate="mal", is_starter=True, order=0),
            Module(kind=PKS, substrate="mmal", order=1),
        ],
    )
    prediction = assembly.assemble(spec)
    assert prediction.variants["hydrolysis"] == canon("CC(=O)C(C)C(=O)O")


def test_pks_methyltransferase_methylates_the_alpha_carbon():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[
            Module(kind=PKS, substrate="mal", is_starter=True, order=0),
            Module(kind=PKS, substrate="mal", n_methylation=True, order=1),
        ],
    )
    prediction = assembly.assemble(spec)
    assert prediction.variants["hydrolysis"] == canon("CC(=O)C(C)C(=O)O")


def test_hybrid_pks_nrps_chain():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[
            nrps_module("glycine", order=0),
            Module(kind=PKS, substrate="mal", order=1),
        ],
    )
    prediction = assembly.assemble(spec)
    assert not prediction.failed
    # The amide survives the polyketide extension.
    assert "N" in prediction.variants["hydrolysis"]


def test_empty_spec_fails_cleanly():
    prediction = assembly.assemble(ClusterSpec(bgc_id="T", modules=[]))
    assert prediction.failed
    assert "no modules" in prediction.notes


def test_variants_cover_several_release_modes():
    spec = ClusterSpec(
        bgc_id="T",
        modules=[nrps_module("alanine", order=0), nrps_module("alanine", order=1)],
    )
    prediction = assembly.assemble(spec)
    assert "hydrolysis" in prediction.variants
    assert "macrolactam" in prediction.variants


def test_dehydratase_without_ketoreductase_is_ignored():
    module = Module(kind=PKS, substrate="mal", dh=True)
    assert module.dh is False
    assert module.reduction_state == "keto"


@pytest.mark.parametrize(
    "kr,dh,er,expected",
    [
        (False, False, False, "keto"),
        (True, False, False, "hydroxyl"),
        (True, True, False, "ene"),
        (True, True, True, "alkyl"),
    ],
)
def test_reduction_state_mapping(kr, dh, er, expected):
    assert Module(kind=PKS, kr=kr, dh=dh, er=er).reduction_state == expected
