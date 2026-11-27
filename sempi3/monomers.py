"""Monomer library.

MIBiG stores a SMILES string with every curated adenylation-domain substrate,
so for NRPS modules the library is mostly a fallback. It is still needed for

* PKS extender and starter units, which MIBiG does not annotate at all, and
* sources that report a substrate by name only, such as antiSMASH.
"""

from __future__ import annotations

from typing import (
    Dict,
    Optional,
)

from rdkit import Chem

from .chem import (
    ALPHA_AMINO_ACID,
    ChemError,
    OMEGA_AMINO_ACID,
)

#: Starter unit laid down by a PKS loading module, as a free acid.
#: After decarboxylation a malonyl-CoA derived starter is an acetyl group.
PKS_STARTERS: Dict[str, str] = {
    "mal": "CC(=O)O",  # acetyl
    "mmal": "CCC(=O)O",  # propionyl
    "emal": "CCCC(=O)O",  # butyryl
    "mxmal": "COCC(=O)O",  # methoxyacetyl
    "acetyl": "CC(=O)O",
    "propionyl": "CCC(=O)O",
    "benzoyl": "O=C(O)c1ccccc1",
    "isobutyryl": "CC(C)C(=O)O",
}

#: Branch added to the alpha carbon by a PKS extender unit. ``None`` means the
#: alpha carbon stays unsubstituted, as it does for malonyl-CoA.
PKS_EXTENDER_BRANCHES: Dict[str, Optional[str]] = {
    "mal": None,  # malonyl-CoA
    "mmal": "C",  # methylmalonyl-CoA
    "emal": "CC",  # ethylmalonyl-CoA
    "mxmal": "OC",  # methoxymalonyl-ACP
    "allylmal": "CC=C",
    "isobutmal": "C(C)C",
    "chlorethmal": "CCl",
}

#: Proteinogenic and common non-proteinogenic NRPS substrates, as free acids.
#: Only used when the source gives a name without a structure.
AMINO_ACIDS: Dict[str, str] = {
    "alanine": "C[C@@H](C(=O)O)N",
    "arginine": "N[C@@H](CCCNC(N)=N)C(=O)O",
    "asparagine": "N[C@@H](CC(N)=O)C(=O)O",
    "aspartic acid": "N[C@@H](CC(=O)O)C(=O)O",
    "cysteine": "N[C@@H](CS)C(=O)O",
    "glutamic acid": "N[C@@H](CCC(=O)O)C(=O)O",
    "glutamine": "N[C@@H](CCC(N)=O)C(=O)O",
    "glycine": "NCC(=O)O",
    "histidine": "N[C@@H](Cc1cnc[nH]1)C(=O)O",
    "isoleucine": "CC[C@H](C)[C@H](N)C(=O)O",
    "leucine": "CC(C)C[C@H](N)C(=O)O",
    "lysine": "NCCCC[C@H](N)C(=O)O",
    "methionine": "CSCC[C@H](N)C(=O)O",
    "phenylalanine": "N[C@@H](Cc1ccccc1)C(=O)O",
    "proline": "O=C(O)[C@@H]1CCCN1",
    "serine": "N[C@@H](CO)C(=O)O",
    "threonine": "C[C@@H](O)[C@H](N)C(=O)O",
    "tryptophan": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
    "tyrosine": "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "valine": "CC(C)[C@H](N)C(=O)O",
    # Frequent non-proteinogenic substrates.
    "ornithine": "NCCC[C@H](N)C(=O)O",
    "2,4-diaminobutyric acid": "NCC[C@H](N)C(=O)O",
    "2,3-diaminopropionic acid": "NC[C@H](N)C(=O)O",
    "4-hydroxyphenylglycine": "N[C@H](C(=O)O)c1ccc(O)cc1",
    "3,5-dihydroxyphenylglycine": "N[C@H](C(=O)O)c1cc(O)cc(O)c1",
    "phenylglycine": "N[C@H](C(=O)O)c1ccccc1",
    "2,3-dihydroxybenzoic acid": "O=C(O)c1cccc(O)c1O",
    "salicylic acid": "O=C(O)c1ccccc1O",
    "pipecolic acid": "O=C(O)[C@@H]1CCCCN1",
    "beta-alanine": "NCCC(=O)O",
    "2-aminobutyric acid": "CC[C@H](N)C(=O)O",
    "2-aminoadipic acid": "N[C@@H](CCCC(=O)O)C(=O)O",
    "2-aminoisobutyric acid": "CC(C)(N)C(=O)O",
    "isovaline": "CC[C@](C)(N)C(=O)O",
    "citrulline": "NC(=O)NCCC[C@H](N)C(=O)O",
    "kynurenine": "N[C@@H](CC(=O)c1ccccc1N)C(=O)O",
    # Substrates named by antiSMASH that are not in the list above.
    "4-hydroxyproline": "O[C@H]1CN[C@@H](C(=O)O)C1",
    "beta-hydroxytyrosine": "N[C@@H]([C@H](O)c1ccc(O)cc1)C(=O)O",
    "beta-lysine": "NC[C@@H](N)CCC(=O)O",
    "n5-hydroxyornithine": "ON(O)CCC[C@H](N)C(=O)O",
    "n5-formyl-n5-hydroxyornithine": "O=CN(O)CCC[C@H](N)C(=O)O",
    "n5-acetyl-n5-hydroxyornithine": "CC(=O)N(O)CCC[C@H](N)C(=O)O",
    "allo-threonine": "C[C@H](O)[C@H](N)C(=O)O",
    "2-hydroxyisovaleric acid": "CC(C)[C@@H](O)C(=O)O",
    "4-methylproline": "C[C@H]1CN[C@@H](C(=O)O)C1",
    "4-butenyl-4-methylthreonine": "C/C=C/C[C@@H](C)[C@H](O)[C@H](N)C(=O)O",
    "glycolic acid": "OCC(=O)O",
    "quinoxaline-2-carboxylic acid": "O=C(O)c1cnc2ccccc2n1",
    "anthranilic acid": "Nc1ccccc1C(=O)O",
}

#: Short names used by antiSMASH and older MIBiG releases.
ALIASES: Dict[str, str] = {
    "ala": "alanine",
    "arg": "arginine",
    "asn": "asparagine",
    "asp": "aspartic acid",
    "cys": "cysteine",
    "gln": "glutamine",
    "glu": "glutamic acid",
    "gly": "glycine",
    "his": "histidine",
    "ile": "isoleucine",
    "leu": "leucine",
    "lys": "lysine",
    "met": "methionine",
    "phe": "phenylalanine",
    "pro": "proline",
    "ser": "serine",
    "thr": "threonine",
    "trp": "tryptophan",
    "tyr": "tyrosine",
    "val": "valine",
    "orn": "ornithine",
    "dab": "2,4-diaminobutyric acid",
    "dpr": "2,3-diaminopropionic acid",
    "hpg": "4-hydroxyphenylglycine",
    "dhpg": "3,5-dihydroxyphenylglycine",
    "phg": "phenylglycine",
    "dhb": "2,3-dihydroxybenzoic acid",
    "sal": "salicylic acid",
    "pip": "pipecolic acid",
    "b-ala": "beta-alanine",
    "bala": "beta-alanine",
    "abu": "2-aminobutyric acid",
    "aad": "2-aminoadipic acid",
    "aib": "2-aminoisobutyric acid",
    "iva": "isovaline",
    "cit": "citrulline",
    "dhpg": "3,5-dihydroxyphenylglycine",
    "hpr": "4-hydroxyproline",
    "boh-tyr": "beta-hydroxytyrosine",
    "r-beta-hydroxytyrosine": "beta-hydroxytyrosine",
    "blys": "beta-lysine",
    "oh-orn": "n5-hydroxyornithine",
    "fo-oh-orn": "n5-formyl-n5-hydroxyornithine",
    "ac-oh-orn": "n5-acetyl-n5-hydroxyornithine",
    "athr": "allo-threonine",
    "dioh-bz": "2,3-dihydroxybenzoic acid",
    "hiv": "2-hydroxyisovaleric acid",
    "d-hiv": "2-hydroxyisovaleric acid",
    "ival": "isovaline",
    "bmt": "4-butenyl-4-methylthreonine",
    "me-pro": "4-methylproline",
    "4me-pro": "4-methylproline",
    "ala-b": "beta-alanine",
    "gly-oh": "glycolic acid",
    "qa": "quinoxaline-2-carboxylic acid",
    "abz": "anthranilic acid",
}


def resolve(name: str) -> Optional[str]:
    """Return the free-acid SMILES for a substrate name, or ``None``."""
    if not name:
        return None
    key = name.strip().lower()
    key = ALIASES.get(key, key)
    if key in AMINO_ACIDS:
        return AMINO_ACIDS[key]
    # Names such as "D-alanine" or "N-methyl-leucine" carry the modification in
    # the name; the assembly engine applies those separately, so strip them.
    for prefix in ("d-", "l-", "(2s)-", "(2r)-"):
        if key.startswith(prefix) and key[len(prefix) :] in AMINO_ACIDS:
            return AMINO_ACIDS[key[len(prefix) :]]
    return None


def starter_smiles(substrate: Optional[str], default: str = "mal") -> str:
    """SMILES of a PKS starter unit."""
    key = (substrate or default).strip().lower()
    return PKS_STARTERS.get(key, PKS_STARTERS[default])


def extender_branch(substrate: Optional[str], default: Optional[str] = None) -> Optional[str]:
    """Branch SMILES added to the alpha carbon by a PKS extender unit."""
    if not substrate:
        return default
    return PKS_EXTENDER_BRANCHES.get(substrate.strip().lower(), default)


def invert_alpha_carbon(smiles: str) -> str:
    """Invert the configuration of a monomer's alpha carbon.

    Used for epimerisation domains, which turn an L-residue into a D-residue.
    Monomers without a stereocentre, such as glycine, come back unchanged.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemError(f"cannot parse monomer SMILES: {smiles!r}")
    for pattern in (ALPHA_AMINO_ACID, OMEGA_AMINO_ACID):
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue
        alpha = mol.GetAtomWithIdx(matches[0][1])
        tag = alpha.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            alpha.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            alpha.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        else:
            return smiles
        return Chem.MolToSmiles(mol)
    return smiles
