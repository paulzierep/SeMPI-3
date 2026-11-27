"""RDKit primitives for growing an assembly-line product.

The chain is grown as a real molecule, not as concatenated SMILES text, so
that exotic monomers (MIBiG stores a SMILES string for every curated
substrate) can be used without writing a fragment by hand for each one.

The growing chain always carries a free carboxylic acid at its active end,
standing in for the thioester that is bound to the carrier protein. Atom map
numbers tag the atoms the next reaction needs, because RDKit atom indices
shift whenever an atom is removed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import (
    List,
    Optional,
    Tuple,
)

from rdkit import Chem
from rdkit.Chem import AllChem  # noqa: F401  (import keeps RDKit warnings quiet)

#: Atom map numbers used to tag the atoms that take part in the next bond.
TAG_CARBOXYL_C = 1  # carbonyl carbon of the active (thioester-mimicking) acid
TAG_CARBOXYL_O = 2  # the -OH oxygen of that acid, the leaving group
TAG_NTERM = 3  # backbone nitrogen of the first monomer, for macrolactamisation

#: Backbone of an alpha-amino acid: free amine, alpha carbon, carboxylic acid.
#: The nitrogen must not already be an amide, so that a side-chain amide such
#: as the one in asparagine is never mistaken for the backbone.
ALPHA_AMINO_ACID = Chem.MolFromSmarts("[NX3;H2,H1;!$([NX3][CX3]=[OX1])][CX4][CX3](=[OX1])[OX2H1]")

#: Fallback for beta- and longer-chain amino acids (e.g. beta-alanine).
OMEGA_AMINO_ACID = Chem.MolFromSmarts("[NX3;H2,H1;!$([NX3][CX3]=[OX1])][CX4][CX4][CX3](=[OX1])[OX2H1]")

#: Any carboxylic acid, used for starter units without a nitrogen.
CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")


class ChemError(ValueError):
    """Raised when a monomer or a chain cannot be handled."""


@dataclass
class Chain:
    """A growing assembly-line intermediate."""

    mol: Chem.Mol
    #: Reduction states applied so far, for reporting.
    history: List[str]

    def copy(self) -> "Chain":
        return Chain(Chem.Mol(self.mol), list(self.history))


def _tagged_atom(mol: Chem.Mol, tag: int) -> Optional[Chem.Atom]:
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == tag:
            return atom
    return None


def _tagged_index(mol: Chem.Mol, tag: int) -> int:
    atom = _tagged_atom(mol, tag)
    if atom is None:
        raise ChemError(f"chain has lost its tag {tag}")
    return atom.GetIdx()


def _clear_tags(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def parse_monomer(smiles: str) -> Chem.Mol:
    """Parse a monomer SMILES, keeping explicit hydrogens off the graph."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ChemError(f"cannot parse monomer SMILES: {smiles!r}")
    return mol


def find_backbone(mol: Chem.Mol) -> Tuple[Optional[int], int, int]:
    """Locate the backbone of a monomer.

    Returns a tuple ``(n_idx, carboxyl_c_idx, hydroxyl_o_idx)``. ``n_idx`` is
    ``None`` for monomers without an amine, such as the aryl acid starters used
    by some NRPS systems.

    Raises:
        ChemError: if the monomer has no carboxylic acid at all.
    """
    for pattern in (ALPHA_AMINO_ACID, OMEGA_AMINO_ACID):
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            match = matches[0]
            # The SMARTS end with ...C(=O)[OH]: carbon, carbonyl O, hydroxyl O.
            return match[0], match[-3], match[-1]
    acid_matches = mol.GetSubstructMatches(CARBOXYLIC_ACID)
    if acid_matches:
        c_idx, _o_double, o_single = acid_matches[0]
        return None, c_idx, o_single
    raise ChemError(f"monomer has no carboxylic acid: {Chem.MolToSmiles(mol)}")


def start_chain(monomer_smiles: str) -> Chain:
    """Begin a chain from a loading-module substrate."""
    mol = parse_monomer(monomer_smiles)
    n_idx, c_idx, o_idx = find_backbone(mol)
    rw = Chem.RWMol(mol)
    rw.GetAtomWithIdx(c_idx).SetAtomMapNum(TAG_CARBOXYL_C)
    rw.GetAtomWithIdx(o_idx).SetAtomMapNum(TAG_CARBOXYL_O)
    if n_idx is not None:
        rw.GetAtomWithIdx(n_idx).SetAtomMapNum(TAG_NTERM)
    return Chain(rw.GetMol(), history=["start"])


def _combine(chain_mol: Chem.Mol, monomer: Chem.Mol) -> Tuple[Chem.RWMol, int]:
    """Combine chain and monomer, returning the editable mol and the index shift."""
    offset = chain_mol.GetNumAtoms()
    combined = Chem.RWMol(Chem.CombineMols(chain_mol, monomer))
    return combined, offset


def _drop_leaving_hydroxyl(rw: Chem.RWMol) -> None:
    """Remove the -OH that leaves when the active acid forms a new bond."""
    o_idx = _tagged_index(rw, TAG_CARBOXYL_O)
    rw.RemoveAtom(o_idx)


def append_amino_acid(chain: Chain, monomer_smiles: str, n_methylation: bool = False) -> Chain:
    """Form an amide bond between the active acid and a new amino acid.

    Args:
        chain: The growing chain.
        monomer_smiles: Free form of the monomer, as stored by MIBiG.
        n_methylation: Add a methyl group to the new backbone nitrogen.
    """
    monomer = parse_monomer(monomer_smiles)
    n_idx, c_idx, o_idx = find_backbone(monomer)
    if n_idx is None:
        raise ChemError(f"monomer cannot extend a chain, it has no amine: {monomer_smiles}")

    rw, offset = _combine(chain.mol, monomer)
    new_n = offset + n_idx
    # Tag the incoming acid before any atom is removed, so the tags survive
    # the index shift caused by removing the leaving hydroxyl.
    rw.GetAtomWithIdx(offset + c_idx).SetAtomMapNum(TAG_CARBOXYL_C + 100)
    rw.GetAtomWithIdx(offset + o_idx).SetAtomMapNum(TAG_CARBOXYL_O + 100)

    old_c = _tagged_index(rw, TAG_CARBOXYL_C)
    rw.AddBond(old_c, new_n, Chem.BondType.SINGLE)
    if n_methylation:
        methyl = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(new_n, methyl, Chem.BondType.SINGLE)

    # The old acid has reacted: drop its hydroxyl and hand the tags over.
    rw.GetAtomWithIdx(old_c).SetAtomMapNum(0)
    _drop_leaving_hydroxyl(rw)
    _promote_incoming_tags(rw)
    mol = _finalise(rw)
    return Chain(mol, chain.history + ["amide"])


def _promote_incoming_tags(rw: Chem.RWMol) -> None:
    """Turn the +100 tags of the incoming monomer into the active tags."""
    for atom in rw.GetAtoms():
        tag = atom.GetAtomMapNum()
        if tag == TAG_CARBOXYL_C + 100:
            atom.SetAtomMapNum(TAG_CARBOXYL_C)
        elif tag == TAG_CARBOXYL_O + 100:
            atom.SetAtomMapNum(TAG_CARBOXYL_O)
        elif tag == TAG_NTERM + 100:
            # Only the first monomer keeps an N-terminus tag.
            atom.SetAtomMapNum(0)


def _finalise(rw: Chem.RWMol) -> Chem.Mol:
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # pragma: no cover - defensive
        raise ChemError(f"assembled chain does not sanitise: {exc}") from exc
    return mol


#: Side chain that a heterocyclisation domain closes onto: the thiol of
#: cysteine (giving a thiazoline) or the hydroxyl of serine and threonine
#: (giving an oxazoline).
CYCLISABLE_SIDE_CHAIN = Chem.MolFromSmarts("[NX3;H2,H1][CX4][CX4][SX2H1,OX2H1]")

#: Temporary tags used while a heterocycle is closed.
TAG_RING_CARBON = 8
TAG_RING_NITROGEN = 9
TAG_LEAVING_OXYGEN = 10


def append_heterocycle(chain: Chain, monomer_smiles: str) -> Chain:
    """Add a residue and close it into a thiazoline or oxazoline ring.

    A heterocyclisation (Cy) domain forms the peptide bond and then
    dehydrates it: the side-chain thiol or hydroxyl attacks the preceding
    carbonyl, water leaves, and a five-membered ring containing a C=N bond is
    left behind. Cysteine gives a thiazoline, serine and threonine an
    oxazoline.

    Raises:
        ChemError: if the monomer has no side chain that can close the ring,
            or the chain does not end in a carbonyl.
    """
    monomer = parse_monomer(monomer_smiles)
    matches = monomer.GetSubstructMatches(CYCLISABLE_SIDE_CHAIN)
    if not matches:
        raise ChemError(f"monomer cannot form a heterocycle: {monomer_smiles}")
    n_idx, _alpha_idx, _beta_idx, side_idx = matches[0]
    _amine, c_idx, o_idx = find_backbone(monomer)

    rw, offset = _combine(chain.mol, monomer)
    rw.GetAtomWithIdx(offset + c_idx).SetAtomMapNum(TAG_CARBOXYL_C + 100)
    rw.GetAtomWithIdx(offset + o_idx).SetAtomMapNum(TAG_CARBOXYL_O + 100)
    rw.GetAtomWithIdx(offset + n_idx).SetAtomMapNum(TAG_RING_NITROGEN)

    side_atom = rw.GetAtomWithIdx(offset + side_idx)
    side_atom.SetNumExplicitHs(0)
    side_atom.SetNoImplicit(True)
    side_atom.SetAtomMapNum(TAG_CYCLISATION_PARTNER)

    ring_carbon = _tagged_index(rw, TAG_CARBOXYL_C)
    rw.AddBond(ring_carbon, offset + n_idx, Chem.BondType.SINGLE)

    carbonyl_oxygen = None
    for neighbour in rw.GetAtomWithIdx(ring_carbon).GetNeighbors():
        bond = rw.GetBondBetweenAtoms(ring_carbon, neighbour.GetIdx())
        if neighbour.GetSymbol() == "O" and bond.GetBondType() == Chem.BondType.DOUBLE:
            carbonyl_oxygen = neighbour.GetIdx()
            break
    if carbonyl_oxygen is None:
        raise ChemError("the chain does not end in a carbonyl to cyclise")

    # Everything that must survive the two removals is tagged first.
    rw.GetAtomWithIdx(ring_carbon).SetAtomMapNum(TAG_RING_CARBON)
    rw.GetAtomWithIdx(carbonyl_oxygen).SetAtomMapNum(TAG_LEAVING_OXYGEN)

    _drop_leaving_hydroxyl(rw)
    rw.RemoveAtom(_tagged_index(rw, TAG_LEAVING_OXYGEN))

    ring_carbon = _tagged_index(rw, TAG_RING_CARBON)
    ring_nitrogen = _tagged_index(rw, TAG_RING_NITROGEN)
    ring_side = _tagged_index(rw, TAG_CYCLISATION_PARTNER)
    rw.GetBondBetweenAtoms(ring_carbon, ring_nitrogen).SetBondType(Chem.BondType.DOUBLE)
    rw.AddBond(ring_carbon, ring_side, Chem.BondType.SINGLE)
    for tag in (TAG_RING_CARBON, TAG_RING_NITROGEN, TAG_CYCLISATION_PARTNER):
        atom = _tagged_atom(rw, tag)
        if atom is not None:
            atom.SetAtomMapNum(0)

    _promote_incoming_tags(rw)
    return Chain(_finalise(rw), chain.history + ["heterocycle"])


def append_pks_unit(
    chain: Chain,
    reduction_state: str,
    alpha_substituent: Optional[str] = None,
) -> Chain:
    """Extend the chain by one polyketide unit.

    The new unit contributes its alpha carbon and a new carboxyl (the new
    thioester). The carbonyl carbon of the previous unit becomes the beta
    position and is processed by this module's reductive domains.

    Args:
        chain: The growing chain.
        reduction_state: ``keto``, ``hydroxyl``, ``ene`` or ``alkyl``.
        alpha_substituent: SMILES of the branch on the alpha carbon, e.g.
            ``"C"`` for methylmalonyl-CoA. ``None`` for malonyl-CoA.
    """
    if reduction_state not in ("keto", "hydroxyl", "ene", "alkyl"):
        raise ChemError(f"unknown reduction state: {reduction_state!r}")

    rw = Chem.RWMol(chain.mol)
    beta_c = _tagged_index(rw, TAG_CARBOXYL_C)

    # Build the new two-carbon unit: alpha carbon and the new acid.
    alpha = rw.AddAtom(Chem.Atom(6))
    new_c = rw.AddAtom(Chem.Atom(6))
    new_o_double = rw.AddAtom(Chem.Atom(8))
    new_o_single = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(alpha, new_c, Chem.BondType.SINGLE)
    rw.AddBond(new_c, new_o_double, Chem.BondType.DOUBLE)
    rw.AddBond(new_c, new_o_single, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(new_o_single).SetNumExplicitHs(1)

    # Tag the new acid now, so the tags survive the index shift caused by
    # removing the leaving hydroxyl below.
    rw.GetAtomWithIdx(new_c).SetAtomMapNum(TAG_CARBOXYL_C + 100)
    rw.GetAtomWithIdx(new_o_single).SetAtomMapNum(TAG_CARBOXYL_O + 100)

    if alpha_substituent:
        branch = parse_monomer(alpha_substituent)
        rw_branch_offset = rw.GetNumAtoms()
        rw.InsertMol(branch)
        rw.AddBond(alpha, rw_branch_offset, Chem.BondType.SINGLE)

    # Attach the new unit to the old carbonyl carbon, which is now the beta
    # carbon, and set that carbon's oxidation level.
    bond_order = Chem.BondType.DOUBLE if reduction_state == "ene" else Chem.BondType.SINGLE
    rw.AddBond(beta_c, alpha, bond_order)
    _set_beta_state(rw, beta_c, reduction_state)

    # Hand the active tags to the new acid.
    _tagged_atom(rw, TAG_CARBOXYL_C).SetAtomMapNum(0)
    _drop_leaving_hydroxyl(rw)
    _promote_incoming_tags(rw)
    mol = _finalise(rw)
    return Chain(mol, chain.history + [f"pks:{reduction_state}"])


def _set_beta_state(rw: Chem.RWMol, beta_idx: int, reduction_state: str) -> None:
    """Turn the former thioester carbonyl into keto, hydroxyl, ene or alkyl."""
    beta = rw.GetAtomWithIdx(beta_idx)
    carbonyl_o = None
    for neighbour in beta.GetNeighbors():
        bond = rw.GetBondBetweenAtoms(beta_idx, neighbour.GetIdx())
        if neighbour.GetSymbol() == "O" and bond.GetBondType() == Chem.BondType.DOUBLE:
            carbonyl_o = neighbour.GetIdx()
            break
    if carbonyl_o is None:
        # A nitrogen-bearing beta position (hybrid NRPS/PKS) keeps its amide.
        return

    if reduction_state == "keto":
        return
    if reduction_state == "hydroxyl":
        rw.GetBondBetweenAtoms(beta_idx, carbonyl_o).SetBondType(Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(carbonyl_o).SetNumExplicitHs(1)
        return
    # ene and alkyl both lose the oxygen entirely.
    rw.RemoveAtom(carbonyl_o)


def terminate(chain: Chain, mode: str, lactone_policy: str = "largest") -> Chem.Mol:
    """Release the finished chain from the assembly line.

    Args:
        chain: The completed chain.
        mode: ``hydrolysis`` (free acid), ``reduction`` (aldehyde),
            ``macrolactam`` (head-to-tail amide) or ``macrolactone``
            (ester to the last free hydroxyl on the backbone).
    """
    rw = Chem.RWMol(chain.mol)
    if mode == "hydrolysis":
        pass
    elif mode == "reduction":
        o_idx = _tagged_index(rw, TAG_CARBOXYL_O)
        rw.RemoveAtom(o_idx)
    elif mode in ("macrolactam", "macrolactone"):
        _cyclise(rw, mode, policy=lactone_policy)
    else:
        raise ChemError(f"unknown release mode: {mode!r}")
    mol = _finalise(rw)
    mol = Chem.Mol(mol)
    _clear_tags(mol)
    Chem.SanitizeMol(mol)
    return mol


#: Temporary tag for the atom a macrocyclisation closes onto.
TAG_CYCLISATION_PARTNER = 7


def _cyclise(rw: Chem.RWMol, mode: str, policy: str = "largest") -> None:
    rw.UpdatePropertyCache(strict=False)
    c_idx = _tagged_index(rw, TAG_CARBOXYL_C)
    if mode == "macrolactam":
        partner = _tagged_atom(rw, TAG_NTERM)
        if partner is None:
            raise ChemError("no N-terminus to cyclise onto")
    else:
        partner_idx = _last_free_hydroxyl(rw, exclude=c_idx, policy=policy)
        if partner_idx is None:
            raise ChemError("no free hydroxyl to cyclise onto")
        partner = rw.GetAtomWithIdx(partner_idx)
        partner.SetNumExplicitHs(0)
        partner.SetNoImplicit(True)
    # Tag the partner so the bond can be made after the leaving group is gone
    # and the indices have shifted.
    partner.SetAtomMapNum(TAG_CYCLISATION_PARTNER)
    _drop_leaving_hydroxyl(rw)
    c_idx = _tagged_index(rw, TAG_CARBOXYL_C)
    partner_idx = _tagged_index(rw, TAG_CYCLISATION_PARTNER)
    rw.AddBond(c_idx, partner_idx, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(partner_idx).SetAtomMapNum(0)


#: Macrocycles in MIBiG products run from about 12 to 29 atoms, with a median
#: of 17. A lactone is therefore placed to give a ring in that range rather
#: than onto whichever hydroxyl happens to sit nearest the acid.
TARGET_RING_SIZE = 16


def _free_hydroxyls(rw: Chem.RWMol, exclude: int) -> list:
    rw.UpdatePropertyCache(strict=False)
    found = []
    for atom in rw.GetAtoms():
        if atom.GetSymbol() != "O" or atom.GetTotalNumHs() != 1:
            continue
        if atom.GetAtomMapNum() == TAG_CARBOXYL_O:
            continue
        neighbours = list(atom.GetNeighbors())
        if len(neighbours) != 1 or neighbours[0].GetIdx() == exclude:
            continue
        found.append(atom.GetIdx())
    return found


def _ring_size(rw: Chem.RWMol, start: int, end: int) -> Optional[int]:
    """Atoms in the ring that closing a bond between two atoms would make."""
    from collections import deque

    queue = deque([(start, 1)])
    seen = {start}
    while queue:
        index, distance = queue.popleft()
        if index == end:
            return distance
        for neighbour in rw.GetAtomWithIdx(index).GetNeighbors():
            if neighbour.GetIdx() not in seen:
                seen.add(neighbour.GetIdx())
                queue.append((neighbour.GetIdx(), distance + 1))
    return None


def _last_free_hydroxyl(
    rw: Chem.RWMol, exclude: int, policy: str = "largest", target: int = TARGET_RING_SIZE
) -> Optional[int]:
    """Pick the hydroxyl a macrolactone should close onto.

    ``nearest`` takes the hydroxyl furthest along the atom ordering, which is
    the one closest to the acid and gives the smallest ring. ``largest`` takes
    the one giving the biggest ring. ``target`` aims for a ring of about
    :data:`TARGET_RING_SIZE` atoms, which is where real macrolactones sit.
    """
    candidates = _free_hydroxyls(rw, exclude)
    if not candidates:
        return None
    if policy == "nearest":
        return candidates[-1]
    sized = [(index, _ring_size(rw, exclude, index)) for index in candidates]
    sized = [(index, size) for index, size in sized if size]
    if not sized:
        return candidates[-1]
    if policy == "largest":
        return max(sized, key=lambda item: item[1])[0]
    return min(sized, key=lambda item: abs(item[1] - target))[0]


def to_smiles(mol: Chem.Mol) -> str:
    """Canonical SMILES with the working tags removed."""
    clean = Chem.Mol(mol)
    _clear_tags(clean)
    return Chem.MolToSmiles(clean)
