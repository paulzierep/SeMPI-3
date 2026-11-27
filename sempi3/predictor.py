"""Predicting a structure from one antiSMASH GenBank file.

This is the entry point for using SeMPI-3 downstream of antiSMASH on a real
genome: give it an annotated cluster, get back a structure, the assembly line
it came from, and antiSMASH's own prediction for comparison.

The substrate reference shipped in ``sempi3/data/substrate_signatures.tsv``
holds one binding-pocket signature per row, from MIBiG and PARAS. It is only
consulted for domains antiSMASH left unnamed.
"""

from __future__ import annotations

from dataclasses import (
    dataclass,
    field,
)
from functools import lru_cache
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)

from . import (
    antismash_io,
    assembly,
    stachelhaus as st,
    tailoring,
)
from .spec import NRPS

SIGNATURE_TABLE = Path(__file__).parent / "data" / "substrate_signatures.tsv"

#: How much of the pocket has to agree before a substrate is believed.
MIN_IDENTITY = 0.6


@lru_cache(maxsize=1)
def _reference() -> List[Tuple[str, str, int, str]]:
    """Signatures with their substrate, support, and that substrate's SMILES."""
    rows: List[Tuple[str, str, int, str]] = []
    if not SIGNATURE_TABLE.exists():  # pragma: no cover - packaged with the module
        return rows
    with open(SIGNATURE_TABLE) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append((parts[0], parts[1], int(parts[2]), parts[3] if len(parts) > 3 else ""))
    return rows


@lru_cache(maxsize=1)
def _structures_by_name() -> Dict[str, str]:
    """Substrate name -> SMILES, from the packaged reference."""
    table: Dict[str, str] = {}
    for _signature, label, _support, smiles in _reference():
        if smiles and label not in table:
            table[label] = smiles
    return table


def structure_for(name: str) -> Optional[str]:
    """The structure of a named substrate, from the monomer table or the reference.

    A predicted substrate with no structure cannot be built into the chain, so
    both sources are consulted before giving up.
    """
    from . import monomers, substrates as substrate_names

    if not name:
        return None
    direct = monomers.resolve(name)
    if direct:
        return direct
    normalised = substrate_names.normalise(name) or name.strip().lower()
    return _structures_by_name().get(normalised)


def predict_substrate(sequence: str) -> Tuple[Optional[str], float, str]:
    """Name the substrate of one adenylation domain, with its confidence.

    Confidence is the share of pocket residues shared with the closest known
    domain, so an exact pocket gives 1.0. The structure is returned alongside
    the name, since a name without one cannot be built into the chain.
    """
    signature = st.extract_signature(sequence)
    if signature.count(st.GAP) > 3:
        return None, 0.0, ""
    best_label: Optional[str] = None
    best_score = 0.0
    best_support = 0
    best_smiles = ""
    for known, label, support, smiles in _reference():
        matches = sum(1 for a, b in zip(signature, known) if a == b and a != st.GAP)
        score = matches / len(signature)
        if score > best_score or (score == best_score and support > best_support):
            best_score, best_label, best_support, best_smiles = score, label, support, smiles
    if best_score < MIN_IDENTITY:
        return None, best_score, ""
    return best_label, best_score, best_smiles


@dataclass
class ClusterPrediction:
    """What SeMPI-3 says about one cluster."""

    bgc_id: str
    smiles: Optional[str]
    release: Optional[str]
    products: List[str] = field(default_factory=list)
    modules: List[Dict[str, object]] = field(default_factory=list)
    variants: Dict[str, str] = field(default_factory=dict)
    tailoring_genes: Dict[str, int] = field(default_factory=dict)
    antismash_smiles: Optional[str] = None
    filled_substrates: int = 0
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        return {
            "bgc_id": self.bgc_id,
            "smiles": self.smiles,
            "release": self.release,
            "products": self.products,
            "modules": self.modules,
            "variants": self.variants,
            "tailoring_genes": self.tailoring_genes,
            "antismash_smiles": self.antismash_smiles,
            "filled_substrates": self.filled_substrates,
            "notes": self.notes,
        }


def predict_all(
    gbk_path: Path,
    fill_substrates: bool = True,
    apply_tailoring: bool = True,
) -> List[ClusterPrediction]:
    """Predict a product for every region in an antiSMASH GenBank file.

    A genome-wide antiSMASH run puts many regions in one file; each is an
    assembly line of its own and gets its own prediction.
    """
    return [
        _predict_record(record, fill_substrates=fill_substrates, apply_tailoring=apply_tailoring)
        for record in antismash_io.load_regions(Path(gbk_path))
    ]


def predict(
    gbk_path: Path,
    fill_substrates: bool = True,
    apply_tailoring: bool = True,
) -> ClusterPrediction:
    """Predict the product of a single-cluster antiSMASH file.

    Args:
        gbk_path: A GenBank file written by antiSMASH.
        fill_substrates: Name the adenylation domains antiSMASH left unnamed.
        apply_tailoring: Decorate the backbone from the cluster's enzymes.
    """
    record = antismash_io.load(Path(gbk_path))
    return _predict_record(record, fill_substrates=fill_substrates, apply_tailoring=apply_tailoring)


def _predict_record(record, fill_substrates: bool, apply_tailoring: bool) -> ClusterPrediction:
    spec = record.spec

    filled = 0
    if fill_substrates:
        for module in spec.modules:
            if module.kind != NRPS or module.substrate or not module.sequence:
                continue
            substrate, confidence, smiles = predict_substrate(module.sequence)
            if substrate:
                # Prefer the structure: a name the monomer table does not know
                # would be dropped when the chain is built.
                module.substrate = f"smiles:{smiles}" if smiles else substrate
                module.confidence = confidence
                filled += 1

    prediction = assembly.assemble(spec, skip_unknown_substrates=True)
    smiles = prediction.smiles
    if apply_tailoring and smiles and record.tailoring_genes:
        genes = tailoring.TailoringGenes(counts=record.tailoring_genes)
        decorated = tailoring.apply(smiles, genes, enabled=list(tailoring.DEFAULT_REACTIONS))
        if decorated:
            smiles = decorated

    modules = []
    for module in spec.modules:
        modules.append(
            {
                "order": module.order,
                "kind": module.kind,
                "gene": module.gene,
                "substrate": (module.substrate or "").replace("smiles:", "") or None,
                "confidence": round(module.confidence, 3) if module.confidence else None,
                "reduction_state": module.reduction_state if module.kind != NRPS else None,
                "epimerization": module.epimerization,
                "heterocyclization": module.heterocyclization,
            }
        )

    return ClusterPrediction(
        bgc_id=record.bgc_id,
        smiles=smiles,
        release=prediction.release,
        products=record.products,
        modules=modules,
        variants=prediction.variants,
        tailoring_genes=record.tailoring_genes,
        antismash_smiles=record.baseline,
        filled_substrates=filled,
        notes=prediction.notes,
    )
