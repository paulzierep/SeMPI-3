"""Encoding a gene cluster as a vector.

The assembly line is only part of what a cluster says about its product. The
annotation also carries the enzymes around the line, how the genes are laid
out, and what kind of cluster antiSMASH thinks it is. This module turns all of
that into one numeric vector, so clusters can be compared to each other and
used by a model.

The groups of features are:

``domains``
    How many of each assembly-line domain the cluster has.
``substrates``
    How many modules take each substrate.
``genes``
    Counts by antiSMASH's gene role: biosynthetic, regulatory, transport,
    resistance.
``families``
    Counts of the enzyme families in the ``sec_met_domain`` annotation, which
    is where tailoring enzymes such as P450s and fatty-acid ligases show up.
``layout``
    Size of the cluster, how many genes, how tightly they are packed, how often
    the strand changes: a crude read on operon structure.
"""

from __future__ import annotations

import math
import re
from dataclasses import (
    dataclass,
    field,
)
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
)

import numpy as np

from .spec import (
    ClusterSpec,
    NRPS,
)

#: Gene roles antiSMASH assigns.
GENE_KINDS = ("biosynthetic", "biosynthetic-additional", "regulatory", "transport", "resistance", "other")

#: Candidate-cluster kinds.
CLUSTER_KINDS = ("single", "neighbouring", "chemical_hybrid", "interleaved")


@dataclass
class ClusterFeatures:
    """One cluster's features, kept named so they can be inspected."""

    bgc_id: str
    values: Dict[str, float] = field(default_factory=dict)

    def vector(self, names: Sequence[str]) -> np.ndarray:
        return np.array([self.values.get(name, 0.0) for name in names], dtype=float)


def _domain_counts(spec: ClusterSpec) -> Dict[str, float]:
    counts: Dict[str, float] = {
        "domain:modules": float(spec.n_modules),
        "domain:nrps_modules": float(sum(1 for m in spec.modules if m.kind == NRPS)),
        "domain:pks_modules": float(sum(1 for m in spec.modules if m.kind != NRPS)),
        "domain:kr": float(sum(1 for m in spec.modules if m.kr)),
        "domain:dh": float(sum(1 for m in spec.modules if m.dh)),
        "domain:er": float(sum(1 for m in spec.modules if m.er)),
        "domain:epimerase": float(sum(1 for m in spec.modules if m.epimerization)),
        "domain:methyltransferase": float(sum(1 for m in spec.modules if m.n_methylation)),
        "domain:heterocyclization": float(sum(1 for m in spec.modules if m.heterocyclization)),
        "domain:non_elongating": float(sum(1 for m in spec.modules if m.non_elongating)),
        "domain:starter": float(sum(1 for m in spec.modules if m.is_starter)),
        "domain:unassigned_substrate": float(
            sum(1 for m in spec.modules if m.kind == NRPS and not m.substrate)
        ),
    }
    if spec.release:
        counts[f"release:{spec.release}"] = 1.0
    return counts


def _substrate_counts(spec: ClusterSpec) -> Dict[str, float]:
    counts: Dict[str, float] = {}
    for module in spec.modules:
        if not module.substrate:
            continue
        name = module.substrate.replace("smiles:", "")[:40].lower()
        key = f"substrate:{name}"
        counts[key] = counts.get(key, 0.0) + 1.0
    return counts


def _gene_and_family_counts(record) -> Dict[str, float]:
    counts: Dict[str, float] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        for kind in feature.qualifiers.get("gene_kind", []):
            counts[f"gene:{kind}"] = counts.get(f"gene:{kind}", 0.0) + 1.0
        for entry in feature.qualifiers.get("sec_met_domain", []):
            family = entry.split("(")[0].strip().split()[0]
            key = f"family:{family}"
            counts[key] = counts.get(key, 0.0) + 1.0
    return counts


def _layout_features(record) -> Dict[str, float]:
    """Size, gene density, packing and strand structure."""
    coding = [f for f in record.features if f.type == "CDS"]
    coding.sort(key=lambda f: int(f.location.start))
    length = max(len(record.seq), 1)

    values: Dict[str, float] = {
        "layout:log_length": math.log10(length),
        "layout:n_genes": float(len(coding)),
        "layout:genes_per_kb": len(coding) / (length / 1000.0),
    }
    if len(coding) < 2:
        return values

    gaps = [int(b.location.start) - int(a.location.end) for a, b in zip(coding, coding[1:])]
    strands = [int(f.location.strand or 0) for f in coding]
    switches = sum(1 for a, b in zip(strands, strands[1:]) if a != b)
    values.update(
        {
            "layout:median_gap": float(np.median(gaps)),
            "layout:overlapping_share": float(sum(1 for g in gaps if g < 0) / len(gaps)),
            # Genes closer than 50 bp are usually in one operon.
            "layout:operon_share": float(sum(1 for g in gaps if g < 50) / len(gaps)),
            "layout:strand_switch_rate": switches / max(len(strands) - 1, 1),
            "layout:coding_fraction": float(
                sum(int(f.location.end) - int(f.location.start) for f in coding) / length
            ),
        }
    )
    for feature in record.features:
        if feature.type == "cand_cluster":
            for kind in feature.qualifiers.get("kind", []):
                if kind in CLUSTER_KINDS:
                    values[f"cluster_kind:{kind}"] = 1.0
            for edge in feature.qualifiers.get("contig_edge", []):
                values["cluster:contig_edge"] = 1.0 if edge == "True" else 0.0
    return values


def build(record, spec: ClusterSpec, bgc_id: str) -> ClusterFeatures:
    """Encode one cluster."""
    values: Dict[str, float] = {}
    values.update(_domain_counts(spec))
    values.update(_substrate_counts(spec))
    values.update(_gene_and_family_counts(record))
    values.update(_layout_features(record))
    for product in getattr(record, "products", []) or []:
        values[f"product:{product}"] = 1.0
    return ClusterFeatures(bgc_id=bgc_id, values=values)


def feature_names(features: Sequence[ClusterFeatures], min_count: int = 3) -> List[str]:
    """Names shared by enough clusters to be worth keeping."""
    counts: Dict[str, int] = {}
    for item in features:
        for name in item.values:
            counts[name] = counts.get(name, 0) + 1
    return sorted(name for name, count in counts.items() if count >= min_count)


def matrix(features: Sequence[ClusterFeatures], names: Sequence[str], normalise: bool = True) -> np.ndarray:
    """Stack the features into a matrix, optionally scaled to unit length.

    Scaling matters for similarity: without it a big cluster would look like
    every other big cluster.
    """
    data = np.vstack([item.vector(names) for item in features])
    if not normalise:
        return data
    # Counts are compressed so that one very common family cannot dominate.
    data = np.log1p(np.clip(data, 0, None))
    norms = np.linalg.norm(data, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return data / norms


def nearest(
    query: np.ndarray,
    reference: np.ndarray,
    exclude: Optional[int] = None,
    top: int = 1,
) -> List[Tuple[int, float]]:
    """Indices and cosine similarities of the closest reference clusters."""
    scores = reference @ query
    if exclude is not None and 0 <= exclude < len(scores):
        scores[exclude] = -np.inf
    order = np.argsort(scores)[::-1][:top]
    return [(int(index), float(scores[index])) for index in order]
