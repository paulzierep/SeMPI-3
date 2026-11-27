#!/usr/bin/env python
"""Can a cluster be encoded well enough to find its relatives?

Every cluster is turned into a vector (`sempi3.features`) and matched against
all the others. The structure of the closest cluster is then taken as the
prediction, and scored exactly like the assembled backbone.

Two things come out of this:

* a retrieval arm for the benchmark, which any real predictor should beat on
  novel clusters and will lose to on clusters that have a close relative in
  the database;
* a measure of how much signal the cluster encoding carries at all, which is
  what decides whether these features are worth feeding to a model.

Similar clusters in MIBiG are often the same compound family, so the result is
reported twice: over everything, and with the nearest neighbour forced to be
below a similarity cut, which is the harder and more honest setting.

Usage:
    python scripts/experiment_retrieval.py --groups domain,substrate,gene,family,layout
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import (
    Dict,
    List,
)

import numpy as np
from Bio import SeqIO

from sempi3 import (
    antismash_io,
    compare,
    features,
    mibig,
)

JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
OUT = Path("results/retrieval_experiment.json")

GROUPS = ("domain", "substrate", "release", "gene", "family", "layout", "cluster_kind", "cluster", "product")


def collect(limit: int | None = None) -> tuple:
    """Features and known structures for every benchmark cluster."""
    items: List[features.ClusterFeatures] = []
    structures: Dict[str, List[str]] = {}
    names: Dict[str, List[str]] = {}

    paths = sorted(JSON_DIR.glob("*.json"))
    if limit:
        paths = paths[:limit]
    for json_path in paths:
        accession = json_path.stem
        gbk = ANTISMASH_DIR / f"{accession}.gbk"
        if not gbk.exists():
            continue
        entry = mibig.load_entry(json_path)
        if not entry.is_pks_nrps or not entry.structures:
            continue
        try:
            parsed = antismash_io.load(gbk)
            record = next(SeqIO.parse(str(gbk), "genbank"))
        except Exception:  # noqa: BLE001
            continue
        if not parsed.spec.modules:
            continue
        record.products = parsed.products  # type: ignore[attr-defined]
        items.append(features.build(record, parsed.spec, accession))
        structures[accession] = entry.structures
        names[accession] = entry.compound_names
    return items, structures, names


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--groups", type=str, default=",".join(GROUPS))
    parser.add_argument("--cut", type=float, default=0.9, help="similarity above which a neighbour is skipped")
    args = parser.parse_args()

    keep = tuple(g.strip() for g in args.groups.split(",") if g.strip())
    items, structures, names = collect(args.limit)
    print(f"clusters encoded: {len(items)}")

    all_names = features.feature_names(items)
    selected = [n for n in all_names if n.split(":")[0] in keep]
    print(f"features: {len(selected)} of {len(all_names)} (groups: {', '.join(keep)})")

    data = features.matrix(items, selected)
    scores_all: List[float] = []
    scores_far: List[float] = []
    similarity_of_match: List[float] = []

    for index, item in enumerate(items):
        neighbours = features.nearest(data[index], data, exclude=index, top=25)
        if not neighbours:
            continue
        best_index, similarity = neighbours[0]
        similarity_of_match.append(similarity)
        prediction = structures[items[best_index].bgc_id][0]
        score = compare.best_match(
            {"retrieved": prediction}, structures[item.bgc_id], names[item.bgc_id], timeout=3
        )
        if score:
            scores_all.append(score.tanimoto)

        # The harder setting: ignore neighbours that look nearly identical.
        for candidate_index, candidate_similarity in neighbours:
            if candidate_similarity < args.cut:
                far_prediction = structures[items[candidate_index].bgc_id][0]
                far_score = compare.best_match(
                    {"retrieved": far_prediction}, structures[item.bgc_id], names[item.bgc_id], timeout=3
                )
                if far_score:
                    scores_far.append(far_score.tanimoto)
                break

    report = {
        "clusters": len(items),
        "features": len(selected),
        "groups": list(keep),
        "mean_tanimoto_nearest": round(float(np.mean(scores_all)), 4) if scores_all else 0.0,
        "median_tanimoto_nearest": round(float(np.median(scores_all)), 4) if scores_all else 0.0,
        "mean_tanimoto_below_cut": round(float(np.mean(scores_far)), 4) if scores_far else 0.0,
        "cut": args.cut,
        "mean_similarity_of_match": round(float(np.mean(similarity_of_match)), 4)
        if similarity_of_match
        else 0.0,
        "scored": len(scores_all),
    }
    print(json.dumps(report, indent=2))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
