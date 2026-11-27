#!/usr/bin/env python
"""Which properties of a product can be read off its gene cluster?

Rather than asking a model for the whole molecule, this asks it for one
property at a time: does the product carry a sugar, how many nitrogens does
it hold, is there a macrocycle, roughly how large is it. Each is a separate
question with its own answer and its own measurable accuracy.

That matters for two reasons. A property the cluster predicts well is one
SeMPI-3 can use to choose between the structures it has assembled. A property
it predicts badly marks a limit worth knowing about rather than guessing at.

Everything is out of fold. Binary targets are scored against always answering
with the majority, counts against always answering with the mean.

Usage:
    python scripts/experiment_structure_features.py --folds 5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
)

import numpy as np
from Bio import SeqIO
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
)

from sempi3 import (
    antismash_io,
    features,
    mibig,
    structure_features as sf,
)

JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
OUT = Path("results/structure_features_experiment.json")

#: Feature groups, so each can be switched off to see what it was worth.
GROUPS = ("domain", "substrate", "release", "gene", "family", "layout", "cluster_kind", "cluster", "product")


def collect(limit: Optional[int] = None):
    """Cluster features paired with the measured properties of the product."""
    items: List[features.ClusterFeatures] = []
    profiles: List[sf.StructureProfile] = []

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
        profile = sf.describe(entry.structures[0])
        if profile is None:
            continue
        record.products = parsed.products  # type: ignore[attr-defined]
        items.append(features.build(record, parsed.spec, accession))
        profiles.append(profile)
    return items, profiles


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--without",
        type=str,
        default="",
        help="Comma-separated feature groups to leave out, to see what they were worth.",
    )
    args = parser.parse_args()

    items, profiles = collect(args.limit)
    print(f"clusters: {len(items)}")

    drop = {g.strip() for g in args.without.split(",") if g.strip()}
    names = [n for n in features.feature_names(items) if n.split(":")[0] not in drop]
    matrix = features.matrix(items, names, normalise=False)
    print(f"features: {len(names)}" + (f" (without {', '.join(sorted(drop))})" if drop else ""))

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(items))
    fold_of = np.empty(len(items), dtype=int)
    for position, index in enumerate(order):
        fold_of[index] = position % args.folds

    report: Dict[str, Dict[str, float]] = {}
    for target in sf.TARGETS:
        truth = np.array([p.get(target) for p in profiles], dtype=float)
        binary = target in sf.BINARY_TARGETS
        predicted = np.zeros(len(truth), dtype=float)

        for fold in range(args.folds):
            test = fold_of == fold
            train = ~test
            if test.sum() == 0 or train.sum() == 0:
                continue
            if binary:
                model = RandomForestClassifier(
                    n_estimators=300,
                    min_samples_leaf=2,
                    class_weight="balanced",
                    random_state=args.seed,
                    n_jobs=-1,
                )
            else:
                model = RandomForestRegressor(
                    n_estimators=300, min_samples_leaf=2, random_state=args.seed, n_jobs=-1
                )
            model.fit(matrix[train], truth[train])
            predicted[test] = model.predict(matrix[test])

        if binary:
            accuracy = float((predicted.round() == truth).mean())
            majority = float(max(truth.mean(), 1 - truth.mean()))
            report[target] = {
                "kind": "binary",
                "accuracy": round(accuracy, 4),
                "majority": round(majority, 4),
                "gain": round(accuracy - majority, 4),
                "positive_rate": round(float(truth.mean()), 4),
            }
        else:
            residual = float(((truth - predicted) ** 2).mean())
            variance = float(((truth - truth.mean()) ** 2).mean())
            r2 = 1.0 - residual / variance if variance else 0.0
            correlation = float(np.corrcoef(truth, predicted)[0, 1]) if variance else 0.0
            report[target] = {
                "kind": "count",
                "r2": round(r2, 4),
                "correlation": round(correlation, 4),
                "mean": round(float(truth.mean()), 2),
                "mean_absolute_error": round(float(np.abs(truth - predicted).mean()), 3),
            }

    print(f"\n{'property':22s} {'kind':7s} {'score':>8s} {'baseline':>9s} {'gain':>8s}")
    for target, values in sorted(report.items(), key=lambda kv: -(kv[1].get("gain", kv[1].get("r2", 0)))):
        if values["kind"] == "binary":
            print(
                f"{target:22s} {'binary':7s} {values['accuracy']:8.3f} "
                f"{values['majority']:9.3f} {values['gain']:+8.3f}"
            )
        else:
            print(
                f"{target:22s} {'count':7s} {values['r2']:8.3f} "
                f"{'0.000':>9s} {values['r2']:+8.3f}   (r={values['correlation']:.2f}, "
                f"mean {values['mean']}, error {values['mean_absolute_error']})"
            )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"features": len(names), "dropped": sorted(drop), "targets": report}, indent=2) + "\n")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
