#!/usr/bin/env python
"""Can a model pick the chain-release mechanism better than the domain rule?

The label is read off the known structure: a macrolactone, a macrolactam, or
neither. The features are what a user has for an unknown cluster: the module
composition and the tailoring domains antiSMASH found.

Predictions are made out of fold, and written to a lookup the benchmark can
use, so the structure benchmark never scores a release that was chosen with
knowledge of that cluster.

Usage:
    python scripts/experiment_release.py --folds 5
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import (
    Dict,
    List,
    Tuple,
)

import numpy as np
from sklearn.ensemble import RandomForestClassifier

from sempi3 import datasets

RELEASE = Path("data/interim/release.jsonl")
OUT = Path("results/release_prediction.json")

#: Feature columns, in a fixed order.
NUMERIC = (
    "n_modules",
    "n_nrps",
    "n_pks",
    "n_kr",
    "n_dh",
    "n_er",
    "n_epimerase",
    "n_methyltransferase",
)


def featurise(row: dict) -> List[float]:
    classes = {c.upper() for c in row.get("classes", [])}
    features = [float(row[name]) for name in NUMERIC]
    features += [
        1.0 if row["has_thioesterase"] else 0.0,
        1.0 if row["has_reductase"] else 0.0,
        1.0 if "NRPS" in classes else 0.0,
        1.0 if "PKS" in classes else 0.0,
        1.0 if classes == {"NRPS"} else 0.0,
        1.0 if classes == {"PKS"} else 0.0,
        float(row["n_nrps"]) / max(row["n_modules"], 1),
    ]
    return features


def domain_rule(row: dict) -> str:
    """The rule SeMPI-3 uses today, for comparison."""
    if row["has_reductase"]:
        return "linear"
    if row["n_pks"] == 0:
        return "macrolactam"
    return "macrolactone"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rows = datasets.read_jsonl(RELEASE)
    print(f"clusters: {len(rows)}   labels: {Counter(r['label'] for r in rows).most_common()}")

    features = np.array([featurise(r) for r in rows])
    labels = np.array([r["label"] for r in rows])
    identifiers = [r["bgc_id"] for r in rows]

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(rows))
    fold_of = np.empty(len(rows), dtype=int)
    for position, index in enumerate(order):
        fold_of[index] = position % args.folds

    predictions: Dict[str, str] = {}
    correct = 0
    for fold in range(args.folds):
        test = fold_of == fold
        train = ~test
        forest = RandomForestClassifier(
            n_estimators=400, min_samples_leaf=2, class_weight="balanced", random_state=args.seed, n_jobs=-1
        )
        forest.fit(features[train], labels[train])
        predicted = forest.predict(features[test])
        for identifier, prediction, truth in zip(np.array(identifiers)[test], predicted, labels[test]):
            predictions[identifier] = prediction
            correct += int(prediction == truth)

    rule_correct = sum(1 for row in rows if domain_rule(row) == row["label"])
    majority = Counter(labels).most_common(1)[0][0]
    majority_correct = sum(1 for label in labels if label == majority)

    report = {
        "clusters": len(rows),
        "majority_class": {"label": majority, "accuracy": round(majority_correct / len(rows), 4)},
        "domain_rule": {"accuracy": round(rule_correct / len(rows), 4)},
        "random_forest": {"accuracy": round(correct / len(rows), 4)},
        "predictions": predictions,
    }
    print(f"  always '{majority}':   {report['majority_class']['accuracy']:.3f}")
    print(f"  current domain rule: {report['domain_rule']['accuracy']:.3f}")
    print(f"  random forest (OOF): {report['random_forest']['accuracy']:.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=0, sort_keys=True) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
