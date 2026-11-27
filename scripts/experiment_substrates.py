#!/usr/bin/env python
"""Does each extra ingredient improve A-domain substrate prediction?

Every variant is evaluated the same way: five folds, grouped by gene cluster,
so no cluster is ever in both the training and the test side. The test set is
always the MIBiG-derived domains, so the numbers are comparable down the
table, and antiSMASH's own call is scored on exactly the same domains.

Added training data is filtered against the test fold by gene accession, so a
cluster under test never contributes training rows through another database.

Usage:
    python scripts/experiment_substrates.py [--folds 5]
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)

import numpy as np

from sempi3 import (
    datasets,
    stachelhaus as st,
    substrates,
)

PARAS_TSV = Path("data/external/paras/dataset_paras.tsv")
PRIETO = Path("data/interim/prieto_adomains.jsonl")
ADOMAINS = Path("data/interim/adomains.jsonl")


def load_prieto() -> List[dict]:
    """A domains from the SeMPI 2.0 supplement, via Prieto et al."""
    if not PRIETO.exists():
        return []
    rows = []
    for line in PRIETO.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        label = substrates.normalise(row.get("substrate"))
        if not label or len(row.get("signature", "")) != len(st.EXTENDED_POSITIONS):
            continue
        rows.append(
            {
                "accessions": [row["accession"]],
                "signature": row["signature"],
                "code": row["code"],
                "substrate": label,
            }
        )
    return rows


def load_paras() -> List[dict]:
    """PARAS rows as {accession, index, signature, code, substrate}."""
    if not PARAS_TSV.exists():
        return []
    rows = []
    with open(PARAS_TSV) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            signature = row["extended_signature"].strip().upper()
            label = substrates.normalise(row["substrate_specificity"])
            if not label or len(signature) != len(st.EXTENDED_POSITIONS):
                continue
            accessions = []
            for part in row["domain_id"].split("|"):
                accession = part.rsplit(".A", 1)[0] if ".A" in part else part
                accessions.append(accession)
            rows.append(
                {
                    "accessions": accessions,
                    "signature": signature,
                    "code": st.code_from_signature(signature),
                    "substrate": label,
                }
            )
    return rows


def nearest_identity_predict(
    train_features: np.ndarray,
    train_labels: np.ndarray,
    test_features: np.ndarray,
    min_identity: float,
) -> List[Optional[str]]:
    """Nearest-signature vote, abstaining below ``min_identity``."""
    table: Dict[str, List[str]] = {}
    for feature, label in zip(train_features, train_labels):
        table.setdefault(feature, []).append(label)
    known = list(table.items())

    predictions: List[Optional[str]] = []
    for feature in test_features:
        exact = table.get(feature)
        if exact:
            predictions.append(Counter(exact).most_common(1)[0][0])
            continue
        best_score = 0.0
        best_labels: List[str] = []
        for candidate, labels in known:
            matches = sum(1 for a, b in zip(feature, candidate) if a == b and a != st.GAP)
            score = matches / len(feature)
            if score > best_score:
                best_score, best_labels = score, list(labels)
            elif score == best_score:
                best_labels.extend(labels)
        if best_score < min_identity or not best_labels:
            predictions.append(None)
        else:
            predictions.append(Counter(best_labels).most_common(1)[0][0])
    return predictions


def forest_predict(
    train_features: np.ndarray,
    train_labels: np.ndarray,
    test_features: np.ndarray,
) -> List[str]:
    """Random forest over one-hot encoded signature columns.

    This is the encoding SeMPI 2.0 and PARAS both use: every position in the
    signature is a categorical feature.
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import OneHotEncoder

    encoder = OneHotEncoder(handle_unknown="ignore")
    train_matrix = encoder.fit_transform([list(f) for f in train_features])
    test_matrix = encoder.transform([list(f) for f in test_features])
    forest = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=0)
    forest.fit(train_matrix, train_labels)
    return list(forest.predict(test_matrix))


def evaluate(
    records: List[dict],
    feature_key: str,
    folds: int,
    paras_rows: Optional[List[dict]] = None,
    min_identity: float = 0.0,
    antismash_fallback: bool = False,
    seed: int = 0,
    model: str = "nearest",
    antismash_first: bool = False,
    strict: bool = False,
) -> Tuple[int, int, int, int]:
    """Return (n, correct, abstained, extra_training_rows_used)."""
    features = np.array([r[feature_key] for r in records])
    labels = np.array([r["label"] for r in records])
    groups = np.array([r["bgc_id"] for r in records])
    calls = np.array([r["antismash_call"] or "" for r in records])

    unique = np.array(sorted(set(groups)))
    rng = np.random.default_rng(seed)
    rng.shuffle(unique)
    fold_of = np.array([{g: i % folds for i, g in enumerate(unique)}[g] for g in groups])

    correct = abstained = extra_used = 0
    for fold in range(folds):
        test = fold_of == fold
        train = ~test
        train_features = list(features[train])
        train_labels = list(labels[train])

        if paras_rows:
            # Genes of the clusters under test, so added data cannot leak them.
            blocked = set()
            blocked_signatures = set()
            for record in np.array(records)[test]:
                blocked.add(record["gene"])
                if strict:
                    blocked_signatures.add(record["signature"])
            for row in paras_rows:
                if any(accession in blocked for accession in row["accessions"]):
                    continue
                if strict and row["signature"] in blocked_signatures:
                    continue
                train_features.append(row["signature"] if feature_key == "signature" else row["code"])
                train_labels.append(row["substrate"])
                extra_used += 1

        if model == "forest":
            predictions = forest_predict(
                np.array(train_features), np.array(train_labels), features[test]
            )
        else:
            predictions = nearest_identity_predict(
                np.array(train_features), np.array(train_labels), features[test], min_identity
            )
        for prediction, truth, call in zip(predictions, labels[test], calls[test]):
            # antiSMASH leads where it commits; the model fills its gaps.
            if antismash_first and call:
                prediction = call
            elif prediction is None and antismash_fallback and call:
                prediction = call
            if prediction is None:
                abstained += 1
            elif prediction == truth:
                correct += 1
    return len(records), correct, abstained, extra_used // max(folds, 1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    raw = datasets.read_jsonl(ADOMAINS)
    records = []
    for row in raw:
        label = substrates.normalise(row.get("substrate"))
        sequence = row.get("sequence") or ""
        if not label or len(sequence) < 50:
            continue
        records.append(
            {
                "bgc_id": row["bgc_id"],
                "gene": row["gene"],
                "label": label,
                "code": st.extract_code(sequence),
                "signature": st.extract_signature(sequence),
                "antismash_call": substrates.normalise(row.get("antismash_call")),
            }
        )
    paras = load_paras()
    print(f"MIBiG domains: {len(records)}   PARAS rows: {len(paras)}\n")

    n = len(records)
    antismash_correct = sum(1 for r in records if r["antismash_call"] == r["label"])
    antismash_abstained = sum(1 for r in records if not r["antismash_call"])
    results = [
        {
            "variant": "antiSMASH (baseline)",
            "n": n,
            "correct": antismash_correct,
            "abstained": antismash_abstained,
        }
    ]

    prieto = load_prieto()
    both = paras + prieto
    print(f"Prieto rows: {len(prieto)}\n")
    variants = [
        ("Stachelhaus 9-aa, MIBiG only", "code", None, "nearest", False, False),
        ("Extended 34-aa, MIBiG only", "signature", None, "nearest", False, False),
        ("Extended 34-aa, + PARAS", "signature", paras, "nearest", False, False),
        ("Extended 34-aa, + PARAS + Prieto", "signature", both, "nearest", False, False),
        ("Random forest, MIBiG only", "signature", None, "forest", False, False),
        ("Random forest, + PARAS", "signature", paras, "forest", False, False),
        ("Random forest, + PARAS + Prieto", "signature", both, "forest", False, False),
        ("Random forest, + both, no identical signatures", "signature", both, "forest", False, True),
        ("antiSMASH leads, nearest fills gaps (+both)", "signature", both, "nearest", True, False),
        ("antiSMASH leads, forest fills gaps (+both)", "signature", both, "forest", True, False),
        ("antiSMASH leads, forest fills gaps (strict)", "signature", both, "forest", True, True),
    ]
    for name, key, extra, model, as_first, strict in variants:
        total, correct, abstained, extra_used = evaluate(
            records, key, args.folds, extra, 0.0, False, model=model, antismash_first=as_first, strict=strict
        )
        results.append(
            {
                "variant": name,
                "n": total,
                "correct": correct,
                "abstained": abstained,
                "extra_training_rows": extra_used,
            }
        )

    print(f"{'variant':46s} {'accuracy':>9s} {'when called':>12s} {'abstained':>10s}")
    for row in results:
        accuracy = row["correct"] / row["n"]
        called = row["n"] - row["abstained"]
        when_called = row["correct"] / called if called else 0.0
        row["accuracy"] = round(accuracy, 4)
        row["accuracy_when_called"] = round(when_called, 4)
        print(f"{row['variant']:46s} {accuracy:9.3f} {when_called:12.3f} {row['abstained']:10d}")

    out = Path("results/substrate_experiment.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
