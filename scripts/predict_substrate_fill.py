#!/usr/bin/env python
"""Predict substrates for the A domains antiSMASH leaves unassigned.

antiSMASH reports ``X`` for roughly one adenylation domain in six. Each one
costs SeMPI-3 a residue of the backbone. This script fills those gaps with the
binding-pocket model, and writes the result as a lookup the benchmark can use.

Predictions are made **out of fold**: clusters are split into groups, and the
model that predicts a cluster's domains is never trained on that cluster, nor
on rows from other databases that belong to its genes. Without that, the
benchmark would be scoring a model that has already seen the answer.

Usage:
    python scripts/predict_substrate_fill.py --folds 5
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import (
    Dict,
    List,
)

import numpy as np

from sempi3 import (
    antismash_io,
    datasets,
    hmm_align,
    mibig,
    stachelhaus as st,
    substrates,
)


def _signature(sequence: str) -> str:
    """Pocket signature, by profile alignment when HMMER is available."""
    found, _method = hmm_align.signature(sequence)
    return found

ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ADOMAINS = Path("data/interim/adomains.jsonl")
OUT = Path("results/substrate_fill.json")


def _load_labelled() -> List[dict]:
    """MIBiG-labelled domains plus the PARAS dataset, in one shape."""
    rows: List[dict] = []
    for raw in datasets.read_jsonl(ADOMAINS):
        label = substrates.normalise(raw.get("substrate"))
        sequence = raw.get("sequence") or ""
        if not label or len(sequence) < 50:
            continue
        rows.append(
            {
                "signature": _signature(sequence),
                "label": label,
                "bgc_id": raw["bgc_id"],
                "genes": {raw["gene"]},
            }
        )

    import csv

    paras_path = Path("data/external/paras/dataset_paras.tsv")
    if paras_path.exists():
        with open(paras_path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                signature = row["extended_signature"].strip().upper()
                label = substrates.normalise(row["substrate_specificity"])
                if not label or len(signature) != len(st.EXTENDED_POSITIONS):
                    continue
                genes = {
                    part.rsplit(".A", 1)[0] if ".A" in part else part
                    for part in row["domain_id"].split("|")
                }
                rows.append({"signature": signature, "label": label, "bgc_id": None, "genes": genes})
    return rows


def _predict(table: Dict[str, List[str]], signature: str) -> str | None:
    exact = table.get(signature)
    if exact:
        return Counter(exact).most_common(1)[0][0]
    best_score = 0.0
    best: List[str] = []
    for candidate, labels in table.items():
        matches = sum(1 for a, b in zip(signature, candidate) if a == b and a != st.GAP)
        score = matches / len(signature)
        if score > best_score:
            best_score, best = score, list(labels)
        elif score == best_score:
            best.extend(labels)
    return Counter(best).most_common(1)[0][0] if best else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    # Clusters in scope: PKS/NRPS with a known structure and antiSMASH output.
    targets = []
    for json_path in sorted(JSON_DIR.glob("*.json")):
        accession = json_path.stem
        if not (ANTISMASH_DIR / f"{accession}.gbk").exists():
            continue
        entry = mibig.load_entry(json_path)
        if entry.is_pks_nrps and entry.structures:
            targets.append(accession)
    print(f"clusters in scope: {len(targets)}")

    labelled = _load_labelled()
    print(f"labelled domains available: {len(labelled)}")

    rng = np.random.default_rng(args.seed)
    order = np.array(targets)
    rng.shuffle(order)
    fold_of = {accession: index % args.folds for index, accession in enumerate(order)}

    filled: Dict[str, Dict[str, str]] = {}
    gaps = 0
    for fold in range(args.folds):
        members = [a for a in targets if fold_of[a] == fold]
        # Collect the gaps in this fold and the genes they sit on.
        pending = []
        blocked_genes = set()
        for accession in members:
            try:
                record = antismash_io.load(ANTISMASH_DIR / f"{accession}.gbk")
            except Exception:  # noqa: BLE001
                continue
            for module in record.spec.modules:
                if module.gene:
                    blocked_genes.add(module.gene)
                if module.kind == "nrps" and not module.substrate and module.sequence:
                    pending.append((accession, module.gene, module.domain_index, module.sequence))
        if not pending:
            continue

        table: Dict[str, List[str]] = {}
        for row in labelled:
            if row["bgc_id"] in members:
                continue
            if row["genes"] & blocked_genes:
                continue
            if row["signature"].count(st.GAP) > 3:
                continue
            table.setdefault(row["signature"], []).append(row["label"])

        for accession, gene, index, sequence in pending:
            gaps += 1
            signature = _signature(sequence)
            if signature.count(st.GAP) > 3:
                continue
            prediction = _predict(table, signature)
            if prediction:
                filled.setdefault(accession, {})[f"{gene}:{index}"] = prediction
        print(f"  fold {fold}: {len(pending)} gaps, training table {len(table)} signatures")

    total_filled = sum(len(v) for v in filled.values())
    print(f"gaps found: {gaps}; filled: {total_filled} across {len(filled)} clusters")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(filled, indent=0, sort_keys=True) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
