#!/usr/bin/env python
"""Choose between assembled structures using predicted product properties.

SeMPI-3 can build several structures for one cluster: a chain released four
different ways, with or without each decoration. Until now one was chosen by
a rule. Here the choice is made by asking the cluster what its product should
look like — how heavy, how many nitrogens, how many rings, a sugar or not —
and keeping the candidate that fits that description best.

The chemistry still comes from the assembly line. The model only ranks what
the chemistry produced, which is why every candidate is a valid molecule.

Property models are trained out of fold, so a cluster is never ranked by a
model that has seen it.

Usage:
    python scripts/experiment_rank.py --folds 5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)

import numpy as np
from Bio import SeqIO
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
)

from sempi3 import (
    antismash_io,
    assembly,
    compare,
    features,
    mibig,
    structure_features as sf,
    tailoring,
)

JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
FILL = Path("results/substrate_fill.json")
OUT = Path("results/rank_experiment.json")

#: Properties used for ranking, with the weight each carries.
#:
#: Weights follow how well each is predicted (see
#: ``scripts/experiment_structure_features.py``): size and nitrogen content
#: are predicted well and count for more than ring counts, which are noisier.
WEIGHTS: Dict[str, float] = {
    "n_heavy_atoms": 1.0,
    "n_nitrogen": 1.0,
    "n_amide": 1.0,
    "n_oxygen": 0.7,
    "n_rings": 0.5,
    "has_macrocycle": 0.5,
    "has_sugar": 0.4,
    "has_aromatic_ring": 0.4,
}


def collect(limit: Optional[int] = None):
    """Clusters with their features, candidates and known structures."""
    data = []
    fill: Dict[str, Dict[str, str]] = {}
    if FILL.exists():
        fill = json.loads(FILL.read_text())

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

        for module in parsed.spec.modules:
            if module.substrate or module.domain_index is None:
                continue
            predicted = fill.get(accession, {}).get(f"{module.gene}:{module.domain_index}")
            if predicted:
                module.substrate = predicted

        prediction = assembly.assemble(parsed.spec, skip_unknown_substrates=True)
        if not prediction.variants:
            continue

        # Candidates: every release the chain allows, bare and decorated.
        candidates: List[Tuple[str, str]] = []
        genes = tailoring.TailoringGenes(counts=parsed.tailoring_genes)
        for release, smiles in prediction.variants.items():
            candidates.append((release, smiles))
            if parsed.tailoring_genes:
                decorated = tailoring.apply(smiles, genes, enabled=list(tailoring.DEFAULT_REACTIONS))
                if decorated and decorated != smiles:
                    candidates.append((f"{release}+tailored", decorated))

        profile = sf.describe(entry.structures[0])
        if profile is None:
            continue
        record.products = parsed.products  # type: ignore[attr-defined]
        data.append(
            {
                "bgc_id": accession,
                "features": features.build(record, parsed.spec, accession),
                "profile": profile,
                "candidates": candidates,
                "rule_choice": prediction.release,
                "rule_smiles": prediction.smiles,
                "structures": entry.structures,
                "names": entry.compound_names,
            }
        )
    return data


def score_candidate(candidate_profile: sf.StructureProfile, wanted: Dict[str, float]) -> float:
    """How far a candidate sits from the predicted description. Lower is better."""
    total = 0.0
    for name, weight in WEIGHTS.items():
        if name not in wanted:
            continue
        predicted = wanted[name]
        actual = candidate_profile.get(name)
        if name in sf.BINARY_TARGETS:
            total += weight * abs(predicted - actual)
        else:
            scale = max(abs(predicted), 1.0)
            total += weight * min(abs(predicted - actual) / scale, 2.0)
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    data = collect(args.limit)
    print(f"clusters with candidates: {len(data)}")
    print(f"mean candidates per cluster: {np.mean([len(d['candidates']) for d in data]):.1f}")

    names = features.feature_names([d["features"] for d in data])
    matrix = features.matrix([d["features"] for d in data], names, normalise=False)

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(data))
    fold_of = np.empty(len(data), dtype=int)
    for position, index in enumerate(order):
        fold_of[index] = position % args.folds

    # Out-of-fold predictions of every property used for ranking.
    wanted: List[Dict[str, float]] = [{} for _ in data]
    for target in WEIGHTS:
        truth = np.array([d["profile"].get(target) for d in data], dtype=float)
        for fold in range(args.folds):
            test = fold_of == fold
            train = ~test
            if test.sum() == 0 or train.sum() == 0:
                continue
            if target in sf.BINARY_TARGETS:
                model = RandomForestClassifier(
                    n_estimators=300, min_samples_leaf=2, class_weight="balanced",
                    random_state=args.seed, n_jobs=-1,
                )
            else:
                model = RandomForestRegressor(
                    n_estimators=300, min_samples_leaf=2, random_state=args.seed, n_jobs=-1
                )
            model.fit(matrix[train], truth[train])
            for index, value in zip(np.where(test)[0], model.predict(matrix[test])):
                wanted[int(index)][target] = float(value)

    rule_scores: List[float] = []
    ranked_scores: List[float] = []
    oracle_scores: List[float] = []
    changed = 0

    for index, entry in enumerate(data):
        best_candidate = None
        best_distance = float("inf")
        for release, smiles in entry["candidates"]:
            profile = sf.describe(smiles)
            if profile is None:
                continue
            distance = score_candidate(profile, wanted[index])
            if distance < best_distance:
                best_distance, best_candidate = distance, (release, smiles)
        if best_candidate is None:
            continue

        rule = compare.best_match(
            {entry["rule_choice"] or "x": entry["rule_smiles"]},
            entry["structures"], entry["names"], timeout=2,
        )
        ranked = compare.best_match(
            {best_candidate[0]: best_candidate[1]}, entry["structures"], entry["names"], timeout=2
        )
        oracle = compare.best_match(
            dict(entry["candidates"]), entry["structures"], entry["names"], timeout=2
        )
        if rule and ranked and oracle:
            rule_scores.append(rule.tanimoto)
            ranked_scores.append(ranked.tanimoto)
            oracle_scores.append(oracle.tanimoto)
            if best_candidate[1] != entry["rule_smiles"]:
                changed += 1

    report = {
        "clusters": len(rule_scores),
        "changed": changed,
        "mean_tanimoto_rule": round(float(np.mean(rule_scores)), 4),
        "mean_tanimoto_ranked": round(float(np.mean(ranked_scores)), 4),
        "mean_tanimoto_oracle": round(float(np.mean(oracle_scores)), 4),
    }
    print(f"\nchoosing among {report['clusters']} clusters ({changed} choices differ from the rule):")
    print(f"  rule          : {report['mean_tanimoto_rule']:.4f}")
    print(f"  ranked        : {report['mean_tanimoto_ranked']:.4f}")
    print(f"  oracle        : {report['mean_tanimoto_oracle']:.4f}  (best candidate, not reachable)")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
