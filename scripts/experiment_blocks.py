#!/usr/bin/env python
"""Learn to emit building blocks from the BGC annotation.

The idea is DeepSeMS's — a token sequence read off the cluster, a token
sequence written out for the product — but the output tokens are building
blocks rather than SMILES characters, so the chemistry engine can turn any
generated sequence into a valid molecule.

Supervision comes from MIBiG's curated modules: for clusters that carry both
an antiSMASH annotation and a hand-curated assembly line, the model learns to
map one onto the other. Our reading of antiSMASH agrees with the curation on
only about two polyketide modules in three, so there is something to learn.

Everything is evaluated out of fold, grouped by cluster, and the report
includes how similar each prediction is to its nearest *training* structure,
so a model that merely recalls a neighbour is visible as one.

Usage:
    python scripts/experiment_blocks.py --folds 5
"""

from __future__ import annotations

import argparse
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
from sklearn.ensemble import RandomForestClassifier

from sempi3 import (
    antismash_io,
    assembly,
    blocks as blocks_module,
    compare,
    mibig,
)
from sempi3.spec import PKS

JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
OUT = Path("results/blocks_experiment.json")


def gene_aliases(record) -> Dict[str, str]:
    """Every name a CDS answers to, mapped to its locus tag."""
    aliases: Dict[str, str] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        names = []
        for key in ("locus_tag", "gene", "protein_id"):
            names.extend(feature.qualifiers.get(key, []))
        if not names:
            continue
        canonical = feature.qualifiers.get("locus_tag", names)[0]
        for name in names:
            aliases[name] = canonical
    return aliases


def module_features(spec, index: int) -> List[float]:
    """What the model sees about one module and its neighbours."""
    module = spec.modules[index]
    total = max(len(spec.modules), 1)

    def encode(mod) -> List[float]:
        if mod is None:
            return [0.0] * 8
        return [
            1.0 if mod.kind == PKS else 0.0,
            float(mod.kr),
            float(mod.dh),
            float(mod.er),
            float(mod.kr_active),
            float(mod.non_elongating),
            float(mod.alpha_methyl),
            float(mod.is_starter),
        ]

    previous = spec.modules[index - 1] if index > 0 else None
    following = spec.modules[index + 1] if index + 1 < len(spec.modules) else None

    features = encode(module) + encode(previous) + encode(following)
    features += [
        index / total,
        float(index),
        float(total - index - 1),
        float(total),
        float(sum(1 for m in spec.modules if m.kind == PKS)),
        float(sum(1 for m in spec.modules if m.kind != PKS)),
        1.0 if module.reduction_override else 0.0,
        float(module.substrate == "mmal") if module.substrate else 0.0,
        float(module.substrate == "mal") if module.substrate else 0.0,
    ]
    return features


def collect(limit: Optional[int] = None) -> Tuple[List[dict], Dict[str, dict]]:
    """Aligned (antiSMASH module -> curated block) pairs, plus cluster data."""
    from Bio import SeqIO

    rows: List[dict] = []
    clusters: Dict[str, dict] = {}
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
        clusters[accession] = {
            "spec": parsed.spec,
            "structures": entry.structures,
            "names": entry.compound_names,
            "tailoring": parsed.tailoring_genes,
        }

        curated = mibig.to_spec(entry)
        if not curated.modules:
            continue
        aliases = gene_aliases(record)

        # Align by gene and by order within that gene, for each kind.
        curated_by_gene: Dict[Tuple[str, str], List] = {}
        for module in curated.modules:
            gene = aliases.get(module.gene or "", module.gene or "")
            curated_by_gene.setdefault((gene, module.kind), []).append(module)
        seen: Dict[Tuple[str, str], int] = {}
        for index, module in enumerate(parsed.spec.modules):
            key = (module.gene or "", module.kind)
            position = seen.get(key, 0)
            seen[key] = position + 1
            candidates = curated_by_gene.get(key)
            if not candidates or position >= len(candidates):
                continue
            target = candidates[position]
            if module.kind != PKS:
                continue
            rows.append(
                {
                    "bgc_id": accession,
                    "index": index,
                    "features": module_features(parsed.spec, index),
                    "target_reduction": target.reduction_state,
                    "rule_reduction": module.reduction_state,
                }
            )
    return rows, clusters


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rows, clusters = collect(args.limit)
    print(f"aligned polyketide modules: {len(rows)} from {len({r['bgc_id'] for r in rows})} clusters")
    print(f"clusters available for structure scoring: {len(clusters)}")
    if not rows:
        print("nothing to train on")
        return 1

    features = np.array([r["features"] for r in rows], dtype=float)
    targets = np.array([r["target_reduction"] for r in rows])
    rule = np.array([r["rule_reduction"] for r in rows])
    groups = np.array([r["bgc_id"] for r in rows])

    unique = np.array(sorted(set(groups)))
    rng = np.random.default_rng(args.seed)
    rng.shuffle(unique)
    fold_of = np.array([{g: i % args.folds for i, g in enumerate(unique)}[g] for g in groups])

    predictions = np.empty(len(rows), dtype=object)
    for fold in range(args.folds):
        test = fold_of == fold
        train = ~test
        if test.sum() == 0 or train.sum() == 0:
            continue
        forest = RandomForestClassifier(
            n_estimators=400, min_samples_leaf=2, class_weight="balanced", random_state=args.seed, n_jobs=-1
        )
        forest.fit(features[train], targets[train])
        predictions[test] = forest.predict(features[test])

    rule_accuracy = float((rule == targets).mean())
    model_accuracy = float((predictions == targets).mean())
    print(f"\nreduction state, against MIBiG curation:")
    print(f"  domain rule : {rule_accuracy:.3f}")
    print(f"  learned     : {model_accuracy:.3f}")
    print(f"  majority    : {Counter(targets).most_common(1)[0][1] / len(targets):.3f}")

    # Apply the learned states and score the structures they produce.
    by_cluster: Dict[str, Dict[int, str]] = {}
    for row, prediction in zip(rows, predictions):
        if prediction:
            by_cluster.setdefault(row["bgc_id"], {})[row["index"]] = str(prediction)

    changed_structures = 0
    scored_rule: List[float] = []
    scored_model: List[float] = []
    novelty: List[float] = []
    training_structures = [data["structures"][0] for data in clusters.values()]

    for accession, overrides in by_cluster.items():
        data = clusters.get(accession)
        if not data:
            continue
        spec = data["spec"]
        base = assembly.assemble(spec, skip_unknown_substrates=True)
        base_score = compare.best_match(
            {base.release or "x": base.smiles}, data["structures"], data["names"], timeout=3
        )
        for index, state in overrides.items():
            spec.modules[index].reduction_override = state
        tuned = assembly.assemble(spec, skip_unknown_substrates=True)
        tuned_score = compare.best_match(
            {tuned.release or "x": tuned.smiles}, data["structures"], data["names"], timeout=3
        )
        if base_score and tuned_score:
            scored_rule.append(base_score.tanimoto)
            scored_model.append(tuned_score.tanimoto)
            if base.smiles != tuned.smiles:
                changed_structures += 1
        # How close is the prediction to the nearest known structure? A model
        # that recalls a neighbour shows up here.
        if tuned.smiles:
            best = 0.0
            predicted_mol = compare.parse(tuned.smiles)
            if predicted_mol is not None:
                for other, reference in zip(clusters, training_structures):
                    if other == accession:
                        continue
                    reference_mol = compare.parse(reference)
                    if reference_mol is None:
                        continue
                    best = max(best, compare.tanimoto(predicted_mol, reference_mol))
                novelty.append(best)

    report = {
        "aligned_modules": len(rows),
        "clusters_trained_on": len({r["bgc_id"] for r in rows}),
        "reduction_accuracy": {
            "domain_rule": round(rule_accuracy, 4),
            "learned": round(model_accuracy, 4),
        },
        "structures_compared": len(scored_rule),
        "structures_changed_by_model": changed_structures,
        "mean_tanimoto_rule": round(float(np.mean(scored_rule)), 4) if scored_rule else 0.0,
        "mean_tanimoto_learned": round(float(np.mean(scored_model)), 4) if scored_model else 0.0,
        "mean_similarity_to_nearest_known": round(float(np.mean(novelty)), 4) if novelty else 0.0,
    }
    print("\nstructures built from the learned blocks:")
    print(f"  domain rule : {report['mean_tanimoto_rule']:.4f}")
    print(f"  learned     : {report['mean_tanimoto_learned']:.4f}")
    print(f"  mean similarity of a prediction to the nearest known structure: "
          f"{report['mean_similarity_to_nearest_known']:.3f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
