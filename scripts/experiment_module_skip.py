#!/usr/bin/env python
"""Learn which detected modules do not belong in the assembly line.

Reading antiSMASH's domains gives more modules than MIBiG's curators record
for 228 of 513 clusters, by a median of two. That matters more than anything
else measured so far: when the module count matches the curation the mean
Tanimoto is 0.304, and when it does not it is 0.212.

So the model is given a decision per module — emit a building block, or skip
it. Labels come from whether a module could be aligned to a curated one.
Everything is out of fold and grouped by cluster, and the report says how
close the predictions land to the nearest known structure, so recall of a
neighbour would be visible.

Usage:
    python scripts/experiment_module_skip.py --folds 5
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
    compare,
    mibig,
)
from sempi3.spec import PKS

JSON_DIR = Path("data/raw/json/mibig_json_4.0")
ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")
OUT = Path("results/module_skip_experiment.json")


def gene_aliases(record) -> Dict[str, str]:
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


def module_features(spec, index: int, contig_edge: float) -> List[float]:
    """Everything available about a module that bears on whether it is real."""
    module = spec.modules[index]
    total = max(len(spec.modules), 1)
    genes = [m.gene for m in spec.modules]
    same_gene = sum(1 for g in genes if g == module.gene)

    return [
        1.0 if module.kind == PKS else 0.0,
        float(module.has_carrier),
        float(module.has_condensation),
        float(module.domain_count),
        float(module.kr),
        float(module.dh),
        float(module.er),
        float(module.is_starter),
        float(module.non_elongating),
        float(bool(module.substrate)),
        float(module.epimerization),
        float(module.n_methylation),
        index / total,
        float(index == 0),
        float(index == total - 1),
        float(total),
        float(same_gene),
        # Is this the first or last module on its gene?
        float(index == 0 or genes[index - 1] != module.gene),
        float(index == total - 1 or genes[index + 1] != module.gene),
        contig_edge,
        float(sum(1 for m in spec.modules if m.has_carrier)),
        float(total - sum(1 for m in spec.modules if m.has_carrier)),
    ]


def collect(limit: Optional[int] = None) -> Tuple[List[dict], Dict[str, dict]]:
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

        contig_edge = 0.0
        for feature in record.features:
            if feature.type == "cand_cluster":
                contig_edge = 1.0 if "True" in feature.qualifiers.get("contig_edge", []) else 0.0
        clusters[accession] = {
            "spec": parsed.spec,
            "structures": entry.structures,
            "names": entry.compound_names,
            "contig_edge": contig_edge,
        }

        curated = mibig.to_spec(entry)
        if not curated.modules:
            continue
        aliases = gene_aliases(record)
        curated_by_gene: Dict[Tuple[str, str], int] = Counter()
        for module in curated.modules:
            gene = aliases.get(module.gene or "", module.gene or "")
            curated_by_gene[(gene, module.kind)] += 1

        seen: Dict[Tuple[str, str], int] = Counter()
        for index, module in enumerate(parsed.spec.modules):
            key = (module.gene or "", module.kind)
            position = seen[key]
            seen[key] += 1
            # Modules beyond what the curation records for this gene are the
            # ones to learn to drop.
            keep = position < curated_by_gene.get(key, 0)
            rows.append(
                {
                    "bgc_id": accession,
                    "index": index,
                    "features": module_features(parsed.spec, index, contig_edge),
                    "keep": bool(keep),
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
    labelled = {r["bgc_id"] for r in rows}
    print(f"modules labelled: {len(rows)} from {len(labelled)} clusters")
    print(f"kept by curation: {sum(1 for r in rows if r['keep'])}, dropped: {sum(1 for r in rows if not r['keep'])}")

    features = np.array([r["features"] for r in rows], dtype=float)
    keep = np.array([r["keep"] for r in rows])
    groups = np.array([r["bgc_id"] for r in rows])

    unique = np.array(sorted(labelled))
    rng = np.random.default_rng(args.seed)
    rng.shuffle(unique)
    fold_of = np.array([{g: i % args.folds for i, g in enumerate(unique)}[g] for g in groups])

    predicted = np.zeros(len(rows), dtype=bool)
    for fold in range(args.folds):
        test = fold_of == fold
        train = ~test
        if test.sum() == 0 or train.sum() == 0:
            continue
        forest = RandomForestClassifier(
            n_estimators=400, min_samples_leaf=2, class_weight="balanced", random_state=args.seed, n_jobs=-1
        )
        forest.fit(features[train], keep[train])
        predicted[test] = forest.predict(features[test])

    accuracy = float((predicted == keep).mean())
    baseline = float(keep.mean())
    print(f"\nkeep/skip decision: learned {accuracy:.3f}, always-keep {baseline:.3f}")

    drops: Dict[str, List[int]] = {}
    for row, decision in zip(rows, predicted):
        if not decision:
            drops.setdefault(row["bgc_id"], []).append(row["index"])

    before: List[float] = []
    after: List[float] = []
    novelty: List[float] = []
    changed = 0
    known = [data["structures"][0] for data in clusters.values()]

    for accession, data in clusters.items():
        if accession not in labelled:
            continue
        spec = data["spec"]
        base = assembly.assemble(spec, skip_unknown_substrates=True)
        base_score = compare.best_match(
            {base.release or "x": base.smiles}, data["structures"], data["names"], timeout=3
        )
        skipped = set(drops.get(accession, []))
        kept_modules = [m for i, m in enumerate(spec.modules) if i not in skipped]
        if not kept_modules:
            kept_modules = spec.modules
        trimmed = assembly.assemble(
            type(spec)(bgc_id=spec.bgc_id, modules=kept_modules, release=spec.release, source=spec.source),
            skip_unknown_substrates=True,
        )
        trimmed_score = compare.best_match(
            {trimmed.release or "x": trimmed.smiles}, data["structures"], data["names"], timeout=3
        )
        if base_score and trimmed_score:
            before.append(base_score.tanimoto)
            after.append(trimmed_score.tanimoto)
            if base.smiles != trimmed.smiles:
                changed += 1
        if trimmed.smiles:
            predicted_mol = compare.parse(trimmed.smiles)
            if predicted_mol is not None:
                best = 0.0
                for other, reference in zip(clusters, known):
                    if other == accession:
                        continue
                    reference_mol = compare.parse(reference)
                    if reference_mol is not None:
                        best = max(best, compare.tanimoto(predicted_mol, reference_mol))
                novelty.append(best)

    report = {
        "modules": len(rows),
        "clusters": len(labelled),
        "keep_accuracy": round(accuracy, 4),
        "always_keep_accuracy": round(baseline, 4),
        "structures_compared": len(before),
        "structures_changed": changed,
        "mean_tanimoto_before": round(float(np.mean(before)), 4) if before else 0.0,
        "mean_tanimoto_after": round(float(np.mean(after)), 4) if after else 0.0,
        "mean_similarity_to_nearest_known": round(float(np.mean(novelty)), 4) if novelty else 0.0,
    }
    print("\nstructures after skipping the modules the model rejects:")
    print(f"  before : {report['mean_tanimoto_before']:.4f}")
    print(f"  after  : {report['mean_tanimoto_after']:.4f}  ({changed} structures changed)")
    print(f"  mean similarity to nearest known structure: {report['mean_similarity_to_nearest_known']:.3f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
