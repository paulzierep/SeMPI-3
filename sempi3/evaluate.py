"""Benchmark harness.

Three predictions are scored against the known MIBiG structure for every
cluster:

``antismash``
    antiSMASH's own structure prediction, read from its GenBank output. The
    baseline SeMPI-3 has to beat.

``sempi3_antismash``
    SeMPI-3 assembled from the antiSMASH domain annotation. This is the fair
    comparison: same input information, different assembly rules.

``sempi3_antismash_oracle``
    The same, but scored on the best of its release variants instead of the one
    it committed to. Not a prediction a user could make: an upper bound that
    shows what better release-mechanism prediction would be worth.

``sempi3_mibig``
    SeMPI-3 assembled from MIBiG's curated modules. Curated substrates are
    partly inferred from the known structure, so this arm is an upper bound on
    the assembly chemistry, not a de-novo prediction.

Every arm except the oracle commits to a single structure per cluster, which
is what a user gets.
"""

from __future__ import annotations

import csv
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import (
    dataclass,
    field,
)
from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Optional,
)

from . import (
    antismash_io,
    assembly,
    compare,
    mibig,
    tailoring,
)

#: The prediction arms, in report order.
ARMS = ("antismash", "sempi3_antismash", "sempi3_antismash_oracle", "sempi3_mibig")

FIELDS = [
    "bgc_id",
    "classes",
    "compound",
    "true_atoms",
    "n_modules_antismash",
    "n_modules_mibig",
    "status",
] + [
    f"{arm}_{metric}"
    for arm in ARMS
    for metric in ("tanimoto", "coverage_pred", "coverage_true", "atoms", "release")
] + [
    "true_smiles",
    "notes",
]


@dataclass
class Result:
    """One benchmark row: all arms for one cluster."""

    bgc_id: str
    classes: List[str] = field(default_factory=list)
    status: str = "ok"
    true_smiles: Optional[str] = None
    compound: Optional[str] = None
    true_atoms: int = 0
    n_modules_antismash: int = 0
    n_modules_mibig: int = 0
    scores: Dict[str, Optional[compare.Similarity]] = field(default_factory=dict)
    smiles: Dict[str, Optional[str]] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_row(self) -> Dict[str, object]:
        row: Dict[str, object] = {
            "bgc_id": self.bgc_id,
            "classes": "|".join(self.classes),
            "compound": self.compound,
            "true_atoms": self.true_atoms,
            "n_modules_antismash": self.n_modules_antismash,
            "n_modules_mibig": self.n_modules_mibig,
            "status": self.status,
            "true_smiles": self.true_smiles,
            "notes": "; ".join(self.notes[:4]),
        }
        for arm in ARMS:
            score = self.scores.get(arm)
            row[f"{arm}_tanimoto"] = round(score.tanimoto, 4) if score else None
            row[f"{arm}_coverage_pred"] = round(score.mcs_coverage_pred, 4) if score else None
            row[f"{arm}_coverage_true"] = round(score.mcs_coverage_true, 4) if score else None
            row[f"{arm}_atoms"] = score.pred_atoms if score else None
            row[f"{arm}_release"] = score.release if score else None
            row[f"{arm}_smiles"] = self.smiles.get(arm)
        return row


def apply_substrate_fill(spec, fill: Dict[str, str]) -> int:
    """Give unassigned NRPS modules a predicted substrate.

    ``fill`` maps ``"<gene>:<domain index>"`` to a substrate name. Only modules
    antiSMASH left empty are touched, so its own calls always win.
    """
    applied = 0
    for module in spec.modules:
        if module.substrate or module.domain_index is None:
            continue
        predicted = fill.get(f"{module.gene}:{module.domain_index}")
        if predicted:
            module.substrate = predicted
            applied += 1
    return applied


def evaluate_entry(
    accession: str,
    json_dir: Path,
    antismash_dir: Path,
    mcs_timeout: int = 5,
    substrate_fill: Optional[Dict[str, Dict[str, str]]] = None,
    release_prediction: Optional[Dict[str, str]] = None,
    tailoring_enabled: Optional[List[str]] = None,
) -> Optional[Result]:
    """Score every arm for one cluster, or return ``None`` if out of scope."""
    json_path = Path(json_dir) / f"{accession}.json"
    if not json_path.exists():
        return None
    entry = mibig.load_entry(json_path)
    if not entry.is_pks_nrps or not entry.structures:
        return None

    result = Result(
        bgc_id=accession,
        classes=entry.classes,
        true_smiles=entry.structures[0],
    )
    true = compare.parse(entry.structures[0])
    result.true_atoms = true.GetNumHeavyAtoms() if true else 0

    # Arm 1 and 2: antiSMASH output.
    gbk_path = Path(antismash_dir) / f"{accession}.gbk"
    if gbk_path.exists():
        try:
            record = antismash_io.load(gbk_path)
        except Exception as exc:  # noqa: BLE001
            result.notes.append(f"antismash parse failed: {exc}")
            record = None
        if record is not None:
            result.n_modules_antismash = record.spec.n_modules
            if record.baseline_smiles:
                result.smiles["antismash"] = record.baseline
                result.scores["antismash"] = compare.best_match(
                    {f"antismash_{i}": s for i, s in enumerate(record.baseline_smiles)},
                    entry.structures,
                    entry.compound_names,
                    timeout=mcs_timeout,
                )
            if substrate_fill:
                applied = apply_substrate_fill(record.spec, substrate_fill.get(accession, {}))
                if applied:
                    result.notes.append(f"filled {applied} substrate(s) with the pocket model")
            prediction = assembly.assemble(
                record.spec,
                skip_unknown_substrates=True,
                release_override=(release_prediction or {}).get(accession),
            )
            result.notes.extend(prediction.notes[:2])
            if prediction.variants and tailoring_enabled and record.tailoring_genes:
                genes = tailoring.TailoringGenes(counts=record.tailoring_genes)
                decorated = tailoring.apply(prediction.smiles, genes, enabled=tailoring_enabled)
                if decorated and decorated != prediction.smiles:
                    prediction.variants[prediction.release] = decorated
                    prediction.smiles = decorated
                    result.notes.append("tailoring applied")
            if prediction.variants:
                # Committed prediction: the single structure a user would get.
                result.smiles["sempi3_antismash"] = prediction.smiles
                result.scores["sempi3_antismash"] = compare.best_match(
                    {prediction.release: prediction.smiles},
                    entry.structures,
                    entry.compound_names,
                    timeout=mcs_timeout,
                )
                # Oracle over release variants: an upper bound, not a prediction.
                result.smiles["sempi3_antismash_oracle"] = prediction.smiles
                result.scores["sempi3_antismash_oracle"] = compare.best_match(
                    prediction.variants, entry.structures, entry.compound_names, timeout=mcs_timeout
                )
    else:
        result.notes.append("no antismash output")

    # Arm 3: MIBiG curated modules.
    spec = mibig.to_spec(entry)
    result.n_modules_mibig = spec.n_modules
    if spec.modules:
        prediction = assembly.assemble(spec, skip_unknown_substrates=True)
        if prediction.variants:
            result.smiles["sempi3_mibig"] = prediction.smiles
            result.scores["sempi3_mibig"] = compare.best_match(
                {prediction.release: prediction.smiles},
                entry.structures,
                entry.compound_names,
                timeout=mcs_timeout,
            )

    for arm in ARMS:
        score = result.scores.get(arm)
        if score and score.compound:
            result.compound = score.compound
            break
    if not any(result.scores.get(arm) for arm in ARMS):
        result.status = "no_prediction"
    return result


def _worker(args):
    accession, json_dir, antismash_dir, timeout, fill_path, release_path, tailoring_enabled = args
    fill = None
    if fill_path:
        with open(fill_path) as handle:
            fill = json.load(handle)
    releases = None
    if release_path:
        with open(release_path) as handle:
            releases = json.load(handle).get("predictions", {})
    try:
        return evaluate_entry(
            accession,
            Path(json_dir),
            Path(antismash_dir),
            mcs_timeout=timeout,
            substrate_fill=fill,
            release_prediction=releases,
            tailoring_enabled=tailoring_enabled,
        )
    except Exception as exc:  # noqa: BLE001 - one bad entry must not stop the run
        return Result(bgc_id=accession, status=f"error: {exc}")


def run(
    json_dir: Path,
    antismash_dir: Path,
    limit: Optional[int] = None,
    workers: int = 4,
    mcs_timeout: int = 5,
    substrate_fill_path: Optional[Path] = None,
    release_prediction_path: Optional[Path] = None,
    tailoring_enabled: Optional[List[str]] = None,
) -> List[Result]:
    """Evaluate every MIBiG entry that has a structure and is PKS/NRPS."""
    accessions = sorted(p.stem for p in Path(json_dir).glob("*.json"))
    if limit:
        accessions = accessions[:limit]
    fill_path = str(substrate_fill_path) if substrate_fill_path else None
    release_path = str(release_prediction_path) if release_prediction_path else None
    jobs = [
        (a, str(json_dir), str(antismash_dir), mcs_timeout, fill_path, release_path, tailoring_enabled)
        for a in accessions
    ]
    results: List[Result] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(_worker, jobs, chunksize=4):
            if result is not None:
                results.append(result)
    return results


def write_csv(results: Iterable[Result], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(FIELDS)
    for arm in ARMS:
        key = f"{arm}_smiles"
        if key not in fields:
            fields.append(key)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for result in results:
            writer.writerow(result.to_row())


def _median(values: List[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def summarise(results: List[Result]) -> Dict[str, object]:
    """Headline numbers, overall and per arm."""
    summary: Dict[str, object] = {"clusters": len(results)}

    per_arm: Dict[str, Dict[str, object]] = {}
    for arm in ARMS:
        scores = [r.scores[arm] for r in results if r.scores.get(arm)]
        tanimotos = [s.tanimoto for s in scores]  # type: ignore[union-attr]
        coverage = [s.mcs_coverage_pred for s in scores]  # type: ignore[union-attr]
        per_arm[arm] = {
            "predicted": len(scores),
            "tanimoto_mean": round(sum(tanimotos) / len(tanimotos), 4) if tanimotos else 0.0,
            "tanimoto_median": round(_median(tanimotos), 4),
            "tanimoto_ge_0.4": sum(1 for t in tanimotos if t >= 0.4),
            "tanimoto_ge_0.6": sum(1 for t in tanimotos if t >= 0.6),
            "coverage_pred_mean": round(sum(coverage) / len(coverage), 4) if coverage else 0.0,
        }
    summary["arms"] = per_arm

    # Head-to-head on clusters where both arms made a prediction.
    for challenger in ("sempi3_antismash", "sempi3_antismash_oracle", "sempi3_mibig"):
        pairs = [
            (r.scores["antismash"].tanimoto, r.scores[challenger].tanimoto)  # type: ignore[union-attr]
            for r in results
            if r.scores.get("antismash") and r.scores.get(challenger)
        ]
        wins = sum(1 for base, mine in pairs if mine > base + 1e-6)
        losses = sum(1 for base, mine in pairs if mine < base - 1e-6)
        summary[f"head_to_head_{challenger}"] = {
            "compared": len(pairs),
            "wins": wins,
            "ties": len(pairs) - wins - losses,
            "losses": losses,
            "mean_delta": round(sum(mine - base for base, mine in pairs) / len(pairs), 4) if pairs else 0.0,
        }

    by_class: Dict[str, Dict[str, float]] = {}
    for label, predicate in (
        ("NRPS_only", lambda r: r.classes == ["NRPS"]),
        ("PKS_only", lambda r: r.classes == ["PKS"]),
        ("hybrid", lambda r: set(r.classes) >= {"PKS", "NRPS"}),
    ):
        subset = [r for r in results if predicate(r)]
        entry: Dict[str, float] = {"clusters": len(subset)}
        for arm in ARMS:
            values = [r.scores[arm].tanimoto for r in subset if r.scores.get(arm)]  # type: ignore[union-attr]
            entry[arm] = round(sum(values) / len(values), 4) if values else 0.0
        by_class[label] = entry
    summary["by_class"] = by_class
    return summary


def write_summary(summary: Dict[str, object], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
