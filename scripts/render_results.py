#!/usr/bin/env python
"""Render benchmark numbers into the README between the RESULTS markers.

Usage:
    python scripts/render_results.py results/benchmark.summary.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

START = "<!-- RESULTS:START -->"
END = "<!-- RESULTS:END -->"

ARM_LABELS = {
    "antismash": "antiSMASH (baseline)",
    "sempi3_antismash": "SeMPI-3 from antiSMASH",
    "sempi3_mibig": "SeMPI-3 from MIBiG modules",
}


def render(summary: dict) -> str:
    arms = summary["arms"]
    lines = [
        f"Benchmark set: **{summary['clusters']} clusters** "
        "(MIBiG 4.0 PKS/NRPS entries with a known structure).",
        "",
        "| Arm | Predicted | Tanimoto mean | Tanimoto median | ≥0.4 | ≥0.6 | Coverage of prediction |",
        "|---|---|---|---|---|---|---|",
    ]
    for arm, label in ARM_LABELS.items():
        data = arms.get(arm)
        if not data:
            continue
        lines.append(
            f"| {label} | {data['predicted']} | {data['tanimoto_mean']:.3f} | "
            f"{data['tanimoto_median']:.3f} | {data['tanimoto_ge_0.4']} | "
            f"{data['tanimoto_ge_0.6']} | {data['coverage_pred_mean']:.3f} |"
        )

    lines += ["", "**Head to head against antiSMASH**, on clusters where both made a prediction:", ""]
    lines += ["| Arm | Compared | Wins | Ties | Losses | Mean Tanimoto difference |", "|---|---|---|---|---|---|"]
    for arm in ("sempi3_antismash", "sempi3_mibig"):
        head = summary.get(f"head_to_head_{arm}")
        if not head:
            continue
        lines.append(
            f"| {ARM_LABELS[arm]} | {head['compared']} | {head['wins']} | {head['ties']} | "
            f"{head['losses']} | {head['mean_delta']:+.3f} |"
        )

    by_class = summary.get("by_class", {})
    if by_class:
        lines += ["", "**Mean Tanimoto by cluster class:**", ""]
        lines += ["| Class | Clusters | antiSMASH | SeMPI-3 (antiSMASH input) | SeMPI-3 (MIBiG modules) |"]
        lines += ["|---|---|---|---|---|"]
        for label, data in by_class.items():
            lines.append(
                f"| {label} | {int(data['clusters'])} | {data['antismash']:.3f} | "
                f"{data['sempi3_antismash']:.3f} | {data['sempi3_mibig']:.3f} |"
            )

    return "\n".join(lines)


def main() -> int:
    summary_path = Path(sys.argv[1] if len(sys.argv) > 1 else "results/benchmark.summary.json")
    readme_path = Path(sys.argv[2] if len(sys.argv) > 2 else "README.md")
    summary = json.loads(summary_path.read_text())
    text = readme_path.read_text()
    if START not in text or END not in text:
        print(f"markers not found in {readme_path}", file=sys.stderr)
        return 1
    head, rest = text.split(START, 1)
    _old, tail = rest.split(END, 1)
    readme_path.write_text(f"{head}{START}\n{render(summary)}\n{END}{tail}")
    print(f"updated {readme_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
