#!/usr/bin/env python
"""Render the headline benchmark numbers into the README.

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
    "sempi3_antismash": "**SeMPI-3**",
    "sempi3_mibig": "SeMPI-3 with curated substrates",
}


def render(summary: dict) -> str:
    """A compact table for the README; the detail lives in docs/benchmark.md."""
    arms = summary["arms"]
    head = summary.get("head_to_head_sempi3_antismash", {})
    by_class = summary.get("by_class", {})

    lines = [
        f"Over the **{summary['clusters']} MIBiG 4.0 PKS/NRPS clusters** with a known structure:",
        "",
        "| | Clusters predicted | Tanimoto | Scoring ≥0.4 | Structure match |",
        "|---|---|---|---|---|",
    ]
    for arm, label in ARM_LABELS.items():
        data = arms.get(arm)
        if not data:
            continue
        lines.append(
            f"| {label} | {data['predicted']} | {data['tanimoto_mean']:.3f} | "
            f"{data['tanimoto_ge_0.4']} | {data['coverage_pred_mean']:.3f} |"
        )

    if head:
        lines += [
            "",
            f"Head to head on the {head['compared']} clusters where both predict: "
            f"**{head['wins']} wins, {head['losses']} losses** "
            f"({head['mean_delta']:+.3f} mean Tanimoto).",
        ]

    if by_class:
        pieces = []
        for label, data in by_class.items():
            name = label.replace("_only", "").replace("_", " ")
            pieces.append(f"{name} {data['sempi3_antismash']:.3f} vs {data['antismash']:.3f}")
        lines += ["", "By class (SeMPI-3 vs antiSMASH): " + ", ".join(pieces) + "."]

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
