"""Benchmark figures.

Three figures, each rendered for a light and a dark page so the README can
serve whichever the reader's theme asks for:

``quality``   how many clusters reach a given score, per arm
``head_to_head``  SeMPI-3 against antiSMASH, cluster by cluster
``by_class``  mean score per cluster class

Colours come from a palette validated for colour-vision deficiency in both
modes. Series are also labelled directly, so identity never rests on colour.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

#: Arm -> label, in report order.
ARM_LABELS: Dict[str, str] = {
    "antismash": "antiSMASH",
    "sempi3_antismash": "SeMPI-3",
    "sempi3_mibig": "SeMPI-3 (MIBiG modules)",
}

#: Short labels for direct annotation, where space is tight.
ARM_SHORT: Dict[str, str] = {
    "antismash": "antiSMASH",
    "sempi3_antismash": "SeMPI-3",
    "sempi3_mibig": "SeMPI-3 (MIBiG)",
}


@dataclass(frozen=True)
class Theme:
    """Colours for one page mode."""

    name: str
    surface: str
    text_primary: str
    text_secondary: str
    grid: str
    reference: str
    series: Dict[str, str]

    @property
    def transparent(self) -> bool:
        return False


LIGHT = Theme(
    name="light",
    surface="#ffffff",
    text_primary="#191919",
    text_secondary="#666666",
    grid="#e4e4e4",
    reference="#b0b0b0",
    series={
        "sempi3_antismash": "#257de0",
        "antismash": "#f47230",
        "sempi3_mibig": "#00b894",
    },
)

DARK = Theme(
    name="dark",
    surface="#191919",
    text_primary="#efefef",
    text_secondary="#9c9c9c",
    grid="#333333",
    reference="#5a5a5a",
    series={
        "sempi3_antismash": "#338ced",
        "antismash": "#e76624",
        "sempi3_mibig": "#00a888",
    },
)

THEMES = (LIGHT, DARK)


def _style_axes(ax, theme: Theme, grid_axis: str = "both") -> None:
    """Recessive grid and axes: the data should be the loudest thing."""
    ax.set_facecolor(theme.surface)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme.grid)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=theme.text_secondary, labelsize=9, length=3, width=1.0)
    ax.grid(True, axis=grid_axis, color=theme.grid, linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)


def _title(ax, text: str, subtitle: Optional[str], theme: Theme) -> None:
    """Title above the axes, with the subtitle on its own line below it."""
    ax.text(
        0.0,
        1.13 if subtitle else 1.04,
        text,
        transform=ax.transAxes,
        color=theme.text_primary,
        fontsize=12,
        ha="left",
        va="bottom",
    )
    if subtitle:
        ax.text(
            0.0,
            1.03,
            subtitle,
            transform=ax.transAxes,
            color=theme.text_secondary,
            fontsize=9,
            ha="left",
            va="bottom",
        )


def _figure_legend(fig, handles, labels, theme: Theme, ncol: int) -> None:
    """One legend for the whole figure, below the panels and out of the data."""
    legend = fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=9,
        loc="lower center",
        ncol=ncol,
        bbox_to_anchor=(0.5, -0.01),
        handlelength=1.6,
        columnspacing=2.2,
    )
    for text in legend.get_texts():
        text.set_color(theme.text_primary)


def _caption(fig, text: str, theme: Theme) -> None:
    fig.text(0.5, 0.005, text, color=theme.text_secondary, fontsize=8, ha="center", va="bottom")


def _label_on_curve(ax, x, y, target: float, text: str, theme: Theme) -> None:
    """Put a direct label where a survival curve crosses ``target``.

    Each series gets its own target height, so the labels never collide the
    way end-of-curve labels do once the curves converge.
    """
    if x.size == 0:
        return
    index = int(np.argmin(np.abs(y - target)))
    ax.annotate(
        text,
        xy=(x[index], y[index]),
        xytext=(7, 5),
        textcoords="offset points",
        color=theme.text_secondary,
        fontsize=8,
        va="bottom",
        ha="left",
    )


def _survival(values: Sequence[float]) -> tuple:
    """Fraction of clusters scoring at least x, for every observed x."""
    clean = np.sort(np.asarray([v for v in values if v == v], dtype=float))
    if clean.size == 0:
        return np.array([0.0]), np.array([0.0])
    fraction = 1.0 - np.arange(clean.size) / clean.size
    return clean, fraction


def plot_quality(df: pd.DataFrame, theme: Theme, path: Path) -> Path:
    """How many clusters reach a given score, for each arm."""
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), facecolor=theme.surface)
    label_heights = {"antismash": 0.30, "sempi3_antismash": 0.52, "sempi3_mibig": 0.74}
    panels = [
        ("tanimoto", "Tanimoto similarity", "Whole-molecule fingerprint similarity"),
        ("coverage_pred", "Structure match", "Share of predicted atoms found in the product"),
        (
            "strict_coverage_pred",
            "Structure match, strict",
            "The same, with bond orders required to agree",
        ),
    ]
    for ax, (metric, title, subtitle) in zip(axes, panels):
        _style_axes(ax, theme)
        for arm, label in ARM_LABELS.items():
            column = f"{arm}_{metric}"
            if column not in df:
                continue
            x, y = _survival(df[column].tolist())
            ax.plot(x, y, color=theme.series[arm], linewidth=2.0, label=label, solid_capstyle="round")
            # Direct label where the curve crosses its own target height, so
            # identity does not rest on colour and labels never collide.
            _label_on_curve(ax, x, y, label_heights[arm], ARM_SHORT[arm], theme)
        _title(ax, title, subtitle, theme)
        ax.set_xlabel("score", color=theme.text_secondary, fontsize=9)
        ax.set_ylabel("clusters scoring at least x", color=theme.text_secondary, fontsize=9)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(0, 1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.11, 1, 0.90])
    _figure_legend(fig, handles, labels, theme, ncol=3)
    return _save(fig, path, theme)


def plot_head_to_head(df: pd.DataFrame, theme: Theme, path: Path) -> Path:
    """Per-cluster comparison of SeMPI-3 against antiSMASH, split by class."""
    groups = [
        ("NRPS", df[df["classes"] == "NRPS"]),
        ("PKS", df[df["classes"] == "PKS"]),
        ("PKS/NRPS hybrid", df[df["classes"].astype(str).str.contains(r"PKS\|NRPS|NRPS\|PKS", regex=True)]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.4), facecolor=theme.surface, sharex=True, sharey=True)
    colour = theme.series["sempi3_antismash"]

    for ax, (label, subset) in zip(axes, groups):
        _style_axes(ax, theme)
        pair = subset.dropna(subset=["antismash_tanimoto", "sempi3_antismash_tanimoto"])
        base = pair["antismash_tanimoto"].to_numpy()
        mine = pair["sempi3_antismash_tanimoto"].to_numpy()
        ax.plot([0, 1], [0, 1], color=theme.reference, linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
        ax.scatter(
            base,
            mine,
            s=18,
            color=colour,
            alpha=0.55,
            linewidths=0.5,
            edgecolors=theme.surface,
            zorder=2,
        )
        wins = int((mine > base + 1e-6).sum())
        losses = int((mine < base - 1e-6).sum())
        _title(ax, label, f"{len(pair)} clusters", theme)
        ax.text(
            0.04,
            0.95,
            f"SeMPI-3 better: {wins}\nantiSMASH better: {losses}",
            transform=ax.transAxes,
            color=theme.text_primary,
            fontsize=9,
            va="top",
            ha="left",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("antiSMASH Tanimoto", color=theme.text_secondary, fontsize=9)
    axes[0].set_ylabel("SeMPI-3 Tanimoto", color=theme.text_secondary, fontsize=9)
    fig.tight_layout(rect=[0, 0.06, 1, 0.88])
    _caption(fig, "Each point is one cluster; above the dashed line SeMPI-3 is closer to the known structure.", theme)
    return _save(fig, path, theme)


def plot_by_class(df: pd.DataFrame, theme: Theme, path: Path) -> Path:
    """Mean Tanimoto per cluster class, one bar group per class."""
    groups = [
        ("NRPS", df[df["classes"] == "NRPS"]),
        ("PKS", df[df["classes"] == "PKS"]),
        ("Hybrid", df[df["classes"].astype(str).str.contains(r"PKS\|NRPS|NRPS\|PKS", regex=True)]),
        ("All", df),
    ]
    arms = list(ARM_LABELS)
    fig, ax = plt.subplots(figsize=(8.8, 4.8), facecolor=theme.surface)
    _style_axes(ax, theme, grid_axis="y")

    width = 0.26
    positions = np.arange(len(groups))
    for index, arm in enumerate(arms):
        column = f"{arm}_tanimoto"
        means = [subset[column].mean() if column in subset and len(subset) else np.nan for _, subset in groups]
        offset = (index - (len(arms) - 1) / 2) * width
        bars = ax.bar(
            positions + offset,
            means,
            width=width * 0.92,  # a 2px-equivalent gap between neighbouring bars
            color=theme.series[arm],
            label=ARM_LABELS[arm],
            linewidth=0,
        )
        for rect, value in zip(bars, means):
            if value != value:
                continue
            ax.annotate(
                f"{value:.2f}",
                xy=(rect.get_x() + rect.get_width() / 2, value),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                color=theme.text_secondary,
            )

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{label}\n{len(subset)} clusters" for label, subset in groups],
        color=theme.text_primary,
        fontsize=9,
    )
    ax.set_ylabel("mean Tanimoto", color=theme.text_secondary, fontsize=9)
    _title(ax, "Backbone similarity by cluster class", "Higher is closer to the known structure", theme)
    handles, labels = ax.get_legend_handles_labels()
    ax.set_ylim(0, max(0.5, ax.get_ylim()[1]))
    fig.tight_layout(rect=[0, 0.10, 1, 0.88])
    _figure_legend(fig, handles, labels, theme, ncol=3)
    return _save(fig, path, theme)


def plot_properties(report: Dict[str, dict], theme: Theme, path: Path) -> Path:
    """How well the cluster alone predicts each property of its product."""
    targets = report.get("targets", {})
    counts = [(k, v) for k, v in targets.items() if v.get("kind") == "count"]
    binaries = [(k, v) for k, v in targets.items() if v.get("kind") == "binary"]
    counts.sort(key=lambda kv: kv[1]["r2"])
    binaries.sort(key=lambda kv: kv[1]["gain"])

    fig, axes = plt.subplots(
        1, 2, figsize=(11.5, 4.6), facecolor=theme.surface, gridspec_kw={"width_ratios": [1.3, 1]}
    )

    ax = axes[0]
    _style_axes(ax, theme, grid_axis="x")
    labels = [k.replace("n_", "").replace("_", " ") for k, _ in counts]
    values = [v["r2"] for _, v in counts]
    positions = np.arange(len(labels))
    ax.barh(positions, values, color=theme.series["sempi3_antismash"], height=0.62, linewidth=0)
    for y, value in zip(positions, values):
        ax.annotate(
            f"{value:.2f}", xy=(value, y), xytext=(4, 0), textcoords="offset points",
            va="center", fontsize=8, color=theme.text_secondary,
        )
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, color=theme.text_primary, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("variance explained", color=theme.text_secondary, fontsize=9)
    _title(ax, "How much of a count the cluster explains", "1.0 would be exact", theme)

    ax = axes[1]
    _style_axes(ax, theme, grid_axis="x")
    labels = [k.replace("has_", "").replace("_", " ") for k, _ in binaries]
    positions = np.arange(len(labels))
    base = [v["majority"] for _, v in binaries]
    got = [v["accuracy"] for _, v in binaries]
    ax.barh(positions + 0.19, got, height=0.34, color=theme.series["sempi3_antismash"],
            label="predicted from the cluster", linewidth=0)
    ax.barh(positions - 0.19, base, height=0.34, color=theme.series["antismash"],
            label="always answer with the majority", linewidth=0)
    for y, value in zip(positions + 0.19, got):
        ax.annotate(f"{value:.2f}", xy=(value, y), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8, color=theme.text_secondary)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, color=theme.text_primary, fontsize=9)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("accuracy", color=theme.text_secondary, fontsize=9)
    _title(ax, "Does the product carry it?", "yes or no, out of fold", theme)

    handles, legend_labels = ax.get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.10, 1, 0.88])
    _figure_legend(fig, handles, legend_labels, theme, ncol=2)
    return _save(fig, path, theme)


def _save(fig, path: Path, theme: Theme) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, facecolor=theme.surface, bbox_inches="tight")
    plt.close(fig)
    return path


def render_all(csv_path: Path, out_dir: Path, properties: Optional[Path] = None) -> List[Path]:
    """Render every figure in both modes from a benchmark table."""
    import json

    df = pd.read_csv(csv_path)
    property_report = None
    if properties and Path(properties).exists():
        property_report = json.loads(Path(properties).read_text())
    written: List[Path] = []
    builders = (
        ("quality", plot_quality),
        ("head_to_head", plot_head_to_head),
        ("by_class", plot_by_class),
    )
    for theme in THEMES:
        suffix = "" if theme.name == "light" else "-dark"
        for name, builder in builders:
            written.append(builder(df, theme, Path(out_dir) / f"{name}{suffix}.png"))
        if property_report:
            written.append(
                plot_properties(property_report, theme, Path(out_dir) / f"properties{suffix}.png")
            )
    return written
