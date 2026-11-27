"""Command line interface."""

from __future__ import annotations

import json
from pathlib import Path

import click

from . import (
    assembly,
    evaluate,
    mibig,
    predictor,
)
from .spec import ClusterSpec

DEFAULT_JSON_DIR = Path("data/raw/json/mibig_json_4.0")
DEFAULT_ANTISMASH_DIR = Path("data/raw/antismash/mibig_antismash_4.0_gbk_as8b1")


@click.group()
@click.version_option(package_name="sempi3")
def cli() -> None:
    """SeMPI-3: predict the backbone of a modular PKS/NRPS product."""


@cli.command()
@click.argument("bgc_id")
@click.option(
    "--json-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_JSON_DIR,
    show_default=True,
    help="Directory of MIBiG JSON entries.",
)
@click.option("--show-variants", is_flag=True, help="Print every release variant.")
def predict(bgc_id: str, json_dir: Path, show_variants: bool) -> None:
    """Predict the backbone of one MIBiG cluster."""
    entry = mibig.load_entry(json_dir / f"{bgc_id}.json")
    spec = mibig.to_spec(entry)
    prediction = assembly.assemble(spec)

    click.echo(f"{bgc_id}  classes={'|'.join(entry.classes)}  modules={spec.n_modules}")
    for module in spec.modules:
        substrate = (module.substrate or "-").replace("smiles:", "")
        flags = "".join(
            [
                "KR" if module.kr else "",
                "DH" if module.dh else "",
                "ER" if module.er else "",
                "E" if module.epimerization else "",
                "MT" if module.n_methylation else "",
            ]
        )
        click.echo(f"  {module.order:>2} {module.kind:<4} {substrate[:40]:<40} {flags}")
    click.echo(f"predicted: {prediction.smiles}")
    if show_variants:
        for release, smiles in prediction.variants.items():
            click.echo(f"  {release}: {smiles}")
    for structure, name in zip(entry.structures, entry.compound_names):
        click.echo(f"known ({name}): {structure}")
    for note in prediction.notes:
        click.echo(f"note: {note}")


@cli.command("run")
@click.argument("gbk", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["text", "json", "smiles", "tabular"]),
    default="text",
    help="tabular is one row per region, for a workflow.",
)
@click.option("--no-fill", is_flag=True, help="Do not name the substrates antiSMASH left unnamed.")
@click.option("--no-tailoring", is_flag=True, help="Report the bare backbone.")
@click.option("--out", type=click.Path(path_type=Path), default=None, help="Write to a file.")
def run_on_file(gbk: Path, output_format: str, no_fill: bool, no_tailoring: bool, out: Path) -> None:
    """Predict the product of an antiSMASH-annotated cluster.

    This is the way to use SeMPI-3 downstream of antiSMASH: point it at a
    GenBank file antiSMASH wrote, for any genome.
    """
    results = predictor.predict_all(gbk, fill_substrates=not no_fill, apply_tailoring=not no_tailoring)
    if not results:
        click.echo("no PKS/NRPS assembly line found in this file", err=True)
        raise SystemExit(0)

    if output_format == "json":
        text = json.dumps([r.to_dict() for r in results], indent=2)
    elif output_format == "smiles":
        text = "\n".join(f"{r.smiles}\t{r.bgc_id}" for r in results)
    elif output_format == "tabular":
        header = [
            "region",
            "products",
            "modules",
            "substrates_predicted",
            "release",
            "smiles",
            "antismash_smiles",
        ]
        rows = ["\t".join(header)]
        for r in results:
            rows.append(
                "\t".join(
                    [
                        r.bgc_id,
                        "|".join(r.products) or "-",
                        str(len(r.modules)),
                        str(r.filled_substrates),
                        r.release or "-",
                        r.smiles or "-",
                        r.antismash_smiles or "-",
                    ]
                )
            )
        text = "\n".join(rows)
    else:
        text = "\n\n".join(_render(r) for r in results)

    if out:
        out.write_text(text + "\n")
        click.echo(f"wrote {out}")
    else:
        click.echo(text)


def _render(result) -> str:
    """Human-readable report for one region."""
    if True:
        lines = [
            f"{result.bgc_id}  products={'|'.join(result.products) or '-'}  modules={len(result.modules)}",
        ]
        for module in result.modules:
            substrate = str(module["substrate"] or "-")[:34]
            extra = module["reduction_state"] or ""
            marks = "".join(
                ["E" if module["epimerization"] else "", "Cy" if module["heterocyclization"] else ""]
            )
            confidence = f" ({module['confidence']})" if module["confidence"] else ""
            lines.append(
                f"  {module['order']:>2} {module['kind']:<4} {substrate:<34}{confidence:<8} {extra} {marks}"
            )
        if result.filled_substrates:
            lines.append(f"substrates named by the pocket model: {result.filled_substrates}")
        if result.tailoring_genes:
            lines.append(f"tailoring enzymes: {result.tailoring_genes}")
        lines.append(f"release: {result.release}")
        lines.append(f"SeMPI-3:   {result.smiles}")
        if result.antismash_smiles:
            lines.append(f"antiSMASH: {result.antismash_smiles}")
        for note in result.notes[:4]:
            lines.append(f"note: {note}")
        return "\n".join(lines)


@cli.command("build-reference")
@click.option(
    "--json-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_JSON_DIR,
    show_default=True,
    help="MIBiG entries, any release.",
)
@click.option(
    "--antismash-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_ANTISMASH_DIR,
    show_default=True,
)
@click.option(
    "--paras",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/external/paras/dataset_paras.tsv"),
    show_default=True,
    help="Optional extra set of curated domains.",
)
@click.option("--out", type=click.Path(path_type=Path), default=None, help="Where to write the table.")
def build_reference(json_dir: Path, antismash_dir: Path, paras: Path, out: Path) -> None:
    """Rebuild the substrate reference from the data currently on disk.

    Run this after a new MIBiG release: fetch the data, then rebuild, and the
    packaged table grows with it. Nothing else has to change.
    """
    from . import reference

    path = reference.build(json_dir, antismash_dir, paras if paras and paras.exists() else None, out)
    click.echo(f"wrote {path}")


@cli.command("spec")
@click.argument("bgc_id")
@click.option(
    "--json-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_JSON_DIR,
    show_default=True,
)
def show_spec(bgc_id: str, json_dir: Path) -> None:
    """Print the parsed assembly line as JSON."""
    entry = mibig.load_entry(json_dir / f"{bgc_id}.json")
    click.echo(json.dumps(mibig.to_spec(entry).to_dict(), indent=2))


@cli.command("assemble")
@click.argument("spec_json", type=click.Path(exists=True, dir_okay=False, path_type=Path))
def assemble_spec(spec_json: Path) -> None:
    """Assemble a structure from a saved cluster specification."""
    with open(spec_json) as handle:
        spec = ClusterSpec.from_dict(json.load(handle))
    click.echo(json.dumps(assembly.assemble(spec).to_dict(), indent=2))


@cli.command("benchmark")
@click.option(
    "--json-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_JSON_DIR,
    show_default=True,
)
@click.option(
    "--antismash-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_ANTISMASH_DIR,
    show_default=True,
)
@click.option("--out", type=click.Path(path_type=Path), default=Path("results/benchmark.csv"), show_default=True)
@click.option("--limit", type=int, default=None, help="Only evaluate the first N entries.")
@click.option("--workers", type=int, default=4, show_default=True)
@click.option("--mcs-timeout", type=int, default=5, show_default=True)
def benchmark(json_dir: Path, antismash_dir: Path, out: Path, limit: int, workers: int, mcs_timeout: int) -> None:
    """Evaluate every arm over MIBiG and write a results table."""
    results = evaluate.run(
        json_dir, antismash_dir, limit=limit, workers=workers, mcs_timeout=mcs_timeout
    )
    evaluate.write_csv(results, out)
    summary = evaluate.summarise(results)
    evaluate.write_summary(summary, out.with_suffix(".summary.json"))
    click.echo(json.dumps(summary, indent=2))
    click.echo(f"wrote {out}")


def main() -> None:
    cli()


if __name__ == "__main__":
    main()
