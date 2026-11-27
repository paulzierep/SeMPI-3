"""Rebuilding the substrate reference.

The table in ``sempi3/data/substrate_signatures.tsv`` maps a binding-pocket
signature to the substrate the domains carrying it activate. It is the only
thing SeMPI-3 needs at prediction time beyond the cluster itself.

It is a derived file: every row comes from a curated domain in MIBiG or in
another collection such as PARAS. When a new release of either appears, fetch
it and rebuild — no code changes, no retraining, and the file stays small
enough to ship with the package.
"""

from __future__ import annotations

import csv
from collections import (
    Counter,
    defaultdict,
)
from pathlib import Path
from typing import (
    Dict,
    Optional,
)

from . import (
    datasets,
    stachelhaus as st,
    substrates,
)

DEFAULT_OUTPUT = Path(__file__).parent / "data" / "substrate_signatures.tsv"


def build_training_table(
    json_dir: Path,
    antismash_dir: Path,
    out: Optional[Path] = None,
) -> Path:
    """Write the table the ranking models are fitted from.

    One row per cluster: the features read off its gene cluster, and the
    properties measured on its known product. Nothing is fitted here, so the
    file stays readable and carries no dependency on a scikit-learn version.
    """
    from Bio import SeqIO

    from . import (
        antismash_io,
        features,
        mibig,
        ranking,
        structure_features as sf,
    )

    items = []
    profiles = []
    for json_path in sorted(Path(json_dir).glob("*.json")):
        accession = json_path.stem
        gbk = Path(antismash_dir) / f"{accession}.gbk"
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
        profile = sf.describe(entry.structures[0])
        if profile is None:
            continue
        record.products = parsed.products
        items.append(features.build(record, parsed.spec, accession))
        profiles.append(profile)

    if not items:
        raise RuntimeError("no clusters found to train on")

    names = features.feature_names(items)
    matrix = features.matrix(items, names, normalise=False)
    targets = {name: [p.get(name) for p in profiles] for name in ranking.WEIGHTS}
    path = ranking.write_training_table(names, matrix, targets, out)
    print(f"{len(items)} clusters, {len(names)} features, {len(targets)} properties")
    return path

#: A signature with more gaps than this did not align well enough to trust.
MAX_GAPS = 3


def _read_fasta(path: Path) -> Dict[str, str]:
    entries: Dict[str, str] = {}
    name: Optional[str] = None
    parts: list = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                entries[name] = "".join(parts)
            name, parts = line[1:].strip(), []
        else:
            parts.append(line.strip())
    if name:
        entries[name] = "".join(parts)
    return entries


def read_substrate_smiles(directory: Path) -> Dict[str, str]:
    """Structures for substrate names, from a PARAS release.

    Our own monomer table covers the common residues; PARAS names several
    hundred more, and without a structure a predicted substrate cannot be
    built into the chain at all.
    """
    path = Path(directory) / "smiles.tsv"
    table: Dict[str, str] = {}
    if not path.exists():
        return table
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index == 0:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[1].strip():
                continue
            name = substrates.normalise(parts[0]) or parts[0].strip().lower()
            table[name] = parts[1].strip()
    return table


def read_paras(directory: Path) -> Dict[str, str]:
    """Signatures and substrates from a PARAS database directory.

    Accepts either the two database files (3,653 domains) or the condensed
    table (2,873), whichever is present.
    """
    directory = Path(directory)
    signatures_path = directory / "extended_signatures.fasta"
    mapping_path = directory / "domain_substrate_mapping.txt"
    rows: Dict[str, str] = {}

    if signatures_path.exists() and mapping_path.exists():
        signatures = _read_fasta(signatures_path)
        with open(mapping_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                identifier = row.get("domain_id", "")
                signature = signatures.get(identifier, "")
                label = substrates.normalise(row.get("specificity"))
                if label and len(signature) == len(st.EXTENDED_POSITIONS):
                    rows[identifier] = f"{signature}\t{label}"
        return rows

    condensed = directory / "dataset_paras.tsv"
    if condensed.exists():
        with open(condensed) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                signature = row.get("extended_signature", "").strip().upper()
                label = substrates.normalise(row.get("substrate_specificity"))
                if label and len(signature) == len(st.EXTENDED_POSITIONS):
                    rows[row.get("domain_id", signature)] = f"{signature}\t{label}"
    return rows


def build(
    json_dir: Path,
    antismash_dir: Path,
    paras_dir: Optional[Path] = None,
    out: Optional[Path] = None,
) -> Path:
    """Write a fresh signature table from the data on disk.

    Args:
        json_dir: MIBiG entries, any release.
        antismash_dir: antiSMASH output for those entries, which supplies the
            domain sequences MIBiG itself does not carry.
        paras_dir: Optional directory holding a PARAS release, which gives
            signatures directly.
        out: Where to write; the packaged table by default.

    Returns:
        The path written.
    """
    out = Path(out) if out else DEFAULT_OUTPUT
    table: Dict[str, Counter] = defaultdict(Counter)

    records = datasets.build_adomain_dataset(Path(json_dir), Path(antismash_dir))
    usable = [
        r
        for r in records
        if substrates.normalise(r.substrate) and r.sequence and len(r.sequence) >= 50
    ]

    # Signatures in the table have to be read the same way as the queries that
    # will be matched against them, or the two will not line up.
    from . import hmm_align

    if hmm_align.available():
        signatures = hmm_align.extract_signatures([r.sequence for r in usable])
        signatures = [
            s if s.count(st.GAP) <= MAX_GAPS else st.extract_signature(r.sequence)
            for s, r in zip(signatures, usable)
        ]
        print(f"signatures read with the profile HMM ({len(usable)} domains)")
    else:
        signatures = [st.extract_signature(r.sequence) for r in usable]
        print("hmmalign not found: signatures read by pairwise alignment")

    from_mibig = 0
    for record, signature in zip(usable, signatures):
        label = substrates.normalise(record.substrate)
        if signature.count(st.GAP) > MAX_GAPS:
            continue
        table[signature][label] += 1
        from_mibig += 1

    from_extra = 0
    if paras_dir and Path(paras_dir).exists():
        for entry in read_paras(paras_dir).values():
            signature, label = entry.split("\t")
            table[signature][label] += 1
            from_extra += 1

    # Every label needs a structure, or a domain carrying it contributes
    # nothing to the chain. Our own monomers come first, then PARAS's table.
    from . import monomers

    extra_smiles = read_substrate_smiles(paras_dir) if paras_dir else {}
    out.parent.mkdir(parents=True, exist_ok=True)
    without_structure = 0
    with open(out, "w") as handle:
        handle.write("signature\tsubstrate\tsupport\tsmiles\n")
        for signature in sorted(table):
            counts = table[signature]
            label, _ = counts.most_common(1)[0]
            smiles = monomers.resolve(label) or extra_smiles.get(label, "")
            if not smiles:
                without_structure += 1
            handle.write(f"{signature}\t{label}\t{sum(counts.values())}\t{smiles}\n")

    print(
        f"{len(table)} signatures from {from_mibig} MIBiG domains and {from_extra} extra rows; "
        f"{without_structure} still have no structure"
    )
    return out
