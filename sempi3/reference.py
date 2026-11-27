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

#: A signature with more gaps than this did not align well enough to trust.
MAX_GAPS = 3


def build(
    json_dir: Path,
    antismash_dir: Path,
    paras_tsv: Optional[Path] = None,
    out: Optional[Path] = None,
) -> Path:
    """Write a fresh signature table from the data on disk.

    Args:
        json_dir: MIBiG entries, any release.
        antismash_dir: antiSMASH output for those entries, which supplies the
            domain sequences MIBiG itself does not carry.
        paras_tsv: Optional extra collection in the PARAS format, which gives
            signatures directly.
        out: Where to write; the packaged table by default.

    Returns:
        The path written.
    """
    out = Path(out) if out else DEFAULT_OUTPUT
    table: Dict[str, Counter] = defaultdict(Counter)

    records = datasets.build_adomain_dataset(Path(json_dir), Path(antismash_dir))
    from_mibig = 0
    for record in records:
        label = substrates.normalise(record.substrate)
        if not label or not record.sequence or len(record.sequence) < 50:
            continue
        signature = st.extract_signature(record.sequence)
        if signature.count(st.GAP) > MAX_GAPS:
            continue
        table[signature][label] += 1
        from_mibig += 1

    from_extra = 0
    if paras_tsv and Path(paras_tsv).exists():
        with open(paras_tsv) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                signature = row.get("extended_signature", "").strip().upper()
                label = substrates.normalise(row.get("substrate_specificity"))
                if not label or len(signature) != len(st.EXTENDED_POSITIONS):
                    continue
                table[signature][label] += 1
                from_extra += 1

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as handle:
        handle.write("signature\tsubstrate\tsupport\n")
        for signature in sorted(table):
            counts = table[signature]
            label, _ = counts.most_common(1)[0]
            handle.write(f"{signature}\t{label}\t{sum(counts.values())}\n")

    print(f"{len(table)} signatures from {from_mibig} MIBiG domains and {from_extra} extra rows")
    return out
