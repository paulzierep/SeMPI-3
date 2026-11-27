#!/usr/bin/env python
"""Turn the SeMPI 2.0 supplementary A-domain table into training rows.

``A_domains.csv`` (Table S2 of Zierep et al., Metabolites 2020) lists 2,868
adenylation domains. The 1,546 rows taken from Prieto et al. name a UniProt
protein and the position of the domain within it, but carry no sequence, so
the proteins are fetched once from UniProt and the domains are located by
aligning the GrsA A domain along each protein.

The result is written as the same kind of row the other datasets use:
a binding-pocket signature and a substrate name.

Usage:
    python scripts/build_prieto_dataset.py --csv ~/Downloads/.../A_domains.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)

from Bio import Align
from Bio.Align import substitution_matrices

from sempi3 import (
    stachelhaus as st,
    substrates,
)

CACHE = Path("data/external/prieto")
FASTA = CACHE / "sequences.fasta"
OUT = Path("data/interim/prieto_adomains.jsonl")

#: UniProt serves this many accessions per request comfortably.
BATCH = 100

#: The adenylation domain of the reference, used to find domains in a protein.
REFERENCE_DOMAIN_END = 560

#: An alignment scoring below this is not another A domain.
MIN_SCORE = 200

#: One-letter and short substrate codes used in the supplementary table.
SPECIFICITY_CODES: Dict[str, str] = {
    "A": "alanine",
    "R": "arginine",
    "N": "asparagine",
    "D": "aspartic acid",
    "C": "cysteine",
    "E": "glutamic acid",
    "Q": "glutamine",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "L": "leucine",
    "K": "lysine",
    "M": "methionine",
    "F": "phenylalanine",
    "P": "proline",
    "S": "serine",
    "T": "threonine",
    "W": "tryptophan",
    "Y": "tyrosine",
    "V": "valine",
}


def fetch_sequences(accessions: List[str]) -> Dict[str, str]:
    """Download the proteins from UniProt, caching the FASTA on disk."""
    CACHE.mkdir(parents=True, exist_ok=True)
    if not FASTA.exists():
        chunks = []
        for start in range(0, len(accessions), BATCH):
            batch = accessions[start : start + BATCH]
            url = "https://rest.uniprot.org/uniprotkb/accessions?" + urllib.parse.urlencode(
                {"accessions": ",".join(batch), "format": "fasta"}
            )
            request = urllib.request.Request(url, headers={"User-Agent": "sempi3-dev/0.1"})
            with urllib.request.urlopen(request, timeout=120) as response:
                chunks.append(response.read().decode())
            print(f"  fetched {min(start + BATCH, len(accessions))}/{len(accessions)}")
            time.sleep(0.3)
        FASTA.write_text("".join(chunks))

    sequences: Dict[str, str] = {}
    name: Optional[str] = None
    parts: List[str] = []
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(parts)
            fields = line[1:].split("|")
            name = fields[1] if len(fields) > 2 else line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if name:
        sequences[name] = "".join(parts)
    return sequences


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def find_domains(protein: str, reference: str, aligner, needed: int = 12) -> List[Tuple[int, int]]:
    """Spans of every adenylation domain in a protein, in order.

    Each hit is masked out before searching again, so a protein carrying
    several modules yields one span per module. The search stops once
    ``needed`` domains are in hand, since the table says how many are wanted.
    """
    cleaned = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in protein.upper())
    spans: List[Tuple[int, int]] = []
    working = cleaned
    for _ in range(max(needed, 1)):
        try:
            alignment = aligner.align(reference, working)[0]
        except (ValueError, IndexError):
            break
        if alignment.score < MIN_SCORE:
            break
        blocks = alignment.aligned[1]
        if len(blocks) == 0:
            break
        start = int(blocks[0][0])
        end = int(blocks[-1][1])
        if end - start < 200:
            break
        spans.append((start, end))
        if len(spans) >= needed:
            break
        working = working[:start] + "X" * (end - start) + working[end:]
    return sorted(spans)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, required=True, help="A_domains.csv from the supplement")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    prieto = [r for r in rows if r["Source"].startswith("Prieto") and r["Specificity"].strip()]
    print(f"Prieto rows with a substrate: {len(prieto)}")

    wanted: Dict[str, List[dict]] = {}
    for row in prieto:
        accession, _, index = row["Name"].rpartition("_")
        if not accession or not index.isdigit():
            continue
        wanted.setdefault(accession, []).append({"index": int(index), "specificity": row["Specificity"].strip()})
    accessions = sorted(wanted)
    if args.limit:
        accessions = accessions[: args.limit]
    print(f"UniProt proteins to resolve: {len(accessions)}")

    sequences = fetch_sequences(accessions)
    print(f"sequences available: {len(sequences)}")

    reference = st.reference_sequence()
    domain_reference = reference[:REFERENCE_DOMAIN_END]
    aligner = _aligner()

    records = []
    missing_protein = missing_domain = unresolved_label = 0
    for position, accession in enumerate(accessions, start=1):
        protein = sequences.get(accession)
        if not protein:
            missing_protein += 1
            continue
        # The table counts domains from zero, so the highest index plus one
        # is how many spans have to be found.
        needed = max(entry["index"] for entry in wanted[accession]) + 1
        spans = find_domains(protein, domain_reference, aligner, needed=needed)
        for entry in sorted(wanted[accession], key=lambda e: e["index"]):
            index = entry["index"]
            if index < 0 or index >= len(spans):
                missing_domain += 1
                continue
            code = entry["specificity"]
            label = substrates.normalise(SPECIFICITY_CODES.get(code.upper(), code))
            if not label:
                unresolved_label += 1
                continue
            start, end = spans[index]
            signature = st.extract_signature(protein[start:end])
            if signature.count(st.GAP) > 3:
                continue
            records.append(
                {
                    "source": "prieto",
                    "accession": accession,
                    "index": index,
                    "signature": signature,
                    "code": st.code_from_signature(signature),
                    "substrate": label,
                }
            )
        if position % 100 == 0:
            print(f"  processed {position}/{len(accessions)} proteins, {len(records)} domains")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    print(
        f"wrote {OUT}: {len(records)} domains "
        f"(missing protein {missing_protein}, missing domain {missing_domain}, "
        f"unresolved label {unresolved_label})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
