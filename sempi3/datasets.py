"""Training sets built from MIBiG and its antiSMASH annotation.

Two supervised problems are set up here:

``adenylation domains``
    Protein sequence of an A domain -> the substrate it activates. MIBiG
    curates the substrate; antiSMASH supplies the domain sequence and its own
    call, which is the baseline to beat. MIBiG gives the gene but almost never
    the domain's coordinates, so domains are aligned by gene and by their order
    within that gene.

``chain release``
    Cluster -> how the product leaves the assembly line. The label is read off
    the known structure: a macrolactone, a macrolactam, or neither.

Both datasets carry the cluster accession, so a benchmark can split by cluster
and never train and test on the same gene cluster.
"""

from __future__ import annotations

import json
from dataclasses import (
    asdict,
    dataclass,
)
from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Optional,
)

from Bio import SeqIO
from rdkit import Chem

from . import mibig
from .spec import NRPS

#: Smallest ring that counts as a macrocycle rather than a normal heterocycle.
MACROCYCLE_MIN_SIZE = 8

#: Ester inside a ring: the bond a macrolactonising thioesterase makes.
LACTONE = Chem.MolFromSmarts("[C;R](=O)[O;R][C;R]")

#: Amide inside a ring: the bond a macrolactamising thioesterase makes.
LACTAM = Chem.MolFromSmarts("[C;R](=O)[N;R]")


@dataclass
class AdomainRecord:
    """One adenylation domain with its curated substrate."""

    bgc_id: str
    gene: str
    index: int
    sequence: str
    substrate: str
    substrate_smiles: Optional[str]
    antismash_call: Optional[str]

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class ReleaseRecord:
    """One cluster with the release mechanism read off its known structure."""

    bgc_id: str
    label: str
    classes: List[str]
    n_modules: int
    n_nrps: int
    n_pks: int
    has_thioesterase: bool
    has_reductase: bool
    n_kr: int
    n_dh: int
    n_er: int
    n_epimerase: int
    n_methyltransferase: int
    true_smiles: str

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def _domain_sequences(gbk_path: Path, domain_name: str) -> Dict[str, List[dict]]:
    """A domain's sequences per gene, in order along the protein.

    MIBiG names a gene by its protein id (``AAQ82561.1``) while antiSMASH tags
    domains with the locus tag (``fscA``). The CDS features carry both, so an
    alias map is built first and every name a CDS answers to becomes a key.
    """
    per_gene: Dict[str, List[dict]] = {}
    aliases: Dict[str, str] = {}
    for record in SeqIO.parse(str(gbk_path), "genbank"):
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
        for feature in record.features:
            if feature.type != "aSDomain":
                continue
            if feature.qualifiers.get("aSDomain", [""])[0] != domain_name:
                continue
            gene = feature.qualifiers.get("locus_tag", ["?"])[0]
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue
            call = None
            for entry in feature.qualifiers.get("specificity", []):
                if entry.lower().startswith("substrate consensus:"):
                    call = entry.split(":", 1)[1].strip()
                    break
            per_gene.setdefault(gene, []).append(
                {
                    "start": int(feature.location.start),
                    "sequence": translation,
                    "call": call,
                }
            )
    for entries in per_gene.values():
        entries.sort(key=lambda e: e["start"])
    # Make every alias of a gene reach the same domain list.
    for alias, canonical in aliases.items():
        if canonical in per_gene and alias not in per_gene:
            per_gene[alias] = per_gene[canonical]
    return per_gene


def build_adomain_dataset(json_dir: Path, antismash_dir: Path) -> List[AdomainRecord]:
    """Pair curated substrates with the antiSMASH domain sequences."""
    records: List[AdomainRecord] = []
    for json_path in sorted(Path(json_dir).glob("*.json")):
        accession = json_path.stem
        gbk_path = Path(antismash_dir) / f"{accession}.gbk"
        if not gbk_path.exists():
            continue
        entry = mibig.load_entry(json_path)
        spec_modules = entry.raw.get("biosynthesis", {}).get("modules", [])
        if not spec_modules:
            continue
        try:
            per_gene = _domain_sequences(gbk_path, "AMP-binding")
        except Exception:  # noqa: BLE001 - a malformed record must not stop the build
            continue
        if not per_gene:
            continue

        seen: Dict[str, int] = {}
        for raw in spec_modules:
            domain = raw.get("a_domain") or {}
            substrates = domain.get("substrates") or []
            gene = domain.get("gene") or (raw.get("genes") or [None])[0]
            if not substrates or not gene:
                continue
            index = seen.get(gene, 0)
            seen[gene] = index + 1
            candidates = per_gene.get(gene)
            if not candidates or index >= len(candidates):
                continue
            first = substrates[0]
            records.append(
                AdomainRecord(
                    bgc_id=accession,
                    gene=gene,
                    index=index,
                    sequence=candidates[index]["sequence"],
                    substrate=first.get("name", "").strip().lower(),
                    substrate_smiles=first.get("structure"),
                    antismash_call=candidates[index]["call"],
                )
            )
    return records


def release_label(smiles: str) -> str:
    """Read the release mechanism off a known structure.

    ``macrolactone`` and ``macrolactam`` require the ester or amide to sit in a
    ring of at least :data:`MACROCYCLE_MIN_SIZE` atoms, so that an ordinary
    five-membered lactam such as a proline residue is not mistaken for a
    macrocyclisation. Everything else is ``linear``.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unknown"
    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings() if len(ring) >= MACROCYCLE_MIN_SIZE]
    if not rings:
        return "linear"
    # Only the ring members of the match are tested: the carbonyl oxygen hangs
    # off the ring and would never be inside it.
    for pattern, label, ring_atom_indices in (
        (LACTONE, "macrolactone", (0, 2, 3)),
        (LACTAM, "macrolactam", (0, 2)),
    ):
        for match in mol.GetSubstructMatches(pattern):
            atoms = {match[i] for i in ring_atom_indices}
            if any(atoms <= ring for ring in rings):
                return label
    return "linear"


def build_release_dataset(json_dir: Path, antismash_dir: Path) -> List[ReleaseRecord]:
    """Label every benchmark cluster with its release mechanism and features."""
    from . import antismash_io  # imported here to keep module import cheap

    records: List[ReleaseRecord] = []
    for json_path in sorted(Path(json_dir).glob("*.json")):
        accession = json_path.stem
        gbk_path = Path(antismash_dir) / f"{accession}.gbk"
        if not gbk_path.exists():
            continue
        entry = mibig.load_entry(json_path)
        if not entry.is_pks_nrps or not entry.structures:
            continue
        try:
            parsed = antismash_io.load(gbk_path)
        except Exception:  # noqa: BLE001
            continue
        spec = parsed.spec
        if not spec.modules:
            continue
        label = release_label(entry.structures[0])
        if label == "unknown":
            continue
        records.append(
            ReleaseRecord(
                bgc_id=accession,
                label=label,
                classes=entry.classes,
                n_modules=spec.n_modules,
                n_nrps=sum(1 for m in spec.modules if m.kind == NRPS),
                n_pks=sum(1 for m in spec.modules if m.kind != NRPS),
                has_thioesterase=spec.release == "hydrolysis",
                has_reductase=spec.release == "reduction",
                n_kr=sum(1 for m in spec.modules if m.kr),
                n_dh=sum(1 for m in spec.modules if m.dh),
                n_er=sum(1 for m in spec.modules if m.er),
                n_epimerase=sum(1 for m in spec.modules if m.epimerization),
                n_methyltransferase=sum(1 for m in spec.modules if m.n_methylation),
                true_smiles=entry.structures[0],
            )
        )
    return records


def write_jsonl(records: Iterable, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for record in records:
            handle.write(json.dumps(record.to_dict()) + "\n")
    return path


def read_jsonl(path: Path) -> List[dict]:
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]
