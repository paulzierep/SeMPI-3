"""Reading antiSMASH GenBank output.

antiSMASH is not a dependency of SeMPI-3. Its *output* is the primary input:
the annotated GenBank files carry the domain architecture, the substrate
predictions, and antiSMASH's own structure prediction, which serves as the
baseline to beat.

MIBiG publishes antiSMASH 8 output for every entry
(``mibig_antismash_4.0_gbk_as8b1.tar.bz2``), so the whole benchmark runs
without ever installing antiSMASH.
"""

from __future__ import annotations

import re
from dataclasses import (
    dataclass,
    field,
)
from pathlib import Path
from typing import (
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
)

from Bio import SeqIO

from .spec import (
    ClusterSpec,
    Module,
    NRPS,
    PKS,
    RELEASE_HYDROLYSIS,
    RELEASE_REDUCTION,
)

#: Role each antiSMASH domain plays when modules are cut out of a protein.
DOMAIN_ROLES: Dict[str, str] = {
    "PKS_KS": "ks",
    "PKS_AT": "at",
    "Trans-AT_docking": "trans_at",
    "PKS_KR": "kr",
    "PKS_DH": "dh",
    "PKS_DH2": "dh",
    "PKS_DHt": "dh",
    "PKS_ER": "er",
    "ACP": "carrier",
    "PKS_PP": "carrier",
    "PP-binding": "carrier",
    "PCP": "carrier",
    "Condensation": "c",
    "Heterocyclization": "cy",
    "AMP-binding": "a",
    "A-OX": "a",
    "Epimerization": "e",
    "MT": "mt",
    "cMT": "mt",
    "nMT": "mt",
    "oMT": "omt",
    "Thioesterase": "te",
    "TD": "td",
    "TIGR01720": "ignore",
}

#: antiSMASH substrate names that mean "no prediction".
UNKNOWN_SUBSTRATES = {"x", "pk", "nrp", "?", ""}

_CONSENSUS = re.compile(r"substrate consensus:\s*(.+)", re.IGNORECASE)
_KR_ACTIVITY = re.compile(r"KR activity:\s*(\w+)", re.IGNORECASE)
_TRANSATOR = re.compile(r"transATor:\s*(.+)", re.IGNORECASE)

#: transATor calls mapped onto the oxidation level of the beta carbon.
#:
#: Trans-AT assembly lines carry no acyltransferase in the module, so the
#: usual KR/DH/ER reading does not apply. antiSMASH instead reports what the
#: ketosynthase clade implies, and that is used directly.
TRANSATOR_REDUCTION = {
    "keto": "keto",
    "beta-oh": "hydroxyl",
    "beta-d-oh": "hydroxyl",
    "beta-l-oh": "hydroxyl",
    "beta-d-ome": "hydroxyl",
    "alphabeta-oh": "hydroxyl",
    "alphame-betaoh": "hydroxyl",
    "alphame-beta-l-oh": "hydroxyl",
    "alphame-beta-d-oh": "hydroxyl",
    "non-elongating-beta-oh": "hydroxyl",
    "db": "ene",
    "edb": "ene",
    "zdb": "ene",
    "shdb": "ene",
    "red/shdb": "ene",
    "alphame-db": "ene",
    "alphame-edb": "ene",
    "alphame-zdb": "ene",
    "alpha-d-me-shdb": "ene",
    "non-elongating-db": "ene",
    "non-elongating-alphame-edb": "ene",
    "red": "alkyl",
}


@dataclass
class AntismashRecord:
    """What SeMPI-3 reads out of one antiSMASH GenBank file."""

    bgc_id: str
    spec: ClusterSpec
    baseline_smiles: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    tailoring_genes: Dict[str, int] = field(default_factory=dict)

    @property
    def baseline(self) -> Optional[str]:
        """antiSMASH's own structure prediction, if it made one."""
        return self.baseline_smiles[0] if self.baseline_smiles else None


def _clean_smiles(raw: str) -> str:
    """Strip the whitespace GenBank line wrapping leaves inside a qualifier."""
    return re.sub(r"\s+", "", raw)


def _consensus_substrate(qualifiers: Dict[str, List[str]]) -> Optional[str]:
    for entry in qualifiers.get("specificity", []):
        match = _CONSENSUS.search(entry)
        if match:
            value = match.group(1).strip()
            if value.lower() in UNKNOWN_SUBSTRATES:
                return None
            return value
    return None


def _is_active(qualifiers: Dict[str, List[str]]) -> bool:
    """False when antiSMASH reports the domain as inactive.

    About one ketoreductase in four is called inactive; such a domain leaves
    the beta-keto group alone.
    """
    for entry in qualifiers.get("specificity", []):
        match = _KR_ACTIVITY.search(entry)
        if match:
            return match.group(1).strip().lower() != "inactive"
    return True


def _transator_call(qualifiers: Dict[str, List[str]]) -> Optional[str]:
    for entry in qualifiers.get("specificity", []):
        match = _TRANSATOR.search(entry)
        if match:
            value = match.group(1).strip().lower()
            return None if value.startswith("(") else value
    return None


@dataclass
class _Domain:
    role: str
    name: str
    locus_tag: str
    start: int
    substrate: Optional[str] = None
    sequence: Optional[str] = None
    index_in_gene: Optional[int] = None
    active: bool = True
    transator: Optional[str] = None


def _within(feature, bounds: Optional[Tuple[int, int]]) -> bool:
    """True when a feature lies inside the region being read."""
    if bounds is None:
        return True
    start, end = bounds
    return int(feature.location.start) >= start and int(feature.location.end) <= end


def _collect_domains(record, bounds: Optional[Tuple[int, int]] = None) -> List[_Domain]:
    """Domains of one record, in assembly-line order.

    Order is gene order along the contig, then position within the protein,
    which is the order antiSMASH itself assumes for modular systems. With
    ``bounds`` only the domains of one region are taken, so a genome-wide
    annotation does not collapse into a single assembly line.
    """
    domains: List[_Domain] = []
    for feature in record.features:
        if feature.type != "aSDomain":
            continue
        if not _within(feature, bounds):
            continue
        if feature.qualifiers.get("aSTool", [""])[0] != "nrps_pks_domains":
            continue
        name = feature.qualifiers.get("aSDomain", [""])[0]
        role = DOMAIN_ROLES.get(name)
        if role in (None, "ignore"):
            continue
        domains.append(
            _Domain(
                role=role,
                name=name,
                locus_tag=feature.qualifiers.get("locus_tag", ["?"])[0],
                start=int(feature.location.start),
                substrate=_consensus_substrate(feature.qualifiers),
                sequence=feature.qualifiers.get("translation", [None])[0],
                active=_is_active(feature.qualifiers),
                transator=_transator_call(feature.qualifiers),
            )
        )
    domains.sort(key=lambda d: d.start)
    # Number the substrate-selecting domains within each gene, so a module can
    # be matched to the same domain in another data source.
    counters: Dict[str, int] = {}
    for domain in domains:
        if domain.role != "a":
            continue
        index = counters.get(domain.locus_tag, 0)
        domain.index_in_gene = index
        counters[domain.locus_tag] = index + 1
    return domains


class _ModuleBuilder:
    """Cuts a flat list of domains into modules."""

    def __init__(self) -> None:
        self.modules: List[Module] = []
        self._current: Optional[Module] = None
        self._has_core: Dict[str, bool] = {}
        self.release: Optional[str] = None
        # A Cy domain sits before the A domain of the module it acts on.
        self._pending_heterocycle = False

    def _flush(self) -> None:
        if self._current is not None:
            self._current.order = len(self.modules)
            self.modules.append(self._current)
        self._current = None
        self._has_core = {}

    def _start(self, kind: str, domain: _Domain) -> None:
        self._flush()
        self._current = Module(kind=kind, gene=domain.locus_tag)

    def add(self, domain: _Domain) -> None:
        role = domain.role
        if self._current is not None:
            self._current.domain_count += 1
            if role == "carrier":
                self._current.has_carrier = True
            elif role in ("c", "cy", "ks"):
                self._current.has_condensation = True
        if role == "ks":
            self._start(PKS, domain)
            self._current.domain_count = 1  # type: ignore[union-attr]
            self._current.has_condensation = True  # type: ignore[union-attr]
            self._has_core["ks"] = True
            if domain.transator:
                call = domain.transator
                self._current.non_elongating = call.startswith("non-elongating")  # type: ignore[union-attr]
                self._current.alpha_methyl = "alphame" in call  # type: ignore[union-attr]
                self._current.reduction_override = TRANSATOR_REDUCTION.get(call)  # type: ignore[union-attr]
            return
        if role == "cy":
            self._pending_heterocycle = True
            return
        if role == "a":
            if self._current is None or self._current.kind == PKS or self._has_core.get("a"):
                self._start(NRPS, domain)
            if self._pending_heterocycle:
                self._current.heterocyclization = True  # type: ignore[union-attr]
                self._pending_heterocycle = False
            if self._current.domain_count == 0:  # type: ignore[union-attr]
                self._current.domain_count = 1  # type: ignore[union-attr]
            self._has_core["a"] = True
            self._current.substrate = domain.substrate  # type: ignore[union-attr]
            self._current.domain_index = domain.index_in_gene  # type: ignore[union-attr]
            self._current.sequence = domain.sequence  # type: ignore[union-attr]
            return
        if role == "at":
            if self._current is None:
                # An AT without a preceding KS starts a PKS loading module.
                self._start(PKS, domain)
                self._current.is_starter = True  # type: ignore[union-attr]
            if self._current is not None and self._current.kind == PKS:
                self._current.substrate = domain.substrate
                self._has_core["at"] = True
            return
        if self._current is None:
            return
        if role == "kr":
            self._current.kr = True
            self._current.kr_active = domain.active
        elif role == "dh":
            self._current.dh = True
        elif role == "er":
            self._current.er = True
        elif role == "e":
            self._current.epimerization = True
        elif role == "mt":
            self._current.n_methylation = True
        elif role == "te":
            self.release = self.release or RELEASE_HYDROLYSIS
        elif role == "td":
            self.release = self.release or RELEASE_REDUCTION

    def finish(self) -> List[Module]:
        self._flush()
        return self.modules


def parse_record(record, bgc_id: str, bounds: Optional[Tuple[int, int]] = None) -> AntismashRecord:
    """Turn one parsed GenBank record into a spec plus the antiSMASH baseline.

    ``bounds`` restricts the read to one region of the record.
    """
    baseline: List[str] = []
    products: List[str] = []
    for feature in record.features:
        if not _within(feature, bounds):
            continue
        if feature.type == "cand_cluster":
            for raw in feature.qualifiers.get("SMILES", []):
                cleaned = _clean_smiles(raw)
                if cleaned and cleaned not in baseline:
                    baseline.append(cleaned)
        if feature.type == "region":
            products.extend(feature.qualifiers.get("product", []))

    from .tailoring import detect as detect_tailoring

    if bounds is None:
        tailoring_genes = detect_tailoring(record).to_dict()
    else:
        inside = type("_Slice", (), {"features": [f for f in record.features if _within(f, bounds)]})()
        tailoring_genes = detect_tailoring(inside).to_dict()

    builder = _ModuleBuilder()
    for domain in _collect_domains(record, bounds):
        builder.add(domain)
    modules = builder.finish()

    # Modules without a carrier protein or a core domain are hmm artefacts;
    # keep the ones that can actually extend a chain.
    notes: List[str] = []
    kept = [m for m in modules if m.kind == NRPS or True]
    if len(kept) != len(modules):  # pragma: no cover - defensive
        notes.append(f"dropped {len(modules) - len(kept)} incomplete modules")

    spec = ClusterSpec(
        bgc_id=bgc_id,
        modules=kept,
        release=builder.release,
        source="antismash",
        classes=sorted(set(products)),
        notes=notes,
    )
    return AntismashRecord(
        bgc_id=bgc_id,
        spec=spec,
        baseline_smiles=baseline,
        products=products,
        tailoring_genes=tailoring_genes,
    )


def load_regions(path: Path) -> List[AntismashRecord]:
    """Every biosynthetic region in a GenBank file, as its own cluster.

    antiSMASH writes one file per region for a single cluster, but a
    genome-wide run produces files holding many regions. Treating those as one
    assembly line would be wrong, so each region is read separately. A file
    with no region features is read whole, which covers MIBiG's own entries.
    """
    path = Path(path)
    found: List[AntismashRecord] = []
    for record in SeqIO.parse(str(path), "genbank"):
        regions = [f for f in record.features if f.type == "region"]
        if not regions:
            parsed = parse_record(record, bgc_id=record.id or path.stem)
            if parsed.spec.modules:
                found.append(parsed)
            continue
        for feature in regions:
            numbers = feature.qualifiers.get("region_number", [])
            label = f"{record.id or path.stem}_region{numbers[0] if numbers else len(found) + 1}"
            bounds = (int(feature.location.start), int(feature.location.end))
            parsed = parse_record(record, bgc_id=label, bounds=bounds)
            if parsed.spec.modules:
                found.append(parsed)
    return found


def load(path: Path) -> AntismashRecord:
    """Read one antiSMASH GenBank file."""
    path = Path(path)
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"no records in {path}")
    merged = parse_record(records[0], bgc_id=path.stem)
    # A cluster split over several records keeps its later modules; this is
    # rare in MIBiG but costs nothing to support.
    for extra in records[1:]:
        more = parse_record(extra, bgc_id=path.stem)
        offset = len(merged.spec.modules)
        for module in more.spec.modules:
            module.order += offset
            merged.spec.modules.append(module)
        merged.baseline_smiles.extend(s for s in more.baseline_smiles if s not in merged.baseline_smiles)
    return merged
