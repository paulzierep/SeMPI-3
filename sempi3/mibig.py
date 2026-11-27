"""Reading MIBiG 4.0 entries.

MIBiG supplies two things SeMPI-3 needs:

* ``biosynthesis.modules`` — curated assembly-line modules, with a SMILES
  string for every adenylation-domain substrate. These are turned into a
  :class:`~sempi3.spec.ClusterSpec`.
* ``compounds[].structure`` — the SMILES of the known product, used as ground
  truth.

Note on what this measures: the curated substrates are an *annotation*, partly
inferred from the known structure. Assembling from them therefore measures the
assembly chemistry, not de-novo prediction from a genome. Section "Evaluation
design" of the README spells this out.
"""

from __future__ import annotations

import json
from dataclasses import (
    dataclass,
    field,
)
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
)

from .spec import (
    ClusterSpec,
    Module,
    NRPS,
    PKS,
    RELEASE_HYDROLYSIS,
)

#: MIBiG module types that belong to an NRPS assembly line.
NRPS_MODULE_TYPES = {"nrps-type1", "nrps-type6"}

#: MIBiG module types that belong to a PKS assembly line.
PKS_MODULE_TYPES = {
    "pks-modular",
    "pks-modular-starter",
    "pks-trans-at",
    "pks-trans-at-starter",
    "pks-iterative",
}

#: Modification domains mapped onto module flags.
KETOREDUCTASE = "ketoreductase"
DEHYDRATASE = "dehydratase"
ENOYLREDUCTASE = "enoylreductase"
METHYLTRANSFERASE = "methyltransferase"
EPIMERASE = "epimerase"
THIOESTERASE = "thioesterase"
THIOREDUCTASE = "thioreductase"


@dataclass
class MibigEntry:
    """One MIBiG cluster: its curated modules and its known products."""

    accession: str
    classes: List[str]
    compounds: List[Dict[str, Any]] = field(default_factory=list)
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def structures(self) -> List[str]:
        """SMILES of every compound with a structure."""
        return [c["structure"] for c in self.compounds if c.get("structure")]

    @property
    def compound_names(self) -> List[str]:
        return [c.get("name", "") for c in self.compounds if c.get("structure")]

    @property
    def is_pks_nrps(self) -> bool:
        return bool({"PKS", "NRPS"} & set(self.classes))


def load_entry(path: Path) -> MibigEntry:
    """Read a single MIBiG JSON file."""
    with open(path) as handle:
        data = json.load(handle)
    classes = [c.get("class", "") for c in data.get("biosynthesis", {}).get("classes", [])]
    return MibigEntry(
        accession=data["accession"],
        classes=classes,
        compounds=data.get("compounds", []),
        raw=data,
    )


def iter_entries(directory: Path) -> Iterator[MibigEntry]:
    """Read every MIBiG JSON file in ``directory``, sorted by accession."""
    for path in sorted(Path(directory).glob("*.json")):
        yield load_entry(path)


def _modification_types(module: Dict[str, Any]) -> List[str]:
    out = []
    for entry in module.get("modification_domains") or []:
        if isinstance(entry, dict):
            if entry.get("active") is False:
                continue
            out.append(entry.get("type", ""))
        else:
            out.append(str(entry))
    return out


def _substrate(module: Dict[str, Any], notes: List[str], order: int) -> Optional[str]:
    """Substrate of an adenylation domain, preferring the curated structure."""
    domain = module.get("a_domain") or {}
    substrates = domain.get("substrates") or []
    if not substrates:
        return None
    if len(substrates) > 1:
        names = ", ".join(s.get("name", "?") for s in substrates)
        notes.append(f"module {order}: {len(substrates)} candidate substrates ({names}), used the first")
    first = substrates[0]
    if first.get("structure"):
        return "smiles:" + first["structure"]
    return first.get("name")


def to_spec(entry: MibigEntry) -> ClusterSpec:
    """Build a :class:`ClusterSpec` from a MIBiG entry's curated modules."""
    notes: List[str] = []
    modules: List[Module] = []
    release: Optional[str] = None

    raw_modules = entry.raw.get("biosynthesis", {}).get("modules", [])
    for order, raw in enumerate(raw_modules):
        mod_type = raw.get("type", "")
        if mod_type in NRPS_MODULE_TYPES:
            kind = NRPS
        elif mod_type in PKS_MODULE_TYPES:
            kind = PKS
        else:
            notes.append(f"module {order}: unsupported type {mod_type!r}, skipped")
            continue
        if raw.get("active") is False:
            notes.append(f"module {order}: inactive, skipped")
            continue

        modifications = _modification_types(raw)
        if THIOESTERASE in modifications:
            release = release or RELEASE_HYDROLYSIS
        genes = raw.get("genes") or []
        modules.append(
            Module(
                kind=kind,
                substrate=_substrate(raw, notes, order) if kind == NRPS else None,
                is_starter=mod_type.endswith("-starter"),
                kr=KETOREDUCTASE in modifications,
                dh=DEHYDRATASE in modifications,
                er=ENOYLREDUCTASE in modifications,
                epimerization=EPIMERASE in modifications,
                n_methylation=METHYLTRANSFERASE in modifications,
                gene=genes[0] if genes else None,
                order=order,
            )
        )

    return ClusterSpec(
        bgc_id=entry.accession,
        modules=modules,
        release=release,
        source="mibig",
        classes=entry.classes,
        notes=notes,
    )
