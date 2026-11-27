"""Data model for biosynthetic gene cluster assembly lines.

A :class:`ClusterSpec` is the common intermediate representation of SeMPI-3.
Every input source (antiSMASH output, MIBiG curated annotation, a hand written
test case) is parsed into this model, and the assembly engine only ever sees
this model. That keeps the chemistry independent of the file formats.
"""

from __future__ import annotations

from dataclasses import (
    dataclass,
    field,
)
from typing import (
    Any,
    Dict,
    List,
    Optional,
)

#: Module kinds understood by the assembly engine.
PKS = "pks"
NRPS = "nrps"

#: Chain release mechanisms.
RELEASE_HYDROLYSIS = "hydrolysis"
RELEASE_REDUCTION = "reduction"
RELEASE_MACROLACTAM = "macrolactam"
RELEASE_MACROLACTONE = "macrolactone"


@dataclass
class Module:
    """One elongation (or loading) module of an assembly line.

    Attributes:
        kind: :data:`PKS` or :data:`NRPS`.
        substrate: Substrate name as reported by the source, e.g. ``"mmal"``
            for methylmalonyl-CoA or ``"ala"`` for alanine. Resolved against
            :mod:`sempi3.monomers`; ``None`` means the source gave no
            prediction.
        is_starter: True for the loading module, which does not elongate an
            existing chain.
        kr, dh, er: Reductive loop domains present in a PKS module. They act on
            the beta-keto group created by *this* module's condensation.
        epimerization: An E domain follows the condensation, giving a
            D-configured residue (NRPS only).
        heterocyclization: A Cy domain replaces the condensation domain, so the
            residue is closed into a thiazoline or oxazoline ring.
        n_methylation: An nMT domain methylates the backbone nitrogen (NRPS).
        gene: Locus tag or protein id the module was found on, for reporting.
        order: Position in the assembly line, starting at 0.
        confidence: Optional score from the prediction source.
        domain_index: Index of the substrate-selecting domain within its gene,
            counting from zero. Together with ``gene`` it identifies the domain
            across data sources.
        sequence: Amino acid sequence of that domain, when the source gives it.
    """

    kind: str
    substrate: Optional[str] = None
    is_starter: bool = False
    kr: bool = False
    dh: bool = False
    er: bool = False
    epimerization: bool = False
    heterocyclization: bool = False
    n_methylation: bool = False
    gene: Optional[str] = None
    order: int = 0
    confidence: Optional[float] = None
    domain_index: Optional[int] = None
    sequence: Optional[str] = None

    def __post_init__(self) -> None:
        if self.kind not in (PKS, NRPS):
            raise ValueError(f"unknown module kind: {self.kind!r}")
        # A dehydratase without a ketoreductase cannot act, and an enoyl
        # reductase needs both. Sources occasionally report such domains;
        # normalising here keeps the assembly rules simple.
        if self.dh and not self.kr:
            self.dh = False
        if self.er and not self.dh:
            self.er = False

    @property
    def reduction_state(self) -> str:
        """Fate of the beta-keto group: ``keto``, ``hydroxyl``, ``ene`` or ``alkyl``."""
        if not self.kr:
            return "keto"
        if not self.dh:
            return "hydroxyl"
        if not self.er:
            return "ene"
        return "alkyl"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "kind": self.kind,
            "substrate": self.substrate,
            "is_starter": self.is_starter,
            "kr": self.kr,
            "dh": self.dh,
            "er": self.er,
            "epimerization": self.epimerization,
            "heterocyclization": self.heterocyclization,
            "n_methylation": self.n_methylation,
            "gene": self.gene,
            "order": self.order,
            "confidence": self.confidence,
            "domain_index": self.domain_index,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Module":
        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        return cls(**{k: v for k, v in data.items() if k in known})


@dataclass
class ClusterSpec:
    """An ordered assembly line extracted from one gene cluster.

    Attributes:
        bgc_id: Identifier of the cluster, e.g. a MIBiG accession.
        modules: Elongation modules in biosynthetic order.
        release: Chain release mechanism, if one could be inferred.
        source: Where the specification came from, e.g. ``"antismash"`` or
            ``"mibig"``. Recorded so benchmark results stay traceable.
        classes: Biosynthetic classes reported by the source.
        notes: Free-form warnings collected while parsing.
    """

    bgc_id: str
    modules: List[Module] = field(default_factory=list)
    release: Optional[str] = None
    source: str = "unknown"
    classes: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def n_modules(self) -> int:
        return len(self.modules)

    @property
    def hybrid(self) -> bool:
        kinds = {m.kind for m in self.modules}
        return kinds == {PKS, NRPS}

    @property
    def complete(self) -> bool:
        """True when every module has a substrate prediction."""
        return bool(self.modules) and all(m.substrate for m in self.modules)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "bgc_id": self.bgc_id,
            "source": self.source,
            "classes": self.classes,
            "release": self.release,
            "notes": self.notes,
            "modules": [m.to_dict() for m in self.modules],
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ClusterSpec":
        return cls(
            bgc_id=data["bgc_id"],
            modules=[Module.from_dict(m) for m in data.get("modules", [])],
            release=data.get("release"),
            source=data.get("source", "unknown"),
            classes=list(data.get("classes", [])),
            notes=list(data.get("notes", [])),
        )
