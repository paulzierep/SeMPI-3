"""Building blocks: the token language between a model and the chemistry.

A product is written as a short sequence of tokens — one per module, then the
release, then any decorations:

    PK:mmal:hydroxyl  PK:mal:ene  AA:leucine:D  REL:macrolactone  TAIL:glycosyl

A model generates that sequence; :func:`blocks_to_spec` turns it back into an
assembly line and the usual engine builds the molecule. Generating blocks
rather than SMILES means every output is a chemically valid structure, and the
model only has to learn the biology.

It is still generative: a sequence can be produced that no cluster in the
training set has. What it cannot do is invent a scaffold the assembly rules
do not cover, which is the price of guaranteed validity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
)

from .spec import (
    ClusterSpec,
    Module,
    NRPS,
    PKS,
)

#: Oxidation levels a polyketide module can leave behind.
REDUCTION_STATES = ("keto", "hydroxyl", "ene", "alkyl")

#: Extender units a polyketide module can use.
EXTENDERS = ("mal", "mmal", "emal", "mxmal")

#: Prefixes of the token language.
TOKEN_PKS = "PK"
TOKEN_NRPS = "AA"
TOKEN_RELEASE = "REL"
TOKEN_TAILORING = "TAIL"


@dataclass
class Block:
    """One module's contribution, as a token."""

    kind: str
    substrate: Optional[str] = None
    reduction: Optional[str] = None
    epimerized: bool = False
    n_methyl: bool = False
    heterocycle: bool = False
    non_elongating: bool = False
    is_starter: bool = False

    def token(self) -> str:
        if self.kind == PKS:
            parts = [TOKEN_PKS, self.substrate or "mal", self.reduction or "keto"]
            if self.non_elongating:
                parts.append("ne")
            if self.is_starter:
                parts.append("start")
        else:
            parts = [TOKEN_NRPS, (self.substrate or "unknown").replace(" ", "_")]
            if self.epimerized:
                parts.append("D")
            if self.n_methyl:
                parts.append("NMe")
            if self.heterocycle:
                parts.append("cy")
        return ":".join(parts)

    @classmethod
    def from_token(cls, token: str) -> "Block":
        parts = token.split(":")
        head = parts[0]
        flags = set(parts[3:]) if head == TOKEN_PKS else set(parts[2:])
        if head == TOKEN_PKS:
            return cls(
                kind=PKS,
                substrate=parts[1] if len(parts) > 1 else "mal",
                reduction=parts[2] if len(parts) > 2 else "keto",
                non_elongating="ne" in flags,
                is_starter="start" in flags,
            )
        return cls(
            kind=NRPS,
            substrate=parts[1].replace("_", " ") if len(parts) > 1 else None,
            epimerized="D" in flags,
            n_methyl="NMe" in flags,
            heterocycle="cy" in flags,
        )

    def to_module(self, order: int) -> Module:
        """Turn the block back into a module the assembly engine understands."""
        if self.kind == PKS:
            return Module(
                kind=PKS,
                substrate=self.substrate,
                reduction_override=self.reduction,
                non_elongating=self.non_elongating,
                is_starter=self.is_starter,
                order=order,
            )
        return Module(
            kind=NRPS,
            substrate=self.substrate,
            epimerization=self.epimerized,
            n_methylation=self.n_methyl,
            heterocyclization=self.heterocycle,
            order=order,
        )


def module_to_block(module: Module) -> Block:
    """Describe a parsed module as a block."""
    if module.kind == PKS:
        return Block(
            kind=PKS,
            substrate=(module.substrate or "mal").lower(),
            reduction=module.reduction_state,
            non_elongating=module.non_elongating,
            is_starter=module.is_starter,
        )
    return Block(
        kind=NRPS,
        substrate=(module.substrate or "unknown").replace("smiles:", ""),
        epimerized=module.epimerization,
        n_methyl=module.n_methylation,
        heterocycle=module.heterocyclization,
    )


def spec_to_blocks(spec: ClusterSpec) -> List[Block]:
    return [module_to_block(module) for module in spec.modules]


def spec_to_tokens(spec: ClusterSpec, tailoring: Optional[Sequence[str]] = None) -> List[str]:
    """The full token sequence for a cluster, release and decorations included."""
    tokens = [block.token() for block in spec_to_blocks(spec)]
    if spec.release:
        tokens.append(f"{TOKEN_RELEASE}:{spec.release}")
    for name in tailoring or []:
        tokens.append(f"{TOKEN_TAILORING}:{name}")
    return tokens


def blocks_to_spec(
    blocks: Sequence[Block],
    bgc_id: str,
    release: Optional[str] = None,
    source: str = "generated",
) -> ClusterSpec:
    """Build an assembly line out of generated blocks."""
    modules = [block.to_module(order) for order, block in enumerate(blocks)]
    return ClusterSpec(bgc_id=bgc_id, modules=modules, release=release, source=source)


def tokens_to_spec(tokens: Sequence[str], bgc_id: str) -> Tuple[ClusterSpec, List[str]]:
    """Parse a generated token sequence into a spec and its decorations."""
    blocks: List[Block] = []
    release: Optional[str] = None
    tailoring: List[str] = []
    for token in tokens:
        if token.startswith(f"{TOKEN_RELEASE}:"):
            release = token.split(":", 1)[1]
        elif token.startswith(f"{TOKEN_TAILORING}:"):
            tailoring.append(token.split(":", 1)[1])
        elif token.startswith((f"{TOKEN_PKS}:", f"{TOKEN_NRPS}:")):
            blocks.append(Block.from_token(token))
    return blocks_to_spec(blocks, bgc_id, release=release), tailoring


def vocabulary(specs: Sequence[ClusterSpec]) -> Dict[str, int]:
    """Token counts over a set of clusters, for inspecting the language."""
    counts: Dict[str, int] = {}
    for spec in specs:
        for token in spec_to_tokens(spec):
            counts[token] = counts.get(token, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: -item[1]))
