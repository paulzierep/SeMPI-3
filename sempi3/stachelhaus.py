"""The Stachelhaus code: the binding-pocket residues of an adenylation domain.

Substrate specificity is set by a handful of residues lining the binding
pocket, not by the domain as a whole. Reading those residues off and comparing
them is the classical way to predict what an A domain activates
(Stachelhaus et al. 1999; Challis et al. 2000), and it is what SeMPI and
NRPSpredictor are built on.

The positions are defined on GrsA (gramicidin S synthetase A), whose A domain
was crystallised with phenylalanine. A query domain is aligned to that
reference and the residues facing the pocket are read at the aligned
positions.

The reference sequence is taken from MIBiG entry BGC0002122 (gramicidin S),
locus ``AF333_17650``; at the positions below it reads ``DAWTIAAIC``, the
published code for phenylalanine, which is the check that the numbering is
right.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
)

from Bio import Align
from Bio.Align import substitution_matrices

#: Pocket positions in GrsA numbering, 1-based.
POCKET_POSITIONS = (235, 236, 239, 278, 299, 301, 322, 330, 331)

#: The code the reference itself must yield, as a self-check.
REFERENCE_CODE = "DAWTIAAIC"

#: The wider binding pocket, 34 positions in GrsA numbering.
#:
#: This is the "extended signature" that NRPSpredictor2 and PARAS use. The
#: positions were recovered from the PARAS dataset: for 400 A domains present
#: in both MIBiG and PARAS, each signature column was matched against this
#: reference by alignment, and the position that reproduced the column most
#: often was taken. Every Stachelhaus position falls inside the resulting set,
#: which is the check that the recovery worked.
EXTENDED_POSITIONS = (
    210, 213, 214, 230, 234, 235, 236, 237, 238, 239, 240, 243,
    278, 279, 299, 300, 301, 302, 303, 320, 321, 322, 323, 324,
    325, 326, 327, 328, 329, 330, 331, 332, 333, 334,
)

#: Indices inside a PARAS extended signature that hold the Stachelhaus code.
SIGNATURE_CODE_INDICES = (5, 6, 9, 12, 14, 16, 21, 29, 30)

#: Character used where the alignment leaves a position uncovered.
GAP = "-"

_REFERENCE_PATH = Path(__file__).parent / "data" / "grsa_reference.txt"


@lru_cache(maxsize=1)
def reference_sequence() -> str:
    """The GrsA reference sequence shipped with the package."""
    return _REFERENCE_PATH.read_text().strip()


@lru_cache(maxsize=1)
def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def code_from_signature(signature: str) -> str:
    """Pull the Stachelhaus code out of a 34-residue extended signature."""
    if not signature or len(signature) != len(EXTENDED_POSITIONS):
        return GAP * len(POCKET_POSITIONS)
    return "".join(signature[i] for i in SIGNATURE_CODE_INDICES)


def _residues_at(sequence: str, reference: Optional[str], positions: Sequence[int]) -> str:
    reference = reference or reference_sequence()
    sequence = (sequence or "").strip().upper().replace("*", "")
    if len(sequence) < 50:
        return GAP * len(positions)
    cleaned = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in sequence)
    try:
        alignment = _aligner().align(reference, cleaned)[0]
    except (ValueError, IndexError):  # pragma: no cover - defensive
        return GAP * len(positions)
    mapping: Dict[int, int] = {}
    for ref_block, query_block in zip(alignment.aligned[0], alignment.aligned[1]):
        ref_start, ref_end = ref_block
        query_start, _query_end = query_block
        for offset in range(ref_end - ref_start):
            mapping[ref_start + offset] = query_start + offset
    out = []
    for position in positions:
        index = mapping.get(position - 1)
        out.append(cleaned[index] if index is not None else GAP)
    return "".join(out)


def extract_signature(sequence: str, reference: Optional[str] = None) -> str:
    """Read the 34-residue extended binding-pocket signature."""
    return _residues_at(sequence, reference, EXTENDED_POSITIONS)


def extract_code(sequence: str, reference: Optional[str] = None) -> str:
    """Read the pocket residues out of one A-domain sequence.

    Args:
        sequence: Amino acid sequence of the domain, which may be a fragment of
            the full protein.
        reference: Override the GrsA reference, for testing.

    Returns:
        One character per pocket position, in order, with ``-`` where the
        alignment does not cover the position.
    """
    reference = reference or reference_sequence()
    sequence = (sequence or "").strip().upper().replace("*", "")
    if len(sequence) < 50:
        return GAP * len(POCKET_POSITIONS)
    # Unknown residues upset the substitution matrix; map them to the
    # matrix's own wildcard.
    cleaned = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in sequence)
    try:
        alignment = _aligner().align(reference, cleaned)[0]
    except (ValueError, IndexError):  # pragma: no cover - defensive
        return GAP * len(POCKET_POSITIONS)

    # Map every reference index that the alignment covers to a query index.
    mapping: Dict[int, int] = {}
    for ref_block, query_block in zip(alignment.aligned[0], alignment.aligned[1]):
        ref_start, ref_end = ref_block
        query_start, _query_end = query_block
        for offset in range(ref_end - ref_start):
            mapping[ref_start + offset] = query_start + offset

    code = []
    for position in POCKET_POSITIONS:
        query_index = mapping.get(position - 1)
        code.append(cleaned[query_index] if query_index is not None else GAP)
    return "".join(code)


def self_check() -> bool:
    """True when the reference reproduces the published phenylalanine code."""
    reference = reference_sequence()
    direct = "".join(reference[p - 1] for p in POCKET_POSITIONS)
    return direct == REFERENCE_CODE


@dataclass
class CodeMatch:
    """A predicted substrate and how close the matching code was."""

    substrate: Optional[str]
    identity: float
    neighbours: int

    def to_dict(self) -> Dict[str, object]:
        return {
            "substrate": self.substrate,
            "identity": round(self.identity, 3),
            "neighbours": self.neighbours,
        }


def _identity(left: str, right: str) -> float:
    positions = [(a, b) for a, b in zip(left, right) if a != GAP and b != GAP]
    if not positions:
        return 0.0
    return sum(1 for a, b in positions if a == b) / len(positions)


class StachelhausPredictor:
    """Nearest-code predictor, the classical Stachelhaus approach.

    Attributes:
        min_identity: Codes less similar than this are not trusted; the
            predictor abstains rather than inventing a residue, because a wrong
            monomer costs the backbone more than a missing one.
    """

    def __init__(self, min_identity: float = 0.75) -> None:
        self.min_identity = min_identity
        self.table: Dict[str, List[str]] = {}

    def fit(self, codes: Sequence[str], labels: Sequence[str]) -> "StachelhausPredictor":
        self.table = {}
        for code, label in zip(codes, labels):
            if not code or code.count(GAP) > 3:
                continue
            self.table.setdefault(code, []).append(label)
        return self

    def predict_one(self, code: str) -> CodeMatch:
        if not code or code.count(GAP) > 3 or not self.table:
            return CodeMatch(None, 0.0, 0)
        exact = self.table.get(code)
        if exact:
            return CodeMatch(_vote(exact), 1.0, len(exact))

        best_identity = 0.0
        best_labels: List[str] = []
        for known, labels in self.table.items():
            identity = _identity(code, known)
            if identity > best_identity:
                best_identity = identity
                best_labels = list(labels)
            elif identity == best_identity:
                best_labels.extend(labels)
        if best_identity < self.min_identity:
            return CodeMatch(None, best_identity, 0)
        return CodeMatch(_vote(best_labels), best_identity, len(best_labels))

    def predict(self, codes: Sequence[str]) -> List[CodeMatch]:
        return [self.predict_one(code) for code in codes]


def _vote(labels: Sequence[str]) -> Optional[str]:
    if not labels:
        return None
    counts: Dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return max(counts.items(), key=lambda item: item[1])[0]
