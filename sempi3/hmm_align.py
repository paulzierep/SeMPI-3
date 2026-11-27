"""Reading binding-pocket residues with a profile HMM.

:mod:`sempi3.stachelhaus` aligns a query adenylation domain to one reference
protein. That works, but a pairwise alignment to a single sequence drifts on
domains far from the reference, which is exactly where antiSMASH also fails
and where a substrate call is most needed.

Aligning to a profile instead is steadier: against PARAS's own signatures, the
profile route agrees on 95.6% of residues where the pairwise route manages
91.2%. The profile is the Pfam AMP-binding model that PARAS ships (CC BY 4.0),
and the alignment is done by ``hmmalign`` from HMMER.

HMMER is optional. When it is missing the pairwise route is used, so the
package still works with nothing but RDKit and Biopython.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
)

from . import stachelhaus as st

#: The profile, shipped with the package.
HMM_FILE = Path(__file__).parent / "data" / "AMP-binding.hmm"

#: Minimum length before a sequence is worth aligning.
MIN_LENGTH = 50


@lru_cache(maxsize=1)
def hmmalign_binary() -> Optional[str]:
    """Path to ``hmmalign``, or ``None`` when HMMER is not installed."""
    return shutil.which("hmmalign")


@lru_cache(maxsize=1)
def available() -> bool:
    """True when profile alignment can be used."""
    return bool(hmmalign_binary()) and HMM_FILE.exists()


def _match_states(a2m: str) -> Dict[int, str]:
    """Residues of one A2M row, keyed by the match state they occupy.

    In A2M a lowercase letter is an insertion, which belongs to no match
    state; an uppercase letter or a gap advances the match-state counter.
    Match states are shared by every row, which is what makes residues from
    different sequences comparable.
    """
    residues: Dict[int, str] = {}
    state = 0
    for character in a2m:
        if character.islower() or character == ".":
            continue
        state += 1
        if character != "-":
            residues[state] = character.upper()
    return residues


def _run_hmmalign(sequences: Dict[str, str]) -> Dict[str, str]:
    """Align sequences to the profile, returning one A2M row per name."""
    binary = hmmalign_binary()
    if not binary:
        raise RuntimeError("hmmalign is not installed")

    with tempfile.TemporaryDirectory() as directory:
        fasta = Path(directory) / "input.fasta"
        with open(fasta, "w") as handle:
            for name, sequence in sequences.items():
                handle.write(f">{name}\n{sequence}\n")
        result = subprocess.run(
            [binary, "--amino", "--outformat", "A2M", str(HMM_FILE), str(fasta)],
            capture_output=True,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"hmmalign failed: {result.stderr.strip()[:200]}")

    rows: Dict[str, str] = {}
    name: Optional[str] = None
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            rows[name] = ""
        elif name is not None:
            rows[name] += line.strip()
    return rows


def extract_signatures(sequences: Sequence[str]) -> List[str]:
    """Binding-pocket signatures for many domains, in one alignment.

    The reference is aligned alongside the queries, so the pocket positions —
    which are defined on the reference — can be followed into every query by
    the match state they land on.

    Sequences too short to align, or that the profile does not cover, come
    back as all gaps, and the caller can fall back to the pairwise route.
    """
    reference = st.reference_sequence()
    usable = {
        f"q{index}": sequence
        for index, sequence in enumerate(sequences)
        if sequence and len(sequence) >= MIN_LENGTH
    }
    blank = st.GAP * len(st.EXTENDED_POSITIONS)
    if not usable:
        return [blank] * len(sequences)

    rows = _run_hmmalign({"reference": reference, **usable})
    reference_row = rows.get("reference")
    if not reference_row:
        return [blank] * len(sequences)

    # Which match state does each pocket position of the reference occupy?
    state_of_residue: Dict[int, int] = {}
    state = 0
    residue = 0
    for character in reference_row:
        if character.islower() or character == ".":
            residue += 1
            continue
        state += 1
        if character != "-":
            residue += 1
            state_of_residue[residue] = state
    target_states = [state_of_residue.get(position) for position in st.EXTENDED_POSITIONS]

    out: List[str] = []
    for index in range(len(sequences)):
        row = rows.get(f"q{index}")
        if not row:
            out.append(blank)
            continue
        residues = _match_states(row)
        out.append(
            "".join(
                residues.get(target, st.GAP) if target else st.GAP for target in target_states
            )
        )
    return out


def extract_signature(sequence: str) -> str:
    """Signature of one domain, by profile alignment."""
    return extract_signatures([sequence])[0]


def signature(sequence: str, prefer_hmm: bool = True) -> Tuple[str, str]:
    """Best available signature, with the method that produced it.

    Returns ``(signature, method)`` where method is ``"hmm"`` or
    ``"pairwise"``. A profile alignment that comes back mostly gaps is
    discarded in favour of the pairwise one.
    """
    if prefer_hmm and available():
        try:
            found = extract_signature(sequence)
        except (RuntimeError, OSError):
            found = st.GAP * len(st.EXTENDED_POSITIONS)
        if found.count(st.GAP) <= 3:
            return found, "hmm"
    return st.extract_signature(sequence), "pairwise"
