"""Rule-based assembly of a backbone structure from a :class:`ClusterSpec`.

The rules follow textbook assembly-line logic:

* The first module primes the chain. An NRPS loading module installs its own
  amino acid; a PKS loading module installs a starter unit (acetyl unless the
  source says otherwise).
* Every further NRPS module forms an amide bond; every further PKS module adds
  a two-carbon unit and sets the oxidation level of the beta carbon from its
  reductive domains.
* The chain is released by hydrolysis, reduction, macrolactamisation or
  macrolactonisation.

What comes out is the *backbone*: tailoring enzymes acting after release
(glycosylation, halogenation, oxidative crosslinking) are not modelled, so the
prediction is deliberately an under-estimate of the final natural product.
"""

from __future__ import annotations

from dataclasses import (
    dataclass,
    field,
)
from typing import (
    Dict,
    List,
    Optional,
)

from . import monomers
from .chem import (
    Chain,
    ChemError,
    append_amino_acid,
    append_heterocycle,
    append_pks_unit,
    start_chain,
    terminate,
    to_smiles,
)
from .spec import (
    ClusterSpec,
    Module,
    NRPS,
    PKS,
    RELEASE_HYDROLYSIS,
    RELEASE_MACROLACTAM,
    RELEASE_MACROLACTONE,
    RELEASE_REDUCTION,
)


@dataclass
class Prediction:
    """The outcome of assembling one cluster.

    Attributes:
        bgc_id: Cluster the prediction belongs to.
        smiles: Backbone for the release mechanism that was chosen.
        release: The release mechanism SeMPI-3 committed to.
        variants: SMILES for every release mechanism that could be built,
            keyed by mechanism. Reported so a benchmark can separate errors in
            the backbone from errors in the release step.
        n_modules: Number of modules that were used.
        notes: Warnings, e.g. modules skipped for want of a substrate.
        failed: True when no structure could be built at all.
    """

    bgc_id: str
    smiles: Optional[str] = None
    variants: Dict[str, str] = field(default_factory=dict)
    n_modules: int = 0
    release: Optional[str] = None
    notes: List[str] = field(default_factory=list)
    failed: bool = False

    def to_dict(self) -> Dict[str, object]:
        return {
            "bgc_id": self.bgc_id,
            "smiles": self.smiles,
            "release": self.release,
            "variants": self.variants,
            "n_modules": self.n_modules,
            "notes": self.notes,
            "failed": self.failed,
        }


class AssemblyError(RuntimeError):
    """Raised when a cluster cannot be assembled."""


def _substrate_smiles(module: Module) -> Optional[str]:
    """SMILES of a module's substrate, from a structure or a name."""
    substrate = module.substrate
    if not substrate:
        return None
    if substrate.startswith("smiles:"):
        return substrate[len("smiles:") :]
    return monomers.resolve(substrate)


def choose_release(spec: ClusterSpec, variants: Dict[str, str]) -> str:
    """Pick the one release mechanism to report for this cluster.

    Reporting the best of several variants flatters a predictor: it is an
    oracle that a real user does not have. SeMPI-3 therefore commits to one
    mechanism, chosen from the domains:

    * a terminal reductase releases an aldehyde;
    * a thioesterase on an all-NRPS line usually closes a macrolactam, and on
      a line containing PKS modules a macrolactone;
    * without either domain the chain is hydrolysed off.

    The oracle is still computed by the benchmark, as an upper bound.
    """
    if spec.release == RELEASE_REDUCTION and RELEASE_REDUCTION in variants:
        return RELEASE_REDUCTION

    kinds = {module.kind for module in spec.modules}
    if kinds == {NRPS}:
        preference = [RELEASE_MACROLACTAM, RELEASE_MACROLACTONE, RELEASE_HYDROLYSIS]
    else:
        preference = [RELEASE_MACROLACTONE, RELEASE_MACROLACTAM, RELEASE_HYDROLYSIS]

    for mode in preference:
        if mode in variants:
            return mode
    return next(iter(variants))


def _release_modes(spec: ClusterSpec, chain: Chain) -> List[str]:
    """Release mechanisms worth building for this cluster."""
    if spec.release:
        modes = [spec.release]
    else:
        modes = [RELEASE_HYDROLYSIS]
    # Macrocyclisation is extremely common in both families, and the source
    # rarely says which one applies, so both are offered as variants.
    for extra in (RELEASE_MACROLACTAM, RELEASE_MACROLACTONE, RELEASE_HYDROLYSIS, RELEASE_REDUCTION):
        if extra not in modes:
            modes.append(extra)
    return modes


#: Release label from a structure-trained model -> the variant to report.
RELEASE_LABELS = {
    "linear": RELEASE_HYDROLYSIS,
    "macrolactone": RELEASE_MACROLACTONE,
    "macrolactam": RELEASE_MACROLACTAM,
}


def assemble(
    spec: ClusterSpec,
    assume_acetyl_starter: bool = True,
    skip_unknown_substrates: bool = False,
    release_override: Optional[str] = None,
) -> Prediction:
    """Build the backbone encoded by ``spec``.

    Args:
        spec: The assembly line to build.
        assume_acetyl_starter: When the first module is a PKS elongation
            module, prime the chain with acetyl, as most assembly lines do.
        skip_unknown_substrates: Carry on when an NRPS module has no substrate,
            leaving a gap, instead of giving up on the cluster.
        release_override: Release mechanism chosen by a model, either a variant
            name or a label such as ``"macrolactone"``. Falls back to the
            domain rule when the named variant could not be built.
    """
    prediction = Prediction(bgc_id=spec.bgc_id, notes=list(spec.notes))
    if not spec.modules:
        prediction.failed = True
        prediction.notes.append("no modules")
        return prediction

    chain: Optional[Chain] = None
    used = 0
    for module in spec.modules:
        try:
            chain, consumed = _apply_module(
                chain,
                module,
                assume_acetyl_starter=assume_acetyl_starter,
                skip_unknown_substrates=skip_unknown_substrates,
                notes=prediction.notes,
            )
        except ChemError as exc:
            prediction.notes.append(f"module {module.order}: {exc}")
            prediction.failed = True
            return prediction
        used += int(consumed)

    if chain is None:
        prediction.failed = True
        prediction.notes.append("no module could start a chain")
        return prediction

    prediction.n_modules = used
    for mode in _release_modes(spec, chain):
        try:
            prediction.variants[mode] = to_smiles(terminate(chain, mode))
        except (ChemError, Exception) as exc:  # noqa: BLE001 - RDKit raises broadly
            prediction.notes.append(f"release {mode}: {exc}")
    if not prediction.variants:
        prediction.failed = True
        return prediction
    override = RELEASE_LABELS.get(release_override or "", release_override)
    if override and override in prediction.variants:
        prediction.release = override
    else:
        prediction.release = choose_release(spec, prediction.variants)
    prediction.smiles = prediction.variants[prediction.release]
    return prediction


def _apply_module(
    chain: Optional[Chain],
    module: Module,
    assume_acetyl_starter: bool,
    skip_unknown_substrates: bool,
    notes: List[str],
):
    """Apply one module, returning the new chain and whether it was used."""
    substrate = _substrate_smiles(module)

    if module.kind == NRPS:
        if not substrate:
            if skip_unknown_substrates:
                notes.append(f"module {module.order}: no substrate, skipped")
                return chain, False
            raise ChemError(f"no substrate for NRPS module {module.order}")
        if module.epimerization:
            substrate = monomers.invert_alpha_carbon(substrate)
        if chain is None:
            if module.n_methylation:
                notes.append(f"module {module.order}: N-methylation on the starter is not modelled")
            return start_chain(substrate), True
        if module.heterocyclization:
            # A Cy domain only closes a ring when the residue has a thiol or
            # hydroxyl side chain; otherwise it just condenses.
            try:
                return append_heterocycle(chain, substrate), True
            except ChemError as exc:
                notes.append(f"module {module.order}: no heterocycle ({exc})")
        return append_amino_acid(chain, substrate, n_methylation=module.n_methylation), True

    # PKS module.
    branch = monomers.extender_branch(module.substrate)
    if branch is None and module.n_methylation:
        # A C-methyltransferase in a PKS module methylates the alpha carbon.
        branch = "C"
    if chain is None:
        starter = monomers.starter_smiles(module.substrate if module.is_starter else None)
        chain = start_chain(starter)
        if module.is_starter:
            return chain, True
        if not assume_acetyl_starter:
            raise ChemError("first PKS module is not a loading module")
        notes.append("assumed an acetyl starter unit")
        return append_pks_unit(chain, module.reduction_state, alpha_substituent=branch), True
    if module.is_starter:
        notes.append(f"module {module.order}: loading module after the chain started, treated as elongation")
    return append_pks_unit(chain, module.reduction_state, alpha_substituent=branch), True
