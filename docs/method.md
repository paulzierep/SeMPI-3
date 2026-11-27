# How SeMPI-3 works

[← back to the README](../README.md)

## How it works

```mermaid
flowchart TD
    MIBIG[("MIBiG 4.0<br/>3,013 curated clusters")]
    AS[("antiSMASH 8 output<br/>published by MIBiG")]
    PARAS[("PARAS<br/>3,653 curated A domains")]

    AS --> PARSE["Read aSDomain features<br/>order by gene, then by position"]
    PARSE --> MOD["Cut into modules<br/>a ketosynthase starts a PKS module<br/>an adenylation domain starts an NRPS module"]

    MOD --> SCOPE{"A modular<br/>assembly line?"}
    SCOPE -- "no: iterative or type III" --> FLAG(["Reported, flagged,<br/>not built"])
    SCOPE -- yes --> Q{"Did antiSMASH name<br/>the substrate?"}
    Q -- "yes, 83%" --> SUB["Use its call"]
    Q -- "no, 17%" --> POCKET["Binding-pocket lookup<br/>34 residues read by profile HMM<br/>matched to 4,562 curated signatures"]
    PARAS --> POCKET
    MIBIG --> POCKET
    POCKET --> SUB

    SUB --> BUILD["Grow the chain in RDKit<br/>amide bonds, polyketide units,<br/>KR/DH/ER set the beta carbon,<br/>Cy domains close thiazolines"]
    BUILD --> REL["Release the chain<br/>hydrolysis, reduction,<br/>macrolactam or macrolactone"]
    REL --> TAIL["Decorate from the cluster's enzymes<br/>sugar, methyl, halogen"]
    TAIL --> RANK["Rank the candidates against<br/>the product the cluster predicts:<br/>size, nitrogens, rings, sugar"]
    RANK --> PRED(["Predicted structure"])

    PRED --> SCORE["Score four ways:<br/>fingerprint similarity,<br/>common substructure loose and strict,<br/>containment"]
    MIBIG --> TRUTH(["Known structure"])
    TRUTH --> SCORE
    ASBASE(["antiSMASH's own SMILES"]) --> SCORE
    AS --> ASBASE

    classDef data fill:#257de0,stroke:#1a5fa8,color:#ffffff
    classDef step fill:#f6f8fa,stroke:#8b949e,color:#191919
    classDef out fill:#00b894,stroke:#00866a,color:#ffffff
    class MIBIG,AS,PARAS data
    class PARSE,MOD,SUB,POCKET,BUILD,REL,TAIL,SCORE step
    class RANK step
    class PRED,TRUTH,ASBASE,FLAG out
```

The numbered walk-through:

1. **Parse** (`antismash_io.py`) — read `aSDomain` features from the antiSMASH
   GenBank, order them by gene and position in the protein, and cut them into
   modules: a new PKS module at every ketosynthase, a new NRPS module at every
   adenylation domain. Substrates come from the `substrate consensus`
   qualifier, loading domains from their own call. A cluster whose product
   does not come off an assembly line at all — an iterative or type III
   synthase — is flagged rather than built, because its module count says
   nothing about the size of its product.
2. **Assemble** (`assembly.py`, `chem.py`) — grow the chain as a real molecule:
   - NRPS modules form an amide bond with their amino acid;
   - heterocyclisation (Cy) domains close the residue into a thiazoline
     (cysteine) or oxazoline (serine, threonine) ring;
   - PKS modules add a two-carbon unit, and the module's KR/DH/ER domains set
     the oxidation level of the beta carbon (keto → hydroxyl → ene → alkyl);
   - epimerase domains invert the alpha carbon, methyltransferases add a methyl;
   - the chain is released by hydrolysis, reduction, macrolactamisation or
     macrolactonisation.
   Tailoring enzymes found in the cluster then decorate the backbone
   (`tailoring.py`): a sugar, a halogen or a methyl group, one each at most.
3. **Fill the gaps** (`hmm_align.py`, `predictor.py`) — antiSMASH reports no
   substrate for about one adenylation domain in six. Those domains are given
   one by reading their 34 binding-pocket residues — against a profile HMM
   where HMMER is available, otherwise by pairwise alignment — and matching
   them to 4,562 curated signatures from MIBiG and PARAS. antiSMASH's own
   calls are never overridden.
4. **Choose** (`ranking.py`) — several structures can be built from one
   cluster. Rather than picking by rule, the cluster is asked what its product
   should look like, and the candidate closest to that description is kept.
5. **Compare** (`compare.py`) — score the prediction against the known
   structure four ways, because no single number says whether a backbone is
   right:

   | Measure | What it answers |
   |---|---|
   | `tanimoto` | Morgan-fingerprint similarity to the whole product |
   | `coverage_pred` | share of the prediction found in the product, bond orders ignored |
   | `strict_coverage_pred` | the same with bond orders required to agree |
   | `substructure` | is the prediction literally contained in the product |

   The strict measure exists because the loose one is blind to exactly the
   differences an assembly line decides. A saturated chain scored against an
   unsaturated one gives a loose coverage of 1.00 and a strict coverage of
   0.75; Tanimoto meanwhile falls from 1.0 to 0.375, which overstates a
   one-bond difference. Reading all four together is more honest than
   trusting any one of them.

The result is a **backbone**. Post-assembly tailoring — glycosylation,
halogenation, oxidative crosslinking, prenylation — is not modelled, so a
prediction is deliberately smaller than the finished natural product. The
`coverage_pred` metric is the one that should be high for a backbone
predictor: it asks whether everything predicted is really present in the
product.

## Background

SeMPI-3 follows [SeMPI](https://academic.oup.com/nar/article/44/W1/W510/2499349)
and [SeMPI 2.0](https://doi.org/10.3390/metabo11010013), which predicted
PKS/NRPS products and screened them against natural-product databases.

It is built on data from other people's work:

* **MIBiG 4.0** — curated gene clusters and the structures they make.
* **antiSMASH 8** — the domain annotation every prediction starts from, and
  the baseline every result is measured against. MIBiG publishes its output
  for every entry, so antiSMASH itself is not a dependency.
* **PARAS/PARASECT** ([Terlouw et al.](https://www.biorxiv.org/content/10.1101/2025.01.08.631717v2.full),
  CC BY 4.0) — 3,653 curated adenylation-domain substrates, their binding-pocket
  signatures, structures for several hundred unusual monomers, and the Pfam
  AMP-binding profile used to read a pocket. Their data is used; their code,
  which is AGPL, is not.
