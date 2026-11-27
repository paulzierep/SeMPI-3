# SeMPI-3

Predict the structure of a secondary metabolite from its biosynthetic gene
cluster. SeMPI-3 reads an **antiSMASH** annotation, cuts the assembly line into
modules, and builds the product with **RDKit**.

Scope: modular type-I **PKS**, **NRPS** and their hybrids.

## Install

```bash
git clone https://github.com/paulzierep/SeMPI-3.git
cd SeMPI-3
conda env create -f environment.yml
conda activate sempi3
pip install -e .
```

Nothing else is needed to predict — the binding-pocket table ships in the
package. antiSMASH is **not** a dependency; SeMPI-3 reads its *output*. HMMER
is optional but recommended (`conda install -c bioconda hmmer`), which lets the
binding pocket be read against a profile instead of by pairwise alignment.

## Run

Run antiSMASH as you normally would, then point SeMPI-3 at the GenBank file it
wrote. A genome-wide run holds many regions; each is predicted separately.

```bash
sempi3 run region001.gbk                     # report: structure, modules, antiSMASH's own call
sempi3 run region001.gbk --format tabular    # one row per region, for a workflow
sempi3 run region001.gbk --format smiles     # just the SMILES
sempi3 run region001.gbk --format json --out prediction.json
sempi3 run region001.gbk --no-tailoring      # bare backbone, undecorated
```

```text
BGC0000290_region1  products=NRPS|T3PKS  modules=7
   0 nrps Dhpg
   1 nrps Tyr                        E        <- epimerase: a D residue
   2 nrps Hpg                        E
   ...
tailoring enzymes: {'halogenase': 2, 'hydroxylase': 4}
release: macrolactam
SeMPI-3:   O=C1N[C@@H](c2ccc(O)cc2Cl)C(=O)N[C@@H](c2cc(O)cc(O)c2)...
antiSMASH: NC(c1ccc(O)cc1)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(c1cc(O)cc(O)c1)...
```

Ours closes the macrolactam, adds the chlorine the halogenases imply, and keeps
the stereochemistry the epimerase domains set. antiSMASH leaves the chain
linear and flat. A cluster whose product does not come off an assembly line —
an iterative or type III synthase — is flagged rather than mispredicted.

## Update

A new release of MIBiG or PARAS is absorbed by rebuilding, not retraining:

```bash
./scripts/fetch_data.sh        # ~95 MB download, benchmark and rebuild only
sempi3 build-reference         # rebuild the binding-pocket table
sempi3 train                   # rebuild the table the ranking models fit from
pytest -q && sempi3 benchmark  # check nothing regressed
```

There are no model files. The one fitted component trains from a packaged table
the first time it is needed, so nothing is pinned to a scikit-learn version.

## Results

<!-- RESULTS:START -->
Over the **1628 MIBiG 4.0 PKS/NRPS clusters** with a known structure:

| | Clusters predicted | Tanimoto | Scoring ≥0.4 | Structure match |
|---|---|---|---|---|
| antiSMASH (baseline) | 1218 | 0.146 | 56 | 0.425 |
| **SeMPI-3** | 1406 | 0.204 | 159 | 0.532 |
| SeMPI-3 with curated substrates | 578 | 0.270 | 131 | 0.584 |

Head to head on the 1216 clusters where both predict: **881 wins, 327 losses** (+0.073 mean Tanimoto).

By class (SeMPI-3 vs antiSMASH): NRPS 0.301 vs 0.219, PKS 0.139 vs 0.091, hybrid 0.217 vs 0.136.
<!-- RESULTS:END -->

![Benchmark quality](docs/figures/quality.png)

![Head to head against antiSMASH](docs/figures/head_to_head.png)

Scores are low in absolute terms for every tool, because all of them predict a
*backbone* while the reference is the finished, tailored natural product. The
comparison between them is what counts.

## More

* [How it works](docs/method.md) — the pipeline, step by step
* [Benchmark](docs/benchmark.md) — how it is measured, what was tried, what
  failed, and what is worth doing next

```bash
pytest -q                      # 112 tests, no downloaded data needed
```

Built on [MIBiG 4.0](https://mibig.secondarymetabolites.org/),
[antiSMASH 8](https://antismash.secondarymetabolites.org/) and
[PARAS](https://github.com/BTheDragonMaster/parasect) (data under CC BY 4.0).
Follows [SeMPI](https://academic.oup.com/nar/article/44/W1/W510/2499349) and
[SeMPI 2.0](https://doi.org/10.3390/metabo11010013).

## License

MIT — see [LICENSE](LICENSE).
