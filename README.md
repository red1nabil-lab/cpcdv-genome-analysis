# CpCDV genome analysis

Analysis code for a whole-genome reassessment of strain designations in **chickpea
chlorotic dwarf virus** (CpCDV; species *Mastrevirus cicerparvi*, family
*Geminiviridae*), based on all 347 CpCDV sequences in GenBank.

Everything needed to reproduce the published numbers is here: the scripts, the cached
GenBank records, the alignments and the identity matrices. The pipeline runs offline.

---

## What the analysis does

| Step | Question |
|---|---|
| 01–03 | Assemble the 347 public records and curate country, host, date and strain from GenBank source qualifiers |
| 04 | Rotate every circular genome to the conserved nonanucleotide `TAATATTAC` so that genomes deposited at arbitrary origins can be compared |
| 05–06 | Align, then compute pairwise identity under both gap conventions |
| 07–08 | Verify divergent genomes against other mastrevirus species |
| 10–11 | Nucleotide diversity by region, country, host and strain; sliding-window diversity |
| 12, 17 | Model selection and maximum-likelihood phylogeny (GTR+G+I, 1000 ultrafast bootstraps) |
| 13b | Place partial sequences by Smith–Waterman local alignment |
| 25 | Realign **all 35,511 genome pairs** independently, to follow the genus demarcation procedure exactly |
| 26–27 | Test the published strain labels against the 94% criterion, with threshold-free checks |
| 28 | Rooting-independent monophyly test with node support |
| 29, 37, 39 | Recombination by maximum chi-square; breakpoint independence; event-window sensitivity |
| 30 | Sensitivity of diversity statistics to gap filtering |
| 31, 35 | Threshold violations after removing recombinants, against matched random-removal nulls |
| 36, 38 | Independent recombination evidence: PHI and the Hudson–Kaplan bound |
| 14, 33, 34 | Figures |
| 15 | Supplementary Table S1 (the dataset inventory) |

Headline results: the published strain designations are broadly concordant with
genome-wide identity (adjusted Rand index 0.990); all 28 within-strain threshold
violations fall inside **strain C**; strains C, L, N and O form a single clade and
strain P nests within strain H; recombination is pervasive (PHI p = 0.005,
Hudson–Kaplan ≥ 197 events, 42 localised events) with breakpoints concentrated in the
small intergenic region and Rep exon 1.

## Running it

```bash
pip install -r requirements.txt
export CPCDV_DATA=$PWD/data          # optional; this is the default
python scripts/01_parse.py
python scripts/03_curate.py
...
```

Scripts are numbered in dependency order and each writes its outputs to `CPCDV_DATA`.
Because the intermediate files are committed, later steps can be run without the earlier
ones — for example, the recombination analysis runs directly from the committed
alignment:

```bash
python scripts/29_recombination.py
python scripts/38_incompat_fixed_cols.py
```

**Configuration.** `scripts/cpcdv_paths.py` resolves three locations, all overridable by
environment variable:

| Variable | Default | Used for |
|---|---|---|
| `CPCDV_DATA` | `./data` | inputs and all analysis outputs |
| `CPCDV_FIGURES` | `./figures` | generated figures |
| `MUSCLE_BIN` | `muscle` | MUSCLE v3.8.31 binary (steps 05 and 08 only) |

> **One modification from the code as run.** The scripts originally hardcoded a single
> local data directory. That constant — and only that constant — now resolves through
> `cpcdv_paths.py`. No analysis logic was altered. Steps 02, 07 and 08 contact NCBI;
> every other step runs offline from the committed data.

### Requirements

- Python 3.13 — Biopython, NumPy, SciPy, matplotlib, pandas, openpyxl, python-docx
- R 4.5 — ape, phangorn, ggtree, ggtreeExtra, ggnewscale, patchwork, jsonlite, dplyr
- MUSCLE v3.8.31 (steps 05 and 08 only; later steps use the committed alignment)

## Data

`data/` holds the inputs and the main intermediates (~7 MB):

| File | Contents |
|---|---|
| `CpCDV_all_seq.fasta` | the 347 GenBank sequences analysed |
| `02_genbank_cache.gb` | cached GenBank flat file, so the pipeline reproduces without network access |
| `03_curated_metadata.csv` | curated country, region, host, collection year and strain per accession |
| `04_full_oriented.fasta` | 270 complete genomes rotated to a common origin |
| `05_full_aln.fasta` | MUSCLE alignment (2,829 columns) |
| `06_pid_pw.npy`, `06_pid_gap.npy` | identity matrices from the multiple alignment, both gap conventions |
| `25_truepw_pid_*_A.npy` | identity matrices from independent pairwise alignment of all 35,511 pairs |
| `12_ml_tree.newick` | maximum-likelihood tree with ultrafast bootstrap support |
| `29_recombination_maxchi.csv` | per-genome recombination results, breakpoints and p-values |

Matrices are NumPy arrays; row/column order is `06_accession_order.json` for the
`06_*` matrices and `25_truepw_meta.json` for the `25_*` matrices. **Note that
`06_pid_*.npy` store identities as fractions (0–1) while `25_truepw_*.npy` store
percentages** — the analysis scripts normalise on load.

All sequence data derive from public GenBank records; accession numbers are listed in
`03_curated_metadata.csv`.

## What is not here

The scripts that generate the manuscript text, its bibliography and the results
document are deliberately excluded — this repository is the analysis, not the write-up.

## Citing

A manuscript reporting these analyses is in preparation. Until it appears, please cite
this repository. Issues and corrections are welcome: the pipeline found three errors in
its own earlier versions, and a fourth found by a reviewer, so more are likely.

## Licence

[MIT](LICENSE) for the code. Sequence data are from GenBank and are not subject to this
licence.
