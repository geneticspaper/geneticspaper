# Ashkenazi mtDNA Levantine-Origin Model (Livni-Skorecki + Mitotree)

A reproducible R pipeline that tests whether the four major Ashkenazi maternal
founder lineages -- **K1a1b1a, K1a9, K2a2a, N1b2** -- are of Middle Eastern /
Levantine origin. It re-implements the branching-process / founder-vs-host model
of **Livni & Skorecki (2025)** and integrates a phylogeographic and ancient-DNA
analysis built on an enriched **Mitotree + GenBank** reference set (Maier et al.
2026), directly addressing the opposing "prehistoric European" argument of
**Costa et al. (2013)**. Behar/Costa rows are provenance flags within the
deduplicated dataset, not separate sample tables.

## What it shows

1. **Founder vs host (Module A).** A Galton-Watson lineage-extinction model plus a
   singleton/doubleton detection curve is retained as the Livni-Skorecki model
   framework, but absorption rates are **not** estimated from Behar/Costa sample
   tables. Empirical founder counts come from the enriched deduped dataset.
   A "Euro-Levantine" pool has < 0.4% model probability of producing >=3 of the
   four majors under Livni-Skorecki assumptions; this is a model sensitivity
   result, not a Mitotree measurement.
2. **Phylogeography (Module B, `ape` + `pegas`).** Rarefying to equal sample sizes
   weakens the Costa et al. "solely European nesting" signal for K1a. K2a and
   N1b remain mixed, so this is evidence against a blanket nesting claim rather
   than a stand-alone proof. The four founders are compared against the published
   non-Jewish rarity benchmark.
3. **Ancient DNA & diffusion (Module C, `ape` + `adegenet`).** The oldest ancient
   K1a records in the enriched dataset are Near Eastern (Anatolia ~8750 BCE,
   Levant); historical founder-clade hits are explicitly tabulated; and
   parent-signature frequency gradients are more compatible with a Near Eastern
   source for K1a/K2a than with a purely European source.
4. **Advanced models & Bayesian synthesis (Module D).** Sweeps negative-binomial
   offspring dispersion, tests piecewise historical expansion timing, and builds
   a multi-channel Bayesian origin synthesis over rarity, antiquity, branch age,
   and equal-$n$ nesting channels.
5. **Origin stress tests & host bounds (Module G).** Falsifies European host
   predictions against the AADR v66 ancient mtDNA record, computes 1000 Genomes
   European host frequency bounds and likelihood ratios, and evaluates prior sensitivity.

The full write-up with all figures and tables is in
[`report/ashkenazi_mtdna_origin.md`](report/ashkenazi_mtdna_origin.md) and as a
scientific PDF manuscript at [`report/ashkenazi_mtdna_paper.pdf`](report/ashkenazi_mtdna_paper.pdf).

## Layout

```
mitomdoel/
  scripts/00_extract_data.py           # stream-extract the Mitotree .xlsx -> CSV
  scripts/01_fetch_genbank_founders.py # fetch GenBank founder metadata
  scripts/02_extract_genes2026.py       # extract Genes 2026 supplementary data
  scripts/03_fetch_public_haplogroup_pages.py # fetch YFull / FTDNA public metadata
  scripts/04_extract_aadr.py           # extract AADR ancient mtDNA records
  scripts/compute_1kg_eu_freq.py       # compute 1000 Genomes EUR clade frequencies
  data/reference/*.csv                 # hand-encoded tables and reference data
  data/derived/*.csv                   # extracted Mitotree sheets (gitignored)
  R/utils.R                            # data loaders, region classifier, GW PMF solvers
  R/lineage_data.R                     # lineage definitions and frequency loaders
  R/00_enrich_dedup.R                  # Mitotree+GenBank enrichment and deduplication
  R/01_branching_model.R               # Module A (branching model)
  R/02_phylogeography.R                # Module B (ape, pegas)
  R/03_ancient_and_diffusion.R         # Module C (ape, adegenet)
  R/06_advanced_models.R               # Module D (dispersion, growth timing, Bayesian synthesis)
  R/07_origin_stress_tests.R           # Module G (AADR ancient falsification & 1KG EUR host bounds)
  R/04_report.R                        # assembles the Markdown report
  R/05_paper.R                         # compiles the scientific LaTeX/PDF paper manuscript
  run_all.R                            # runs everything
  report/ashkenazi_mtdna_origin.{md,qmd} # Markdown report and Quarto source
  report/ashkenazi_mtdna_paper.{pdf,tex} # Scientific paper PDF and LaTeX source
  outputs/{figures,tables}/            # generated artifacts
```

## Requirements

- **R >= 4.4** (developed on 4.6.1).
- R packages: `ape`, `pegas`, `adegenet` (everything else is base R).
- **Python 3** (standard library only) for one-off data extractions.
- **pdflatex / LaTeX system** (or `tinytex`) for compiling `report/ashkenazi_mtdna_paper.pdf`.

Install the R packages into a project-local library (keeps the system R clean):

```r
dir.create(".Rlib", showWarnings = FALSE)
install.packages(c("ape", "pegas", "adegenet"),
                 lib = ".Rlib",
                 repos = "https://cloud.r-project.org")
```

`R/utils.R` automatically adds `.Rlib` to the library path if present.

## Run

```bash
# 1. Extract the Mitotree supplementary workbook to CSV (path optional; the
#    default points at the media-1(1).xlsx used during development).
python3 scripts/00_extract_data.py "/path/to/media-1(1).xlsx"

# 2. Fetch additional GenBank founder metadata (optional but recommended)
python3 scripts/01_fetch_genbank_founders.py your.email@example.org

# 3. Compute 1000 Genomes EUR clade frequencies (auto-run by run_all.R if missing)
python3 scripts/compute_1kg_eu_freq.py

# 4. Run the whole pipeline (dedup -> branching -> phylogeography -> ancient DNA -> advanced models -> stress tests -> report & paper)
Rscript run_all.R

# 5. (optional) render the Quarto report if quarto + pandoc are available
quarto render report/ashkenazi_mtdna_origin.qmd
```

Frequency sources can be overridden via environment variables before running `run_all.R`:
- `MTDNA_ASHKENAZI_FREQ=brook|mitotree`
- `MTDNA_NONJEW_FREQ=1kg_eur|livni_skorecki|mitotree`
- `MTDNA_FREQ_COMPARE=1` (writes `D_frequency_source_comparison.csv`)

Outputs land in `outputs/figures/` (PNG), `outputs/tables/` (CSV),
`report/ashkenazi_mtdna_origin.md`, and `report/ashkenazi_mtdna_paper.pdf`.

## Data provenance

- `S1 Public Samples` and `S10 Mitotree Structure` are extracted from the Mitotree
  supplementary workbook (`media-1(1).xlsx`; Maier et al. 2026). The public
  release here is 60,705 samples (57,551 modern, 3,150 ancient), a subset of the
  full ~330,000-sequence tree.
- `data/external/genbank_founder_records.csv` is fetched from GenBank E-utilities
  for K1a1b1a, K1a9, K2a2a and N1b2. Most accessions are already represented in
  Mitotree; deduplication prevents double-counting.
- `data/derived/enriched_deduped_samples.csv` is the empirical analysis dataset:
  Mitotree + GenBank query metadata + Behar/Costa provenance flags, deduped by
  accession/sample.
- `data/reference/` holds tables transcribed from Livni & Skorecki (2025):
  Table 3 (regional parent-clade frequencies), Table 7 (founder frequencies in
  non-Jews) and Table 2 (founder-event parameters). These are benchmarks/model
  inputs, not the empirical sample set.

## Limitations (please read)

Mitotree/GenBank is a cross-study reference set, not a random Ashkenazi
population sample. Modern `Country` records where a carrier lives, not where a
lineage originated, so modern distributions are not used as direct evidence of
origin. Singleton/doubleton absorption rates cannot be estimated from this
reference set. The public sample subset and reconstructed (not
per-sample-aligned) motifs limit resolution. This pipeline was built to *test*
the Levantine hypothesis against the European one; it finds weighted evidence
leaning toward a Levantine/Near Eastern origin, especially for K1a/K1a1b1a and
K2a/K2a2a, but it is not a substitute for peer-reviewed primary analysis.

## References

- Livni J. & Skorecki K. (2025). *Distinguishing between founder and host population mtDNA lineages in the Ashkenazi population.* Human Gene (SSRN 5035272).
- Costa M.D. et al. (2013). *A substantial prehistoric European ancestry amongst Ashkenazi maternal lineages.* Nature Communications 4:2543.
- Behar D.M. et al. (2006). *The matrilineal ancestry of Ashkenazi Jewry.* American Journal of Human Genetics 78:487–497.
- Maier P.A. et al. (2026). *Mitotree: the universal human mitochondrial reference phylogeny.* bioRxiv.
- Waldman S. et al. (2022). *Genome-wide data from medieval German Jews show that the Ashkenazi founder event pre-dated the 14th century.* Cell 185:4703–4716.
- Brace S. et al. (2022). *Genomes from a medieval mass burial show Ashkenazi-associated hereditary diseases pre-date the 12th century.* Current Biology 32:4350–4359.
- Pallares-Vina L. et al. (2026). *Uncovering a medieval pogrom: genetic history of a Jewish community in Catalonia (Spain).* Genes 17(3):358.
- Fernandez E. et al. (2014). *Ancient DNA analysis of 8000 B.C. Near Eastern farmers.* PLoS Genetics 10(6):e1004401.
- Grossman S. et al. (2019). *Rare human mitochondrial HV lineages spread from the Near East and Caucasus during post-LGM and Neolithic expansions.* Scientific Reports 9:12280.
- Brook K.A. (2022). *The Maternal Genetic Lineages of Ashkenazic Jews.* Academic Studies Press.
- Atzmon G. et al. (2010). *Abraham's children in the genome era: major Jewish Diaspora populations comprise distinct genetic clusters with shared Middle Eastern ancestry.* American Journal of Human Genetics 86:850–859.
- Carmi S. et al. (2014). *Sequencing an Ashkenazi reference panel supports population-targeted personal genomics and illuminates Jewish and European origins.* Nature Communications 5:4835.
- Agranat-Tamir L. et al. (2020). *The genomic history of the Bronze Age Southern Levant.* Cell 181:1146–1157.
- 1000 Genomes Project Consortium (2015). *A global reference for human genetic variation.* Nature 526:68–74.
