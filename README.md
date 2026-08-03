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

The full write-up with all figures and tables is
[`report/ashkenazi_mtdna_origin.md`](report/ashkenazi_mtdna_origin.md).

## Layout

```
mitomdoel/
  scripts/00_extract_data.py     # stream-extract the Mitotree .xlsx -> CSV
  scripts/01_fetch_genbank_founders.py # fetch GenBank founder metadata
  data/reference/*.csv           # hand-encoded tables from the papers
  data/derived/*.csv             # extracted Mitotree sheets (gitignored)
  R/utils.R                      # data loaders, region classifier, Galton-Watson
  R/00_enrich_dedup.R            # Mitotree+GenBank enrichment and deduplication
  R/01_branching_model.R         # Module A
  R/02_phylogeography.R          # Module B (ape, pegas)
  R/03_ancient_and_diffusion.R   # Module C (ape, adegenet)
  R/04_report.R                  # assembles the Markdown report
  run_all.R                      # runs everything
  report/ashkenazi_mtdna_origin.{md,qmd}
  outputs/{figures,tables}/      # generated artifacts
```

## Requirements

- **R >= 4.4** (developed on 4.6.1).
- R packages: `ape`, `pegas`, `adegenet` (everything else is base R).
- **Python 3** (standard library only) for the one-off data extraction.

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

# 3. Run the whole pipeline (dedup -> branching model -> phylogeography -> ancient DNA -> report)
Rscript run_all.R

# 4. (optional) render the Quarto report if quarto + pandoc are available
quarto render report/ashkenazi_mtdna_origin.qmd
```

Outputs land in `outputs/figures/` (PNG), `outputs/tables/` (CSV), and
`report/ashkenazi_mtdna_origin.md`.

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

- Livni J. & Skorecki K. (2025). *Distinguishing between founder and host
  population mtDNA lineages in the Ashkenazi population.* Human Gene (SSRN 5035272).
- Costa M.D. et al. (2013). *A substantial prehistoric European ancestry amongst
  Ashkenazi maternal lineages.* Nature Communications 4:2543.
- Behar D.M. et al. (2006). Major Ashkenazi maternal founder analysis.
- Maier P.A. et al. (2026). *Mitotree: the universal human mitochondrial reference
  phylogeny at 10x the resolution.* bioRxiv.
