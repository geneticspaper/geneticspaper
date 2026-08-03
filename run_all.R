#!/usr/bin/env Rscript
# run_all.R -- end-to-end pipeline for the Ashkenazi mtDNA origin model.
#
# Prerequisite: data/derived/*.csv must exist. Regenerate from the Mitotree
# workbook with:  python3 scripts/00_extract_data.py [path/to/media-1(1).xlsx]
#
# Usage:  Rscript run_all.R
#
# Optional frequency-source overrides (defaults: Brook + 1KG EUR subclades):
#   MTDNA_ASHKENAZI_FREQ=brook|mitotree
#   MTDNA_NONJEW_FREQ=1kg_eur|livni_skorecki|mitotree
#   MTDNA_FREQ_COMPARE=1   # also write D_frequency_source_comparison.csv

t0 <- Sys.time()
source(file.path("R", "utils.R"))
init_frequency_sources_from_env()

if (!file.exists(file.path(DERIVED, "s1_public_samples.csv"))) {
  stop("Missing data/derived/s1_public_samples.csv -- run scripts/00_extract_data.py first.")
}

kg_freq <- file.path(REFERENCE, "1kg_eur_clade_freq.csv")
if (!file.exists(kg_freq)) {
  message("Regenerating 1000 Genomes EUR (CEU+GBR+FIN+IBS+TSI) clade frequencies...")
  status <- system2("python3", file.path("scripts", "compute_1kg_eu_freq.py"))
  if (status != 0L || !file.exists(kg_freq)) {
    stop("Failed to build ", kg_freq, " -- run python3 scripts/compute_1kg_eu_freq.py")
  }
}

source(file.path("R", "00_enrich_dedup.R"))
source(file.path("R", "01_branching_model.R"))
source(file.path("R", "02_phylogeography.R"))
source(file.path("R", "03_ancient_and_diffusion.R"))
source(file.path("R", "06_advanced_models.R"))
source(file.path("R", "07_origin_stress_tests.R"))
source(file.path("R", "04_report.R"))
source(file.path("R", "05_paper.R"))

make_enriched_deduped()
run_branching_model()
run_phylogeography()
run_ancient_and_diffusion()
run_advanced_models()
run_origin_stress_tests()
if (nzchar(Sys.getenv("MTDNA_FREQ_COMPARE", unset = ""))) {
  frequency_source_comparison()
}
build_report()
build_paper()

message(sprintf(paste0(
  "\nAll modules complete in %.1f s.\n",
  "  Markdown report: report/ashkenazi_mtdna_origin.md\n",
  "  Scientific PDF:  report/ashkenazi_mtdna_paper.pdf"),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
