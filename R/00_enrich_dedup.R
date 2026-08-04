# 00_enrich_dedup.R
# Build the analysis dataset from Mitotree S1 plus additional GenBank
# E-utilities metadata and explicit Behar/Costa provenance flags.
#
# Important: Behar and Costa rows are already present inside Mitotree S1. This
# step does not append paper tables as extra population samples. It enriches
# provenance, appends GenBank-only records not already represented by accession,
# and removes duplicate accessions/samples before empirical counts are computed.

if (!exists("ROOT")) source(file.path("R", "utils.R"))

norm_acc <- function(x) {
  x <- trimws(ifelse(is.na(x), "", x))
  x <- sub("\\.[0-9]+$", "", x)
  x[x %in% c("", ".")] <- NA
  x
}

make_enriched_deduped <- function() {
  message("== Module 0: enriched deduped Mitotree+GenBank+Behar+Costa dataset ==")
  raw <- load_samples()

  # Mitotree can credit a later compilation instead of the paper that first
  # published a sample. I14741 is one of the Erfurt individuals introduced by
  # Waldman et al. (2022); Akbari et al. (2026) only reuses it. Keep the
  # primary-study attribution in analysis tables and reports.
  i14741 <- raw$Subject == "I14741"
  if (any(i14741)) {
    raw$Study[i14741] <- "Waldman 2022"
    raw$Author[i14741] <- "Waldman"
    raw$Year[i14741] <- "2022"
    raw$PubMed[i14741] <- "36455558"
    raw$URL[i14741] <- "https://doi.org/10.1016/j.cell.2022.11.002"
    raw$Title[i14741] <- paste(
      "Genome-wide data from medieval German Jews show that the Ashkenazi",
      "founder event pre-dated the 14th century"
    )
    raw$Journal[i14741] <- "Cell 185:4703-4716"
  }
  raw$source_dataset <- "Mitotree_S1"
  raw$AccessionBase <- norm_acc(raw$Accession)
  raw$SampleBase <- norm_acc(raw$Sample)
  text <- paste(raw$Study, raw$Title, raw$Journal, raw$StudyLabel, raw$Note, sep = " ")
  raw$in_behar <- grepl("Behar", raw$Study, ignore.case = TRUE)
  raw$in_costa <- grepl("Costa", raw$Study, ignore.case = TRUE)
  raw$in_genes2026 <- FALSE
  raw$in_genbank <- grepl("^[A-Z]{1,4}[0-9]+(\\.[0-9]+)?$", raw$Accession)
  raw$in_jewish_study <- grepl(
    "(?i)ashkenazi|jewish|jewry|\\bjews\\b|\\bjew\\b|sephard|mizrahi|hebrew|sobib",
    text,
    perl = TRUE
  )
  raw$provenance <- paste0(
    "mitotree",
    ifelse(raw$in_genbank, ";genbank_accession", ""),
    ifelse(raw$in_behar, ";behar", ""),
    ifelse(raw$in_costa, ";costa", ""),
    ifelse(raw$in_jewish_study, ";jewish_study", "")
  )
  raw$haplogroup_source <- "Mitotree"

  gb_path <- file.path(ROOT, "data", "external", "genbank_founder_records.csv")
  if (file.exists(gb_path)) {
    gb <- read.csv(gb_path, stringsAsFactors = FALSE, colClasses = "character")
    gb$Accession <- gb$accession_version
    gb$AccessionBase <- norm_acc(gb$accession_version)
    existing <- unique(na.omit(raw$AccessionBase))
    gb_new <- gb[is.na(gb$AccessionBase) | !(gb$AccessionBase %in% existing), ]
    if (nrow(gb_new)) {
      template <- raw[0, ]
      add <- template[rep(1, nrow(gb_new)), , drop = FALSE]
      add[,] <- NA
      add$Subject <- gb_new$uid
      add$Sample <- gb_new$accession_version
      add$SubjectType <- "Modern"
      add$Haplotype <- NA
      add$Haplogroup <- gb_new$query_haplogroup
      # GenBank query assignment is lower confidence than Mitotree's tree call;
      # it is kept separate in haplogroup_source and provenance.
      add$MitotreeHaplogroup <- gb_new$query_haplogroup
      add$Length <- gb_new$length
      add$Study <- "GenBank E-utilities"
      add$Author <- NA
      add$Year <- NA
      add$PubMed <- NA
      add$URL <- paste0("https://www.ncbi.nlm.nih.gov/nuccore/", gb_new$accession_version)
      add$Title <- gb_new$title
      add$Journal <- NA
      add$Accession <- gb_new$accession_version
      add$StudyLabel <- gb_new$uid
      add$Country <- NA
      add$Note <- gb_new$query_term
      add$region <- region_of(add$Country)
      add$source_dataset <- "GenBank_EUtilities"
      add$AccessionBase <- gb_new$AccessionBase
      add$SampleBase <- norm_acc(add$Sample)
      add$in_behar <- FALSE
      add$in_costa <- FALSE
      add$in_genes2026 <- FALSE
      add$in_genbank <- TRUE
      add$in_jewish_study <- FALSE
      add$provenance <- "genbank_eutilities"
      add$haplogroup_source <- "GenBank_query"
      raw <- rbind(raw, add)
    }
  }

  genes_path <- file.path(ROOT, "data", "external", "genes2026_mtdna.csv")
  if (file.exists(genes_path)) {
    genes <- read.csv(genes_path, stringsAsFactors = FALSE, colClasses = "character")
    genes <- genes[nzchar(genes$Sample) & nzchar(genes$Haplogroup), ]
    if (nrow(genes)) {
      template <- raw[0, ]
      add <- template[rep(1, nrow(genes)), , drop = FALSE]
      add[,] <- NA
      add$Subject <- paste0("Genes2026_", genes$Sample)
      add$Sample <- genes$Sample
      add$SubjectType <- "Ancient"
      add$Haplotype <- genes$Haplotype
      add$Haplogroup <- genes$Haplogroup
      add$MitotreeHaplogroup <- trimws(genes$Haplogroup)
      add$Length <- NA
      add$ConfidentBP <- NA
      add$BirthYear <- NA
      add$AgeEstimateMean <- 1348
      add$Study <- "Pallares-Vina 2026"
      add$Author <- "Pallares-Vina"
      add$Year <- "2026"
      add$PubMed <- NA
      add$URL <- "https://doi.org/10.3390/genes17030358"
      add$Title <- "Uncovering a Medieval Pogrom: Genetic History of a Jewish Community in Catalonia (Spain)"
      add$Journal <- "Genes 17(3):358"
      add$Accession <- NA
      add$StudyLabel <- genes$Sample
      add$Country <- "Spain"
      add$Note <- paste0("Roquetes/Tarrega medieval Jewish pogrom; mtDNA quality=", genes$Quality)
      add$region <- region_of(add$Country)
      add$source_dataset <- "Genes_2026_Roquetes"
      add$AccessionBase <- NA
      add$SampleBase <- paste0("genes2026_", genes$Sample)
      add$in_behar <- FALSE
      add$in_costa <- FALSE
      add$in_genes2026 <- TRUE
      add$in_genbank <- FALSE
      add$in_jewish_study <- TRUE
      add$provenance <- "genes2026;medieval_jewish;roquetes_tarrega"
      add$haplogroup_source <- "Genes2026_Haplogrep3"
      raw <- rbind(raw, add)
    }
  }

  # Deduplication priority: exact accession > sample > subject. Prefer Mitotree
  # tree-called rows over GenBank-query-only rows, then fuller metadata.
  raw$dedup_key <- ifelse(!is.na(raw$AccessionBase), paste0("acc:", raw$AccessionBase),
                   ifelse(!is.na(raw$SampleBase), paste0("sample:", raw$SampleBase),
                          paste0("subject:", raw$Subject)))
  raw$priority <- ifelse(raw$source_dataset == "Mitotree_S1", 2L, 1L) +
    ifelse(!is.na(raw$Country), 1L, 0L) +
    ifelse(raw$in_jewish_study, 1L, 0L)
  raw <- raw[order(raw$dedup_key, -raw$priority), ]

  # Merge provenance flags across duplicate keys, then keep the best row.
  flag_sum <- aggregate(cbind(in_behar, in_costa, in_genes2026, in_genbank, in_jewish_study) ~ dedup_key,
                        raw, function(x) any(x, na.rm = TRUE))
  prov_sum <- aggregate(provenance ~ dedup_key, raw, function(x) {
    paste(unique(unlist(strsplit(paste(x, collapse = ";"), ";"))), collapse = ";")
  })
  keep <- !duplicated(raw$dedup_key)
  out <- merge(raw[keep, ], flag_sum, by = "dedup_key", suffixes = c("", "_any"), all.x = TRUE)
  out <- merge(out, prov_sum, by = "dedup_key", suffixes = c("", "_merged"), all.x = TRUE)
  for (nm in c("in_behar", "in_costa", "in_genes2026", "in_genbank", "in_jewish_study")) {
    out[[nm]] <- out[[paste0(nm, "_any")]]
    out[[paste0(nm, "_any")]] <- NULL
  }
  out$provenance <- out$provenance_merged
  out$provenance_merged <- NULL
  out$priority <- NULL

  out_path <- file.path(DERIVED, "enriched_deduped_samples.csv")
  write.csv(out, out_path, row.names = FALSE)

  summary <- data.frame(
    metric = c("mitotree_rows", "rows_after_genbank_append", "deduped_rows",
               "unique_accessions", "behar_flagged", "costa_flagged", "genes2026_flagged",
               "genbank_flagged", "jewish_study_flagged"),
    value = c(nrow(load_samples()), nrow(raw), nrow(out),
              length(unique(na.omit(out$AccessionBase))),
              sum(out$in_behar), sum(out$in_costa), sum(out$in_genes2026), sum(out$in_genbank),
              sum(out$in_jewish_study))
  )
  save_table(summary, "00_enriched_dedup_summary.csv")
  message(sprintf("Wrote %d deduped rows -> %s", nrow(out), out_path))
  print(summary)
  invisible(out)
}

if (sys.nframe() == 0) make_enriched_deduped()
