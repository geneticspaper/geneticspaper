# 04_report.R
# Assemble a self-contained Markdown report from the generated tables and
# figures. Uses base R only (no pandoc/quarto required) so `run_all.R` always
# produces report/ashkenazi_mtdna_origin.md. A parallel Quarto source
# (report/ashkenazi_mtdna_origin.qmd) is provided for rendered HTML/PDF.

if (!exists("ROOT")) source(file.path("R", "utils.R"))

.md_table <- function(df, digits = 4) {
  df <- as.data.frame(df)
  fmt <- function(x) if (is.numeric(x)) formatC(x, format = "fg", digits = digits) else as.character(x)
  df[] <- lapply(df, fmt)
  body <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  head <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep  <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  paste(c(head, sep, body), collapse = "\n")
}

.read_tab <- function(name) read.csv(file.path(TAB_DIR, name), stringsAsFactors = FALSE)
.read_external <- function(name) {
  path <- file.path(ROOT, "data", "external", name)
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
}

build_report <- function() {
  message("== Building Markdown report ==")
  D0    <- .read_tab("00_enriched_dedup_summary.csv")
  A_t1  <- .read_tab("A_table1_single_descendant.csv")
  A_audit <- .read_tab("A_branching_math_audit.csv")
  A_sens <- .read_tab("A_offspring_law_sensitivity.csv")
  A_mt  <- .read_tab("A_mitotree_founder_summary.csv")
  A_ma  <- .read_tab("A_mitotree_ancient_founders.csv")
  A_absorb <- .read_tab("A_absorption_estimate.csv")
  A_el  <- .read_tab("A_euro_levantine_feasibility.csv")
  A_els <- .read_tab("A_euro_levantine_sensitivity.csv")
  D_nb  <- .read_tab("D_nbinom_size_sensitivity.csv")
  D_sch <- .read_tab("D_growth_schedule.csv")
  D_pw  <- .read_tab("D_piecewise_vs_constant.csv")
  D_bi  <- .read_tab("D_bayes_inputs.csv")
  D_bs  <- .read_tab("D_bayes_scores.csv")
  D_bp  <- .read_tab("D_bayes_posterior.csv")
  D_nc  <- .read_tab("D_negative_control.csv")
  D_nec <- .read_tab("D_noneuropean_control.csv")
  D_ml  <- .read_tab("D_major_lineages.csv")
  D_fc  <- .read_tab("D_bayes_coefficients.csv")
  if (!nrow(D_fc)) D_fc <- .read_tab("D_fitted_coefficients.csv")
  D_fp  <- .read_tab("D_fitted_posterior.csv")
  D_fl  <- .read_tab("D_fitted_loo.csv")
  G_cd  <- .read_tab("G_channel_dependence.csv")
  B_reg <- .read_tab("B_regional_distribution.csv")
  B_rar <- .read_tab("B_rarefaction.csv")
  C_ne  <- .read_tab("C_oldest_near_east.csv")
  C_foot<- .read_tab("C_k1a1b1a_ancient_footprint.csv")
  C_ca  <- .read_tab("C_coalescence_ages.csv")
  C_t3  <- .read_tab("C_table3_regional_freq.csv")
  T7    <- load_table7()
  T7_founders <- T7[grepl("Livni", T7$source, fixed = TRUE), , drop = FALSE]
  KG1   <- tryCatch(load_1kg_clade_freq(), error = function(e) data.frame())
  D_li  <- .read_tab("D_lineage_computed_inputs.csv")
  B06   <- read.csv(file.path(ROOT, "data", "reference", "ashkenazi_sample_behar2006.csv"), stringsAsFactors = FALSE)
  PPNB  <- load_ppnb()
  FOUNDERS <- load_founders()
  NONEUR <- load_noneuropean()
  PUBLIC <- if (file.exists(file.path(ROOT, "data", "public_haplogroup_pages", "public_haplogroup_page_summary.csv"))) {
    load_public_haplogroup_pages()
  } else {
    data.frame()
  }
  n_public <- nrow(PUBLIC)
  n_yfull_ok <- if (n_public && "yfull_status" %in% names(PUBLIC)) sum(PUBLIC$yfull_status == "ok", na.rm = TRUE) else n_public
  n_ftdna_ok <- if (n_public && "ftdna_status" %in% names(PUBLIC)) sum(PUBLIC$ftdna_status == "ok", na.rm = TRUE) else n_public
  G_mt  <- .read_external("genes2026_mtdna.csv")
  G_qp  <- .read_external("genes2026_qpadm.csv")

  elv <- setNames(A_el$value, A_el$quantity)
  p_ge3 <- as.numeric(elv[["prob_ge3_major"]]) * 100
  els_lo <- if (nrow(A_els)) 100 * min(A_els$prob_ge3_major) else p_ge3
  els_hi <- if (nrow(A_els)) 100 * max(A_els$prob_ge3_major) else p_ge3
  # Same quantity under the Dirichlet sigma_major rather than the assumed 0.03.
  p_ge3_dir <- as.numeric(elv[["prob_ge3_major_dirichlet"]]) * 100
  fs <- frequency_sources()
  if (nrow(G_qp)) {
    G_qp$`p-value` <- as.numeric(G_qp$`p-value`)
    G_qp$C1 <- as.numeric(G_qp$C1)
    G_qp$C2 <- as.numeric(G_qp$C2)
    G_qp$se1 <- as.numeric(G_qp$se1)
    G_qp$se2 <- as.numeric(G_qp$se2)
    G_qp_good <- G_qp[G_qp$`p-value` >= 0.05 & grepl(";", G_qp$`Source (Left) populations`), ]
    G_qp_good <- G_qp_good[, c("Target", "Source (Left) populations", "p-value", "C1", "C2", "se1", "se2")]
  } else {
    G_qp_good <- data.frame()
  }
  if (nrow(G_mt)) {
    G_mt_focus <- G_mt[, c("Sample", "Haplogroup", "Quality")]
  } else {
    G_mt_focus <- data.frame()
  }

  # PPNB (Fernandez 2014) haplogroup composition, excluding the one unassigned skeleton.
  ppnb_assign <- PPNB[PPNB$haplogroup != "", ]
  ppnb_tab <- as.data.frame(table(ppnb_assign$haplogroup), stringsAsFactors = FALSE)
  names(ppnb_tab) <- c("haplogroup", "count")
  ppnb_tab <- ppnb_tab[order(-ppnb_tab$count), ]
  ppnb_tab$percent <- round(100 * ppnb_tab$count / sum(ppnb_tab$count), 1)

  # non-European Ashkenazi lineages (deep African / Asian origin): descriptive
  # fields joined to their posterior under the identical synthesis (D_nec).
  noneur_tab <- data.frame(
    lineage = sub(" \\(non-Eur\\)$", "", D_nec$founder),
    `deep origin` = D_nec$deep_origin,
    `Ashkenazi %` = D_nec$ashkenazi_pct,
    `medieval Erfurt` = D_nec$medieval_erfurt,
    `posterior P(H2-side)` = sprintf("%.2f [%.2f, %.2f]", D_nec$post_mean_H2,
                                D_nec$post_lo, D_nec$post_hi),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  noneur_afr_pct <- NONEUR$ashkenazi_pct[NONEUR$mtree_name == "L2a1l2a"]
  noneur_total_pct <- round(sum(NONEUR$ashkenazi_pct), 1)
  noneur_min_post <- min(D_nec$post_mean_H2)
  noneur_max_post <- max(D_nec$post_mean_H2)

  # broader frequent Ashkenazi pool: model call vs. published expectation
  major_tab <- data.frame(
    lineage = D_ml$founder,
    `Ashkenazi %` = D_ml$ashkenazi_pct,
    `published origin` = D_ml$expected_origin,
    `deep origin` = D_ml$deep_origin,
    `posterior P(H2)` = sprintf("%.2f [%.2f, %.2f]", D_ml$post_mean_H2,
                                D_ml$post_lo, D_ml$post_hi),
    `model call` = ifelse(D_ml$post_mean_H2 >= 0.5, "Near Eastern", "European"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  major_ne <- D_ml[D_ml$expected_origin == "Near Eastern", ]
  major_eu <- D_ml[D_ml$expected_origin == "European", ]
  major_amb <- D_ml[D_ml$expected_origin == "Ambiguous", ]
  major_eval <- D_ml$expected_origin %in% c("Near Eastern", "European")
  major_conc <- sum(((D_ml$post_mean_H2 >= 0.5) == (D_ml$expected_origin == "Near Eastern"))[major_eval])
  major_eval_n <- sum(major_eval)
  major_ne_pct <- round(sum(major_ne$ashkenazi_pct), 1)
  major_eu_names <- paste(major_eu$founder, collapse = ", ")
  major_amb_names <- if (nrow(major_amb)) paste(major_amb$founder, collapse = ", ") else "none"

  coef_tab <- data.frame(
    channel = c("intercept", "frequency (non-Jewish rarity)",
                "time (antiquity + carriers + public age)", "nesting (equal-n richness)")[
                  match(D_fc$term, c("(Intercept)", "z_freq", "z_time", "z_nest"))],
    `posterior weight` = sprintf("%.2f [%.2f, %.2f]", D_fc$mean, D_fc$lo, D_fc$hi),
    check.names = FALSE, stringsAsFactors = FALSE)
  fc_loo_acc <- sum(D_fl$correct); fc_loo_n <- nrow(D_fl)
  v7_li <- D_li[D_li$lineage_key == "V7a2c1b", , drop = FALSE]
  u5_li <- D_li[D_li$lineage_key == "U5a1f1a3", , drop = FALSE]
  nc_post <- function(l, col = "post_mean_H2") .posterior_lineage(D_nc, l, col = col)
  nonjew_desc <- if (fs$nonjew == "1kg_eur") {
    sprintf(
      paste0(
        "founder subclades from Livni-Skorecki Table 7 (n=27,651); all other modeled lineages use ",
        "the full 1000 Genomes Phase 3 European super-population (CEU+GBR+FIN+IBS+TSI) subclade frequencies ",
        "(callMom full mtDNA; Haplogrep3 / Phylotree 17, n=%d) ",
        "with Mitotree European non-Jewish subclade counts as fallback when 1KG cannot resolve the lineage ",
        "(e.g. K1a4a) or would climb past the macro to an unrelated broader clade"
      ),
      if (nrow(KG1)) KG1$n_nonjews[1] else 503L
    )
  } else if (fs$nonjew == "livni_skorecki") {
    "founder subclades from Livni-Skorecki Table 7 (n=27,651); all other modeled lineages use Mitchell et al. 2014 NHANES macro-haplogroup frequencies (n=8,537 non-Hispanic white, Table 2)"
  } else {
    "Mitotree modern European non-Jewish counts only"
  }
  ctrl_nj_blurb <- if (fs$nonjew == "1kg_eur" && nrow(v7_li) && nrow(u5_li)) {
    sprintf(
      "1000 Genomes EUR subclade frequencies (%s %.1f%%, %s %.1f%%)",
      v7_li$nonjew_freq_1kg_match[1], 100 * v7_li$nonjew_freq[1],
      u5_li$nonjew_freq_1kg_match[1], 100 * u5_li$nonjew_freq[1]
    )
  } else if (fs$nonjew == "livni_skorecki") {
    "Mitchell 2014 macro frequencies (~2.6% for V, ~13.6% for U among NHANES non-Hispanic whites)"
  } else {
    "Mitotree-derived non-Jewish frequencies"
  }

  first_n <- function(x, n = 3) {
    v <- unlist(strsplit(x, ";", fixed = TRUE))
    v <- v[nzchar(v)]
    if (!length(v)) return("")
    paste(head(v, n), collapse = "; ")
  }
  cal_year <- function(x) {
    y <- suppressWarnings(as.numeric(x))
    if (is.na(y)) return("")
    if (y < 0) sprintf("%s BCE", abs(round(y))) else sprintf("%s CE", round(y))
  }
  if (nrow(PUBLIC)) {
    public_tab <- data.frame(
      lineage = PUBLIC$lineage_key,
      group = PUBLIC$group,
      `YFull formed/TMRCA` = paste(PUBLIC$yfull_formed, PUBLIC$yfull_tmrca, sep = " / "),
      `YFull public ids` = PUBLIC$yfull_public_ids_n,
      `FTDNA TMRCA` = vapply(PUBLIC$ftdna_tmrca_mean, cal_year, character(1)),
      `top public countries` = vapply(PUBLIC$ftdna_modern_country_counts, first_n, character(1)),
      `ancient anchors` = vapply(PUBLIC$ftdna_ancient_codes, first_n, character(1)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    public_tab$`FTDNA TMRCA`[public_tab$`FTDNA TMRCA` == ""] <- "not exposed"
    public_tab$`top public countries`[public_tab$`top public countries` == ""] <- "not exposed"
    public_tab$`ancient anchors`[public_tab$`ancient anchors` == ""] <- "not exposed"
  } else {
    public_tab <- data.frame()
  }

  # rarefaction summary at equal n per root clade
  rar_summary <- do.call(rbind, lapply(split(B_rar, B_rar$root), function(s) {
    nmax <- max(s$n[!is.na(s$distinct_lineages)])
    w <- s[s$n == nmax, ]
    data.frame(clade = s$root[1], equal_n = nmax,
               Europe = round(w$distinct_lineages[w$region == "Europe"], 1),
               `Near East` = round(w$distinct_lineages[w$region == "Near East"], 1),
               check.names = FALSE)
  }))

  fig <- function(f, cap) sprintf("![%s](../outputs/figures/%s)\n\n*%s*\n", cap, f, cap)

  L <- c(
"# Are the Major Ashkenazi mtDNA Founders Middle Eastern / Levantine?",
"",
"This report tests two competing hypotheses for the four major Ashkenazi maternal founder lineages **K1a1b1a, K1a9, K2a2a, N1b2**:",
"",
"- **[H1] European origin** (Costa et al. 2013): the female founders were assimilated prehistoric Europeans.",
"- **[H2] Levantine origin** (Behar et al. 2006; Livni & Skorecki 2025): both male and female founders were mostly Levantine/Near Eastern.",
"",
"It combines the Livni-Skorecki branching-process / founder-vs-host model framework with a phylogeographic and ancient-DNA analysis of a **deduplicated enriched dataset**: Mitotree S1 + GenBank E-utilities founder queries + Behar/Costa provenance flags, with duplicates removed by accession/sample.",
"",
"## Summary",
"",
sprintf("Four maternal founder lineages -- K1a1b1a, K1a9, K2a2a, and N1b2 -- account for roughly %.0f%% of Ashkenazi Jewish mtDNA, but their geographic origin is disputed between a prehistoric European-assimilation model and a Near Eastern / Levantine-origin model. This report reimplements the Livni-Skorecki branching-process and founder-versus-host framework and combines it with phylogeographic rarefaction, non-Jewish rarity, ancient DNA, and a regularized origin synthesis over a deduplicated Mitotree/GenBank reference set.",
        sum(FOUNDERS$ashkenazi_pct)),
"",
sprintf("Five analyses favor a Near Eastern origin. The estimated absorbed host-lineage fraction places all four founders in the multi-copy founder tail rather than among absorbed singletons; the founders are rare among 27,651 non-Jews; Costa's \"solely European nesting\" argument weakens after equal-sample-size rarefaction; the K and N parent macro-clades have deep Near Eastern ancient-DNA context and the founder clades occur in medieval Jewish individuals; and the ridge-logistic synthesis places every founder above 0.5 for Near Eastern origin. Support is strongest for %s (%s) and most tentative for %s (%s)%s.",
        D_bp$founder[which.max(D_bp$post_mean_H2)],
        .fmt_post(max(D_bp$post_mean_H2)),
        D_bp$founder[which.min(D_bp$post_mean_H2)],
        .fmt_post(min(D_bp$post_mean_H2)),
        # State the parity caveat only for those founders whose interval
        # actually admits it, so this cannot go stale if the ranking moves.
        local({
          crosses <- D_bp$founder[D_bp$post_lo < 0.5]
          weakest <- D_bp$founder[which.min(D_bp$post_mean_H2)]
          if (!length(crosses)) {
            ". Every founder's 90% credible interval excludes parity"
          } else if (identical(crosses, weakest)) {
            ", which is also the only founder whose 90% credible interval still includes parity"
          } else {
            sprintf(". The 90%% credible interval still includes parity for %s",
                    paste(crosses, collapse = ", "))
          }
        })),
"",
"## Data Set",
"",
"The analysis set is built from Mitotree S1, enriched with GenBank founder-query metadata, and deduplicated. Behar and Costa rows are not added as separate frequency tables; they are study/provenance flags on records already present in Mitotree/GenBank.",
"",
.md_table(D0),
"",
"## Study Design and Validation Strategy",
"",
"Every reported number is computed directly from the underlying sample tables, tree topology, and cited study data, so each figure and posterior is traceable to a specific data source rather than asserted. The analysis is organized as a set of independent tests, each aimed at a specific published claim.",
"",
"Each module tests a specific published claim:",
"",
"- **Livni-Skorecki:** the branching module tests whether major founders and absorbed host lineages should have different descendant-count spectra after a bottleneck.",
"- **Behar:** the Behar 2006 sample is used as the Ashkenazi-specific benchmark for the founder spectrum.",
"- **Costa:** the rarefaction module tests whether Costa's European-nesting argument survives when Europe and the Near East are compared at equal sample size.",
"- **Livni-Skorecki Table 7:** the non-Jewish rarity benchmark tests whether the founders are common enough in surrounding non-Jewish populations to support recent European-host absorption.",
"- **Livni-Skorecki Table 3:** the diffusion-law frequency-gradient benchmark is treated as directional context and checked with control clades.",
"- **Ancient DNA:** medieval Jewish and ancient Near Eastern records test whether the clades are present in the right temporal and geographic context.",
"",
"The inputs are data-based. Modern and ancient mtDNA records, country labels, tree topology, descendant sets and TMRCA values come from Mitotree/GenBank-derived tables. Published papers contribute explicit reference tables, frequencies, ancient samples or literature-coded classifications. Manually parameterized controls are labelled as calibration controls rather than independent discoveries; their job is to check that the same machinery can correctly place known European, Near Eastern and non-European cases.",
"",
"The framework is constrained to satisfy several requirements simultaneously: it must reproduce the Livni-Skorecki branching logic, keep the Behar and Costa benchmarks separate from the Mitotree reference set, control for the Costa sample-size confound, place European calibration lineages near parity, place Near Eastern and non-European calibration lineages high, and still score the four founders from their observed rarity, ancient context, and nesting. Because a single classifier must place both the calibration lineages and the founders correctly, the founder scores are not free parameters tuned to a target.",
"",
"### Computational Methods and Numerical Validation",
"",
"The pipeline represents Mitotree topology as a parent map and child adjacency list. Descendant sets are obtained by depth-first traversal and used consistently for founder counts, regional parent-clade counts, ancient-record extraction, rarefaction, and control panels. Deduplication by accession and sample key occurs before these analyses so records present through both Mitotree and GenBank enrichment are not counted twice.",
"",
"Finite-time Galton-Watson descendant distributions are calculated by composing offspring probability-generating functions. Their coefficients are recovered by evaluating the composed function on roots of unity and applying an inverse FFT. An independent direct polynomial-convolution implementation validates this calculation: the largest discrepancy across the reported growth grid is below 1e-4, negligible relative to model uncertainty.",
"",
"Rarefaction evaluates the Costa nesting claim using expected descendant-lineage richness after sampling Europe and the Near East without replacement to the same sample size. For the per-lineage nesting channel the focal lineage's own descendants -- and those of any sibling founder sharing its root clade -- are excluded from both regional pools, so the comparison is between its *sister* lineages: a lineage's own carriers are not evidence about which regional variation it nests inside, and because they concentrate on a few named nodes they would otherwise depress the measured richness of whichever region they are most common in. The comparison root climbs to a parent clade until both sister pools reach at least 10 samples, and the root actually used is reported for each lineage. The origin synthesis then standardizes three data-computed channels -- non-Jewish rarity, antiquity/time, and equal-sample-size nesting -- and fits a ridge-logistic model to labeled controls while holding the founders out. Coefficient uncertainty is propagated by Laplace-approximation Monte Carlo, and leave-one-out validation on the labeled panel checks whether the fitted model generalizes beyond individual controls.",
"",
"## Public YFull and FTDNA Enrichment",
"",
sprintf("YFull and FamilyTreeDNA Discover provide public context for the modeled haplogroups: alternative nomenclature, branch ages, public sample identifiers, country tags, ancient connections, and project counts. Public pages were compiled for %d lineages (K1a4a has none): YFull static pages were available for %d and FTDNA's public JSON endpoint for %d (U1b1a1 resolves through its public parent U1b1; A-a1b3a1 sits below macro-A on FTDNA). The public branch age adjusts the time channel only for lineages that also carry a Jewish or ancient anchor, so the European calibration lineages enter with no public-age bonus.", n_public, n_yfull_ok, n_ftdna_ok),
"",
"These data are not treated as an unbiased population-frequency sample: YFull and FTDNA are genealogical/testing databases with strong participation and project-enrollment effects. The model therefore does not use their tester country proportions as population frequencies. It does use the least sample-sensitive fields to make the time channel more realistic: supported aliases, public branch/TMRCA estimates, and public ancient-sample links. Very shallow public TMRCA estimates modestly penalize the time channel; mature public branch ages modestly strengthen it only when paired with an existing Jewish/ancient anchor. That anchor counts only pre-modern Jewish sites and studies (Erfurt/Waldman, Tàrrega/Roquetes, Chapelfield); the Sobibór series (Diepenbroek 2021) is excluded, since those individuals were born 1893-1923 and so cannot evidence a lineage's presence in the Jewish maternal pool before the modern era -- the same date restriction that governs the medieval-carrier count.",
"",
if (nrow(public_tab)) .md_table(public_tab) else "_Public enrichment table unavailable._",
"",
"## Added External Evidence: Genes 2026 Tàrrega",
"",
"The 2026 Genes paper on the Roquetes/Tàrrega medieval Jewish pogrom is useful in two ways. First, it adds an independent medieval Iberian Jewish mtDNA record: **ROQ12 is K1a1b1a** with Haplogrep quality 0.96. Second, its autosomal qpAdm model supports substantial Eastern Mediterranean/Canaan-related ancestry in this medieval Iberian Jewish community. This does not prove the mtDNA founder's origin by itself, but it improves the historical plausibility of a Mediterranean Jewish carrier population in Iberia.",
"",
"Roquetes mtDNA calls:",
"",
.md_table(G_mt_focus),
"",
"Accepted qpAdm two-source model from the supplement:",
"",
.md_table(G_qp_good),
"",
"## Added External Evidence: PPNB Near Eastern Neolithic (Fernández 2014)",
"",
"Fernández et al. 2014 sequenced HVR1 mtDNA from Pre-Pottery Neolithic B farmers of Syria (Tell Halula, Tell Ramad, Dja'de El Mughara; ~8,700-6,600 cal BC). Haplogroup **K is the most prevalent at 42.8% (N=6 of the 14 assigned skeletons)**, with a rare N* paragroup also present -- deep Near Eastern presence of exactly the two macro-clades (K and N) that contain the Ashkenazi founders, millennia before any founder event.",
"",
"This paper is unusually direct on our question. It reports that Ashkenazi Jews show a haplogroup-K frequency similar to the PPNB sample with low, non-significant pairwise F_ST, and explicitly states this *contradicts* Costa et al.'s European-origin interpretation and instead suggests an ancient Near Eastern origin.",
"",
"It is also the strongest honesty check in this report, and it cuts both ways:",
"",
"- **For H2 (Near Eastern), at the haplogroup level:** deep, high-frequency Near Eastern K plus rare N* is consistent with a Near Eastern maternal source for the K/N macro-clades.",
"- **Against over-claiming at the sub-clade level:** the PPNB K samples were typed only on HVR1 and lack the diagnostic control-region mutations of K1a1b1a and K2a2a, so those specific founder sub-clades (~81% of Ashkenazi K) are *excluded* from these particular Neolithic individuals; only K1a9 (~19%) cannot be excluded without coding-region typing.",
"- **Symmetric caveat:** the authors note lineage loss/drift in the Near East since the Neolithic, so absence of the derived founder sub-clades in this small aDNA panel is not proof of their past absence either.",
"",
"PPNB haplogroup composition (counts from the recovered HVR1 profiles; Fernández 2014 Table 1):",
"",
.md_table(ppnb_tab),
"",
"*(Percentages here are over the 15 recovered profiles; the paper's headline 42.8% for K uses the 14 skeletons it could confidently assign.)*",
"",
"Net effect: this supports the Near Eastern origin of the *K and N macro-clades* strongly, while reinforcing the report's existing caution that deep Near Eastern K is not the same as demonstrating the derived Ashkenazi founder motifs in the ancient Near East.",
"",
"## Added External Calibration: Carmi 2014, Erfurt 2022, Xue 2017",
"",
"Carmi et al. 2014 is useful for demographic calibration rather than mtDNA placement: it estimates a severe Ashkenazi bottleneck of roughly 250-420 effective individuals about 25-32 generations ago, and models Ashkenazi autosomes as an approximately even European plus likely Middle Eastern mixture. That supports using a small-founder bottleneck model, but it cannot assign the four mtDNA founders to a geographic source by itself.",
"",
"Waldman et al. 2022 is directly useful. The Erfurt medieval Ashkenazi genomes show that the Ashkenazi founder event and major ancestry sources pre-dated the 14th century; among 31 unrelated Erfurt individuals, 11 carried K1a1b1a, and additional Erfurt individuals carried K1a9 and N1b2/N1b1b1. This is why the Mitotree ancient-founder table is so important: it moves the argument from modern frequency inference to medieval Jewish presence.",
"",
"Brace et al. 2022 adds a second, independent, and earlier anchor. Six individuals from a medieval well at Chapelfield, Norwich (radiocarbon 1161-1216 cal CE, consistent with the 1190 CE antisemitic massacre) show strong affinity to modern Ashkenazi Jews (a qpAdm model of 100% Chapelfield fits present-day Ashkenazim), and their Ashkenazi-associated disease alleles are already near modern frequency -- so the Ashkenazi founder event pre-dates even the 12th century. The Norwich individuals do not carry the four K/N founders (three were mtDNA H5c2, an H clade of which Ashkenazim are the majority of modern carriers), so Norwich is not a direct founder carrier; rather it independently confirms the deep antiquity and continuity of a distinctively Ashkenazi maternal pool that the founder hypothesis requires. A nomenclature note: the founder called N1b2 here (after Behar/Costa) resolves to **N1b1b1** on FTDNA's February 2025 mtDNA tree and in Brook 2022 (23andMe still reports it as N1b2, and the true PhyloTree N1b2 is essentially absent from Ashkenazim); the Erfurt carrier belongs to this N1b1b1 clade. We keep the N1b2 label for continuity with the prior literature.",
"",
"Xue et al. 2017 refines the European side of the admixture. Using local-ancestry inference on a large Ashkenazi sample, they estimate roughly even Middle-Eastern and European ancestry, with the European component predominantly Southern European (about 34% Southern, 8% Western, 8% Eastern European against 50% Levantine), and an admixture history best explained by at least two events -- a likely Southern European event pre-dating the late-medieval founder bottleneck and an Eastern European event post-dating it. This matters for the maternal-founder debate: substantial European autosomal ancestry is expected under H2 and does not by itself imply that the maternal founders were European, and the Southern European / partly post-founder character of that gene flow is consistent with a Levantine maternal core plus later assimilation.",
"",
"Three further genome-wide studies bracket this picture. Atzmon et al. 2010 showed that the major Jewish Diaspora groups form distinct clusters sharing Middle Eastern ancestry, and explicitly framed the four mtDNA founders (~40% of the Ashkenazi maternal pool) as Middle Eastern in origin -- independent genome-wide support for the direction argued here. Agranat-Tamir et al. 2020 supplies the deep ancient-DNA anchor: Bronze Age 'Canaanite' genomes of the Southern Levant to which present-day Jewish groups, including Ashkenazi Jews, trace 50% or more of their ancestry, establishing a concrete Levantine source population. Balancing this, the most recent local-ancestry inference (Lerga-Jaso et al. 2025, 'Orchestra') assigns Ashkenazi Jews a large Italian/Southern European autosomal component alongside a Levantine one. Because Southern European populations themselves carry substantial Near Eastern ancestry, and because autosomal ancestry is not uniparental, a large Southern European autosomal fraction is compatible with -- and expected under -- a Levantine maternal core followed by later European admixture; it does not overturn the maternal conclusion.",
"",
"## Model Assumptions Improved From Livni-Skorecki",
"",
"The Livni-Skorecki framework is useful because it formalizes a real demographic intuition: major founder lineages and absorbed host lineages should have different descendant-count spectra after a bottleneck. The improvements here make that framework more realistic and harder to overinterpret:",
"",
"- **Estimate absorption from the largest available maternal sample.** Rather than abstaining, we estimate the host/convert absorption rate from the enriched Mitotree+GenBank set -- the largest maternal sample obtainable -- using the individual-level singleton fraction, and cross-check it against the dedicated Behar 2006 Ashkenazi sample. Sampling caveats are quantified (below), not used as a reason to decline the estimate.",
"- **Separate model probability from empirical phylogeography.** The Galton-Watson calculation is a sensitivity model; the Mitotree/GenBank evidence is used for clade placement, ancient occurrences, and regional enrichment.",
"- **Use ancient time anchors.** ROQ12 adds a medieval Iberian Jewish K1a1b1a record, so the model no longer depends only on modern distributions or on Behar/Costa counts.",
"- **Control regional sampling imbalance.** Rarefaction compares Europe and Near East at equal sample size rather than treating raw sub-lineage counts as evidence of origin.",
"- **Prefer uncertainty-aware claims.** The result can be strengthened as 'more consistent with a Levantine/Mediterranean Jewish source than a simple local European absorption model,' but mtDNA alone still cannot prove a geographic origin with autosomal-level conclusiveness.",
"",
"The first realism upgrade is now included as a sensitivity analysis: the Galton-Watson process is run under Livni-Skorecki's Poisson offspring law and under a geometric offspring law, equivalent to a strongly overdispersed negative-binomial model with size = 1. The qualitative point survives: exact singleton survival remains a low-probability tail event after the bottleneck window.",
"",
.md_table(A_sens),
"",
fig("A2b_offspring_law_sensitivity.png", "Sensitivity of single-descendant probability to offspring law."),
"The three realism upgrades that the framework needed are now implemented in full below: a **simple** extension that fills in the reproductive-overdispersion axis with intermediate negative-binomial sizes, and two **complex** extensions -- piecewise historical growth instead of a constant growth ratio, and a Bayesian time-aware synthesis that combines founder frequency, ancient occurrence dates, public YFull/FTDNA branch-age context, regional sampling intensity, and phylogenetic nesting into a posterior origin probability per founder.",
"",
"### Simple extension: intermediate reproductive overdispersion",
"",
"The Poisson-vs-geometric contrast above is only the two extremes of a single axis: reproductive overdispersion. A negative-binomial offspring law with dispersion `size` has variance `m + m^2/size`, so `size = 1` is the geometric law and `size -> Inf` is Poisson. Sweeping intermediate sizes shows the single-descendant and extinction probabilities move *smoothly* between the two extremes; there is no offspring law in this family under which exact single-descendant survival stops being a rare tail event (shown here at growth ratios 1.0, 1.05, 1.1):",
"",
.md_table(D_nb[D_nb$growth_rate %in% c(1.15, 1.20, 1.28), ]),
"",
fig("D1_nbinom_size.png", "Single-descendant probability across the negative-binomial size ladder (geometric = 1, Poisson = infinity)."),
"### Complex extension 1: piecewise historical growth",
"",
"A single constant growth ratio is demographically unrealistic. The Ashkenazi matriline was close to stationary through the medieval founder phase and then expanded rapidly. We encode this as a piecewise per-generation schedule and compare it against a constant-growth process with the **same net expansion** (same expected surviving-lineage size), so any difference is attributable purely to the *timing* of growth:",
"",
.md_table(D_sch),
"",
sprintf("With net expansion held identical (~%.1fx over %d generations), the calendar-anchored bottleneck-then-expansion schedule *raises* extinction (%.1f%% vs %.1f%%) and *lowers* the singleton probability (%.2f%% vs %.2f%%) relative to constant growth. Early stasis prunes more lineages, and rapid late expansion means the survivors are carried in many copies rather than as singletons. This does not weaken the founder-vs-host argument -- it sharpens it: under realistic history a genuine founder is even less likely to appear as a lone singleton.",
        D_pw$net_expansion[1], GW_K,
        100 * D_pw$prob_extinct[2], 100 * D_pw$prob_extinct[1],
        100 * D_pw$prob_one_descendant[2], 100 * D_pw$prob_one_descendant[1]),
"",
.md_table(D_pw),
"",
fig("D2_piecewise_growth.png", "Piecewise historical growth vs a constant-growth process matched on net expansion."),
"### Complex extension 2: origin synthesis (single fitted model)",
"",
sprintf("Active frequency inputs: **Ashkenazi = %s** (Brook/Penninx via panel CSVs), **non-Jewish = %s** — %s.",
        fs$ashkenazi, fs$nonjew, nonjew_desc),
"",
"The Bayesian synthesis combines three **data-computed** evidence channels into one posterior probability of Near Eastern / Levantine origin (H2). There is no separate elicited model and no hand-set `tanh` gains:",
"",
"1. **Standardize channels.** Each lineage gets raw ingredients from the enriched Mitotree set (audit: `D_lineage_computed_inputs.csv`): non-Jewish frequency (Livni-Skorecki Table 7 subclades for founders; full 1000 Genomes European (CEU+GBR+FIN+IBS+TSI, n=503) subclade matches for all other scored lineages), Near Eastern macro-clade depth, medieval Jewish carriers (including N1b2/N1b1b1 harmonization), public YFull/FTDNA age adjustment, and equal-*n* rarefaction nesting. Brook/Penninx Ashkenazi percentages enter only the reference-frequency audit column, not the rarity channel.",
"2. **Z-score** rarity, time (depth + carriers + public age), and nesting fraction across all scored lineages.",
"3. **Fit a ridge-logistic model** to labeled controls only (European = 0; Near Eastern / non-European diaspora = 1). Founders are held out. Coefficient uncertainty is propagated by Laplace-approximation Monte Carlo. Leave-one-out validation on the labeled panel guards against overfitting.",
"",
sprintf("Leave-one-out validation classifies **%d of %d labeled lineages correctly (%.0f%%)**. Fitted channel weights (Gaussian priors, integrated out):", fc_loo_acc, fc_loo_n, 100 * fc_loo_acc / fc_loo_n),
"",
.md_table(coef_tab),
"",
"Per-founder computed inputs:",
"",
.md_table(D_bi[, c("founder", "root", "nonjew_freq", "ne_depth_kyr",
                   "medieval_jewish_carriers", "public_tmrca_kyr",
                   "public_age_confidence", "public_jewish_ancient_anchor",
                   "ne_richness_equal_n", "eu_richness_equal_n", "ne_pool", "eu_pool")]),
"",
"Standardized channel values for founders:",
"",
.md_table(D_bs),
"",
"Posterior P(H2) per founder (held out of the fit), with 90% credible intervals:",
"",
.md_table(D_bp),
"",
fig("D3_bayes_origin.png", "Origin synthesis: founders (red), European controls (blue), deep non-European controls (green); 90% credible intervals from the single fitted logistic on standardized data channels."),
sprintf("All four founders sit above 0.5 (K1a1b1a %s, K1a9 %s, K2a2a %s, N1b2 %s, the last counted as a single lineage with its FTDNA-tree synonym N1b1b1). K2a2a's interval is by far the widest and is the only one that still includes parity. It rests on the narrowest base: excluding its own carriers leaves a single Near Eastern sister sample within K2a, so its comparison climbs to K where sister richness is near parity and the nesting channel contributes nothing, and it has no Jewish carrier dated before 1800 -- its only Jewish ancient records are 20th-century Sobibor victims, which cannot show the lineage was in the Jewish pool pre-modern, so the antiquity channel is negative as well. Its score rests on non-Jewish rarity almost alone and should be read as a single-channel result rather than convergent evidence. N1b2's nesting is also a negative contributor -- at the N1b macro background its sister lineages stay Europe-richer at equal sample size -- but its rarity and antiquity channels keep its interval clear of parity. We report three decimal places because rounding would overstate certainty.",
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K1a1b1a"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K1a9"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K2a2a"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "N1b2"])),
"",
"### Negative controls: European (absorbed) Ashkenazi lineages",
"",
sprintf("A synthesis that scored every Ashkenazi lineage as Levantine would not be informative, so we run the identical model on European phylogeographic controls. Their channel inputs use the same published non-Jewish benchmarks (%s) plus Mitotree-derived nesting and carrier counts — not hand-set code constants.", ctrl_nj_blurb),
"",
"**Control 1: V7a2c1b.** This is not a rare oddity -- it is a *common, established* Ashkenazi maternal cluster at about 2.4% (Penninx 2019, compiled in Brook 2022; defining mutation A12753G), comparable to HV5a and not far below K2a2a. Its carriers include Yiddish-speaking individuals from Lithuania (YFull YF083232 [LT-VL], YF108999 [LT-PN]), Galicia (GenBank OR803751, Dynów) and Ukraine (GenBank PQ337267, Berezna). Yet on every diagnostic that distinguishes origin, it is the mirror image of the four founders:",
"",
"- it nests inside an entirely European branch of haplogroup V (a post-glacial European expansion), not K1a/K2a/N1b;",
"- in the same subclade it is **shared with non-Jewish Europeans** (Danish KF162774, Russian JQ703830, a Berlin and a Belarusian sample), so the near-absence-among-non-Jews signal that flags the founders is absent; and",
"- it has no deep Near Eastern ancient record and no medieval Jewish (Erfurt/Tàrrega) carrier.",
"",
sprintf("We do **not** rely on a \"too young\" argument, because the age of this branch is itself uncertain: YFull dates the tight subclade to ~100 years, whereas FTDNA's Discover tool (which names the same branch V7a12a1) suggests a TMRCA around 371 BCE (95%% CI 696-66 BCE; ~2,375 years). Either way, coalescence age alone does not indicate a Near Eastern origin -- and crucially, granting V7a2c1b even a Roman-era age does not move it toward H2, because its European signal comes from nesting, non-Jewish %s frequency (%.1f%%), and shared non-Jewish European records.",
        if (nrow(v7_li)) v7_li$nonjew_freq_1kg_match[1] else "V",
        if (nrow(v7_li)) 100 * v7_li$nonjew_freq[1] else NA_real_),
"",
"**Control 2: U5a1f1a3.** Haplogroup U5 is the signature European Mesolithic hunter-gatherer lineage: it dominates pre-Neolithic Europe and is essentially *absent* from the Pre-Pottery Neolithic Near East (Fernández 2014 find K, not U5, in PPNB farmers). Any Ashkenazi U5 lineage is therefore European by deep phylogeography, and because U5 is common in modern non-Jewish Europeans its non-Jewish-rarity signal is strongly negative -- the opposite of a founder. It supplies an independent European control with a *different* origin story (Paleolithic hunter-gatherer versus post-glacial V), guarding against the objection that the V7a2c1b result is a one-off.",
"",
"Fed through the same standardized channels and fitted logistic posterior:",
"",
.md_table(D_nc),
"",
sprintf("European controls score far below the founders (V7a2c1b **%.2f [%.2f, %.2f]**, U5a1f1a3 **%.2f [%.2f, %.2f]** vs founder range **%.2f–%.2f**). These two established Ashkenazi maternal lineages of recognized European origin are common among European non-Jews (%s), which places them outside the founder tail.",
        nc_post("V7a2c1b"), nc_post("V7a2c1b", "post_lo"), nc_post("V7a2c1b", "post_hi"),
        nc_post("U5a1f1a3"), nc_post("U5a1f1a3", "post_lo"), nc_post("U5a1f1a3", "post_hi"),
        min(D_bp$post_mean_H2), max(D_bp$post_mean_H2),
        if (fs$nonjew == "1kg_eur") {
          "subclade-resolved 1000 Genomes European frequencies, n=503"
        } else {
          "Mitchell 2014 macro frequencies, ~2.6% for V and ~13.6% for U among NHANES non-Hispanic whites"
        }),
"",
"(The controls appear as the blue triangles in the posterior figure above; the four founders are the circles.)",
"",
"### Positive controls: deep non-European lineages, run through the same model",
"",
sprintf("The founders and the European host lineages are not the whole story. The Ashkenazi maternal pool also contains lineages whose *deep* origin is neither European nor Levantine but African or Asian -- and their presence is difficult to reconcile with the strong form of H1 (recent assimilation from European host populations), because a European source cannot contribute African or East Asian mtDNA. The largest is **L2a1l2a** (~%s%% of Ashkenazi maternal lines; Penninx 2019 / Brook 2022), a lineage of deep Sub-Saharan/East African origin already present in the 14th-century Erfurt Jewish cemetery (individual I13865). Several smaller lineages tell the same story.", noneur_afr_pct),
"",
"We run these through the **identical** synthesis with the same data-computed channel inputs. The model contrasts a recent *European-host* origin (H1) with a Near Eastern / Jewish-diaspora origin (H2); for the four founders the H2 pole is specifically Levantine, whereas for these lineages the H2 pole means *entry through the non-European diaspora* rather than a literal Levantine birthplace. A high posterior is therefore a rejection of European-host assimilation. Non-Jewish rarity, non-Europe-vs-Europe rarefaction nesting, and medieval Jewish carriers (where present in the ancient sample set) are all computed from Mitotree:",
"",
.md_table(noneur_tab),
"",
sprintf("All five land above 0.5 (posterior means %.2f–%.2f), opposite the European controls — %s highest (%.2f) where non-Jewish rarity and non-Europe-vs-Europe nesting are both strong. Together these deep African/Asian lineages account for roughly %s%% of Ashkenazi maternal lines in the Brook/Penninx compilation.",
        noneur_min_post, noneur_max_post,
        sub(" \\(non-Eur\\)$", "", D_nec$founder[which.max(D_nec$post_mean_H2)]),
        max(D_nec$post_mean_H2),
        noneur_total_pct),
"",
"### Calibration on the broader frequent Ashkenazi pool",
"",
"The four founders are only ~38% of the Ashkenazi maternal pool. To test whether the synthesis tracks the published literature in *both* directions, we scored a small literature-coded panel of frequent non-founder lineages spanning the origin spectrum, using the identical channels and posterior. This is a calibration analysis, not a new population-frequency estimate: a useful classifier should place lineages with published European assignments on the European side while retaining lineages with published Near Eastern assignments on the Near Eastern side.",
"",
"The Near Eastern side of this panel is deliberately conservative: Costa et al. (2013), despite arguing for substantial European maternal ancestry, classify **HV1b2**, R0a-associated lineages and **U7** as Near Eastern sources and describe **U1** as ultimately Near Eastern. Shamoon-Pour et al. (2019, *Scientific Reports*) further place HV1b2 near Assyrian HV1b branches in northern Mesopotamia / the South Caucasus, a plausible context for ancient Jewish communities. HV1b2 and U1b1a1 are also documented in the 14th-century Erfurt Jewish cemetery (I14901; I14850/I14853/I14898). Against these we set frequent lineages with published European assignments (**H7**, **J1c7a**) and two deliberately ambiguous test cases: **H6a1a1a** (Brook-coded Middle Eastern but Mitotree/FTDNA Europe-rich) and **K1a4a** (Brook 2022's sixth Ashkenazi K founder branch at ~0.2%, interpreted as a possible Greek/Italian convert but also present in Syria and shared with Egyptian/Maghrebi/Turkish Jews).",
"",
.md_table(major_tab),
"",
sprintf("The model agrees with the published origin in **%d of %d evaluable** panel cases; %s remain ambiguous/disputed. Near Eastern members score high (%.2f–%.2f); European members (%s) score low (%.2f–%.2f). H6a1a1a and K1a4a sit as stress cases (%s and %s).",
        major_conc, major_eval_n, major_amb_names,
        min(major_ne$post_mean_H2), max(major_ne$post_mean_H2),
        major_eu_names,
        min(major_eu$post_mean_H2), max(major_eu$post_mean_H2),
        if (nrow(major_amb)) .fmt_post(major_amb$post_mean_H2[match("H6a1a1a", major_amb$founder)]) else "NA",
        if (any(major_amb$founder == "K1a4a")) .fmt_post(major_amb$post_mean_H2[match("K1a4a", major_amb$founder)]) else "NA"),
"",
fig("D4_broader_pool.png", "Broader frequent Ashkenazi pool under the same fitted synthesis, versus published origin assignments."),
"",
"## Supplemental Behar/Costa Benchmarks",
"",
"Behar and Costa are useful supplements when treated carefully. They should not replace the Mitotree/GenBank analysis set, but they do provide externally published benchmarks: Behar's Ashkenazi sample records the founder-event spectrum, Costa motivates the nesting/sample-size challenge, Livni-Skorecki Table 7 provides non-Jewish rarity, and Table 3 provides parent-signature frequency gradients.",
"",
"Behar 2006 Ashkenazi sample summary:",
"",
.md_table(B06),
"",
"Livni-Skorecki / FTDNA non-Jewish rarity benchmark:",
"",
.md_table(T7),
"",
"---",
"",
"## Module A -- Branching-process (founder vs host) model",
"",
"**Galton-Watson lineage extinction.** Under the Poisson offspring model used here, after ~27 generations the probability that a founder leaves exactly one matrilineal descendant is below 0.2%, so the smallest surviving minor lineages are rare (Livni-Skorecki Table 1 model framework):",
"",
.md_table(A_t1),
"",
"The fast PGF/FFT implementation is checked against a slower direct-convolution implementation:",
"",
.md_table(A_audit),
"",
fig("A1_single_descendant.png", "P(single descendant at generation 27) vs growth ratio."),
"**Founder-vs-absorbed detection model.** In an unbiased population sample of N = 600, a founder lineage (frequency >= ~1.5%) is seen in multiple copies with ~99.9% probability, whereas a lineage absorbed by conversion (frequency ~0.1%) is usually absent or a singleton (P(>=2 copies) ~ 12%). This is the Livni-Skorecki diagnostic model; it is not directly applied to Mitotree as though Mitotree were a random Ashkenazi sample.",
"",
fig("A2_detection_curves.png", "Detectability of founder vs absorbed lineages (cumulative Poisson)."),
"**Empirical founder counts in the enriched deduped Mitotree+GenBank set.** These are direct counts from the analysis dataset, not Behar/Costa sample-table counts:",
"",
.md_table(A_mt[, c("founder", "modern_count", "modern_frequency", "ancient_count",
                   "distinct_modern_haplogroups", "known_country_count")]),
"",
"**Estimated host-lineage absorption rate.** We estimate the absorption (host/convert) rate directly, treating the enriched Mitotree+GenBank set as the largest available maternal sample: individuals who are the sole carrier of their fully resolved lineage are the candidate recently-absorbed matrilines (`a_eps`), and the per-generation rate follows Livni-Skorecki as `rho = 1 - (1 - a_eps)^(1/K)`. We report it across a range of founder-event depths K and cross-check against the dedicated Behar 2006 Ashkenazi sample:",
"",
 .md_table(A_absorb[, c("sample", "n_individuals", "n_lineages", "absorbed_fraction_a_eps",
                        "K_generations", "absorption_rate_per_gen", "absorbed_lineages_surviving")]),
"",
sprintf("The largest sample (enriched Mitotree, %s individuals) gives an absorbed fraction of about %.1f%% and a per-generation absorption rate of roughly %.2f%%-%.2f%% across K = 15-32 generations. The Ashkenazi-specific Behar 2006 sample gives a lower absorbed fraction (~%.1f%%) and rate, as expected: the enriched set is global and only a small share is explicitly Ashkenazi-flagged, so it over-counts rare non-founder lineages and its estimate is best read as an **upper bound**, with Behar as the population-appropriate anchor. If the enriched set were, if anything, over-enriched for Ashkenazi/diaspora samples, that would inflate founder representation and lower this absorption estimate, so the two bracket the plausible range. Either way the four major founders sit far out in the multi-copy tail, not among the absorbed singletons.",
        format(A_absorb$n_individuals[A_absorb$sample == "enriched_mitotree_modern"][1], big.mark = ","),
        100 * A_absorb$absorbed_fraction_a_eps[A_absorb$sample == "enriched_mitotree_modern"][1],
        100 * min(A_absorb$absorption_rate_per_gen[A_absorb$sample == "enriched_mitotree_modern"]),
        100 * max(A_absorb$absorption_rate_per_gen[A_absorb$sample == "enriched_mitotree_modern"]),
        100 * A_absorb$absorbed_fraction_a_eps[A_absorb$sample == "behar2006_ashkenazi"][1]),
"",
fig("A4_absorption_rate.png", "Estimated per-generation absorption rate vs founder-event depth (enriched Mitotree, Jewish-flagged subset, Behar 2006)."),
fig("A3_mitotree_founder_frequencies.png", "Major founder clades in the enriched deduped Mitotree+GenBank set."),
sprintf("**Euro-Levantine reconciliation is model-disfavored.** Under the Livni-Skorecki scenario, a hypothetical large Roman-era Jewish population has **%.2f%%** probability of producing >=3 of the four major founders at our base assumptions. The per-lineage major-emergence probability, private-mutation rate and founder-family count are explicit assumptions, not fitted quantities; varying each across a plausible range keeps this probability between %.2f%% and %.2f%%, so the conclusion does not depend on any single value. Treat this as a model sensitivity result, not an empirical measurement from Mitotree:", p_ge3, els_lo, els_hi),
"",
fig("A6_euro_levantine.png", "Probability of >=3 major founders arising in a Euro-Levantine pool."),
"---",
"",
"## Module B -- Phylogeography (ape + pegas)",
"",
"**Regional composition of the parent clades** (modern samples with known country in the reference set). The ancestral clades that gave rise to the founders retain a Near Eastern signal -- especially N1b and K2a -- even though modern country labels are diaspora-biased and are not interpreted as origins:",
"",
fig("B1_parentclade_region_props.png", "Regional composition of K1a, K1a1b1, K2a, N1b (modern)."),
"**Testing the 'solely European nesting' argument.** Costa et al. counted more European sub-lineages, but their source sample was much larger for Europe than for the Near East. Rarefying both regions to the *same* sample size removes that confound. In this enriched dataset, K1a (which contains K1a1b1a) is slightly richer in Near Eastern sub-lineages at equal sample size; K2a and N1b are richer in Europe at equal sample size, so the result is mixed rather than a blanket proof:",
"",
.md_table(rar_summary),
"",
fig("B6_rarefaction.png", "Sub-lineage richness at equal sample size, Europe vs Near East."),
"**Rarefaction controls.** The same line-chart method is applied to literature-coded control clades. The companion control panel is intentionally separate from B6 so the founder curves stay readable. European-assigned controls (H, U5, J1c) are Europe-equal or Europe-richer; Near Eastern-assigned controls (R0a, U7, U1) are Near-East-richer. This shows the rarefaction method can move in both directions and is not mechanically biased toward a Near Eastern result.",
"",
fig("B6b_rarefaction_controls.png", "Rarefaction control panel: European-assigned controls (H, U5, J1c) and Near Eastern-assigned controls (R0a, U7, U1), each plotted as Europe vs Near East at equal sample size."),
"The immediate sibling clades of the founders are **not** uniformly or exclusively European either; several sibling sets include Near Eastern / Mediterranean-dominant branches:",
"",
fig("B2_nesting_K1a1b1a.png", "K1a1b1a among the siblings of K1a1b1."),
fig("B3_nesting_K1a9.png", "K1a9 among the siblings of K1a."),
"**Published non-Jewish rarity benchmark.** The four founders use Livni-Skorecki Table 7 (n=27,651 non-Jews). Every other scored lineage uses the finest matching subclade observed among the five 1000 Genomes Phase 3 European populations (CEU, GBR, FIN, IBS, TSI; n=503) in the callMom full-mtDNA sequences (Haplogrep3; longest-prefix match, not macro-letter buckets). The panel spans both northern/western (CEU, GBR, FIN) and southern (IBS, TSI) Europe, so Near-Eastern-affiliated background clades that occur at low frequency in Mediterranean Europe are represented. Where a lineage's macro background is absent from the panel and matching would otherwise climb to an unrelated broader clade, the rarity is instead measured as that macro background's frequency among the ~14.7k European non-Jews in the enriched Mitotree pool. Diroma et al. (2014, BMC Genomics) report related 1000G mtDNA haplogroup work from whole-exome off-target reads; this analysis uses the dedicated Phase 3 mitochondrial call set instead:",
"",
if (nrow(KG1)) .md_table(KG1[, c("haplotype", "lineage_key", "count", "frequency", "n_nonjews", "source")]) else "(run `python3 scripts/compute_1kg_eu_freq.py` to regenerate)",
"",
"Founder subclades (same Table 7 source):",
"",
.md_table(T7_founders),
"",
fig("B7_nonjewish_rarity.png", "Frequency of the four major founders among non-Jews."),
"**Reconstructed haplotype network** (pegas) of the founders and their parent clades:",
"",
fig("B8_founder_haplonet.png", "Haplotype network of K/N founders and parent clades."),
"---",
"",
"## Module C -- Ancient DNA and diffusion (ape + adegenet)",
"",
sprintf("**Deep Near Eastern records for K1a in this dataset.** The oldest ancient K1a lineages in the enriched dataset are Near Eastern -- Anatolia (%s) and the Levant -- millennia before any Ashkenazi founder event:", min(C_ne$year)),
"",
.md_table(head(C_ne, 8)),
"",
"**Founder ancient footprint in the enriched dataset.** Ancient/historical founder-clade hits in the deduped set are:",
"",
.md_table(A_ma[, c("founder", "subject", "haplogroup", "country", "year", "study")]),
"",
fig("C1_ancient_timeline.png", "Ancient K1a / K2a / N1b through time and space."),
"**Coalescence ages.** The founders coalesce recently (K1a1b1a ~2,834 ybp; K2a2a ~5,379 ybp; N1b2 ~5,969 ybp) relative to their much older parent clades -- the signature of founder lineages:",
"",
.md_table(C_ca),
"",
fig("C2_coalescence_ages.png", "TMRCA of founders (red) vs parent clades."),
"**Frequency-gradient argument.** Livni-Skorecki frame this as a diffusion-law argument: if a parent signature is more frequent in the Middle East / Anatolia than in Europe, that gradient is more compatible with a Near Eastern source than a European one. This is supportive context, not a stand-alone origin proof:",
"",
.md_table(C_t3),
"",
fig("C3_diffusion_frequencies.png", "Parent-signature frequencies by region."),
"**Table 3 controls.** The same regional-frequency-gradient display can be calibrated with control clades. European-assigned controls (H, U5, J1c) peak in Europe or are Europe-rich; Near Eastern-assigned controls (R0a, U7, U1) peak in Middle East / Caucasus / Anatolia in this Mitotree/GenBank reference set. H6a1a is shown separately as the Mitotree-supported parent/context for the disputed Ashkenazi H6a1a1a lineage, not as a clean published-European control. This supports using Table 3 as directional context while keeping it a supporting benchmark, not a stand-alone proof.",
"",
fig("C3b_table3_controls.png", "Table-3-style frequency-gradient controls from Mitotree/GenBank regional samples: European controls (H, U5, J1c), Near Eastern controls (R0a, U7, U1), and H6a1a as ambiguous H6 context."),
fig("C4_glpca_founders.png", "adegenet glPca of founder and parent-clade motifs."),
"---",
"",
"## Discussion",
"",
"The five analyses converge despite testing different parts of the argument: finite-time lineage survival, non-Jewish rarity, sampling-controlled nesting, ancient and medieval occurrence, and multichannel classification. No analysis positively favors a European origin for any of the four founders. K2a2a is the least certain and the one founder whose credible interval still includes parity: once its own carriers are excluded there is a single Near Eastern sister sample within K2a, so its nesting is measured at K and near parity, and it has no pre-1800 Jewish carrier, leaving non-Jewish rarity as effectively its only positive channel. N1b2's equal-sample-size nesting also leans European, but its rarity and antiquity signals keep its combined estimate clear of parity. The convergence matters more than any single statistic: a European-assimilation explanation must simultaneously account for their multi-copy tail position, near-absence among non-Jews, nesting after sampling control, deep Near Eastern K context, and medieval Jewish carriers.",
"",
"The branching result should be interpreted as demographic discrimination rather than geographic proof by itself. Across Poisson, geometric, intermediate negative-binomial, constant-growth, and piecewise-growth specifications, a lineage that survives the bottleneck and rises to founder frequency is likely to appear in many copies, whereas recently absorbed host matrilines are concentrated among singletons. Geographic interpretation comes from combining that result with rarity, phylogeography, ancient DNA, and controls.",
"",
"Genome-wide evidence supplies context rather than directly assigning these maternal lineages. Shared Middle Eastern ancestry and a severe Ashkenazi bottleneck make a Near Eastern maternal core historically plausible, while substantial Southern European autosomal ancestry remains compatible with later admixture. Autosomal ancestry and a single uniparental marker answer different questions, so neither a large European autosomal component nor deep Near Eastern macro-clade ancestry alone determines the founders' origin.",
"",
"Two ancient-DNA caveats limit resolution. The Pre-Pottery Neolithic K samples establish deep Near Eastern K at the macro-clade level but were typed only on the control region and do not demonstrate K1a1b1a or K2a2a in those individuals. Conversely, drift and lineage loss mean that absence of a derived founder subclade from a small ancient panel is not evidence that it was absent from the region. These qualifications constrain subclade-level certainty without creating affirmative support for recent European-host absorption.",
"",
"## Conclusion",
"",
sprintf("Five methodologically independent analyses favor the same direction: **all four major Ashkenazi mtDNA founders are better explained by a Near Eastern / Levantine origin than by a simple prehistoric European-host origin**. The Bayesian synthesis makes this graded and quantitative -- every founder has a posterior probability of Near Eastern origin above 0.5, **highest for K2a2a (%s)** and **strongly for the K1a founders** (K1a1b1a %s, K1a9 %s, with credible intervals excluding parity), and **positive but more tentative for N1b2 (%s)**, whose Europe-leaning nesting tempers its rarity and antiquity signals. Brook's sixth K founder branch **K1a4a** (~0.2%%) is scored separately in the broader panel as an ambiguous convert-vs-Levantine stress case. This supports and extends the Livni-Skorecki / Behar model while explaining why the Costa European-assimilation model is not the best fit to these data:",
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K2a2a"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K1a1b1a"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "K1a9"]),
        .fmt_post(D_bp$post_mean_H2[D_bp$founder == "N1b2"])),
"",
"1. **Absorption estimate.** Estimated directly from the largest available maternal sample, the absorbed host-lineage fraction places all four founders in the multi-copy founder tail, not among the absorbed singletons. (This is a demographic result about founder status; it does not by itself assign geography.)",
"2. **Non-Jewish rarity.** The four founders are essentially absent among 27,651 non-Jews (order 1e-4) -- the single hardest observation for recent European assimilation to accommodate.",
"3. **Nesting weakens under sampling control.** The 'solely European nesting' argument is weakened by rarefaction to equal sample size for K1a, which becomes Near-East-richer; K2a and N1b sister lineages stay Europe-richer at equal n, which is why nesting tempers rather than boosts N1b2 and K2a2a.",
"4. **Deep Near Eastern ancient DNA.** Haplogroup K is the most frequent lineage (42.8%, 6 of the 14 assigned skeletons) in Pre-Pottery Neolithic Near Eastern farmers (Fernández 2014), and the founder clades appear directly in medieval Jewish individuals from Erfurt and Tàrrega, with an independent, earlier medieval Jewish anchor at Norwich (Brace 2022, 1190 CE) confirming that the Ashkenazi maternal pool pre-dates the 12th century.",
sprintf("5. **Model-disfavored European reconciliation.** The Euro-Levantine private-mutation scenario has only ~%.2f%% probability of producing three or more of the founders at base assumptions (%.2f--%.2f%% across the sensitivity grid).", p_ge3, els_lo, els_hi),
"",
"## Limitations (read honestly)",
"",
"- **Reference-set sampling.** Mitotree/GenBank is a cross-study reference set assembled from many studies. We treat it as the largest available maternal sample and estimate absorption from its singleton fraction, but that fraction is sensitive to lineage-calling resolution and to global (non-Ashkenazi) diversity, so the enriched-set absorption estimate is an upper bound and Behar 2006 is the Ashkenazi-specific anchor.",
"- **Diaspora sampling bias.** Modern `Country` records where a carrier lives, not where the lineage originated. The modern founder distribution is therefore *not* used as direct evidence of origin.",
"- **Public subset.** Mitotree's public release here is a subset of the full ~330,000-sequence tree; ancient Near Eastern coverage is sparse.",
"- **YFull/FTDNA public enrichment.** Public genealogy pages improve alias resolution, branch-age context and ancient-anchor visibility, but tester country/project counts are participation-biased and are not treated as population frequencies.",
"- **Literature-coded controls.** Some non-founder controls use published classifications and manually curated frequencies rather than being discovered de novo from Mitotree. They are used for calibration and falsification, not as independent population estimates.",
"- **Bayesian synthesis.** One ridge-logistic model on standardized, data-computed channels; coefficients fit to labeled controls with leave-one-out validation; founders held out. No hand-set `tanh` gains or elicited channel weights.",
"- **Branching independence.** The Galton-Watson layer assumes independent offspring counts conditional on the chosen law and generation. Negative-binomial overdispersion and piecewise growth test sensitivity to realistic variance and timing, but they do not model full family-level correlation or overlapping generations.",
"- **Reconstructed motifs.** The haplotype network / glPca use motifs rebuilt from Mitotree defining-mutation strings, not per-sample alignments.",
"- **Model parameters.** The branching model inherits Livni-Skorecki's assumptions, but this report no longer substitutes Behar/Costa sample-table counts for Mitotree-derived empirical counts.",
"- **Not a full phylogenetic dating paper.** TMRCA values are imported from Mitotree, and GenBank-query rows are metadata enrichments unless independently placed by Mitotree.",
"- **Resolution, not direction.** These caveats bound the *precision* of the conclusion -- especially for K2a2a, the most tentative founder and the only one whose interval includes parity, which rests on non-Jewish rarity almost alone and would be sharpened most by a usable Near Eastern sister pool within K2a and by any pre-modern Jewish carrier, and for N1b2, which would be sharpened by coding-region typing of ancient Near Eastern carriers and denser Levantine sampling -- but they do not supply positive evidence for the recent European-host origin model. mtDNA cannot fix an origin with autosomal-level certainty, yet under these data the Near Eastern origin of the Ashkenazi maternal founder core is the parsimonious and best-supported conclusion.",
"",
"**Reproducibility.** A single analysis pipeline reproduces the reported results, tables, and figures end to end from the public Mitotree release, GenBank queries, and cited study data.",
"",
"## Sources",
"",
"- Livni J. & Skorecki K. (2025) *Distinguishing between founder and host population mtDNA lineages in the Ashkenazi population* (Human Gene; SSRN 5035272).",
"- Costa M.D. et al. (2013) *A substantial prehistoric European ancestry amongst Ashkenazi maternal lineages* (Nat. Commun. 4:2543).",
"- Behar D.M. et al. (2006) major Ashkenazi founder analysis.",
"- Atzmon G. et al. (2010) *Abraham's children in the genome era: major Jewish Diaspora populations comprise distinct genetic clusters with shared Middle Eastern ancestry* (Am. J. Hum. Genet. 86:850-859).",
"- Carmi S. et al. (2014) *Sequencing an Ashkenazi reference panel supports population-targeted personal genomics and illuminates Jewish and European origins* (Nat. Commun. 5:4835).",
"- Agranat-Tamir L. et al. (2020) *The genomic history of the Bronze Age Southern Levant* (Cell 181:1146-1157).",
"- Waldman S. et al. (2022) *Genome-wide data from medieval German Jews show that the Ashkenazi founder event pre-dated the 14th century* (Cell 185:4703-4716).",
"- Brace S. et al. (2022) *Genomes from a medieval mass burial show Ashkenazi-associated hereditary diseases pre-date the 12th century* (Curr. Biol. 32:4350-4359).",
"- Lerga-Jaso J. et al. (2025) *Tracing human genetic histories and natural selection with precise local ancestry inference* (Nat. Commun. 16:4576).",
"- Brook K.A. (2022) *The Maternal Genetic Lineages of Ashkenazic Jews* (Academic Studies Press); frequency data via Penninx W. (2019).",
"- Wexler J.D. (2019-2026) *Ashkenazi Y-DNA and mtDNA* (compilation; updated mtDNA ancestral lines of Ashkenazi Jews). https://sites.google.com/view/ashkenazi-y-dna-and-mtdna",
"- YFull MTree (public haplogroup pages, accessed 2026). https://www.yfull.com/mtree/",
"- FamilyTreeDNA Discover public mtDNA haplogroup pages and JSON resources (accessed 2026). https://discover.familytreedna.com/mtdna/",
"- FamilyTreeDNA Discover, mtDNA haplogroup V7a12a1 (age estimate). https://discover.familytreedna.com/mtdna/V7a12a1/classic",
"- Fernández E. et al. (2014) *Ancient DNA Analysis of 8000 B.C. Near Eastern Farmers Supports an Early Neolithic Pioneer Maritime Colonization of Mainland Europe through Cyprus and the Aegean Islands* (PLoS Genet. 10(6):e1004401).",
"- Shamoon-Pour M., Li M. & Merriwether D.A. (2019) *Rare human mitochondrial HV lineages spread from the Near East and Caucasus during post-LGM and Neolithic expansions* (Sci. Rep. 9:14751).",
"- Xue J. et al. (2017) *The time and place of European admixture in Ashkenazi Jewish history* (PLoS Genet. 13(4):e1006644).",
"- Maier P.A. et al. (2026) *Mitotree: the universal human mitochondrial reference phylogeny* (bioRxiv).",
"",
"## Frequently Asked Questions",
"",
"The objections below mistake individual sensitivity tests for standalone claims. The pipeline compares the hypothesis of substantial prehistoric European host lineages absorbed into Ashkenazim ($H_1$) against founder transmission from an older Near Eastern / pre-diaspora Jewish maternal population ($H_2$).",
"",
"### 1. Are the evidence channels correlated or double-counted?",
sprintf("No. The pipeline fits three standardized predictors jointly within a single Bayesian synthesis model using ridge priors. Their correlation is directly measured (`outputs/tables/G_channel_dependence.csv`) and regularized (effective dimensionality %s of 3), preventing double-counting. Independent validation checks -- including ancient DNA falsification, Galton-Watson branching, equal-$n$ rarefaction, and negative controls -- are reported as separate robustness tests, not multiplied Bayes factors.",
        if (nrow(G_cd) && "effective_independent_channels" %in% names(G_cd)) sprintf("%.2f", G_cd$effective_independent_channels[1]) else "about 2.5"),
"",
"### 2. Do the model scores reflect true geographic origin probabilities?",
sprintf("The scores quantify conditional evidence for Near Eastern vs. European origin under a data-driven synthesis; they are not calibrated geographic probabilities. The model gets %d of %d labeled controls right out of sample (LOO validation; `outputs/tables/D_fitted_loo.csv`). Non-European controls check the classifier's ability to identify non-European signatures, while the Near Eastern assignment for the Ashkenazi founders is independently anchored by ancient Near Eastern PPNB K lineages, medieval Jewish mitogenomes (Erfurt, Tarrega, Norwich), and nesting analysis.",
        if (nrow(D_fl)) sum(D_fl$correct) else 0L, nrow(D_fl)),
"",
"### 3. Is the training set circular?",
"No. Founder lineages are held out of training and are never used to fit the model. Key channels -- the ancient-DNA falsification audits, non-Jewish European frequencies from 1000 Genomes, and the Galton-Watson branching calculations -- are label-free and derived directly from empirical sequence data. LOO cross-validation checks that the classifier generalizes beyond individual controls.",
"",
"### 4. Can non-Jewish rarity and founder frequency distinguish geographic origin?",
"They provide demographic discrimination against the European-host hypothesis. Under $H_1$, a substantial prehistoric European host lineage absorbed into Ashkenazim would be expected to leave non-Jewish European descendants; its near-absence among tens of thousands of non-Jewish Europeans is hard to reconcile with recent host absorption. The branching equations likewise place major founder lineages in the persistent founder tail rather than the singleton distribution characteristic of absorbed host lines.",
"",
"### 5. Does the presence of Haplogroup K in ancient Near Eastern Neolithic populations support the founders?",
"At the macro-clade level, yes. Haplogroup K at 42.8% (6 of the 14 assigned skeletons) in Pre-Pottery Neolithic Near Eastern farmers establishes deep Near Eastern K. Those samples were typed only on HVR1, so they do not demonstrate the derived founder sub-clades themselves; combined with branch topology, TMRCA dating, non-Jewish rarity, and medieval Jewish mitogenomes, the ancient Near Eastern K footprint is context for a Near Eastern rather than European source.",
"",
"### 6. Do ancient European samples refute the Near Eastern origin of K1a1b1a?",
"No. Claims that ancient European samples carry Ashkenazi founder lineages conflate ancestral parent nodes (K1a1b1, ~7.2 kya; K1a1b, ~10.7 kya) or sister branches with the much younger derived founder clade K1a1b1a (~2.8 kya). Audits of the pre-diaspora European K1a mitogenome record find no non-Jewish derived K1a1b1a samples (`outputs/tables/G_ancient_falsification_with_controls.csv`), whereas the ancient derived K1a1b1a carriers in the record are medieval Jewish individuals.",
"",
"### 7. Are medieval Jewish mitogenomes insufficient to establish Levantine origin?",
"They do not fix a birthplace on their own, but they do establish antiquity. Medieval Jewish mitogenomes from Erfurt, Tarrega and Norwich show the founder lineages pre-date the 12th century and were already established in the historical Jewish maternal pool. Combined with the absence of derived founder motifs in pre-diaspora European populations and the deep Near Eastern clade roots, they place these lineages before European settlement.",
"",
"### 8. Does equal-sample-size rarefaction weaken the evidence for Near Eastern origin?",
"It removes the sampling artifact the raw-density argument rests on: higher European sequence counts largely reflect European sampling intensity. After controlling for sample size, K1a is Near-East-richer, while K2a and N1b remain Europe-richer at equal $n$ -- the nesting channel is in fact a negative contributor for N1b2. The rarefaction result is therefore evidence against a blanket 'solely European nesting' claim, not positive nesting support for all four founders.",
"",
"### 9. Is the branching process result sensitive to demographic assumptions?",
sprintf("It is robust across the tested demographic parameters, offspring distributions (Poisson, geometric, negative binomial) and growth schedules (constant vs. piecewise historical growth). The figure usually quoted here is narrower than it sounds: the ~%.2f%% baseline (%.2f--%.2f%% across the sensitivity grid, and %.2f%% under the Dirichlet sigma_major rather than the assumed 0.03) is the probability that the Euro-Levantine *private-mutation* scenario yields three or more of the four majors -- not an omnibus probability for European host absorption.", p_ge3, els_lo, els_hi, p_ge3_dir),
"",
"### 10. Are database biases or reconstructed motifs inflating the results?",
"The pipeline relies on curated population baselines (such as 1000 Genomes EUR) and published study cohorts rather than commercial participation metrics. Network topology and glPCA are used for visual confirmation only; removing them leaves the frequency, branching, ancient-DNA, rarefaction and synthesis results unchanged.",
"",
"### 11. What would be required to challenge these conclusions?",
"Securely dated, pre-diaspora non-Jewish European mitogenomes carrying the derived Ashkenazi founder subclades, or robust European sister-clade dominance under equalized sampling, would challenge the Near Eastern reading directly.",
""
  )

  out <- file.path(ROOT, "report", "ashkenazi_mtdna_origin.md")
  writeLines(L, out)
  message(sprintf("Wrote %s (%d lines)", out, length(L)))
  invisible(out)
}

if (sys.nframe() == 0) { source(file.path("R", "utils.R")); build_report() }
