# 07_origin_stress_tests.R
# Module G -- additive stress tests that strengthen the origin inference and
# reduce its dependence on the circular parts of the channel synthesis. Nothing
# here modifies modules A-D or their outputs; it reads their tables and the AADR
# ancient-DNA compilation and writes new, independent evidence.
#
#   G1 (Path 1) Ancient-DNA falsification of the Costa prehistoric-European model.
#               Costa (2013) predicts the founder sub-clades are prehistoric
#               European maternal lineages; if so they should occur among
#               pre-diaspora ancient Europeans. Tested against the AADR v66
#               ancient mtDNA record (Mallick et al. 2024).
#   G2 (Path 2) European host-population frequency bound and likelihood ratio.
#               Turns the rarity observation into a non-circular falsification of
#               the European-host origin: the pooled zero-count European sample
#               bounds any host frequency, and the founder-effect enrichment that
#               Costa's model must then invoke is quantified.
#   G5 (Path 5) De-risking: adversarial (Costa-favouring) prior sensitivity,
#               channel-dependence accounting, and two-sided ancient controls.
#
# (A phylogenetic ancestral-region reconstruction was evaluated and dropped: on
# the modern reference tree it is driven by European over-sampling -- the same
# ascertainment bias the equal-sample-size rarefaction in module B exists to
# correct -- so it cannot serve as independent evidence here.)

if (!exists("load_analysis_samples")) source(file.path("R", "utils.R"))

# Founder -> AADR/phylotree sub-clade name and macro-clade, with Ashkenazi
# frequency (Brook/Penninx) and Livni-Skorecki Table 7 non-Jewish counts.
FOUNDER_STRESS <- data.frame(
  founder  = c("K1a1b1a", "K1a9", "K2a2a", "N1b2"),
  aadr_sub = c("K1a1b1a", "K1a9", "K2a2a", "N1b1b1"),
  macro    = c("K1a", "K1a", "K2a", "N1b"),
  ashk_pct = c(20.70, 6.10, 5.30, 6.08),
  t7_count = c(17, 6, 12, 7),
  t7_n     = 27651L,
  stringsAsFactors = FALSE
)

# European calibration lineages used as two-sided ancient controls: established
# Ashkenazi lineages of recognised European origin. Under Costa they SHOULD have
# pre-diaspora European ancient carriers -- the method must move both ways.
EU_CONTROL_STRESS <- data.frame(
  founder  = c("V (V7a2c1b)", "U5a1 (U5a1f1a3)"),
  aadr_sub = c("V", "U5a1"),
  macro    = c("HV0;V", "U5"),
  # Ancestral labels sharing no prefix with the sub-clade (see .resolves_subclade).
  anc_extra = c("HV0;HV", ""),
  stringsAsFactors = FALSE
)

# Calendar year below which an ancient European carrier would count as
# "pre-diaspora" evidence for a prehistoric-European origin. 500 CE is generous
# to Costa: it admits the entire prehistoric plus Roman record and excludes only
# the established medieval Jewish communities.
PREDIASPORA_YEAR <- 500L
NE_REGIONS_G <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus")

# --------------------------------------------------------------------------
# Shared helpers.
# --------------------------------------------------------------------------
# Clade membership respecting mtDNA name boundaries (a digit right after a
# digit-ending clade denotes a sibling, not a descendant: H1 !-> H10).
.in_clade <- function(hg, clade) {
  hg <- as.character(hg)
  ok <- !is.na(hg) & (hg == clade | startsWith(hg, clade))
  tail_digit <- substr(hg, nchar(clade) + 1L, nchar(clade) + 1L)
  clade_digit <- grepl("[0-9]$", clade)
  ok & !(hg != clade & clade_digit & grepl("^[0-9]$", tail_digit))
}

# Some clades are written in AADR under labels that do not share a prefix with
# the parent clade name, so prefix matching alone would make a sub-clade fall
# outside its own macro pool (V descends from HV0 but is never written "HV0*").
# A macro spec may therefore list several labels separated by ";".
.in_clade_spec <- function(hg, spec) {
  parts <- trimws(strsplit(as.character(spec), ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(rep(FALSE, length(hg)))
  Reduce(`|`, lapply(parts, function(p) .in_clade(hg, p)))
}

# Is a record's haplogroup call resolved deeply enough to reveal sub-clade
# `target` had the individual carried it?
#
# Three cases, for target = K1a1b1a:
#   "K1a1b1a", "K1a1b1a1", "K1a1b1a+16362"  -> a hit, informative
#   "K1a9", "K2a2a", "H1", "U5"             -> diverges at a resolved position,
#                                              a genuine informative non-hit
#   "K", "K1a", "K1a1b1", "K1a+195"         -> ANCESTRAL to the target and so
#                                              unresolved: the individual could
#                                              have been K1a1b1a and the call
#                                              would look identical
# Only the first two belong in a denominator bounding the target's frequency.
# Counting the third group as zero-observations would claim information that
# low-resolution calls do not carry.
.resolves_subclade <- function(hg, target, anc_extra = character(0)) {
  hg <- as.character(hg)
  # Strip extra-mutation / back-mutation / quote-notation suffixes so that
  # "K1a+195" is judged at its K1a resolution.
  base <- sub("[+@'].*$", "", hg)
  hit <- .in_clade(hg, target)
  # `.in_clade` is written for a vector `hg` against a SCALAR `clade`. Calling
  # it with a scalar target against a vector of bases silently collapses its
  # digit-boundary guard (substr returns one element), leaving plain
  # startsWith, so the ancestry test is evaluated per element instead.
  ancestral <- base != target &
    vapply(base, function(b) isTRUE(.in_clade(target, b)), logical(1),
           USE.NAMES = FALSE)
  # Some ancestors share no prefix with the target: HV0 and HV are ancestral to
  # V but are never written "V*", so prefix matching alone would score a record
  # called HV0 as an informative non-hit when a V carrier could look exactly
  # like that. Exact membership only -- a label resolved BELOW one of these
  # (HV0b, HV1a) sits on a different branch and does exclude the target.
  ancestral <- ancestral | (base %in% anc_extra & !hit)
  out <- hit | !ancestral
  out[is.na(hg) | !nzchar(base)] <- FALSE
  out
}

# Self-contained reader for a previously written outputs table (module G reads
# module D's outputs but must not depend on module D being sourced).
.g_read_tab <- function(name) {
  p <- file.path(TAB_DIR, name)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else data.frame()
}

load_aadr_mtdna <- function() {
  p <- file.path(ROOT, "data", "external", "aadr_mtdna.csv")
  if (!file.exists(p)) return(data.frame())
  d <- read.csv(p, stringsAsFactors = FALSE)
  d$date_bp <- suppressWarnings(as.numeric(d$date_bp))
  d$year    <- suppressWarnings(as.numeric(d$year))
  d$region  <- as.character(region_of(d$country))
  d$is_jewish <- as.integer(d$is_jewish) %in% 1L
  # Ancient = has a real pre-present date (drop present-day reference genomes).
  d$ancient <- is.finite(d$date_bp) & d$date_bp > 50
  d
}

# 95% upper confidence bound on a binomial proportion given k successes in n
# (Clopper-Pearson). For k = 0 this is the rule-of-three-style 1 - 0.05^(1/n).
.binom_upper95 <- function(k, n) {
  if (!is.finite(n) || n <= 0) return(NA_real_)
  if (k > n) stop("binomial upper bound: k (", k, ") exceeds n (", n, ")")
  if (k <= 0) return(1 - 0.05^(1 / n))        # = qbeta(0.95, 1, n)
  if (k >= n) return(1)                       # all observed: upper bound is 1
  # One-sided 95% level, matching the k = 0 branch and the documented level.
  stats::qbeta(0.95, k + 1, n - k)
}

# --------------------------------------------------------------------------
# G1 -- ancient-DNA falsification of the prehistoric-European (Costa) model.
# --------------------------------------------------------------------------
ancient_carrier_detail <- function(info = FOUNDER_STRESS) {
  aadr <- load_aadr_mtdna()
  if (!nrow(aadr)) return(data.frame())
  anc <- aadr[aadr$ancient, ]
  do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
    sub <- anc[.in_clade(anc$mt_haplogroup, info$aadr_sub[i]), ]
    if (!nrow(sub)) return(NULL)
    data.frame(
      founder = info$founder[i], aadr_sub = info$aadr_sub[i],
      id = sub$id, mt = sub$mt_haplogroup, year = sub$year,
      country = sub$country, region = sub$region,
      is_jewish = as.integer(sub$is_jewish), group = sub$group_id,
      stringsAsFactors = FALSE
    )
  }))
}

ancient_falsification <- function(info = FOUNDER_STRESS, label = "founder") {
  aadr <- load_aadr_mtdna()
  if (!nrow(aadr)) return(data.frame())
  anc <- aadr[aadr$ancient, ]
  do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
    macro <- anc[.in_clade_spec(anc$mt_haplogroup, info$macro[i]), ]
    sub   <- anc[.in_clade(anc$mt_haplogroup, info$aadr_sub[i]), ]
    eu <- function(d) d[d$region == "Europe" & !d$is_jewish, ]
    ne <- function(d) d[d$region %in% NE_REGIONS_G & !d$is_jewish, ]
    predia <- function(d) d[is.finite(d$year) & d$year < PREDIASPORA_YEAR, ]
    n_macro_eu_predia <- nrow(predia(eu(macro)))
    sub_eu_predia     <- predia(eu(sub))
    n_sub_eu_predia   <- nrow(sub_eu_predia)
    sub_jewish        <- sub[sub$is_jewish, ]
    # Ascertainment set for an UNCONDITIONAL frequency bound: every
    # pre-diaspora non-Jewish ancient European record whose call is resolved
    # deeply enough to have revealed this sub-clade. Records called only at an
    # ancestral level (K, K1a, K1a1b1) carry no information about it and are
    # excluded; records in other clades are informative non-hits and are kept.
    eu_predia_all <- predia(eu(anc))
    n_eu_predia_resolvable <- sum(
      .resolves_subclade(eu_predia_all$mt_haplogroup, info$aadr_sub[i],
                         if ("anc_extra" %in% names(info) && nzchar(info$anc_extra[i]))
                           trimws(strsplit(info$anc_extra[i], ";", fixed = TRUE)[[1]])
                         else character(0))
    )
    # 95% upper bound on the founder's share of the pre-diaspora European macro
    # pool (the frequency Costa's prehistoric-European lineage could have had and
    # still escaped detection).
    host_share_hi <- .binom_upper95(n_sub_eu_predia, n_macro_eu_predia)
    data.frame(
      kind = label, founder = info$founder[i], aadr_sub = info$aadr_sub[i],
      macro = info$macro[i],
      n_sub_ancient = nrow(sub),
      n_sub_jewish = nrow(sub_jewish),
      pct_sub_jewish = if (nrow(sub)) round(100 * nrow(sub_jewish) / nrow(sub), 1) else NA_real_,
      n_macro_eu_predia = n_macro_eu_predia,
      n_eu_predia_resolvable = n_eu_predia_resolvable,
      n_sub_eu_predia_nonjewish = n_sub_eu_predia,
      n_sub_ne_ancient = nrow(ne(sub)),
      oldest_sub_eu_nonjewish_year = if (n_sub_eu_predia) min(sub_eu_predia$year) else NA_real_,
      oldest_sub_ne_year = if (nrow(ne(sub))) min(ne(sub)$year, na.rm = TRUE) else NA_real_,
      host_share_upper95 = host_share_hi,
      # p-value that a lineage as common in the European host as in Ashkenazim
      # (Costa's model with no drift) would go unseen: (1 - f_ashk)^n_macro.
      # (1 - f)^n with f a POPULATION frequency, so n must be the
      # unconditional resolution-sufficient sample, not the macro-clade count.
      p_costa_no_drift = if ("ashk_pct" %in% names(info) && n_eu_predia_resolvable > 0)
        (1 - info$ashk_pct[i] / 100)^n_eu_predia_resolvable else NA_real_,
      # The probability underflows to exactly 0 in double precision at these n
      # (K1a1b1a is 10^-760), so report the log as well as the raw value.
      log10_p_costa_no_drift = if ("ashk_pct" %in% names(info) && n_eu_predia_resolvable > 0)
        n_eu_predia_resolvable * log10(1 - info$ashk_pct[i] / 100) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

# --------------------------------------------------------------------------
# G2 -- European host-frequency bound and host-model likelihood ratio.
# --------------------------------------------------------------------------
# Pools the zero-count European non-Jewish observations (1000G EUR + pre-diaspora
# ancient European macro carriers) to bound any European host frequency, then
# quantifies the founder-effect enrichment Costa's European-host model must
# invoke and a Bayes factor for a Near-Eastern vs European-host origin.
host_model_likelihood <- function(info = FOUNDER_STRESS) {
  fals <- ancient_falsification(info)
  eur_n <- tryCatch({
    kg <- load_1kg_clade_freq(); if (nrow(kg)) kg$n_nonjews[1] else 503L
  }, error = function(e) 503L)
  do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
    fr <- fals[fals$founder == info$founder[i], ]
    # Ancient contribution: pre-diaspora non-Jewish ancient Europeans whose call
    # could have revealed the sub-clade. Previously this was the macro-clade
    # count, which is a WITHIN-macro denominator and so not on the same scale as
    # the unconditional 1000G panel it is pooled with.
    anc_eu_n <- if (nrow(fr)) fr$n_eu_predia_resolvable else 0L
    k_eu <- if (nrow(fr)) fr$n_sub_eu_predia_nonjewish else 0L
    # Reported separately as well as pooled: pooling a present-day panel with a
    # record spanning millennia assumes the frequency was constant over that
    # span, which is an assumption, not a measurement.
    f_hi_modern  <- .binom_upper95(0L, eur_n)
    f_hi_ancient <- .binom_upper95(k_eu, anc_eu_n)
    n_eu <- eur_n + anc_eu_n
    f_hi_eu <- .binom_upper95(k_eu, n_eu)
    f_ashk <- info$ashk_pct[i] / 100
    # Broad non-Jewish benchmark (Livni-Skorecki Table 7; blends NE non-Jews).
    f_t7 <- info$t7_count[i] / info$t7_n[i]
    # Enrichment the founder-effect must supply under a European host origin.
    enrichment <- f_ashk / f_hi_eu
    # Bayes factor (Near-Eastern vs European-host) for the European zero count.
    # H2 (NE origin): expected European host frequency ~ 0 (use f_t7 as a small
    # incidental-presence rate). H1 (European host): the lineage was a
    # substantial prehistoric European maternal lineage; use a conservative
    # host-frequency prior anchored at 1% (Costa: "substantial" ancestry).
    f_h1 <- 0.01
    loglik <- function(f) k_eu * log(f + 1e-12) + (n_eu - k_eu) * log(1 - f + 1e-12)
    log10_bf <- (loglik(f_t7) - loglik(f_h1)) / log(10)
    data.frame(
      founder = info$founder[i],
      eur_panel_n = eur_n, ancient_eu_predia_n = anc_eu_n,
      european_nonjewish_n = n_eu, european_nonjewish_obs = k_eu,
      european_host_freq_upper95 = f_hi_eu,
      host_freq_upper95_modern_1kg = f_hi_modern,
      host_freq_upper95_ancient = f_hi_ancient,
      ashkenazi_freq = f_ashk,
      broad_nonjewish_freq_t7 = f_t7,
      founder_effect_enrichment_needed = enrichment,
      log10_BF_NEorigin_vs_Europeanhost = log10_bf,
      stringsAsFactors = FALSE
    )
  }))
}

# --------------------------------------------------------------------------
# G5 -- de-risking: adversarial priors and channel dependence.
# --------------------------------------------------------------------------
# How strong a Costa-favouring prior is needed to overturn each founder call.
adversarial_prior_sensitivity <- function() {
  post <- .g_read_tab("D_fitted_posterior.csv")
  if (!nrow(post)) post <- .g_read_tab("D_bayes_posterior.csv")
  fdf <- if ("group" %in% names(post)) post[post$group == "founder", ] else post
  id <- if ("lineage" %in% names(fdf)) fdf$lineage else fdf$founder
  p0 <- fdf$post_mean_H2
  logit <- function(p) log(pmin(pmax(p, 1e-9), 1 - 1e-9) /
                             (1 - pmin(pmax(p, 1e-9), 1 - 1e-9)))
  eta <- logit(p0)
  # Adversarial prior expressed as European:NE prior odds; delta = -log(odds).
  # Applied to the reported posterior mean (the only quantity in the summary
  # table), so these are shifted means, not the means of the shifted model.
  at <- function(odds) plogis(eta - log(odds))
  data.frame(
    # The baseline is the fitted model's posterior mean. The logistic has a
    # free, fitted intercept, so it is NOT a flat-prior (0.5-centred) baseline.
    founder = id, posterior_fitted = round(p0, 3),
    post_prior_2to1_european = round(at(2), 3),
    post_prior_4to1_european = round(at(4), 3),
    post_prior_9to1_european = round(at(9), 3),
    # European:NE prior odds at which the posterior crosses 0.5.
    tipping_prior_odds_european = round(exp(eta), 1),
    stringsAsFactors = FALSE
  )
}

# The three channels are not independent votes; report their dependence.
channel_dependence <- function() {
  inp <- .g_read_tab("D_fitted_inputs.csv")
  cols <- c("z_freq", "z_time", "z_nest")
  if (!nrow(inp) || !all(cols %in% names(inp))) return(data.frame())
  Z <- as.matrix(inp[, cols])
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  C <- stats::cor(Z)
  ev <- eigen(C, symmetric = TRUE)$values
  eff_dim <- sum(ev)^2 / sum(ev^2)   # participation ratio (effective # channels)
  data.frame(
    cor_freq_time = C["z_freq", "z_time"],
    cor_freq_nest = C["z_freq", "z_nest"],
    cor_time_nest = C["z_time", "z_nest"],
    effective_independent_channels = round(eff_dim, 2),
    nominal_channels = length(cols),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# Runner.
# --------------------------------------------------------------------------
run_origin_stress_tests <- function() {
  message("== Module G: ancient-DNA falsification, host likelihood, de-risking ==")
  aadr <- load_aadr_mtdna()
  if (!nrow(aadr)) {
    message("  AADR ancient mtDNA not available (data/external/aadr_mtdna.csv); ",
            "run python3 scripts/04_extract_aadr.py. Skipping module G.")
    return(invisible(NULL))
  }
  message(sprintf("  AADR mtDNA records: %d (%d ancient)",
                  nrow(aadr), sum(aadr$ancient)))

  # ---- G1 ancient falsification (founders + two-sided European controls) ----
  detail <- ancient_carrier_detail()
  save_table(detail, "G_ancient_carrier_detail.csv")
  fals <- ancient_falsification(FOUNDER_STRESS, "founder")
  fals_ctrl <- ancient_falsification(EU_CONTROL_STRESS, "eu_control")
  fals_all <- rbind(fals[, intersect(names(fals), names(fals_ctrl))],
                    fals_ctrl[, intersect(names(fals), names(fals_ctrl))])
  save_table(fals, "G_ancient_falsification.csv")
  save_table(fals_all, "G_ancient_falsification_with_controls.csv")
  for (i in seq_len(nrow(fals))) {
    message(sprintf(
      "  %-8s ancient sub-clade carriers=%d (%.0f%% Jewish); pre-diaspora non-Jewish European=%d of %d macro; host share <%.3f%% (95%%)",
      fals$founder[i], fals$n_sub_ancient[i], fals$pct_sub_jewish[i],
      fals$n_sub_eu_predia_nonjewish[i], fals$n_macro_eu_predia[i],
      100 * fals$host_share_upper95[i]))
  }

  # ---- G2 host-model likelihood / enrichment ----
  hml <- host_model_likelihood()
  save_table(hml, "G_host_likelihood.csv")
  for (i in seq_len(nrow(hml))) {
    message(sprintf(
      "  %-8s European host freq <%.3f%% (95%%); Ashkenazi %.1f%% requires %.0fx founder-effect enrichment; log10 BF(NE/EU-host)=%.1f",
      hml$founder[i], 100 * hml$european_host_freq_upper95[i],
      100 * hml$ashkenazi_freq[i], hml$founder_effect_enrichment_needed[i],
      hml$log10_BF_NEorigin_vs_Europeanhost[i]))
  }

  # ---- G5 de-risking ----
  adv <- adversarial_prior_sensitivity()
  save_table(adv, "G_adversarial_prior.csv")
  print(adv)
  chd <- channel_dependence()
  if (nrow(chd)) { save_table(chd, "G_channel_dependence.csv"); print(chd) }

  # ---- Figures ----
  .stress_figures(detail, fals, fals_ctrl, hml, adv)

  invisible(list(detail = detail, falsification = fals, host = hml,
                 adversarial = adv, channels = chd))
}

.stress_figures <- function(detail, fals, fals_ctrl, hml, adv) {
  # G1: ancient carriers of each founder sub-clade through time, by context.
  if (nrow(detail)) {
    open_png("G1_ancient_founder_timeline.png", width = 1250, height = 820)
    op <- par(mar = c(4.5, 7, 3.2, 1))
    fl <- unique(FOUNDER_STRESS$founder)
    detail$fy <- match(detail$founder, fl)
    col <- ifelse(detail$is_jewish == 1, "#c1121f",
           ifelse(detail$region == "Europe", "#264653",
           ifelse(detail$region %in% NE_REGIONS_G, "#f3a712", "#adb5bd")))
    plot(detail$year, jitter(detail$fy, 0.5), pch = 19, col = col, cex = 1.4,
         yaxt = "n", ylab = "", xlab = "calendar year (negative = BCE)",
         ylim = c(0.5, length(fl) + 0.5),
         main = "Ancient carriers of the Ashkenazi founder sub-clades (AADR v66)")
    axis(2, at = seq_along(fl), labels = fl, las = 1)
    abline(v = c(0, PREDIASPORA_YEAR), lty = c(3, 2), col = c("grey60", "#c1121f"))
    text(PREDIASPORA_YEAR, length(fl) + 0.35, "500 CE", col = "#c1121f", cex = 0.8, pos = 4)
    legend("topleft", bty = "n", pch = 19,
           col = c("#c1121f", "#264653", "#f3a712", "#adb5bd"),
           legend = c("Jewish", "Europe (non-Jewish)", "Near East", "other"))
    par(op); dev.off()
  }
  # G2: European host-frequency upper bound vs Ashkenazi frequency.
  if (nrow(hml)) {
    open_png("G2_host_frequency_bound.png", width = 1150, height = 760)
    op <- par(mar = c(4.5, 7, 3.2, 1)); y <- seq_len(nrow(hml))
    xr <- range(c(hml$european_host_freq_upper95, hml$ashkenazi_freq)) * 100
    plot(NA, xlim = c(min(xr) * 0.5, max(xr) * 1.3), ylim = c(0.5, nrow(hml) + 0.5),
         log = "x", yaxt = "n", ylab = "", xlab = "frequency (%) [log scale]",
         main = "European host-frequency bound vs Ashkenazi frequency")
    axis(2, at = y, labels = hml$founder, las = 1)
    points(hml$european_host_freq_upper95 * 100, y, pch = "|", cex = 1.6, col = "#264653")
    points(hml$ashkenazi_freq * 100, y, pch = 19, cex = 1.6, col = "#c1121f")
    segments(hml$european_host_freq_upper95 * 100, y, hml$ashkenazi_freq * 100, y,
             col = "grey60", lty = 3)
    legend("bottomright", bty = "n", pch = c(124, 19), col = c("#264653", "#c1121f"),
           legend = c("European host freq upper 95% (0 observed)", "Ashkenazi frequency"))
    par(op); dev.off()
  }
  # G3 figure omitted (values reported in table); G5 adversarial prior curve.
  if (nrow(adv)) {
    open_png("G5_adversarial_prior.png", width = 1150, height = 760)
    op <- par(mar = c(4.5, 4.8, 3.2, 1))
    odds <- 10^seq(-1, 2.2, length.out = 60)   # NE-favouring < 1 < European-favouring
    plot(NA, xlim = range(odds), ylim = c(0, 1), log = "x",
         xlab = "prior odds European : Near Eastern",
         ylab = "posterior P(Near Eastern origin)",
         main = "Robustness to a Costa-favouring prior")
    cols <- c("#c1121f", "#e5793a", "#264653", "#2a9d8f")
    logit <- function(p) log(p / (1 - p))
    for (i in seq_len(nrow(adv))) {
      e <- logit(pmin(pmax(adv$posterior_fitted[i], 1e-6), 1 - 1e-6))
      lines(odds, plogis(e - log(odds)), col = cols[(i - 1) %% 4 + 1], lwd = 2)
    }
    abline(h = 0.5, lty = 2, col = "grey50"); abline(v = 1, lty = 3, col = "grey60")
    legend("bottomleft", bty = "n", lwd = 2, col = cols[seq_len(nrow(adv))],
           legend = adv$founder)
    par(op); dev.off()
  }
}

if (sys.nframe() == 0) { source(file.path("R", "utils.R")); run_origin_stress_tests() }
