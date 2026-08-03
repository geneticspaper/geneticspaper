# 06_advanced_models.R
# Module D -- realism upgrades to the Livni-Skorecki branching framework and a
# Bayesian, time-aware synthesis of the origin evidence.
#
# This module implements the three extensions that the base report previously
# only flagged as "further work":
#
#   SIMPLE  (D1) Intermediate negative-binomial offspring sizes. The base model
#               contrasts only Poisson (size -> Inf) and geometric (size = 1)
#               reproduction. Here we sweep intermediate dispersions so the
#               reader can see the whole family, not just the two extremes.
#
#   COMPLEX (D2) Piecewise historical growth. Instead of a single constant
#               growth ratio m, we let m_t vary across generations to mimic a
#               near-stationary medieval founder phase followed by rapid modern
#               expansion, and show that the *timing* of growth changes lineage
#               survival even when the net expansion is held fixed.
#
#   COMPLEX (D3) A Bayesian, time-aware origin model. For each founder we
#               combine founder frequency, ancient occurrence dates, public
#               YFull/FTDNA branch-age context, regional sampling intensity, and
#               phylogenetic nesting into a posterior probability of a Near
#               Eastern / Levantine origin, with credible intervals from Monte
#               Carlo over elicited priors. This is an explicitly *subjective*
#               Bayesian evidence synthesis (priors on channel weights), not a
#               likelihood fit to a random population sample, and is reported as
#               such.

if (!exists("gw_descendant_pmf")) source(file.path("R", "utils.R"))
if (!exists("compute_lineage_inputs")) source(file.path("R", "lineage_data.R"))

.adv_read_tab <- function(name) {
  p <- file.path(TAB_DIR, name)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else data.frame()
}

# Founder -> macro (super) clade used by the ancient / nesting / rarefaction
# tables (K1a for the two K1a founders, K2a for K2a2a, N1b for N1b2).
FOUNDER_ROOT <- c(K1a1b1a = "K1a", K1a9 = "K1a", K2a2a = "K2a", N1b2 = "N1b")

# Public YFull/FTDNA enrichment is useful for branch age, alias and ancient-anchor
# context, but not as a population-frequency source. These helpers convert the
# cached public-page summary into modest time-channel adjustments.
parse_ybp <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  out[!nzchar(x)] <- NA_real_
  out
}

public_haplogroup_context <- function(lineage_key) {
  empty <- data.frame(
    public_tmrca_kyr = NA_real_, public_age_confidence = 0,
    public_jewish_ancient_anchor = 0L, public_yfull_ids = NA_real_,
    public_ftdna_placements = NA_real_, stringsAsFactors = FALSE
  )
  pub <- tryCatch(load_public_haplogroup_pages(), error = function(e) data.frame())
  if (!nrow(pub)) return(empty)
  row <- pub[pub$lineage_key == lineage_key, ]
  if (!nrow(row)) return(empty)
  row <- row[1, ]

  yfull_kyr <- parse_ybp(row$yfull_tmrca) / 1000
  ftdna_year <- suppressWarnings(as.numeric(row$ftdna_tmrca_mean))
  ftdna_kyr <- if (is.na(ftdna_year)) NA_real_ else (1950 - ftdna_year) / 1000
  tmrca <- stats::median(c(yfull_kyr, ftdna_kyr), na.rm = TRUE)
  if (!is.finite(tmrca)) tmrca <- NA_real_

  yfull_ids <- suppressWarnings(as.numeric(row$yfull_public_ids_n))
  placements <- suppressWarnings(as.numeric(row$ftdna_placements_n))
  support_n <- sum(c(yfull_ids, placements), na.rm = TRUE)
  conf <- if (support_n > 0) min(1, log10(1 + support_n) / 2.5) else 0

  # Keep this narrow: Jewish/medieval anchors improve time context, while broad
  # ancient relatives from unrelated prehistoric contexts do not automatically
  # count as evidence for H2.
  # Match only unambiguous Jewish ancient studies/sites in the FTDNA ancient list.
  # Bare Reich-lab sample-ID prefixes (I138.../I148...) are NOT used: they match
  # generic non-Jewish ancient samples (e.g. HV1a's I148xx) and would spuriously
  # anchor European lineages.
  anchor_text <- paste(row$ftdna_ancient_codes, row$ftdna_ancient_studies, sep = ";")
  jewish_anchor <- grepl("Waldman|Pallar|Tarrega|Roquetes|ROQ|Sobibor|Diepenbroek|Erfurt|Chapelfield",
                         anchor_text, ignore.case = TRUE)

  data.frame(
    public_tmrca_kyr = tmrca,
    public_age_confidence = conf,
    public_jewish_ancient_anchor = as.integer(jewish_anchor),
    public_yfull_ids = yfull_ids,
    public_ftdna_placements = placements,
    stringsAsFactors = FALSE
  )
}

public_time_adjustment <- function(public_tmrca_kyr, public_age_confidence,
                                   public_jewish_ancient_anchor,
                                   medieval_jewish_carriers) {
  if (is.na(public_tmrca_kyr) || public_age_confidence <= 0) return(0)
  recent_penalty <- if (public_tmrca_kyr < 1) -0.25 * public_age_confidence else 0
  anchored <- medieval_jewish_carriers > 0 || public_jewish_ancient_anchor > 0
  mature_anchor_bonus <- if (anchored && public_tmrca_kyr >= 1.5) 0.12 * public_age_confidence else 0
  recent_penalty + mature_anchor_bonus
}

# --------------------------------------------------------------------------
# D1. SIMPLE -- intermediate negative-binomial offspring sizes.
# --------------------------------------------------------------------------
# Poisson is the size -> Inf limit and geometric is size = 1. We sweep a ladder
# of dispersions to show the single-descendant and extinction probabilities
# move smoothly between the two extremes already reported.
nbinom_size_sensitivity <- function(k = 15L,
                                     growth = c(1.0, 1.025, 1.05, 1.075, 1.1),
                                     sizes = c(0.5, 1, 2, 5, 20)) {
  rows <- list()
  add <- function(label, model, size, m) {
    pmf <- gw_descendant_pmf(m, k = k, model = model, size = size)
    data.frame(
      offspring_law = label,
      nb_size = if (is.finite(size)) size else Inf,
      growth_rate = m,
      var_mean_ratio = if (model == "poisson") 1 else 1 + m / size,
      prob_extinct = pmf[1L],
      prob_one_descendant = pmf[2L],
      prob_ge_three_descendants = sum(pmf[4:length(pmf)]),
      stringsAsFactors = FALSE
    )
  }
  for (m in growth) {
    for (s in sizes) rows[[length(rows) + 1L]] <- add(sprintf("NB(size=%g)", s), "nbinom", s, m)
    rows[[length(rows) + 1L]] <- add("Poisson (size=Inf)", "poisson", Inf, m)
  }
  out <- do.call(rbind, rows)
  out[order(out$growth_rate, out$nb_size), ]
}

# --------------------------------------------------------------------------
# D2. COMPLEX -- piecewise (time-varying) historical growth.
# --------------------------------------------------------------------------
# A single constant growth ratio is demographically unrealistic: the Ashkenazi
# matriline was near-stationary through the medieval founder phase and then
# expanded rapidly. We encode a piecewise schedule and compare it against a
# constant-growth process with the *same* net expansion (same E[Z_k]), so any
# difference is attributable purely to the timing of growth.
default_growth_schedule <- function(k = 15L) {
  # Three phases over k generations: founder bottleneck (near-stationary),
  # steady growth, then rapid modern expansion.
  phase <- cut(seq_len(k), breaks = c(0, round(k / 3), round(2 * k / 3), k),
               labels = c("bottleneck", "steady", "expansion"))
  m_t <- c(bottleneck = 0.98, steady = 1.10, expansion = 1.28)[as.character(phase)]
  data.frame(generation = seq_len(k), phase = as.character(phase),
             growth_ratio = as.numeric(m_t), stringsAsFactors = FALSE)
}

piecewise_growth_model <- function(k = 15L, schedule = default_growth_schedule(k),
                                   model = "poisson") {
  m_t <- schedule$growth_ratio
  # Constant-growth comparator with identical net expansion prod(m_t).
  m_const <- prod(m_t)^(1 / k)
  summ <- function(label, m) {
    pmf <- gw_descendant_pmf(m, k = k, model = model)
    data.frame(
      model = label,
      net_expansion = prod(if (length(m) == 1) rep(m, k) else m),
      mean_growth_ratio = if (length(m) == 1) m else exp(mean(log(m))),
      prob_extinct = pmf[1L],
      prob_one_descendant = pmf[2L],
      prob_ge_three_descendants = sum(pmf[4:length(pmf)]),
      stringsAsFactors = FALSE
    )
  }
  comparison <- rbind(
    summ("constant (matched net growth)", m_const),
    summ("piecewise (bottleneck->expansion)", m_t)
  )
  list(schedule = schedule, m_const = m_const, comparison = comparison)
}

# --------------------------------------------------------------------------
# D3. COMPLEX -- unified origin synthesis (single fitted logistic model).
# --------------------------------------------------------------------------
# Raw channels (non-Jewish rarity, antiquity/carriers/public age, equal-n nesting)
# are computed from enriched Mitotree data; Ashkenazi lineage frequencies come from
# Brook/Wexler (Penninx 2019 basis). Channels are z-scored; a ridge-logistic model
# is fit to labeled controls with founders held out.

# Assemble per-founder inputs from the enriched Mitotree sample set (no hand-set
# CSV numerics; medieval carriers are counted from ancient records only).
bayesian_origin_inputs <- function(founders = names(FOUNDER_ROOT)) {
  comp <- compute_lineage_inputs_batch(founders, nesting_modes = "ne_eu",
                                       rarefy_roots = founders)
  do.call(rbind, lapply(seq_len(nrow(comp)), function(i) {
    f <- comp$lineage_key[i]
    pc <- public_haplogroup_context(f)
    data.frame(
      founder = f, root = comp$macro_root[i],
      ref_freq = comp$ref_freq[i],
      nonjew_freq = comp$nonjew_freq[i],
      oldest_ne_year = if (is.na(comp$ne_depth_kyr[i])) NA_real_
        else 1950 - comp$ne_depth_kyr[i] * 1000,
      ne_depth_kyr = comp$ne_depth_kyr[i],
      medieval_jewish_carriers = comp$medieval_jewish_carriers[i],
      ne_richness_equal_n = comp$ne_richness[i],
      eu_richness_equal_n = comp$eu_richness[i],
      ne_pool = comp$ne_pool[i], eu_pool = comp$eu_pool[i],
      mitotree_node = comp$mitotree_node[i],
      rarefy_root = comp$rarefy_root[i],
      public_tmrca_kyr = pc$public_tmrca_kyr,
      public_age_confidence = pc$public_age_confidence,
      public_jewish_ancient_anchor = pc$public_jewish_ancient_anchor,
      public_yfull_ids = pc$public_yfull_ids,
      public_ftdna_placements = pc$public_ftdna_placements,
      stringsAsFactors = FALSE
    )
  }))
}

# --------------------------------------------------------------------------
# D3. Origin synthesis -- standardized data channels + fitted logistic model.
# --------------------------------------------------------------------------
# All numeric inputs come from lineage_data.R. Channels are z-scored, then a
# ridge-logistic model is fit to labeled controls (European = 0, Near Eastern /
# non-European diaspora = 1). Founders are held out. No hand-set tanh gains or
# elicited channel weights.

.zscore <- function(x) .zscore_by(x, x)

# Standardize `x` using the location/scale of a reference vector `ref`. Splitting
# the reference out (rather than always self-standardizing) lets leave-one-out
# recompute the channel with the held-out row excluded from mean/sd.
.zscore_by <- function(x, ref) {
  mu <- mean(ref, na.rm = TRUE); s <- stats::sd(ref, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  z <- (x - mu) / s; z[is.na(z)] <- 0; z
}

# Rank-based (van der Waerden) normal scores for the non-Jewish rarity channel.
# The founder rarities come from Livni-Skorecki Table 7 (n=27,651) while the
# controls/panel come from 1000 Genomes EUR (n=503) and Mitchell/Mitotree
# fallbacks -- very different sampling frames and clade resolutions. A raw
# z-score of log-rarity lets those absolute-scale differences distort the
# channel; ranking each lineage within the scored set and mapping to normal
# quantiles keeps only the comparable information (relative rarity ordering:
# founders rarer than European controls) and is invariant to the scale/precision
# of the instrument that produced each frequency. `ref` is the ranking pool so
# leave-one-out can exclude the held-out row.
.rank_normal_by <- function(x, ref) {
  refok <- ref[!is.na(ref)]
  n <- length(refok)
  if (n < 2L) return(rep(0, length(x)))
  fr <- vapply(x, function(v) if (is.na(v)) NA_real_ else sum(refok <= v) / (n + 1),
               numeric(1))
  z <- stats::qnorm(fr); z[!is.finite(z)] <- 0; z
}
.rank_normal <- function(x) .rank_normal_by(x, x)

.syn_rows <- function(syn, groups) {
  syn$posterior[syn$posterior$group %in% groups, ]
}

.syn_merge_inputs <- function(syn, groups) {
  merge(syn$posterior[syn$posterior$group %in% groups, ],
        syn$inputs[syn$inputs$group %in% groups, ],
        by = c("lineage", "group", "label"), all.x = TRUE)
}

negative_control_summary <- function(syn) {
  dnec <- load_noneuropean()
  mj <- load_major_lineages()
  out <- .syn_merge_inputs(syn, "eu_control")
  out$founder <- paste0(out$lineage, " (control)")
  out[, c("founder", "lineage", "nonjew_freq", "rarity", "depth_kyr", "carriers",
          "pub_adj", "nest_frac", "z_freq", "z_time", "z_nest",
          "post_mean_H2", "post_lo", "post_med", "post_hi")]
}

noneuropean_summary <- function(syn) {
  meta <- load_noneuropean()
  out <- .syn_merge_inputs(syn, "noneu_control")
  out$founder <- paste0(out$lineage, " (non-Eur)")
  out <- merge(out, meta[, c("mtree_name", "deep_origin", "ashkenazi_pct",
                              "medieval_erfurt")],
               by.x = "lineage", by.y = "mtree_name", all.x = TRUE)
  out$medieval_erfurt <- ifelse(nzchar(out$medieval_erfurt), out$medieval_erfurt, "--")
  out
}

major_lineage_summary <- function(syn) {
  meta <- load_major_lineages()
  out <- .syn_merge_inputs(syn, "panel")
  out <- merge(out, meta[, c("mtree_name", "expected_origin", "ashkenazi_pct",
                              "deep_origin")],
               by.x = "lineage", by.y = "mtree_name", all.x = TRUE)
  out$founder <- out$lineage
  post <- out[, c("founder", "expected_origin", "post_mean_H2",
                  "post_lo", "post_med", "post_hi")]
  post$model_call <- ifelse(post$post_mean_H2 >= 0.5, "Near Eastern", "European")
  post$concordant <- ifelse(post$expected_origin %in% c("Near Eastern", "European"),
                            post$model_call == post$expected_origin, NA)
  list(posterior = post, summary = out)
}

# Legacy name used by report tables.
bayesian_origin_model <- function(p0 = 0.5) {
  syn <- origin_synthesis_model()
  f_inp <- bayesian_origin_inputs()
  f_post <- syn$posterior[syn$posterior$group == "founder", ]
  f_post$founder <- f_post$lineage
  fin <- syn$inputs[syn$inputs$group == "founder", ]
  scores <- data.frame(
    founder = fin$lineage,
    z_freq = fin$z_freq, z_time = fin$z_time, z_nest = fin$z_nest,
    rarity = fin$rarity, nest_frac = fin$nest_frac,
    stringsAsFactors = FALSE
  )
  list(inputs = f_inp, scores = scores, posterior = f_post[, c("founder",
    "post_mean_H2", "post_lo", "post_med", "post_hi")],
    prior_sensitivity = data.frame())
}

negative_control_model <- function(p0 = 0.5) {
  syn <- origin_synthesis_model()
  summary_tbl <- negative_control_summary(syn)
  list(posterior = summary_tbl[, c("founder", "post_mean_H2", "post_lo",
                                   "post_med", "post_hi")],
       summary = summary_tbl)
}

noneuropean_model <- function(p0 = 0.5) {
  syn <- origin_synthesis_model()
  summary_tbl <- noneuropean_summary(syn)
  list(posterior = summary_tbl[, c("founder", "post_mean_H2", "post_lo",
                                   "post_med", "post_hi")],
       summary = summary_tbl)
}

major_lineage_model <- function(p0 = 0.5) {
  major_lineage_summary(origin_synthesis_model())
}

collect_channel_inputs <- function() {
  eps <- 1e-6
  rar <- function(f) log10(1e-3 / (f + eps))
  pub_of <- function(key, carriers) {
    pc <- public_haplogroup_context(key)
    public_time_adjustment(pc$public_tmrca_kyr, pc$public_age_confidence,
                           pc$public_jewish_ancient_anchor, carriers)
  }
  rows <- list()
  add <- function(lineage, group, label, row, key = lineage) {
    nc <- nesting_channels(row)
    carriers <- ifelse(is.na(row$medieval_jewish_carriers), 0, row$medieval_jewish_carriers)
    nest_frac <- nc$nest_frac
    rows[[length(rows) + 1L]] <<- data.frame(
      lineage = lineage, group = group, label = label,
      nonjew_freq = row$nonjew_freq,
      rarity = rar(row$nonjew_freq),
      depth_kyr = ifelse(is.na(row$ne_depth_kyr), 0, row$ne_depth_kyr),
      carriers = pmin(carriers, 4),
      pub_adj = pub_of(key, carriers),
      nest_frac = nest_frac,
      mitotree_node = row$mitotree_node,
      stringsAsFactors = FALSE)
  }

  founder_keys <- names(FOUNDER_ROOT)
  founder_inp <- compute_lineage_inputs_batch(
    founder_keys, nesting_modes = "ne_eu", rarefy_roots = founder_keys
  )
  for (i in seq_len(nrow(founder_inp)))
    add(founder_inp$lineage_key[i], "founder", NA_real_, founder_inp[i, ])

  neg_inp <- compute_lineage_inputs_batch(c("V7a2c1b", "U5a1f1a3"), nesting_modes = "ne_eu")
  for (i in seq_len(nrow(neg_inp)))
    add(neg_inp$lineage_key[i], "eu_control", 0, neg_inp[i, ])

  dnec <- load_noneuropean()
  ne_inp <- compute_lineage_inputs_batch(
    dnec$mtree_name,
    nesting_modes = rep("noneu_eu", nrow(dnec))
  )
  for (i in seq_len(nrow(ne_inp)))
    add(ne_inp$lineage_key[i], "noneu_control", 1, ne_inp[i, ])

  mj <- load_major_lineages()
  lab <- ifelse(mj$expected_origin == "Near Eastern", 1,
                ifelse(mj$expected_origin == "European", 0, NA_real_))
  mj_inp <- compute_lineage_inputs_batch(
    mj$mtree_name,
    nesting_modes = rep("ne_eu", nrow(mj))
  )
  for (i in seq_len(nrow(mj_inp)))
    add(mj_inp$lineage_key[i], "panel", lab[i], mj_inp[i, ])

  out <- do.call(rbind, rows)
  # Rarity uses rank-based normal scores (see .rank_normal) because its inputs
  # come from heterogeneous non-Jewish frequency sources; the antiquity/time and
  # nesting channels come from a single instrument (Mitotree) and keep the plain
  # z-score.
  out$z_freq <- .rank_normal(out$rarity)
  out$z_time <- (.zscore(out$depth_kyr) + .zscore(out$carriers) +
                   .zscore(out$pub_adj)) / 3
  out$z_nest <- .zscore(out$nest_frac)
  out
}

# Penalized (Gaussian-prior) logistic MAP for label ~ z_freq + z_time + z_nest.
# The three channels are scored so that positive values point to a Near Eastern
# origin (H2), so their weights are constrained to be non-negative (lower = 0);
# an unconstrained fit on a small panel can otherwise assign a channel a
# sign that contradicts its construction. The intercept is left free.
.CHANNEL_LOWER <- c(-Inf, 0, 0, 0)

.fit_logistic_map <- function(X, y, prior_sd, lower = .CHANNEL_LOWER) {
  neg_log_post <- function(b) {
    eta <- as.vector(X %*% b)
    ll <- sum(y * eta - log1p(exp(eta)))
    lp <- -0.5 * sum((b / prior_sd)^2)
    -(ll + lp)
  }
  optim(rep(0, ncol(X)), neg_log_post, method = "L-BFGS-B",
        lower = lower, hessian = TRUE)
}

# Full Bayesian (Laplace) logistic fit + coefficient posterior draws.
fit_logistic_bayes <- function(dat, tau = 2, tau0 = 5, M = 40000L, seed = 7L) {
  tr <- dat[!is.na(dat$label), ]
  X  <- cbind(1, tr$z_freq, tr$z_time, tr$z_nest)
  colnames(X) <- c("(Intercept)", "z_freq", "z_time", "z_nest")
  y  <- tr$label
  prior_sd <- c(tau0, tau, tau, tau)
  opt <- .fit_logistic_map(X, y, prior_sd)
  map <- opt$par
  cov <- tryCatch(solve(opt$hessian), error = function(e) diag(prior_sd^2))
  cov <- (cov + t(cov)) / 2
  set.seed(seed)
  ch <- tryCatch(chol(cov),
                 error = function(e) chol(cov + diag(1e-6, ncol(cov))))
  gen <- function(nn) sweep(matrix(rnorm(nn * ncol(X)), nn, ncol(X)) %*% ch, 2, map, "+")
  # Constrained channels (columns 2:4) must stay non-negative. Draw from the
  # Laplace-approximate posterior and reject draws that violate the constraint
  # (a truncated normal); if the MAP sits far enough on the boundary that too
  # few survive, fall back to clamping so the routine always returns M draws.
  cons <- which(is.finite(.CHANNEL_LOWER) & .CHANNEL_LOWER == 0)
  raw <- gen(6L * M)
  ok <- rowSums(raw[, cons, drop = FALSE] < 0) == 0
  draws <- raw[ok, , drop = FALSE]
  if (nrow(draws) >= M) {
    draws <- draws[seq_len(M), , drop = FALSE]
  } else {
    draws <- gen(M)
    for (j in cons) draws[, j] <- pmax(draws[, j], 0)
  }
  colnames(draws) <- colnames(X)
  list(map = setNames(map, colnames(X)), cov = cov, draws = draws,
       prior_sd = prior_sd)
}

# Posterior P(H2) for every lineage (founders included) under the fitted model.
fitted_predict <- function(fit, dat) {
  X <- cbind(1, dat$z_freq, dat$z_time, dat$z_nest)
  p <- plogis(X %*% t(fit$draws))
  data.frame(
    lineage = dat$lineage, group = dat$group, label = dat$label,
    post_mean_H2 = rowMeans(p),
    post_lo = apply(p, 1, quantile, 0.05),
    post_med = apply(p, 1, quantile, 0.5),
    post_hi = apply(p, 1, quantile, 0.95),
    stringsAsFactors = FALSE)
}

fitted_coef_summary <- function(fit) {
  d <- fit$draws
  data.frame(
    term = colnames(d), mean = colMeans(d),
    lo = apply(d, 2, quantile, 0.05), med = apply(d, 2, quantile, 0.5),
    hi = apply(d, 2, quantile, 0.95), stringsAsFactors = FALSE)
}

# Leave-one-out over the labeled rows. Both the standardization AND the
# regression coefficients are refit per fold with the held-out row excluded from
# the mean/sd (and rank) pool, so a lineage never informs its own channel scores.
# This is a fully out-of-sample check rather than one that leaks the held-out
# point through a globally-fixed standardization. Channel transforms mirror
# collect_channel_inputs (rank-normal rarity, z-scored time/nesting).
fitted_loo <- function(dat, tau = 2, tau0 = 5) {
  lab_idx <- which(!is.na(dat$label))
  prior_sd <- c(tau0, tau, tau, tau)
  z_channels <- function(ref) {
    zf <- .rank_normal_by(dat$rarity, ref$rarity)
    zt <- (.zscore_by(dat$depth_kyr, ref$depth_kyr) +
             .zscore_by(dat$carriers, ref$carriers) +
             .zscore_by(dat$pub_adj, ref$pub_adj)) / 3
    zn <- .zscore_by(dat$nest_frac, ref$nest_frac)
    cbind(1, zf, zt, zn)
  }
  p <- vapply(seq_along(lab_idx), function(k) {
    hold <- lab_idx[k]
    Xall <- z_channels(dat[-hold, , drop = FALSE])
    tr_idx <- setdiff(lab_idx, hold)
    b <- .fit_logistic_map(Xall[tr_idx, , drop = FALSE], dat$label[tr_idx], prior_sd)$par
    plogis(sum(Xall[hold, ] * b))
  }, numeric(1))
  tr <- dat[lab_idx, ]
  call <- as.integer(p >= 0.5)
  data.frame(lineage = tr$lineage, group = tr$group, label = tr$label,
             loo_p = round(p, 3), loo_call = call,
             correct = call == tr$label, stringsAsFactors = FALSE)
}

origin_synthesis_model <- function() {
  dat  <- collect_channel_inputs()
  fit  <- fit_logistic_bayes(dat)
  list(inputs = dat, coef = fitted_coef_summary(fit),
       posterior = fitted_predict(fit, dat), loo = fitted_loo(dat), fit = fit)
}

# Backward-compatible alias.
fitted_calibrated_model <- origin_synthesis_model

# --------------------------------------------------------------------------
# Runner -- tables and figures.
# --------------------------------------------------------------------------
run_advanced_models <- function() {
  message("== Module D: advanced branching + Bayesian origin models ==")
  fs <- frequency_sources()
  message(sprintf("Frequency sources: Ashkenazi=%s, non-Jewish=%s",
                  fs$ashkenazi, fs$nonjew))
  save_table(data.frame(
    channel = c("ashkenazi", "non_jewish"),
    source = c(fs$ashkenazi, fs$nonjew),
    stringsAsFactors = FALSE
  ), "D_frequency_sources.csv")

  audit_keys <- c(
    names(FOUNDER_ROOT), "V7a2c1b", "U5a1f1a3",
    load_noneuropean()$mtree_name, load_major_lineages()$mtree_name
  )
  audit_modes <- c(
    rep("ne_eu", length(FOUNDER_ROOT) + 2L),
    rep("noneu_eu", nrow(load_noneuropean())),
    rep("ne_eu", nrow(load_major_lineages()))
  )
  audit_rarefy <- c(
    names(FOUNDER_ROOT),
    rep(NA_character_, 2L + nrow(load_noneuropean()) + nrow(load_major_lineages()))
  )
  audit <- compute_lineage_inputs_batch(
    audit_keys, nesting_modes = audit_modes, rarefy_roots = audit_rarefy
  )
  audit$panel <- c(
    rep("founder", length(FOUNDER_ROOT)),
    rep("eu_control", 2L),
    rep("noneu_control", nrow(load_noneuropean())),
    rep("major_panel", nrow(load_major_lineages()))
  )
  save_table(audit, "D_lineage_computed_inputs.csv")

  # ---- D1. negative-binomial size sensitivity ----
  nb <- nbinom_size_sensitivity()
  save_table(nb, "D_nbinom_size_sensitivity.csv")

  open_png("D1_nbinom_size.png", width = 1150, height = 800)
  op <- par(mar = c(4.5, 4.8, 3, 1))
  sizes <- sort(unique(nb$nb_size))
  pal <- colorRampPalette(c("#c1121f", "#e5793a", "#f3a712", "#2a9d8f", "#264653"))(length(sizes))
  plot(NA, xlim = range(nb$growth_rate),
       ylim = range(nb$prob_one_descendant * 100),
       xlab = "Growth ratio m (daughters / mother)",
       ylab = "P(exactly one descendant at gen 15)  [%]",
       main = "Offspring overdispersion: negative-binomial size sweep")
  for (i in seq_along(sizes)) {
    s <- nb[nb$nb_size == sizes[i], ]
    s <- s[order(s$growth_rate), ]
    lines(s$growth_rate, s$prob_one_descendant * 100, type = "b", pch = 19,
          col = pal[i], lwd = 2)
  }
  legend("topright", bty = "n", lwd = 2, pch = 19, col = pal,
         legend = ifelse(is.finite(sizes), sprintf("NB size = %g", sizes),
                         "Poisson (size=Inf)"))
  grid(); par(op); dev.off()

  # ---- D2. piecewise historical growth ----
  pw <- piecewise_growth_model()
  save_table(pw$schedule, "D_growth_schedule.csv")
  save_table(pw$comparison, "D_piecewise_vs_constant.csv")
  print(pw$comparison)

  open_png("D2_piecewise_growth.png", width = 1200, height = 820)
  op <- par(mar = c(4.5, 4.8, 4, 1))
  sch <- pw$schedule
  plot(sch$generation, sch$growth_ratio, type = "s", lwd = 3, col = "#264653",
       ylim = range(c(sch$growth_ratio, pw$m_const, 1)) + c(-0.03, 0.03),
       xlab = "generation (1 = founder phase)", ylab = "per-generation growth ratio m_t",
       main = "Piecewise historical growth vs matched constant growth")
  abline(h = 1, lty = 3, col = "grey60")
  abline(h = pw$m_const, lty = 2, lwd = 2, col = "#c1121f")
  pc <- pw$comparison
  legend("topleft", bty = "n", cex = 0.95,
         legend = c(
           sprintf("piecewise: P(1 desc)=%.2f%%, P(extinct)=%.1f%%",
                   100 * pc$prob_one_descendant[2], 100 * pc$prob_extinct[2]),
           sprintf("constant m=%.3f: P(1 desc)=%.2f%%, P(extinct)=%.1f%%",
                   pw$m_const, 100 * pc$prob_one_descendant[1], 100 * pc$prob_extinct[1])),
         lty = c(1, 2), lwd = c(3, 2), col = c("#264653", "#c1121f"))
  par(op); dev.off()

  # ---- D3. unified origin synthesis (single fitted model) ----
  syn <- origin_synthesis_model()
  bm_inp <- bayesian_origin_inputs()
  bm_post <- syn$posterior[syn$posterior$group == "founder", ]
  bm_post$founder <- bm_post$lineage
  bm_sc <- syn$inputs[syn$inputs$group == "founder",
                      c("lineage", "rarity", "nest_frac", "z_freq", "z_time", "z_nest")]
  names(bm_sc)[1] <- "founder"
  save_table(bm_inp, "D_bayes_inputs.csv")
  save_table(bm_sc, "D_bayes_scores.csv")
  save_table(bm_post[, c("founder", "post_mean_H2", "post_lo", "post_med", "post_hi")],
             "D_bayes_posterior.csv")
  save_table(syn$coef, "D_bayes_coefficients.csv")
  save_table(syn$inputs[, c("lineage", "group", "label", "rarity", "depth_kyr",
                            "carriers", "pub_adj", "nest_frac",
                            "z_freq", "z_time", "z_nest")],
             "D_fitted_inputs.csv")
  save_table(syn$coef, "D_fitted_coefficients.csv")
  save_table(syn$posterior, "D_fitted_posterior.csv")
  save_table(syn$loo, "D_fitted_loo.csv")
  print(bm_post)

  nc <- negative_control_summary(syn)
  save_table(nc, "D_negative_control.csv")
  for (i in seq_len(nrow(nc))) {
    message(sprintf("Negative control %s posterior P(H2) = %.2f [%.2f, %.2f]",
                    nc$founder[i], nc$post_mean_H2[i],
                    nc$post_lo[i], nc$post_hi[i]))
  }

  ne_ctrl <- noneuropean_summary(syn)
  save_table(ne_ctrl, "D_noneuropean_control.csv")
  for (i in seq_len(nrow(ne_ctrl))) {
    message(sprintf("Non-European control %s posterior P(H2) = %.2f [%.2f, %.2f]",
                    ne_ctrl$founder[i], ne_ctrl$post_mean_H2[i],
                    ne_ctrl$post_lo[i], ne_ctrl$post_hi[i]))
  }

  ml <- major_lineage_summary(syn)
  save_table(ml$summary, "D_major_lineages.csv")
  evaluable <- !is.na(ml$posterior$concordant)
  n_conc <- sum(ml$posterior$concordant[evaluable])
  message(sprintf("Broader-pool concordance with published origin: %d/%d evaluable (%d ambiguous)",
                  n_conc, sum(evaluable), sum(!evaluable)))
  message(sprintf("LOO labeled-panel accuracy: %d/%d (%.0f%%)",
                  sum(syn$loo$correct), nrow(syn$loo),
                  100 * mean(syn$loo$correct)))
  print(syn$coef[, c("term", "mean", "lo", "hi")])

  po <- rbind(
    data.frame(founder = bm_post$founder, post_mean_H2 = bm_post$post_mean_H2,
               post_lo = bm_post$post_lo, post_hi = bm_post$post_hi,
               kind = "founder", stringsAsFactors = FALSE),
    data.frame(founder = nc$founder, post_mean_H2 = nc$post_mean_H2,
               post_lo = nc$post_lo, post_hi = nc$post_hi,
               kind = "eu_control", stringsAsFactors = FALSE),
    data.frame(founder = ne_ctrl$founder, post_mean_H2 = ne_ctrl$post_mean_H2,
               post_lo = ne_ctrl$post_lo, post_hi = ne_ctrl$post_hi,
               kind = "noneu_control", stringsAsFactors = FALSE)
  )
  po <- po[order(po$post_mean_H2), ]

  open_png("D3_bayes_origin.png", width = 1200, height = 960)
  op <- par(mar = c(4.8, 8.5, 3.5, 1))
  y <- seq_len(nrow(po))
  col_map <- c(founder = "#c1121f", eu_control = "#4361ee", noneu_control = "#2a9d8f")
  pch_map <- c(founder = 19L, eu_control = 17L, noneu_control = 15L)
  col_pt <- col_map[po$kind]
  pch_pt <- pch_map[po$kind]
  plot(po$post_mean_H2, y, xlim = c(0, 1), ylim = c(0.5, nrow(po) + 0.5),
       pch = pch_pt, cex = 1.6, col = col_pt, yaxt = "n", ylab = "",
       xlab = "posterior P(H2 side of European-host contrast)",
       main = "Origin synthesis: fitted logistic on standardized data channels")
  segments(po$post_lo, y, po$post_hi, y, col = col_pt, lwd = 3)
  abline(v = 0.5, lty = 2, col = "grey50")
  axis(2, at = y, labels = po$founder, las = 1)
  text(po$post_mean_H2, y + 0.30, sprintf("%s [%.3f, %.3f]",
       .fmt_post(po$post_mean_H2), po$post_lo, po$post_hi), cex = 0.75)
  legend("bottomright", bty = "n", pch = c(19, 17, 15),
         col = c("#c1121f", "#4361ee", "#2a9d8f"),
         legend = c("founder lineage (Levantine H2)",
                    "European negative control",
                    "deep non-European lineage (positive control)"))
  par(op); dev.off()

  # ---- D4. broader frequent Ashkenazi pool ----
  pm <- rbind(
    data.frame(founder = bm_post$founder, expected_origin = "founder (Levantine)",
               post_mean_H2 = bm_post$post_mean_H2, post_lo = bm_post$post_lo,
               post_hi = bm_post$post_hi, stringsAsFactors = FALSE),
    data.frame(founder = ml$posterior$founder, expected_origin = ml$posterior$expected_origin,
               post_mean_H2 = ml$posterior$post_mean_H2, post_lo = ml$posterior$post_lo,
               post_hi = ml$posterior$post_hi, stringsAsFactors = FALSE)
  )
  pm <- pm[order(pm$post_mean_H2), ]
  col_pm <- c(`founder (Levantine)` = "#8d99ae", `Near Eastern` = "#c1121f",
              European = "#4361ee", Ambiguous = "#6c757d")
  pch_pm <- c(`founder (Levantine)` = 19L, `Near Eastern` = 18L,
              European = 18L, Ambiguous = 15L)
  open_png("D4_broader_pool.png", width = 1200, height = 820)
  op <- par(mar = c(4.8, 8.0, 3.5, 1))
  y <- seq_len(nrow(pm))
  cc <- col_pm[pm$expected_origin]; pp <- pch_pm[pm$expected_origin]
  plot(pm$post_mean_H2, y, xlim = c(0, 1), ylim = c(0.5, nrow(pm) + 0.5),
       pch = pp, cex = 1.9, col = cc, yaxt = "n", ylab = "",
       xlab = "posterior P(Near Eastern origin / H2)",
       main = "Broader frequent Ashkenazi pool: model vs. published origin")
  segments(pm$post_lo, y, pm$post_hi, y, col = cc, lwd = 3)
  abline(v = 0.5, lty = 2, col = "grey50")
  axis(2, at = y, labels = pm$founder, las = 1)
  text(pm$post_mean_H2, y + 0.30, sprintf("%.2f", pm$post_mean_H2), cex = 0.8)
  legend("bottomright", bty = "n", pch = c(19, 18, 18, 15),
         col = c("#8d99ae", "#c1121f", "#4361ee", "#6c757d"),
         legend = c("four founders (reference)",
                    "expected Near Eastern (published assignments)",
                    "expected European (published assignments)",
                    "ambiguous/disputed assignment"))
  par(op); dev.off()

  frequency_source_comparison()

  invisible(list(nbinom = nb, piecewise = pw, synthesis = syn,
                 negative_control = nc, noneuropean_control = ne_ctrl,
                 major_lineages = ml))
}

# Compare posteriors under published (Brook + Livni-Skorecki) vs Mitotree-only freqs.
frequency_source_comparison <- function() {
  configs <- list(
    brook_1kg = list(ashkenazi = "brook", nonjew = "1kg_eur"),
    brook_livni = list(ashkenazi = "brook", nonjew = "livni_skorecki"),
    mitotree_only = list(ashkenazi = "mitotree", nonjew = "mitotree")
  )
  rows <- list()
  for (tag in names(configs)) {
    set_frequency_sources(
      ashkenazi = configs[[tag]]$ashkenazi,
      nonjew = configs[[tag]]$nonjew
    )
    syn <- origin_synthesis_model()
    post <- syn$posterior
    post$config <- tag
    rows[[length(rows) + 1L]] <- post[, c("config", "lineage", "group",
                                         "post_mean_H2", "post_lo", "post_hi")]
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  save_table(out, "D_frequency_source_comparison.csv")
  set_frequency_sources(ashkenazi = "brook", nonjew = "1kg_eur")
  out
}

if (sys.nframe() == 0) { source(file.path("R", "utils.R")); run_advanced_models() }
