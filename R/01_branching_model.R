# 01_branching_model.R
# Module A -- re-implementation of the Livni & Skorecki (2025) branching-process
# / founder-vs-host model for Ashkenazi mtDNA lineages.
#
# Reproduces:
#   * Table 1  -- P(single matriline descendant after 15 generations) vs growth.
#   * Fig 4/5  -- founder-vs-absorbed detectability (cumulative Poisson).
#   * Mitotree/GenBank enriched founder counts for K1a1b1a, K1a9, K2a2a, N1b2.
#   * The Livni-Skorecki detectability equations as model framework. Absorption
#     is estimated from the enriched reference set's singleton fraction, with
#     the Behar 2006 Ashkenazi sample as a separate, lineage-specific anchor
#     (see estimate_absorption / A_absorption_estimate.csv).
#   * Section 2.4.4 / Fig 6 -- feasibility of the "Euro-Levantine" scenario
#                      (probability of >=3 major founders ~ 1.5%, or ~5.3% with a
#                       Dirichlet-derived rather than assumed major-frequency term).
#
# All equations are cited to the preprint (ssrn-5035272 / hum.gen.2025.201445).

if (!exists("gw_descendant_pmf")) source(file.path("R", "utils.R"))
if (!exists("build_topology")) source(file.path("R", "02_phylogeography.R"))

# --------------------------------------------------------------------------
# A1. Table 1 -- probability of a single descendant after 15 generations.
#     (Livni-Skorecki eq 1; Galton-Watson lineage extinction.)
# --------------------------------------------------------------------------
# Growth grid for the sensitivity sweep. The old c(1.0 .. 1.1) implied 150
# founder matrilines expanding to ~450-750 women by the present, against an
# actual Ashkenazi maternal population of order 10^6. The grid is centred on
# GW_M_EFFECTIVE (the single demographic constant, see R/utils.R), with 1.10
# retained as a conservative low end and 1.40 as a high one.
GROWTH_GRID <- sort(unique(round(
  c(1.10, 1.15, 1.20, GW_M_EFFECTIVE, 1.28, 1.30, 1.35, 1.40), 4)))

reproduce_table1 <- function(k = GW_K, model = "poisson") {
  growth <- GROWTH_GRID
  p1 <- vapply(growth, function(m) gw_prob_single(m, k = k, model = model), numeric(1))
  data.frame(growth_rate = growth, prob_one_descendant = round(p1, 7))
}

# Slow but transparent independent check of the PGF/FFT Galton-Watson
# implementation. It explicitly composes offspring distributions by truncated
# polynomial convolution. Used as a generated audit table, not for production.
conv_trunc_direct <- function(a, b, jmax) {
  out <- numeric(jmax + 1L)
  for (i in seq_along(a)) for (j in seq_along(b)) {
    k <- (i - 1L) + (j - 1L)
    if (k <= jmax) out[k + 1L] <- out[k + 1L] + a[i] * b[j]
  }
  out
}

# jmax must cover the bulk of Z_k or the truncated convolution disagrees with
# the exact result for reasons of its own grid, not of either method: at
# k = 27 a jmax of 150 produced a spurious 4.8% mismatch in the audit table.
gw_descendant_pmf_direct <- function(m, k = GW_K,
                                     jmax = max(150L, ceiling(40 * m^k))) {
  off <- dpois(0:jmax, lambda = m)
  off <- off / sum(off)
  dist <- numeric(jmax + 1L)
  dist[2L] <- 1 # Z_0 = 1 founder mother.
  for (gen in seq_len(k)) {
    out <- numeric(jmax + 1L)
    power <- numeric(jmax + 1L)
    power[1L] <- 1 # previous PGF to power 0
    for (x in 0:jmax) {
      if (x > 0L) power <- conv_trunc_direct(power, dist, jmax)
      out <- out + off[x + 1L] * power
    }
    dist <- out / sum(out)
  }
  dist
}

branching_math_audit <- function(k = GW_K) {
  growth <- c(1.0, 1.025, 1.05, 1.075, 1.1)
  do.call(rbind, lapply(growth, function(m) {
    fast <- gw_descendant_pmf(m, k = k, jmax = 400L)[2L]
    direct <- gw_descendant_pmf_direct(m, k = k, jmax = 150L)[2L]
    data.frame(
      growth_rate = m,
      pgf_fft_prob_one = fast,
      direct_convolution_prob_one = direct,
      absolute_difference = abs(fast - direct)
    )
  }))
}

offspring_law_sensitivity <- function(k = GW_K) {
  growth <- GROWTH_GRID
  do.call(rbind, lapply(c("poisson", "geometric"), function(model) {
    data.frame(
      offspring_model = model,
      growth_rate = growth,
      prob_extinct = vapply(growth, function(m) gw_descendant_pmf(m, k = k, model = model)[1L], numeric(1)),
      prob_one_descendant = vapply(growth, function(m) gw_prob_single(m, k = k, model = model), numeric(1)),
      prob_ge_three_descendants = vapply(growth, function(m) {
        pmf <- gw_descendant_pmf(m, k = k, model = model)
        sum(pmf[4:length(pmf)])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
}

# --------------------------------------------------------------------------
# A2. Founder-vs-absorbed detectability.
#     A lineage of population frequency nu is expected to appear lambda = N*nu
#     times in a sample of N (eq 20). The probability it is seen in >= j copies
#     follows the cumulative Poisson (eqs 18-19).
# --------------------------------------------------------------------------
prob_at_least <- function(j, lambda) stats::ppois(j - 1, lambda, lower.tail = FALSE)

detection_curves <- function(N = 600L,
                             freqs = c(major = 0.043, minor = 0.015,
                                       thin_minor = 0.0068,
                                       absorbed_G0 = 0.001)) {
  js <- 1:6
  do.call(rbind, lapply(names(freqs), function(nm) {
    nu <- freqs[[nm]]
    lambda <- N * nu
    data.frame(class = nm, freq = nu, N = N, lambda = lambda,
               copies_ge = js,
               prob = vapply(js, prob_at_least, numeric(1), lambda = lambda))
  }))
}

# --------------------------------------------------------------------------
# A3. Mitotree empirical founder summary.
# --------------------------------------------------------------------------
mitotree_founder_summary <- function(founders = c("K1a1b1a", "K1a9", "K2a2a", "N1b2")) {
  samples <- load_analysis_samples()
  struct <- load_structure()
  topo <- build_topology(struct)
  n_modern <- sum(samples$SubjectType == "Modern")
  n_ancient <- sum(samples$SubjectType == "Ancient")
  do.call(rbind, lapply(founders, function(f) {
    desc <- descendants(f, topo)
    modern <- samples[samples$SubjectType == "Modern" &
                        samples$MitotreeHaplogroup %in% desc, ]
    ancient <- samples[samples$SubjectType == "Ancient" &
                         samples$MitotreeHaplogroup %in% desc, ]
    known <- modern[!is.na(modern$Country), ]
    regions <- sort(table(as.character(known$region)), decreasing = TRUE)
    countries <- sort(table(known$Country), decreasing = TRUE)
    data.frame(
      founder = f,
      descendant_nodes = length(desc),
      modern_count = nrow(modern),
      modern_frequency = nrow(modern) / n_modern,
      ancient_count = nrow(ancient),
      distinct_modern_haplogroups = length(unique(modern$MitotreeHaplogroup)),
      known_country_count = nrow(known),
      top_regions = paste(paste(names(regions), as.integer(regions), sep = ":"), collapse = "; "),
      top_countries = paste(paste(names(head(countries, 8)), as.integer(head(countries, 8)), sep = ":"), collapse = "; "),
      n_modern_analysis = n_modern,
      n_ancient_analysis = n_ancient,
      stringsAsFactors = FALSE
    )
  }))
}

mitotree_ancient_founder_rows <- function(founders = c("K1a1b1a", "K1a9", "K2a2a", "N1b2")) {
  samples <- load_analysis_samples()
  struct <- load_structure()
  topo <- build_topology(struct)
  out <- do.call(rbind, lapply(founders, function(f) {
    desc <- descendants(f, topo)
    ancient <- samples[samples$SubjectType == "Ancient" &
                         samples$MitotreeHaplogroup %in% desc, ]
    if (!nrow(ancient)) return(NULL)
    data.frame(
      founder = f,
      subject = ancient$Subject,
      haplogroup = ancient$MitotreeHaplogroup,
      country = ancient$Country,
      year = ancient$AgeEstimateMean,
      study = ancient$Study,
      title = ancient$Title,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out)) data.frame() else out[order(out$founder, out$year), ]
}

# --------------------------------------------------------------------------
# A4. Livni-Skorecki absorption-rate estimation.
#     We DO estimate the host/convert absorption rate. The estimator uses the
#     singleton fraction of a maternal sample: individuals who are the sole
#     representative of their mtDNA lineage are the candidate recently-absorbed
#     (host / convert) matrilines, whereas founder lineages recur in many
#     copies. We estimate it from the largest available maternal sample (the
#     enriched, deduplicated Mitotree + GenBank set) and cross-check against the
#     dedicated Behar 2006 Ashkenazi population sample as an Ashkenazi-specific
#     benchmark. Caveats about sampling are quantified rather than used to
#     abstain from estimation (see report).
# --------------------------------------------------------------------------
# Absorption ratio a_eps = fraction of contemporary individuals whose matriline
# was absorbed (host / convert), estimated from the singleton fraction.
# Per-generation absorption fraction rho solves (1 - a_eps) = (1 - rho)^K, i.e.
# rho = 1 - (1 - a_eps)^(1/K)   (Livni-Skorecki eqs 9-12).
absorption_rate <- function(a_eps, K) 1 - (1 - a_eps)^(1 / K)

# Number of absorbed lineages ever entering the population (eq 7/8), then
# reduced by lineage extinction (surviving fraction 1 - s0).
absorbed_lineages <- function(rho, r, M_current, s0 = 2/3, K = NULL) {
  # Continuous-growth integral of per-generation absorbed mothers:
  #   int_0^K M0 e^{rt} dt = (M_current - M0) / r = (M_current / r) (1 - e^{-rK}).
  # The (1 - e^{-rK}) factor is NOT negligible over a finite window -- dropping
  # it overstates the count. The size of the correction depends on r = log(m):
  # 1.51x at K = 15 under the retired m = 1.075, but only ~1.002x at K = 25
  # under the current m = 1.28. Applied whenever K is supplied.
  raw <- (rho / r) * M_current
  if (!is.null(K) && is.finite(K) && is.finite(r) && r > 0) {
    raw <- raw * (1 - exp(-r * K))
  }
  surviving <- raw * (1 - s0)
  list(raw = raw, surviving = surviving)
}

# How far is a_eps actually determined by absorption, as opposed to by how many
# samples were drawn and how finely the tree labels them?
#
# Rarefying ONE unchanging pool and truncating its clade labels answers that
# directly, and the answer is unflattering: the singleton fraction of the
# enriched set runs from 0.22 at n = 57,531 to ~0.96 at n = 100, and from 0.00
# to 0.22 as the label is truncated from 1 character to full depth. Under Ewens
# sampling E[singleton fraction] = theta/(theta + n - 1), so it is a function of
# sample size by construction. The estimate is therefore only comparable
# between samples of equal n at equal clade resolution, and the enriched set's
# LOW value reflects its large n, not less absorption -- which is the opposite
# of reading it as an upper bound.
absorption_identifiability <- function(haplogroups,
                                       n_grid = c(100L, 565L, 2000L, 10000L),
                                       depth_grid = c(1L, 3L, 4L, 6L),
                                       K = GW_K, reps = 40L, seed = 1L) {
  hg <- haplogroups[!is.na(haplogroups) & nzchar(haplogroups)]
  if (!length(hg)) return(data.frame())
  one <- function(x) sum(table(x) == 1) / length(x)
  old <- if (exists(".Random.seed", .GlobalEnv))
    get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  rows <- list()
  full_n <- length(hg)
  for (n in c(n_grid[n_grid < full_n], full_n)) {
    a <- if (n < full_n) mean(replicate(reps, one(sample(hg, n)))) else one(hg)
    rows[[length(rows) + 1L]] <- data.frame(
      varied = "sample_size", setting = as.character(n),
      n = n, clade_resolution = "full",
      a_eps = round(a, 4), absorption_rate_per_gen = round(absorption_rate(a, K), 5),
      stringsAsFactors = FALSE)
  }
  for (d in depth_grid) {
    h2 <- substr(hg, 1L, d)
    a <- one(h2)
    rows[[length(rows) + 1L]] <- data.frame(
      varied = "clade_resolution", setting = paste0(d, "_chars"),
      n = full_n, clade_resolution = paste0(d, " chars"),
      a_eps = round(a, 4), absorption_rate_per_gen = round(absorption_rate(a, K), 5),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Individual-level singleton fraction: the share of sampled individuals who are
# the only carrier of their (fully resolved) mtDNA lineage. This is the a_eps
# estimator (candidate absorbed host matrilines).
singleton_fraction <- function(haplogroups) {
  hg <- haplogroups[!is.na(haplogroups) & nzchar(haplogroups)]
  if (!length(hg)) return(c(n = 0, n_lineages = 0, a_eps = NA_real_))
  tab <- table(hg)
  c(n = length(hg), n_lineages = length(tab), a_eps = sum(tab == 1) / length(hg))
}

# Estimate the absorption rate from the enriched Mitotree sample (treated as the
# largest available maternal sample) and from Behar 2006 (Ashkenazi benchmark),
# across a grid of founder-event depths K. Growth r, founder size M0, and
# extinction s0 are taken from the founder-event parameter table so the derived
# absorbed-lineage count is internally consistent (no external census guess).
estimate_absorption <- function(K_grid = c(25L, GW_K, 32L)) {
  params <- load_params()
  pv <- setNames(params$selected, params$symbol)
  m  <- as.numeric(pv[["m"]])          # per-generation growth ratio
  M0 <- as.numeric(pv[["M0"]])         # founder families
  s0 <- as.numeric(pv[["s0"]])         # lineage extinction ratio
  r  <- log(m)                          # continuous growth rate

  samples <- load_analysis_samples()
  modern  <- samples[samples$SubjectType == "Modern", ]
  jw <- modern$in_jewish_study %in% c(TRUE, "TRUE")

  behar <- load_behar()
  bv <- setNames(behar$value, behar$quantity)

  defs <- list(
    enriched_mitotree_modern = list(
      stats = singleton_fraction(modern$MitotreeHaplogroup),
      label = "Enriched Mitotree+GenBank modern (largest available sample)"),
    jewish_flagged_modern = list(
      stats = singleton_fraction(modern$MitotreeHaplogroup[jw]),
      label = "Ashkenazi/Jewish-flagged modern subset (small)"),
    behar2006_ashkenazi = list(
      stats = c(n = as.numeric(bv[["sample_size"]]),
                n_lineages = as.numeric(bv[["n_lineages"]]),
                a_eps = as.numeric(bv[["n_singletons"]]) / as.numeric(bv[["sample_size"]])),
      label = "Behar 2006 Ashkenazi population sample (benchmark)")
  )

  do.call(rbind, lapply(names(defs), function(nm) {
    st <- defs[[nm]]$stats
    a_eps <- unname(st[["a_eps"]])
    do.call(rbind, lapply(K_grid, function(K) {
      rho <- absorption_rate(a_eps, K)
      M_current <- M0 * m^K
      absL <- absorbed_lineages(rho, r, M_current, s0, K = K)
      data.frame(
        sample = nm,
        description = defs[[nm]]$label,
        n_individuals = unname(st[["n"]]),
        n_lineages = unname(st[["n_lineages"]]),
        absorbed_fraction_a_eps = round(a_eps, 4),
        K_generations = K,
        absorption_rate_per_gen = round(rho, 5),
        absorbed_lineages_surviving = round(absL$surviving, 1),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# --------------------------------------------------------------------------
# A5. Feasibility of the "Euro-Levantine" reconciliation (section 2.4.4).
#     Could a large (non-founder-effect) Roman-era Jewish population have
#     generated >=3 of the four major Ashkenazi founders by private mutation?
# --------------------------------------------------------------------------
# Per-generation HVS-I (np 16051-16400) mutation probability, from the Soares
# et al. (2009) clock of one substitution per 8,650 years, expressed at this
# project's generation length. The previous hardcoded 4.3e-3 implied a 37-year
# generation while `gens` assumed 25-year ones -- two different generation
# lengths inside one function.
MU_HVS1 <- GEN_YEARS / 8650

# Probability a founder lineage reaches `tau` of the present-day pool. Under
# neutral expansion from F founder lineages the present shares are approximately
# Dirichlet(1,...,1), so P(share >= tau) = (1 - tau)^(F - 1). Deriving it beats
# asserting it: at F = 150, tau = 2% this gives 0.049, not the 0.03 previously
# assumed.
sigma_major_dirichlet <- function(founder_families, tau_major = 0.02) {
  (1 - tau_major)^(founder_families - 1)
}

euro_levantine_feasibility <- function(pool_families = 4000, r = 0.05, gens = 40,
                                       mu_hv = MU_HVS1, founder_families = 150,
                                       sigma_major = 0.03, tau_major = 0.02) {
  families_at_event <- pool_families * exp(r * gens)
  total_mothers     <- (pool_families / r) * (exp(r * gens) - 1)
  n_private_haplos  <- total_mothers * mu_hv
  # Expected NEW HVS-I mutations on ONE present-day matriline over the window.
  #
  # This previously divided every mutation that ever arose by the CURRENT
  # population (n_private_haplos / families_at_event), which is only valid if
  # each mutation leaves exactly one living carrier. In a growing population
  # mutations that arose early leave many descendants, so that ratio understates
  # the per-lineage count by (r*gens)/(1 - exp(-r*gens)) = 2.31x at r = 0.05,
  # gens = 40. A forward simulation of the same process gives 0.176 of present
  # matrilines carrying a within-window mutation against the old formula's
  # 0.074, confirming the correction.
  per_lineage_mut <- mu_hv * gens
  freq_private    <- 1 - exp(-per_lineage_mut)
  expected_in_founders <- founder_families * freq_private
  # sigma_major: probability one such private lineage reaches major (~2%)
  # frequency. Reported both as the previously assumed constant and as the
  # Dirichlet-derived value, because the headline is sensitive to it.
  sigma_derived <- sigma_major_dirichlet(founder_families, tau_major)
  lambda_major  <- expected_in_founders * sigma_major
  lambda_derived <- expected_in_founders * sigma_derived
  list(
    families_at_event = families_at_event,
    total_mothers = total_mothers,
    n_private_haplotypes = n_private_haplos,
    per_lineage_mutations = per_lineage_mut,
    freq_private = freq_private,
    expected_private_in_founders = expected_in_founders,
    sigma_major_assumed = sigma_major,
    sigma_major_dirichlet = sigma_derived,
    lambda_major = lambda_major,
    prob_ge2_major = prob_at_least(2, lambda_major),
    prob_ge3_major = prob_at_least(3, lambda_major),
    lambda_major_dirichlet = lambda_derived,
    prob_ge3_major_dirichlet = prob_at_least(3, lambda_derived)
  )
}

# Sensitivity of the Euro-Levantine feasibility result to its three explicit
# assumptions (per-lineage major-emergence probability, private-mutation rate,
# and number of founder families). Reported so the "model-disfavored" conclusion
# is shown to hold across a plausible range rather than resting on one value,
# replacing the earlier practice of pinning sigma_major to a target percentage.
# `r` and `gens` are now swept too: they drive the per-lineage mutation count
# and so the headline, but were previously held fixed. The sigma_major grid
# spans the Dirichlet-derived value (0.049 at F = 150) as well as the lower
# value the base case assumes.
euro_levantine_sensitivity <- function(
    sigma_grid = c(0.01, 0.02, 0.03, 0.049, 0.06, 0.10),
    mu_grid = c(2.9e-3, MU_HVS1, 4.3e-3, 6.0e-3),
    founders_grid = c(100, 150, 250, 350),
    r_grid = c(0.03, 0.05, 0.07),
    gens_grid = c(30, 40, 50)) {
  base <- list(sigma_major = 0.03, mu_hv = MU_HVS1, founder_families = 150,
               r = 0.05, gens = 40)
  rows <- list()
  vary <- function(param, values) {
    for (v in values) {
      args <- base; args[[param]] <- v
      el <- euro_levantine_feasibility(
        r = args$r, gens = args$gens, mu_hv = args$mu_hv,
        founder_families = args$founder_families,
        sigma_major = args$sigma_major)
      rows[[length(rows) + 1L]] <<- data.frame(
        parameter = param, value = v,
        lambda_major = round(el$lambda_major, 4),
        prob_ge2_major = el$prob_ge2_major,
        prob_ge3_major = el$prob_ge3_major,
        prob_ge3_major_dirichlet_sigma = el$prob_ge3_major_dirichlet,
        stringsAsFactors = FALSE)
    }
  }
  vary("sigma_major", sigma_grid)
  vary("mu_hv", mu_grid)
  vary("founder_families", founders_grid)
  vary("r", r_grid)
  vary("gens", gens_grid)
  do.call(rbind, rows)
}

# --------------------------------------------------------------------------
# Runner -- produce tables and figures.
# --------------------------------------------------------------------------
run_branching_model <- function() {
  message("== Module A: Livni-Skorecki branching model ==")

  t1 <- reproduce_table1()
  save_table(t1, "A_table1_single_descendant.csv")
  print(t1)

  audit <- branching_math_audit()
  save_table(audit, "A_branching_math_audit.csv")

  sens <- offspring_law_sensitivity()
  save_table(sens, "A_offspring_law_sensitivity.csv")

  dc <- detection_curves(N = 600L)
  save_table(dc, "A_detection_curves.csv")

  mt_founders <- mitotree_founder_summary()
  save_table(mt_founders, "A_mitotree_founder_summary.csv")
  print(mt_founders[, c("founder", "modern_count", "modern_frequency",
                        "ancient_count", "distinct_modern_haplogroups")])

  mt_ancient <- mitotree_ancient_founder_rows()
  save_table(mt_ancient, "A_mitotree_ancient_founders.csv")

  # ---- Absorption-rate estimate (largest available sample + Behar benchmark) ----
  absorb <- estimate_absorption()
  save_table(absorb, "A_absorption_estimate.csv")
  # Identifiability check: how much of a_eps is sample size and label depth
  # rather than absorption (see absorption_identifiability).
  enr <- load_analysis_samples()
  absorb_id <- absorption_identifiability(
    enr$MitotreeHaplogroup[enr$SubjectType == "Modern"])
  save_table(absorb_id, "A_absorption_identifiability.csv")
  print(absorb[, c("sample", "n_individuals", "absorbed_fraction_a_eps",
                   "K_generations", "absorption_rate_per_gen")])
  headline <- absorb[absorb$sample == "enriched_mitotree_modern" &
                       absorb$K_generations == 25L, ]
  message(sprintf(
    "Absorption estimate (enriched Mitotree, K=25): a_eps=%.3f, rho=%.4f/gen",
    headline$absorbed_fraction_a_eps, headline$absorption_rate_per_gen))

  el <- euro_levantine_feasibility()
  el_df <- data.frame(quantity = names(el), value = unlist(el))
  save_table(el_df, "A_euro_levantine_feasibility.csv")
  el_sens <- euro_levantine_sensitivity()
  save_table(el_sens, "A_euro_levantine_sensitivity.csv")
  message(sprintf(
    "Euro-Levantine P(>=3 major founders) = %.4f%% at base assumptions; %.3f-%.3f%% across sensitivity grid",
    100 * el$prob_ge3_major,
    100 * min(el_sens$prob_ge3_major), 100 * max(el_sens$prob_ge3_major)))

  # ---- Figure A1: single-descendant probability vs growth ----
  open_png("A1_single_descendant.png")
  op <- par(mar = c(4.5, 4.8, 3, 1))
  plot(t1$growth_rate, t1$prob_one_descendant * 100, type = "b", pch = 19,
       col = "#264653", lwd = 2, xlab = "Growth ratio m (daughters / mother)",
       ylab = sprintf("P(exactly one descendant at gen %d)  [%%]", GW_K),
       main = "Galton-Watson: smallest surviving minor lineage")
  grid(); par(op); dev.off()

  # ---- Figure A2: founder vs absorbed detection curves ----
  open_png("A2_detection_curves.png")
  op <- par(mar = c(4.5, 4.8, 3, 1))
  cls <- unique(dc$class)
  cols <- c(major = "#c1121f", minor = "#e5793a", thin_minor = "#f3a712",
            absorbed_G0 = "#264653")
  plot(NA, xlim = c(1, 6), ylim = c(0, 1), xlab = "seen in >= k copies",
       ylab = "probability", main = "Detectability in a sample of N = 600")
  for (cl in cls) {
    sub <- dc[dc$class == cl, ]
    lines(sub$copies_ge, sub$prob, type = "b", pch = 19, lwd = 2,
          col = cols[[cl]])
  }
  abline(v = 2.5, lty = 3, col = "grey50")
  legend("topright", legend = sprintf("%s (nu=%.3f%%)", cls,
         100 * sapply(cls, function(c) dc$freq[dc$class == c][1])),
         col = cols[cls], lwd = 2, pch = 19, bty = "n", cex = 0.9)
  text(2.5, 0.05, "founders: >=3 copies\nabsorbed: singletons", pos = 4, cex = 0.8)
  par(op); dev.off()

  # ---- Figure A2b: offspring law sensitivity ----
  open_png("A2b_offspring_law_sensitivity.png")
  op <- par(mar = c(4.5, 4.8, 3, 1))
  plot(NA, xlim = range(sens$growth_rate), ylim = range(sens$prob_one_descendant * 100),
       xlab = "Growth ratio m (daughters / mother)",
       ylab = sprintf("P(exactly one descendant at gen %d)  [%%]", GW_K),
       main = "Sensitivity to offspring law")
  cols <- c(poisson = "#264653", geometric = "#c1121f")
  for (model in unique(sens$offspring_model)) {
    sub <- sens[sens$offspring_model == model, ]
    lines(sub$growth_rate, sub$prob_one_descendant * 100, type = "b", pch = 19,
          col = cols[[model]], lwd = 2)
  }
  legend("topright", legend = c("Poisson", "Geometric / NB(size=1)"),
         col = cols, lwd = 2, pch = 19, bty = "n")
  grid(); par(op); dev.off()

  # ---- Figure A3: Mitotree founder frequencies ----
  open_png("A3_mitotree_founder_frequencies.png")
  op <- par(mar = c(4.5, 5.2, 3, 1))
  vals <- mt_founders$modern_frequency * 100
  bp <- barplot(vals, names.arg = mt_founders$founder, col = "#264653",
                ylab = "frequency in Mitotree modern public samples (%)",
                main = "Major Ashkenazi founder clades in Mitotree S1")
  text(bp, vals, labels = sprintf("n=%d\n%.3f%%", mt_founders$modern_count, vals),
       pos = 3, xpd = NA, cex = 0.85)
  par(op); dev.off()

  # ---- Figure A4: absorption rate per generation vs founder-event depth ----
  open_png("A4_absorption_rate.png")
  op <- par(mar = c(4.5, 4.8, 3, 1))
  samp <- unique(absorb$sample)
  cols <- c(enriched_mitotree_modern = "#c1121f",
            jewish_flagged_modern = "#f3a712",
            behar2006_ashkenazi = "#264653")
  plot(NA, xlim = range(absorb$K_generations),
       ylim = c(0, max(absorb$absorption_rate_per_gen * 100, na.rm = TRUE)),
       xlab = "founder-event depth K (generations)",
       ylab = "estimated absorption rate per generation  [%]",
       main = "Host-lineage absorption rate (estimated)")
  for (s in samp) {
    sub <- absorb[absorb$sample == s, ]
    lines(sub$K_generations, sub$absorption_rate_per_gen * 100, type = "b",
          pch = 19, lwd = 2, col = cols[[s]])
  }
  legend("topright", bty = "n", lwd = 2, pch = 19, col = cols[samp],
         legend = c("enriched Mitotree (largest sample)",
                    "Jewish-flagged subset", "Behar 2006 (benchmark)")[match(samp, names(cols))])
  grid(); par(op); dev.off()

  # ---- Figure A6: Euro-Levantine major-founder probability ----
  open_png("A6_euro_levantine.png")
  op <- par(mar = c(4.5, 4.8, 3, 1))
  jj <- 0:6
  pj <- dpois(jj, el$lambda_major)
  bp <- barplot(pj, names.arg = jj, col = ifelse(jj >= 3, "#c1121f", "#adb5bd"),
                xlab = "number of major founders from a Euro-Levantine pool",
                ylab = "probability",
                main = sprintf("Euro-Levantine scenario (lambda=%.2f): P(>=3)=%.2f%%",
                               el$lambda_major, 100 * el$prob_ge3_major))
  par(op); dev.off()

  invisible(list(table1 = t1, sensitivity = sens, detection = dc, mitotree_founders = mt_founders,
                 absorption = absorb, euro_levantine = el))
}

if (sys.nframe() == 0) run_branching_model()
