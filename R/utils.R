# utils.R -- shared helpers for the Ashkenazi mtDNA origin pipeline.
# Base R only (no external deps) so the core pipeline always runs; the
# phylogeography module optionally uses ape/pegas/adegenet when available.

# ---- Paths -----------------------------------------------------------------
proj_root <- function() {
  # Assume scripts are sourced from the project root or R/ directory.
  wd <- getwd()
  if (file.exists(file.path(wd, "data", "derived"))) return(wd)
  up <- normalizePath(file.path(wd, ".."))
  if (file.exists(file.path(up, "data", "derived"))) return(up)
  wd
}

ROOT       <- proj_root()
DERIVED    <- file.path(ROOT, "data", "derived")
REFERENCE  <- file.path(ROOT, "data", "reference")
FIG_DIR    <- file.path(ROOT, "outputs", "figures")
TAB_DIR    <- file.path(ROOT, "outputs", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# Use a project-local library for optional popgen packages if present.
local_lib <- file.path(ROOT, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# ---- Data loaders ----------------------------------------------------------
load_samples <- function() {
  d <- read.csv(file.path(DERIVED, "s1_public_samples.csv"),
                stringsAsFactors = FALSE, colClasses = "character",
                na.strings = character(0))
  # Normalise blank / '.' placeholders to NA for the fields we use.
  blankify <- function(x) { x[x %in% c(".", "", "NA")] <- NA; x }
  d$Country          <- blankify(d$Country)
  d$AgeEstimateMean  <- suppressWarnings(as.numeric(blankify(d$AgeEstimateMean)))
  d$BirthYear        <- suppressWarnings(as.numeric(blankify(d$BirthYear)))
  d$region           <- region_of(d$Country)
  d
}

load_analysis_samples <- function() {
  enriched <- file.path(DERIVED, "enriched_deduped_samples.csv")
  if (file.exists(enriched)) {
    d <- read.csv(enriched, stringsAsFactors = FALSE, colClasses = "character",
                  na.strings = character(0))
    blankify <- function(x) { x[x %in% c(".", "", "NA")] <- NA; x }
    d$Country          <- blankify(d$Country)
    d$AgeEstimateMean  <- suppressWarnings(as.numeric(blankify(d$AgeEstimateMean)))
    d$BirthYear        <- suppressWarnings(as.numeric(blankify(d$BirthYear)))
    d$region           <- region_of(d$Country)
    return(d)
  }
  load_samples()
}

load_structure <- function() {
  read.csv(file.path(DERIVED, "s10_mitotree_structure.csv"),
           stringsAsFactors = FALSE, colClasses = "character",
           na.strings = character(0))
}

load_founders   <- function() read.csv(file.path(REFERENCE, "founders.csv"), stringsAsFactors = FALSE)
load_table3     <- function() read.csv(file.path(REFERENCE, "table3_parent_freq.csv"), stringsAsFactors = FALSE)
load_table7     <- function() read.csv(file.path(REFERENCE, "table7_nonjewish_freq.csv"), stringsAsFactors = FALSE)
load_1kg_clade_freq <- function() {
  read.csv(file.path(REFERENCE, "1kg_eur_clade_freq.csv"), stringsAsFactors = FALSE)
}
load_behar      <- function() read.csv(file.path(REFERENCE, "ashkenazi_sample_behar2006.csv"), stringsAsFactors = FALSE)
load_params     <- function() read.csv(file.path(REFERENCE, "founder_event_params.csv"), stringsAsFactors = FALSE)
load_ppnb       <- function() read.csv(file.path(REFERENCE, "ppnb_fernandez2014.csv"), stringsAsFactors = FALSE)
load_noneuropean <- function() read.csv(file.path(REFERENCE, "noneuropean_ashkenazi_lineages.csv"), stringsAsFactors = FALSE)
load_major_lineages <- function() read.csv(file.path(REFERENCE, "major_ashkenazi_lineages.csv"), stringsAsFactors = FALSE)

# ---- Frequency source configuration ----------------------------------------
# Ashkenazi (within-Ashkenazim) frequencies: Brook/Penninx via founders +
# panel CSVs, or Mitotree modern counts. Non-Jewish rarity: Livni-Skorecki
# Table 7 founders + full 1000 Genomes European (EUR: CEU+GBR+FIN+IBS+TSI,
# n=503) subclades (Haplogrep3), legacy Mitchell macro rows, or Mitotree only.
#
# Environment overrides (optional):
#   MTDNA_ASHKENAZI_FREQ=brook|mitotree
#   MTDNA_NONJEW_FREQ=1kg_eur|livni_skorecki|mitotree
.freq_sources <- new.env(parent = emptyenv())
.freq_sources$ashkenazi <- "brook"
.freq_sources$nonjew <- "1kg_eur"

set_frequency_sources <- function(ashkenazi = NULL, nonjew = NULL) {
  if (!is.null(ashkenazi)) {
    .freq_sources$ashkenazi <- match.arg(ashkenazi, c("brook", "mitotree"))
  }
  if (!is.null(nonjew)) {
    .freq_sources$nonjew <- match.arg(nonjew, c("1kg_eur", "livni_skorecki", "mitotree"))
  }
  invisible(frequency_sources())
}

frequency_sources <- function() {
  list(ashkenazi = .freq_sources$ashkenazi, nonjew = .freq_sources$nonjew)
}

init_frequency_sources_from_env <- function() {
  a <- Sys.getenv("MTDNA_ASHKENAZI_FREQ", unset = "")
  n <- Sys.getenv("MTDNA_NONJEW_FREQ", unset = "")
  if (nzchar(a)) set_frequency_sources(ashkenazi = a)
  if (nzchar(n)) set_frequency_sources(nonjew = n)
  frequency_sources()
}

# Brook/Penninx Ashkenazi %: public_haplogroup_page_targets.csv (all modeled
# lineages), then founders / major / non-European panel CSVs.
ashkenazi_frequency <- function(lineage_key) {
  pct <- function(x) {
    v <- suppressWarnings(as.numeric(gsub("^~\\s*", "", as.character(x))))
    if (!is.finite(v)) return(NA_real_)
    v / 100
  }
  targets <- tryCatch(
    read.csv(file.path(REFERENCE, "public_haplogroup_page_targets.csv"),
             stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(targets) && "ashkenazi_pct" %in% names(targets)) {
    hit <- targets$lineage_key == lineage_key
    if (any(hit)) {
      v <- pct(targets$ashkenazi_pct[hit][1])
      if (is.finite(v)) return(v)
    }
  }
  founders <- load_founders()
  hit <- founders$founder == lineage_key
  if (any(hit)) return(pct(founders$ashkenazi_pct[hit][1]))
  major <- load_major_lineages()
  hit <- major$mtree_name == lineage_key
  if (any(hit)) return(pct(major$ashkenazi_pct[hit][1]))
  ne <- load_noneuropean()
  hit <- ne$mtree_name == lineage_key
  if (any(hit)) return(pct(ne$ashkenazi_pct[hit][1]))
  NA_real_
}
load_public_haplogroup_pages <- function() {
  read.csv(file.path(ROOT, "data", "public_haplogroup_pages", "public_haplogroup_page_summary.csv"),
           stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0))
}

# ---- Geography classifier --------------------------------------------------
# Coarse region assignment used for the founder-lineage phylogeography.
# "Near East" is split into Levant / Arabia-Mesopotamia / Anatolia / Caucasus
# so we can inspect the Levantine signal specifically, then combined as needed.
REGION_LEVELS <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus",
                   "North_Africa", "Europe", "Central_South_Asia", "East_Asia",
                   "Sub_Saharan_Africa", "Oceania", "Americas", "Unknown")

region_of <- function(country) {
  c_levant  <- c("Israel", "Lebanon", "Jordan", "Syrian Arab Republic", "Syria",
                 "Palestinian Territory", "Palestine")
  c_arabmes <- c("Iraq", "Iran", "Kuwait", "Saudi Arabia", "Yemen",
                 "United Arab Emirates", "Oman", "Qatar", "Bahrain")
  c_anatol  <- c("Turkey", "Cyprus")
  c_caucus  <- c("Armenia", "Georgia", "Azerbaijan")
  c_nafr    <- c("Egypt", "Morocco", "Tunisia", "Algeria", "Libya", "Sudan",
                 "Western Sahara")
  regionize <- function(x) {
    if (is.na(x)) return("Unknown")
    # Country labels often carry a parenthetical population annotation
    # ("Israel (Druze)", "Georgia (Republic of Abkhazia)"). The Near-East and
    # North-Africa lists below are exact-match, so strip the annotation before
    # testing them -- otherwise every annotated Levantine/Caucasus/North-African
    # label falls through all five tests and lands in "Unknown".
    xb <- trimws(sub("\\s*\\(.*$", "", x))
    # Republics of the Russian Federation must be tested on the FULL label and
    # before the Europe keyword set, which matches the bare word "Russia":
    # the North Caucasus is part of the Near-East pool, not of Europe, and the
    # Siberian republics are not European at all.
    if (grepl(paste0("Adygea|Chechen|Chechnya|Dagestan|Kabardino|Ossetia|",
                     "Karachay|Ingush|Abkhazia"), x)) return("Caucasus")
    if (grepl("Yakut|Sakha|Buryat|Khakass|\\bTuva\\b|Altai|Evenk|Chukot", x))
      return("East_Asia")
    if (grepl("Antilles", x)) return("Americas")
    if (xb %in% c_levant)  return("Levant")
    if (xb %in% c_arabmes) return("Arabia_Mesopotamia")
    if (xb %in% c_anatol)  return("Anatolia")
    if (xb %in% c_caucus)  return("Caucasus")
    if (xb %in% c_nafr)    return("North_Africa")
    # Europe: match a broad keyword set (includes parenthetical annotations).
    europe_kw <- paste(c("Denmark","United Kingdom","England","Scotland","Wales",
      "Ireland","Iceland","Norway","Sweden","Finland","Estonia","Latvia",
      "Lithuania","Poland","Germany","Netherlands","Belgium","Luxembourg",
      "France","Spain","Portugal","Italy","Malta","Switzerland","Austria",
      "Czech","Slovakia","Hungary","Slovenia","Croatia","Bosnia","Serbia",
      "Montenegro","Kosovo","Albania","North Macedonia","Macedonia","Greece",
      "Bulgaria","Romania","Moldova","Ukraine","Belarus","Russian Federation",
      "Russia","Sardinia","Sicily","Orkney","Canary Islands","Faroe",
      "Channel Islands","Gibraltar","Crimea"),
      collapse = "|")
    if (grepl(europe_kw, x)) return("Europe")
    csa_kw <- paste(c("Pakistan","India","Bangladesh","Sri Lanka","Nepal",
      "Bhutan","Afghanistan","Tajikistan","Uzbekistan","Turkmenistan",
      "Kyrgyzstan","Kazakhstan","Jammu","Kashmir"), collapse = "|")
    if (grepl(csa_kw, x)) return("Central_South_Asia")
    ea_kw <- paste(c("China","Tibet","Japan","Korea","Mongolia","Vietnam",
      "Cambodia","Thailand","Laos","Myanmar","Malaysia","Indonesia",
      "Philippines","Taiwan","Singapore","Brunei"), collapse = "|")
    if (grepl(ea_kw, x)) return("East_Asia")
    oce_kw <- paste(c("Papua","Solomon","Vanuatu","Fiji","Tonga","Samoa",
      "Tuvalu","Tokelau","Niue","Guam","Mariana","Polynesia","Wallis",
      "Cook Islands","Australia","New Zealand","Nauru","Kiribati",
      "Marshall","Micronesia","Palau"), collapse = "|")
    if (grepl(oce_kw, x)) return("Oceania")
    ame_kw <- paste(c("United States","Canada","Mexico","Guatemala","Panama",
      "Colombia","Ecuador","Peru","Brazil","Chile","Argentina","Uruguay",
      "Paraguay","Bolivia","Venezuela","Puerto Rico","Barbados","Cuba",
      "Native American","Greenland","USA","Dominican Republic","Bahamas",
      "Belize","Saint Lucia","Curacao","Guadeloupe","Sint Maarten","Haiti"),
      collapse = "|")
    if (grepl(ame_kw, x)) return("Americas")
    # Default: treat remaining African countries as Sub-Saharan.
    ssa_kw <- paste(c("Madagascar","Zambia","Gambia","Angola","Kenya","Namibia",
      "Nigeria","Botswana","Burkina","South Africa","Ethiopia","Sierra Leone",
      "Tanzania","Comoros","Mozambique","Senegal","Central African","Cameroon",
      "Guinea","Chad","Somalia","Ghana","Mali","Niger","Congo","Uganda",
      "Rwanda","Burundi","Malawi","Zimbabwe","Eritrea","Djibouti","Benin",
      "Togo","Liberia","Mauritania","Gabon","Lesotho","Swaziland"),
      collapse = "|")
    if (grepl(ssa_kw, x)) return("Sub_Saharan_Africa")
    "Unknown"
  }
  out <- vapply(country, regionize, character(1))
  factor(out, levels = REGION_LEVELS)
}

# Broad Near-East grouping (Levant + Arabia/Mesopotamia + Anatolia + Caucasus).
NEAR_EAST_REGIONS <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus")
is_near_east <- function(region) as.character(region) %in% NEAR_EAST_REGIONS

# ---- Galton-Watson branching process (Livni & Skorecki eqs 1-3) ------------
# Offspring probability-generating function f1(s), i.e. the number of daughters
# per mother with mean = growth ratio m. Three offspring laws are supported:
#   * "poisson"   -- neutral model, variance = mean = m (Livni-Skorecki default);
#   * "geometric" -- strongly overdispersed, variance = m + m^2 (= NB size 1);
#   * "nbinom"    -- negative binomial with mean m and dispersion `size`
#                    (variance = m + m^2/size). Poisson is the size -> Inf
#                    limit and geometric is size = 1, so `size` continuously
#                    tunes reproductive overdispersion between the two.
# Returns Pr(offspring = j) for j = 0..jmax.
gw_offspring_pmf <- function(m, jmax = 400L,
                             model = c("poisson", "geometric", "nbinom"),
                             size = 1) {
  model <- match.arg(model)
  j <- 0:jmax
  if (model == "poisson") {
    p <- dpois(j, lambda = m)
  } else if (model == "geometric") {
    # Geometric with mean m: P(k) = (1-q) q^k, mean = q/(1-q) => q = m/(1+m).
    q <- m / (1 + m)
    p <- (1 - q) * q^j
  } else {
    # Negative binomial parameterised by mean m and dispersion `size`.
    p <- dnbinom(j, size = size, mu = m)
  }
  p / sum(p)
}

# Distribution of descendants after k generations, via repeated PGF composition
# evaluated on a complex grid (DFT) for speed and stability.
#
# `m` may be a scalar (constant growth ratio, the Livni-Skorecki assumption) or
# a length-k vector giving a *piecewise / time-varying* per-generation growth
# ratio m_t. In the time-varying case the offspring PGF f_t of generation t is
# composed in order, g_k(z) = f_1(f_2(...f_k(z)...)); the expected surviving
# lineage size after k generations is prod(m_t) regardless of ordering.
# ---- Founder-event timing --------------------------------------------------
# These were previously implicit: k = 15 was a bare default with no generation
# length anywhere in the repo, which placed the founder event at ~1650 CE --
# three centuries AFTER the 14th-century Erfurt carriers the analysis relies on,
# and inconsistent with the K = 25-32 grid used by estimate_absorption().
#
# GEN_YEARS: matrilineal generation interval. Fenner (2005, AJPA 128:415-423)
#   gives a cross-cultural female interval; Tremblay & Vezina (2000) report
#   ~29-30 yr for maternal lines and Wang et al. (2023) ~23-26 yr for recent
#   female generations. 26 yr is a defensible central value.
# FOUNDER_EVENT_YEAR: Carmi et al. (2014) date the Ashkenazi bottleneck to
#   600-800 years before present (25-32 generations); Waldman et al. (2022)
#   require it to pre-date the 14th-century Erfurt community.
GEN_YEARS          <- 26
PRESENT_YEAR       <- 2025
FOUNDER_EVENT_YEAR <- 1325
GW_K <- as.integer(round((PRESENT_YEAR - FOUNDER_EVENT_YEAR) / GEN_YEARS))  # 27

# Net maternal expansion from the founder event to the present. The Ashkenazi
# population grew from order 25,000 around 1325 CE to roughly 8 million by 1900
# (~320x over ~22 generations, m ~ 1.30 for that stretch) and was then flat
# across the 20th century. Spread over all GW_K generations that is an
# EFFECTIVE ratio of 320^(1/27) = 1.238.
#
# Every per-generation growth figure in the pipeline is derived from this one
# constant so they cannot drift apart: `default_growth_schedule` rescales its
# phase ratios to reproduce it, and `founder_event_params.csv` records the same
# effective m. (They previously disagreed by a factor of 23 in net expansion.)
NET_EXPANSION_TARGET <- 320
GW_M_EFFECTIVE <- NET_EXPANSION_TARGET^(1 / GW_K)   # 1.238

# Closed-form offspring PGFs. Evaluating these directly avoids truncating the
# offspring law at `jmax` before composing, which was a second, unnecessary
# approximation on top of the grid.
.gw_offspring_pgf <- function(z, m, model, size) {
  switch(model,
    poisson   = exp(m * (z - 1)),                      # Poisson(m)
    geometric = 1 / (1 + m - m * z),                   # P(j) = (1-q)q^j, q=m/(1+m)
    nbinom    = (1 + (m / size) * (1 - z))^(-size))    # NB(mu = m, size)
}

# Exact P(Z_k = 0), P(Z_k = 1), P(Z_k = 2) with no grid at all.
# For G_k(s) = f_1(f_2(...f_k(s)...)), write A_{k+1}(s) = s and A_j = f_j o A_{j+1}.
# Iterating a = A(0), b = A'(0), c = A''(0) down from j = k to 1 gives the three
# low-order coefficients exactly: p0 = a, p1 = b, p2 = c/2. This is the right
# tool for Table 1 and the detection argument, which only ever read these cells,
# and it is immune to the aliasing described in `gw_descendant_pmf`.
gw_low_probs <- function(m, k = GW_K,
                         model = c("poisson", "geometric", "nbinom"),
                         size = 1) {
  model <- match.arg(model)
  if (length(m) == 1L) m <- rep(m, k)
  if (length(m) != k) stop("`m` must have length 1 or k")
  f1 <- function(s, mu) switch(model,
          poisson   = mu * exp(mu * (s - 1)),
          geometric = mu / (1 + mu - mu * s)^2,
          nbinom    = mu * (1 + (mu / size) * (1 - s))^(-size - 1))
  f2 <- function(s, mu) switch(model,
          poisson   = mu^2 * exp(mu * (s - 1)),
          geometric = 2 * mu^2 / (1 + mu - mu * s)^3,
          nbinom    = mu^2 * (1 + 1 / size) * (1 + (mu / size) * (1 - s))^(-size - 2))
  a <- 0; b <- 1; cc <- 0
  for (j in k:1) {
    mu <- m[j]
    A  <- .gw_offspring_pgf(a, mu, model, size)
    B  <- f1(a, mu) * b
    C  <- f2(a, mu) * b^2 + f1(a, mu) * cc
    a <- A; b <- B; cc <- C
  }
  c(p0 = a, p1 = b, p2 = cc / 2, p_ge3 = 1 - a - b - cc / 2)
}

gw_descendant_pmf <- function(m, k = GW_K, jmax = NULL,
                              model = c("poisson", "geometric", "nbinom"),
                              size = 1) {
  model <- match.arg(model)
  if (length(m) == 1L) m <- rep(m, k)
  if (length(m) != k) stop("`m` must have length 1 or k")
  # The inverse FFT on N points does NOT truncate mass above the grid -- it
  # ALIASES it back onto low j (coeff_j picks up every c_{j + N*l}), i.e. onto
  # exactly the cells this analysis reads. The grid must therefore cover the
  # bulk of Z_k, whose mean is prod(m). At m = 1.075, k = 15 the old fixed
  # jmax = 400 was harmless, but at m = 1.4, k = 25 it inflated P(Z_k = 1) by 81x.
  if (is.null(jmax)) {
    N <- 2^ceiling(log2(max(512, 64 + 40 * prod(m))))
  } else {
    N <- jmax + 1L
    if (prod(m) > N / 20)
      warning(sprintf(
        "gw_descendant_pmf: E[Z_k] = %.0f against grid N = %d; low-order coefficients will be aliased",
        prod(m), N))
  }
  omega <- exp(-2i * pi * (0:(N - 1)) / N)
  # Compose innermost (last generation) first so g_k(z) = f_1(f_2(...f_k(z))).
  gz <- omega
  for (t in rev(seq_len(k))) gz <- .gw_offspring_pgf(gz, m[t], model, size)
  coeffs <- Re(fft(gz, inverse = TRUE) / N)
  coeffs[coeffs < 0] <- 0
  coeffs / sum(coeffs)
}

# Convenience: probability a founder has exactly one matrilineal descendant
# after k generations (reproduces Livni-Skorecki Table 1).
gw_prob_single <- function(m, k = GW_K, ...) {
  # Exact: no grid, so no aliasing (see gw_low_probs).
  unname(gw_low_probs(m, k = k, ...)[["p1"]])
}

# ---- Plot helper -----------------------------------------------------------
open_png <- function(name, width = 1100, height = 800, res = 130) {
  png(file.path(FIG_DIR, name), width = width, height = height, res = res)
}
save_table <- function(df, name) {
  utils::write.csv(df, file.path(TAB_DIR, name), row.names = FALSE)
  invisible(df)
}

.posterior_lineage <- function(tbl, lineage, col = "post_mean_H2",
                               id_cols = c("lineage", "founder")) {
  for (id in id_cols) {
    if (!id %in% names(tbl)) next
    hit <- tbl[[col]][match(lineage, tbl[[id]])]
    if (length(hit) && !is.na(hit)) return(hit)
  }
  NA_real_
}

# Posterior display: three decimals so 0.996 reads as 0.996, not 1.00.
.fmt_post <- function(x, digits = 3L) {
  if (!length(x)) return(character(0))
  fmt <- paste0("%.", digits, "f")
  vapply(x, function(p) {
    if (is.na(p)) return(NA_character_)
    sprintf(fmt, p)
  }, character(1))
}

PALETTE <- c(
  Levant = "#c1121f", Arabia_Mesopotamia = "#e5793a", Anatolia = "#f3a712",
  Caucasus = "#8a5a44", North_Africa = "#dda15e", Europe = "#264653",
  Central_South_Asia = "#2a9d8f", East_Asia = "#457b9d",
  Sub_Saharan_Africa = "#606c38", Oceania = "#6a4c93", Americas = "#9d4edd",
  Unknown = "#adb5bd", `Near East` = "#c1121f")
