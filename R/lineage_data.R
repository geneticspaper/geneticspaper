# lineage_data.R -- compute Bayesian synthesis channel inputs from the enriched
# Mitotree sample set and tree topology. Reference CSVs retain metadata only
# (published origin labels, Brook frequencies, notes); all numeric model inputs
# are derived at runtime from data.

if (!exists("load_analysis_samples")) source(file.path("R", "utils.R"))
if (!exists("build_topology")) source(file.path("R", "02_phylogeography.R"))

.lineage_ctx <- new.env(parent = emptyenv())

.lineage_context <- function() {
  if (!exists("ready", .lineage_ctx)) {
    samples <- load_analysis_samples()
    struct  <- load_structure()
    topo    <- build_topology(struct)
    targets <- read.csv(
      file.path(REFERENCE, "public_haplogroup_page_targets.csv"),
      stringsAsFactors = FALSE
    )
    .lineage_ctx$samples  <- samples
    .lineage_ctx$topo     <- topo
    .lineage_ctx$targets  <- targets
    .lineage_ctx$n_modern <- sum(samples$SubjectType == "Modern")
    .lineage_ctx$ready    <- TRUE
  }
  list(
    samples = .lineage_ctx$samples,
    topo = .lineage_ctx$topo,
    targets = .lineage_ctx$targets,
    n_modern = .lineage_ctx$n_modern
  )
}

# Macro clade used for rarefaction / ancient-depth lookups.
MACRO_ROOT <- c(
  K1a1b1a = "K1a", K1a9 = "K1a", K2a2a = "K2a", N1b2 = "N1b",
  V7a2c1b = "V", U5a1f1a3 = "U5",
  L2a1l2a = "L2a", M1a1b1c = "M1a", M33c3 = "M33", N9a3a1b1 = "N9a3",
  "A-a1b3a1" = "A",
  HV1b2 = "HV1b", R0a = "R0a", U7a5 = "U7", U1b1a1 = "U1",
  H7 = "H", H6a1a1a = "H6", J1c7a = "J1c",
  # Costa (2013) published-origin calibration lineages. These must be listed
  # explicitly: `infer_macro_root`'s fallback returns the lineage itself for any
  # key that is a tree node, so without an entry here the "macro" background
  # collapses onto the lineage and the Mitchell-2014 macro channel drops out.
  # Kept in step with the same map in scripts/compute_1kg_eu_freq.py.
  H1 = "H", H3 = "H", H5 = "H", HV1a = "HV1",
  T2 = "T2", W = "W", I = "I", U4 = "U4", K1a4a = "K1a4"
)

NON_EUROPE_REGIONS <- setdiff(
  REGION_LEVELS,
  c("Europe", "Unknown", "Oceania", "Americas")
)

# Studies / sites that are unambiguously Jewish ancient DNA. This is only a
# supplement to the per-sample `in_jewish_study` flag and must never list a
# general ancient-DNA compilation: Akbari et al. 2026 (the Allen Ancient DNA
# Resource) is such a compendium -- 279 of its samples appear here, none Jewish --
# so counting it as "Jewish" spuriously inflated medieval-carrier counts (e.g.
# T2 = 11 fake carriers) and is deliberately excluded.
JEWISH_ANCIENT_STUDIES <- paste(
  c("Waldman", "Diepenbroek", "Pallares", "Tarrega", "ROQ", "Erfurt",
    "Brace", "Chapelfield"),
  collapse = "|"
)

resolve_mitotree_node <- function(lineage_key) {
  ctx <- .lineage_context()
  topo <- ctx$topo
  row  <- ctx$targets[ctx$targets$lineage_key == lineage_key, , drop = FALSE]
  cands <- unique(c(lineage_key, row$ftdna_name, row$yfull_name))
  cands <- cands[nzchar(cands)]
  hit <- cands[cands %in% topo$nodes]
  if (length(hit)) {
    return(list(node = hit[1], match_method = "exact", candidates = cands))
  }
  best <- NA_character_; best_len <- 0L
  for (c in cands) {
    m <- topo$nodes[startsWith(topo$nodes, c) | startsWith(c, topo$nodes)]
    if (!length(m)) next
    pick <- m[which.max(nchar(m))]
    # Reject overly coarse prefix matches (e.g. macrohaplogroup A for A-a1b3a1).
    if (nchar(pick) < max(4L, nchar(c) - 1L)) next
    if (nchar(pick) > best_len) {
      best <- pick; best_len <- nchar(pick)
    }
  }
  if (!is.na(best)) {
    return(list(node = best, match_method = "prefix", candidates = cands))
  }
  list(node = NA_character_, match_method = "none", candidates = cands)
}

.parent_up <- function(node, topo) {
  if (is.na(node) || !node %in% names(topo$parent)) return(NA_character_)
  p <- topo$parent[[node]]
  if (is.na(p) || p == "." || p == node) NA_character_ else p
}

infer_macro_root <- function(lineage_key, node = NULL) {
  if (lineage_key %in% names(MACRO_ROOT)) return(MACRO_ROOT[[lineage_key]])
  ctx <- .lineage_context()
  topo <- ctx$topo
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  # Walk up to the first ancestor that is a tree node. The membership test must
  # come AFTER a parent step, or a lineage that is itself a node returns itself
  # and no walk ever happens.
  cur <- .parent_up(node, topo)
  for (i in seq_len(8L)) {
    if (is.na(cur)) break
    if (cur %in% topo$nodes) return(cur)
    cur <- .parent_up(cur, topo)
  }
  # Name-prefix fallback, specific -> coarse. Guarded so a key shorter than 3
  # characters does not turn the sequence around and try coarse -> specific.
  hi <- min(4L, nchar(lineage_key))
  if (hi >= 3L) {
    for (L in hi:3L) {
      p <- substr(lineage_key, 1L, L)
      if (p %in% topo$nodes) return(p)
    }
  }
  lineage_key
}

lineage_sample_mask <- function(lineage_key, node = NULL) {
  ctx <- .lineage_context()
  samples <- ctx$samples
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  if (!is.na(node)) {
    desc <- descendants(node, ctx$topo)
    return(samples$MitotreeHaplogroup %in% desc)
  }
  # No tree node: match haplogroup strings that start with the lineage key.
  # Respect mtDNA name boundaries -- a subclade name continues with the opposite
  # token type, so a digit right after a digit-ending key names a sibling, not a
  # descendant (H1 must not capture H10/H11).
  hg <- samples$MitotreeHaplogroup
  ok <- !is.na(hg) & startsWith(hg, lineage_key)
  nxt <- substr(hg, nchar(lineage_key) + 1L, nchar(lineage_key) + 1L)
  if (grepl("[0-9]$", lineage_key)) {
    ok <- ok & !(hg != lineage_key & grepl("^[0-9]$", nxt))
  }
  ok
}

FOUNDER_KEYS <- c("K1a1b1a", "K1a9", "K2a2a", "N1b2")

.noneu_diaspora_keys <- function() {
  tryCatch(load_noneuropean()$mtree_name, error = function(e) character(0))
}

# 1000 Genomes EUR benchmarks apply only where "absent among European non-Jews"
# is meaningful. Deep non-European diaspora lineages (L2a, M1a, …) are correctly
# 0% in EUR but that must not inflate the shared rarity z-score.
.use_1kg_eu_benchmark <- function(lineage_key) {
  if (frequency_sources()$nonjew != "1kg_eur") return(FALSE)
  if (lineage_key %in% FOUNDER_KEYS) return(FALSE)
  if (lineage_key %in% .noneu_diaspora_keys()) return(FALSE)
  TRUE
}

macro_europe_fraction <- function(macro_root) {
  ctx <- .lineage_context()
  topo <- ctx$topo
  if (!macro_root %in% topo$nodes) return(NA_real_)
  desc <- descendants(macro_root, topo)
  mod <- ctx$samples[ctx$samples$SubjectType == "Modern" &
                       ctx$samples$MitotreeHaplogroup %in% desc, ]
  known <- mod[!is.na(mod$region) & as.character(mod$region) != "Unknown", ]
  if (!nrow(known)) return(NA_real_)
  mean(as.character(known$region) == "Europe")
}

# Livni-Skorecki Table 7 founder subclades; Mitchell 2014 macro rows retained in
# table7_nonjewish_freq.csv for legacy livni_skorecki mode only.
.MACRO_TO_TABLE7 <- c(
  K1a = "K", K2a = "K",
  N1b = "N", N9a3 = "N",
  V = "V",
  U5 = "U", U7 = "U", U1 = "U",
  H = "H", H6 = "H",
  J1c = "J",
  HV1b = "HV",
  R0a = "R",
  L2a = "L2",
  M1a = "M", M33 = "M",
  A = "A"
)

.table7_aliases <- function(lineage_key) {
  keys <- lineage_key
  if (lineage_key == "N1b2") keys <- c(keys, "N1b1b1")
  ctx <- .lineage_context()
  row <- ctx$targets[ctx$targets$lineage_key == lineage_key, , drop = FALSE]
  if (nrow(row)) keys <- c(keys, row$ftdna_name, row$yfull_name)
  unique(keys[nzchar(keys)])
}

.table7_subclade_freq <- function(lineage_key) {
  T7 <- tryCatch(load_table7(), error = function(e) data.frame())
  if (!nrow(T7)) return(NA_real_)
  sub <- T7[grepl("Livni", T7$source, fixed = TRUE), , drop = FALSE]
  if (!nrow(sub)) sub <- T7[nchar(T7$haplotype) > 1L, , drop = FALSE]
  keys <- .table7_aliases(lineage_key)
  m <- vapply(sub$haplotype, function(h) any(vapply(keys, function(k)
    grepl(h, k, fixed = TRUE) || grepl(k, h, fixed = TRUE), logical(1))),
    logical(1))
  if (!any(m)) return(NA_real_)
  min(sub$frequency[m], na.rm = TRUE)
}

.table7_macro_freq <- function(macro_root) {
  T7 <- tryCatch(load_table7(), error = function(e) data.frame())
  if (!nrow(T7) || !macro_root %in% names(.MACRO_TO_TABLE7)) return(NA_real_)
  letter <- .MACRO_TO_TABLE7[[macro_root]]
  hit <- T7[T7$haplotype == letter &
              grepl("Mitchell", T7$source, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(hit)) return(NA_real_)
  hit$frequency[1]
}

.kg1_clade_freq_detail <- function(lineage_key) {
  kg <- tryCatch(load_1kg_clade_freq(), error = function(e) data.frame())
  if (!nrow(kg)) return(list(freq = NA_real_, clade = NA_character_, source = NA_character_))
  hit <- kg[kg$lineage_key == lineage_key, , drop = FALSE]
  if (!nrow(hit)) return(list(freq = NA_real_, clade = NA_character_, source = NA_character_))
  n <- hit$n_nonjews[1]
  cnt <- hit$count[1]
  # Unobserved in EUR (count=0) is not literal frequency 0 for log-rarity; use
  # a Jeffreys-style pseudocount so NE panel rows do not dominate z-scoring.
  freq <- if (cnt == 0L && is.finite(n) && n > 0L) {
    0.5 / (n + 1)
  } else {
    hit$frequency[1]
  }
  list(
    freq = freq,
    clade = hit$haplotype[1],
    source = "1kg_eur"
  )
}

.kg1_clade_freq <- function(lineage_key) {
  .kg1_clade_freq_detail(lineage_key)$freq
}

.published_nonjew_frequency <- function(lineage_key, node = NULL) {
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  if (lineage_key %in% FOUNDER_KEYS) {
    return(list(
      freq = .table7_subclade_freq(lineage_key),
      source = "livni_skorecki"
    ))
  }
  if (.use_1kg_eu_benchmark(lineage_key)) {
    kg <- .kg1_clade_freq_detail(lineage_key)
    if (is.finite(kg$freq)) {
      # The 1000G EUR panel (n=503) is the primary rarity benchmark, but it can
      # still be too small to contain some Near-Eastern subclades. When
      # best_clade_match could not find the lineage or its macro background and
      # instead climbed *past* the macro to an unrelated broader clade (e.g.
      # HV1b2 -> HV, counting HV0b/HV6/HV17a that are not HV1b), the reported
      # frequency is spurious and far too high. In that case only -- 1000G blind
      # to the macro -- fall back to the lineage's macro-background frequency
      # among European non-Jews in the much larger enriched Mitotree pool, which
      # resolves the macro directly. (With the wider EUR panel this now fires for
      # fewer lineages -- HV1b, U7a and R0a resolve directly at n=503.)
      macro <- infer_macro_root(lineage_key, node)
      climbed_past_macro <- !is.na(macro) && !is.na(kg$clade) &&
        startsWith(macro, kg$clade) && nchar(kg$clade) < nchar(macro)
      if (climbed_past_macro) {
        mm <- .mitotree_macro_nonjew_freq(lineage_key, node)
        if (is.finite(mm)) {
          return(list(freq = mm, source = "mitotree_macro_background"))
        }
      }
      return(list(freq = kg$freq, source = kg$source))
    }
  }
  macro <- infer_macro_root(lineage_key, node)
  t7_macro <- .table7_macro_freq(macro)
  t7_sub <- .table7_subclade_freq(lineage_key)
  mt <- .nonjew_freq_mitotree(lineage_key, node)
  vals <- c(t7_sub, t7_macro, mt)
  finite <- vals[is.finite(vals)]
  if (!length(finite)) return(list(freq = NA_real_, source = NA_character_))
  best <- min(finite)
  src_map <- c(
    livni_skorecki = is.finite(t7_sub) && t7_sub == best,
    mitchell2014 = is.finite(t7_macro) && t7_macro == best,
    mitotree = is.finite(mt) && mt == best
  )
  list(
    freq = best,
    source = paste(names(src_map)[src_map], collapse = "+")
  )
}

.nonjew_freq_mitotree <- function(lineage_key, node = NULL) {
  ctx <- .lineage_context()
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  mask <- lineage_sample_mask(lineage_key, node)
  modern <- ctx$samples[mask & ctx$samples$SubjectType == "Modern", ]
  if (!nrow(modern)) return(NA_real_)
  nj <- sum(
    !modern$in_jewish_study %in% c(TRUE, "TRUE") &
      as.character(modern$region) == "Europe",
    na.rm = TRUE
  )
  # The numerator counts European non-Jews, so the denominator must be the
  # European non-Jewish modern pool (~14.7k), NOT every modern sample worldwide.
  # Dividing by ctx$n_modern would put this estimator on a different scale from
  # the 1000G and Mitchell-2014 frequencies it is compared and `min()`-ed with
  # in `.nonjew_frequency_detail` -- see `.mitotree_macro_nonjew_freq`, which
  # computes the same kind of quantity on the correct denominator.
  all_modern <- ctx$samples[ctx$samples$SubjectType == "Modern", ]
  denom <- sum(
    !(all_modern$in_jewish_study %in% c(TRUE, "TRUE")) &
      as.character(all_modern$region) == "Europe",
    na.rm = TRUE
  )
  if (!denom) return(NA_real_)
  # Same Jeffreys-style floor the two sibling estimators use, so an absent
  # lineage stays readable on the log-rarity scale instead of being a literal 0.
  if (nj == 0L) return(0.5 / (denom + 1))
  nj / denom
}

# Frequency of a lineage's MACRO background among European non-Jews in the large
# enriched Mitotree modern pool (~14.7k European non-Jewish samples). Used only
# when the 1000G EUR panel cannot resolve the macro and would climb to
# an unrelated broader clade (see .published_nonjew_frequency). Expressed as a
# proper frequency among European non-Jews so it is on the same scale as the
# 1000G values it substitutes for. A Jeffreys floor keeps a macro that is absent
# even here readable on the log-rarity scale rather than a literal zero.
.mitotree_macro_nonjew_freq <- function(lineage_key, node = NULL) {
  ctx <- .lineage_context()
  macro <- infer_macro_root(lineage_key, node)
  if (is.na(macro) || !macro %in% ctx$topo$nodes) return(NA_real_)
  desc <- descendants(macro, ctx$topo)
  mod <- ctx$samples[ctx$samples$SubjectType == "Modern", ]
  eu_nonjew <- !(mod$in_jewish_study %in% c(TRUE, "TRUE")) &
    as.character(mod$region) == "Europe"
  denom <- sum(eu_nonjew, na.rm = TRUE)
  if (!denom) return(NA_real_)
  num <- sum(eu_nonjew & (mod$MitotreeHaplogroup %in% desc), na.rm = TRUE)
  if (num == 0L) return(0.5 / (denom + 1))
  num / denom
}

.nonjew_frequency_detail <- function(lineage_key, node = NULL) {
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  mt <- .nonjew_freq_mitotree(lineage_key, node)
  src <- frequency_sources()$nonjew
  if (src %in% c("livni_skorecki", "1kg_eur")) {
    pub <- .published_nonjew_frequency(lineage_key, node)
    if (is.finite(pub$freq)) return(pub)
    if (is.finite(mt)) return(list(freq = mt, source = "mitotree_fallback"))
    return(list(freq = NA_real_, source = NA_character_))
  }
  list(freq = mt, source = if (is.finite(mt)) "mitotree" else NA_character_)
}

nonjew_frequency <- function(lineage_key, node = NULL) {
  .nonjew_frequency_detail(lineage_key, node)$freq
}

.ref_frequency_detail <- function(lineage_key, node = NULL) {
  ctx <- .lineage_context()
  mask <- lineage_sample_mask(lineage_key, node)
  modern_n <- sum(mask & ctx$samples$SubjectType == "Modern", na.rm = TRUE)
  mt <- if (ctx$n_modern) modern_n / ctx$n_modern else NA_real_
  src <- frequency_sources()$ashkenazi
  if (src == "brook") {
    brook <- ashkenazi_frequency(lineage_key)
    if (is.finite(brook)) return(list(freq = brook, source = "brook"))
    if (is.finite(mt)) return(list(freq = mt, source = "mitotree_fallback"))
    return(list(freq = NA_real_, source = NA_character_))
  }
  list(freq = mt, source = if (is.finite(mt)) "mitotree" else NA_character_)
}

ref_frequency <- function(lineage_key, node = NULL) {
  .ref_frequency_detail(lineage_key, node)$freq
}

# Jewish ancient carriers of a clade in the AADR compilation, used where the
# Mitotree public ancient subset does not carry the label. Prefix matching
# respects mtDNA name boundaries (a digit after a digit-ending clade names a
# sibling, not a descendant).
.count_jewish_ancients_aadr <- function(clades) {
  p <- file.path(ROOT, "data", "external", "aadr_mtdna.csv")
  if (!file.exists(p)) return(0L)
  d <- tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(0L)
  yr_all <- suppressWarnings(as.numeric(d$year))
  d <- d[as.integer(d$is_jewish) %in% 1L &
           (is.na(yr_all) | yr_all <= MEDIEVAL_MAX_YEAR) &
           is.finite(suppressWarnings(as.numeric(d$date_bp))) &
           suppressWarnings(as.numeric(d$date_bp)) > 50, , drop = FALSE]
  if (!nrow(d)) return(0L)
  hg <- as.character(d$mt_haplogroup)
  hit <- rep(FALSE, length(hg))
  for (cl in clades[nzchar(clades)]) {
    ok <- !is.na(hg) & startsWith(hg, cl)
    nxt <- substr(hg, nchar(cl) + 1L, nchar(cl) + 1L)
    if (grepl("[0-9]$", cl)) ok <- ok & !(hg != cl & grepl("^[0-9]$", nxt))
    hit <- hit | ok
  }
  sum(hit)
}

# Latest calendar year that still counts as a historical (pre-modern) anchor.
# The channel exists to show a lineage was already in the Jewish maternal pool
# before the modern era, so 20th-century individuals carry no such evidence:
# without this cutoff K2a2a scored two "medieval" carriers that were in fact
# Sobibor Holocaust victims born in 1913 and 1923 (Diepenbroek 2021).
MEDIEVAL_MAX_YEAR <- 1800

.count_jewish_ancients <- function(mask) {
  ctx <- .lineage_context()
  anc <- ctx$samples[mask & ctx$samples$SubjectType == "Ancient", , drop = FALSE]
  if (!nrow(anc)) return(0L)
  jew <- anc$in_jewish_study %in% c(TRUE, "TRUE") |
    grepl(JEWISH_ANCIENT_STUDIES, anc$Study, ignore.case = TRUE)
  yr <- suppressWarnings(as.numeric(anc$AgeEstimateMean))
  sum(jew & (is.na(yr) | yr <= MEDIEVAL_MAX_YEAR), na.rm = TRUE)
}

# Historical N1b2 (Behar/Costa/23andMe) = N1b1b1 on FTDNA/Brook. Mitotree's ancient
# subset does not always map the Waldman 2022 Erfurt carrier onto the N1b2 label.
.n1b2_alias_keys <- function() {
  ctx <- .lineage_context()
  row <- ctx$targets[ctx$targets$lineage_key == "N1b2", , drop = FALSE]
  unique(c("N1b1b1", row$ftdna_name, row$yfull_name, "N1b2"))
}

medieval_jewish_carriers <- function(lineage_key, node = NULL) {
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  n <- .count_jewish_ancients(lineage_sample_mask(lineage_key, node))
  keys <- if (lineage_key == "N1b2") .n1b2_alias_keys() else lineage_key
  for (alias in setdiff(keys, lineage_key)) {
    alias_node <- resolve_mitotree_node(alias)$node
    if (!is.na(alias_node)) {
      n <- max(n, .count_jewish_ancients(lineage_sample_mask(alias, alias_node)))
    }
  }
  # Mitotree's public ancient subset carries the medieval Jewish individuals
  # unevenly -- it has the Erfurt K1a1b1a carriers but none of the K1a9 or
  # N1b1b1 ones -- so the count would otherwise depend on which lineage happened
  # to be included rather than on the record. Take the fuller of the two
  # sources for every lineage so the channel is measured the same way across
  # them, instead of hard-setting any single lineage's number.
  max(n, .count_jewish_ancients_aadr(keys))
}

lineage_ne_depth_kyr <- function(lineage_key, node = NULL, present_year = 1950L) {
  if (is.null(node)) node <- resolve_mitotree_node(lineage_key)$node
  ne_regions <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus")
  ctx <- .lineage_context()
  if (!is.na(node) && node %in% ctx$topo$nodes) {
    desc <- descendants(node, ctx$topo)
    anc <- ctx$samples[ctx$samples$SubjectType == "Ancient" &
                         ctx$samples$MitotreeHaplogroup %in% desc, ]
    anc <- anc[as.character(anc$region) %in% ne_regions &
                 !is.na(anc$AgeEstimateMean), ]
    if (nrow(anc)) {
      return((present_year - min(anc$AgeEstimateMean, na.rm = TRUE)) / 1000)
    }
  }
  oldest_ne_depth_kyr(infer_macro_root(lineage_key, node), present_year)
}

oldest_ne_depth_kyr <- function(macro_root, present_year = 1950L) {
  C_all <- tryCatch(
    read.csv(file.path(TAB_DIR, "C_ancient_records.csv"), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  ne_regions <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus")
  if (nrow(C_all)) {
    sub <- C_all[C_all$super_clade == macro_root &
                   C_all$region %in% ne_regions &
                   !is.na(C_all$year), ]
    if (nrow(sub)) return((present_year - min(sub$year)) / 1000)
  }
  ctx <- .lineage_context()
  topo <- ctx$topo
  if (!macro_root %in% topo$nodes) return(NA_real_)
  desc <- descendants(macro_root, topo)
  anc <- ctx$samples[ctx$samples$SubjectType == "Ancient" &
                       ctx$samples$MitotreeHaplogroup %in% desc, ]
  anc <- anc[as.character(anc$region) %in% ne_regions & !is.na(anc$AgeEstimateMean), ]
  if (!nrow(anc)) return(NA_real_)
  (present_year - min(anc$AgeEstimateMean, na.rm = TRUE)) / 1000
}

.rarefy_once <- function(root, regions, reps = 100L, seed = 1L,
                         exclude_clade = NULL) {
  ctx <- .lineage_context()
  tryCatch(
    rarefy_sublineages(ctx$samples, ctx$topo, root,
                       regions = regions, reps = reps, seed = seed,
                       exclude_clade = exclude_clade),
    error = function(e) data.frame()
  )
}

# Smallest regional pool a nesting value may rest on. Below this the equal-n
# comparison is being made at an n of a handful of samples, which is noise with
# a decisive-looking sign, so the channel is reported as NA and imputed neutral
# rather than being allowed to dominate the synthesis.
NEST_MIN_POOL <- 10L

# (Unused: the climb is bounded by the declared macro background instead, which
# is the principled frame. Kept only to document why. ) Without any bound, a lineage
# whose macro root equals the lineage itself has its whole root clade removed by
# `exclude_clade` and climbs until something is big enough: W, N9a3a1b1 and
# A-a1b3a1 all ended at node N+8701 with ~13,000 European and ~1,950 Near
# Eastern samples, so three unrelated panel lineages received the same
# macrohaplogroup-N number, which then set the location and scale of z_nest for
# the whole panel. Beyond this bound the comparison is no longer about the
# lineage, so the channel reports NA instead.
NEST_MAX_POOL <- 2000L

.rarefy_at_root <- function(root, mode = c("ne_eu", "noneu_eu"),
                           min_pool = NEST_MIN_POOL, max_up = 8L,
                           exclude_clade = NULL, ceiling_root = NULL) {
  mode <- match.arg(mode)
  ctx <- .lineage_context()
  regions <- if (mode == "ne_eu") {
    list(Europe = "Europe", `Near East` = NEAR_EAST_REGIONS)
  } else {
    list(Europe = "Europe", NonEurope = NON_EUROPE_REGIONS)
  }
  try_root <- root
  best <- NULL
  for (attempt in seq_len(max_up)) {
    r <- .rarefy_once(try_root, regions, exclude_clade = exclude_clade)
    if (is.null(r) || !nrow(r)) {
      # Same macro-background ceiling as below. This branch is the one that
      # fires when excluding the focal clade empties its own root entirely
      # (W, R0a and the other lineages whose macro root IS the lineage), so
      # without the check here they climb away regardless.
      if (!is.null(ceiling_root) && !is.na(ceiling_root) &&
          identical(try_root, ceiling_root)) break
      parent <- .parent_up(try_root, ctx$topo)
      if (is.na(parent)) break
      try_root <- parent
      next
    }
    best <- r
    pools <- tapply(r$pool_size, r$region, max, na.rm = TRUE)
    pos_lab <- if (mode == "ne_eu") "Near East" else "NonEurope"
    pos_pool_val <- if (pos_lab %in% names(pools)) pools[[pos_lab]] else NA_real_
    neg_pool_val <- if ("Europe" %in% names(pools)) pools[["Europe"]] else NA_real_
    if (is.finite(pos_pool_val) && is.finite(neg_pool_val) &&
        pos_pool_val >= min_pool && neg_pool_val >= min_pool) break
    # Never climb above the declared macro background: that is the frame this
    # lineage is defined against, and beyond it the comparison is about some
    # unrelated broader clade. W, N9a3a1b1 and A-a1b3a1 otherwise all ended at
    # node N+8701 and received the same macrohaplogroup-N number.
    if (!is.null(ceiling_root) && !is.na(ceiling_root) &&
        identical(try_root, ceiling_root)) break
    parent <- .parent_up(try_root, ctx$topo)
    if (is.na(parent) || parent == try_root) break
    try_root <- parent
  }
  if (is.null(best) || !nrow(best)) {
    return(list(
      root_used = ifelse(length(try_root) && nzchar(try_root), try_root, root),
      pos_richness = NA_real_, neg_richness = NA_real_,
      pos_se = NA_real_, neg_se = NA_real_,
      pos_pool = NA_real_, neg_pool = NA_real_, ne_richness = NA_real_,
      eu_richness = NA_real_, ne_pool = NA_real_, eu_pool = NA_real_,
      noneu_richness = NA_real_
    ))
  }
  pos_lab <- if (mode == "ne_eu") "Near East" else "NonEurope"
  pools <- tapply(best$pool_size, best$region, max, na.rm = TRUE)
  pos_pool <- if (pos_lab %in% names(pools)) pools[[pos_lab]] else NA_real_
  neg_pool <- if ("Europe" %in% names(pools)) pools[["Europe"]] else NA_real_
  pool_vals <- c(pos_pool, neg_pool)
  pool_vals <- pool_vals[is.finite(pool_vals) & pool_vals > 0]
  nmax <- if (length(pool_vals) >= 2L) min(pool_vals) else if (length(pool_vals)) pool_vals[1] else NA_real_
  w <- best[best$n == nmax, , drop = FALSE]
  getv <- function(rg, col) {
    v <- w[[col]][w$region == rg]
    if (!length(v)) return(NA_real_)
    mean(v, na.rm = TRUE)
  }
  pos_r <- getv(pos_lab, "distinct_lineages")
  neg_r <- getv("Europe", "distinct_lineages")
  pos_se <- if ("se" %in% names(w)) getv(pos_lab, "se") else NA_real_
  neg_se <- if ("se" %in% names(w)) getv("Europe", "se") else NA_real_
  list(
    # Report the root that actually produced `best`, not `try_root`: the climb
    # loop keeps advancing `try_root` after recording `best`, so on a loop that
    # ends by exhausting max_up (or on an empty rarefaction at the next level)
    # `try_root` is one level above the clade these richness numbers came from.
    root_used = as.character(best$root[1]),
    pos_richness = pos_r, neg_richness = neg_r,
    pos_se = pos_se, neg_se = neg_se,
    pos_pool = pos_pool, neg_pool = neg_pool,
    ne_richness = if (mode == "ne_eu") pos_r else NA_real_,
    eu_richness = neg_r,
    ne_pool = if (mode == "ne_eu") pos_pool else NA_real_,
    eu_pool = neg_pool,
    noneu_richness = if (mode == "noneu_eu") pos_r else NA_real_
  )
}

compute_lineage_inputs <- function(lineage_key,
                                   macro_root = NULL,
                                   nesting_mode = c("ne_eu", "noneu_eu"),
                                   rarefy_root = NULL) {
  nesting_mode <- match.arg(nesting_mode)
  res <- resolve_mitotree_node(lineage_key)
  node <- res$node
  if (is.null(macro_root)) macro_root <- infer_macro_root(lineage_key, node)
  rr_root <- macro_root
  if (!is.null(rarefy_root) && length(rarefy_root) == 1L &&
      !is.na(rarefy_root)) {
    rr_root <- rarefy_root
  }
  # The lineage's own carriers are excluded from both pools: nesting asks where
  # its SISTER lineages come from.
  # Exclude the lineage itself AND any sibling founder sharing its root clade:
  # leaving K1a1b1a's carriers in the pool while scoring K1a9 is the same
  # circularity one step removed.
  excl <- unique(c(node, vapply(setdiff(FOUNDER_KEYS, lineage_key),
                                function(k) resolve_mitotree_node(k)$node,
                                character(1))))
  excl <- excl[!is.na(excl)]
  rr <- .rarefy_at_root(rr_root, mode = nesting_mode, exclude_clade = excl,
                        ceiling_root = macro_root)
  macro_eu_frac <- macro_europe_fraction(macro_root)
  sub_nj <- .nonjew_freq_mitotree(lineage_key, node)
  kg_d <- .kg1_clade_freq_detail(lineage_key)
  macro_nj <- .table7_macro_freq(macro_root)
  ref_d <- .ref_frequency_detail(lineage_key, node)
  nj_d <- .nonjew_frequency_detail(lineage_key, node)
  data.frame(
    lineage_key = lineage_key,
    mitotree_node = node,
    match_method = res$match_method,
    macro_root = macro_root,
    macro_europe_fraction = macro_eu_frac,
    rarefy_root = rr$root_used,
    ref_freq = ref_d$freq,
    ref_freq_source = ref_d$source,
    nonjew_freq_subclade = sub_nj,
    nonjew_freq_1kg_clade = kg_d$freq,
    nonjew_freq_1kg_match = kg_d$clade,
    nonjew_freq_mitchell_macro = macro_nj,
    nonjew_freq = nj_d$freq,
    nonjew_freq_source = nj_d$source,
    ne_depth_kyr = lineage_ne_depth_kyr(lineage_key, node),
    medieval_jewish_carriers = medieval_jewish_carriers(lineage_key, node),
    ne_richness = rr$ne_richness,
    eu_richness = rr$eu_richness,
    noneu_richness = rr$noneu_richness,
    ne_pool = rr$ne_pool,
    eu_pool = rr$eu_pool,
    pos_richness = rr$pos_richness,
    neg_richness = rr$neg_richness,
    pos_pool = rr$pos_pool,
    neg_pool = rr$neg_pool,
    pos_se = rr$pos_se,
    neg_se = rr$neg_se,
    modern_count = sum(lineage_sample_mask(lineage_key, node) &
                         .lineage_context()$samples$SubjectType == "Modern"),
    data_complete = !is.na(node) | sum(lineage_sample_mask(lineage_key, node)) > 0,
    stringsAsFactors = FALSE
  )
}

compute_lineage_inputs_batch <- function(lineage_keys,
                                         nesting_modes = NULL,
                                         rarefy_roots = NULL) {
  if (is.null(nesting_modes)) {
    nesting_modes <- rep("ne_eu", length(lineage_keys))
  } else if (length(nesting_modes) == 1L) {
    nesting_modes <- rep(nesting_modes, length(lineage_keys))
  }
  if (is.null(rarefy_roots)) {
    rarefy_roots <- rep(NA_character_, length(lineage_keys))
  } else if (length(rarefy_roots) == 1L) {
    rarefy_roots <- rep(rarefy_roots, length(lineage_keys))
  }
  out <- do.call(rbind, lapply(seq_along(lineage_keys), function(i) {
    compute_lineage_inputs(
      lineage_keys[i],
      nesting_mode = nesting_modes[i],
      rarefy_root = rarefy_roots[i]
    )
  }))
  # (A `founder_strength` min-max rescale of ref_freq used to be added here. It
  # was normalised over whichever batch happened to be passed, so the same
  # lineage got a different value depending on the call, and nothing read it.)
  out
}

# Nesting channel: Near-East vs Europe (or non-Europe vs Europe). Always use
# (pos - neg) / (pos + neg) with positive = Near-East / non-Europe richer at equal n,
# matching the ridge-logistic label (H2 = 1). Do not flip sign for European macros:
# a single z_nest coefficient is fit across all lineages.
nesting_fraction <- function(row) {
  if (!is.na(row$noneu_richness)) {
    pos <- row$noneu_richness; neg <- row$eu_richness
    pos_pool <- row$pos_pool
  } else {
    pos <- row$ne_richness; neg <- row$eu_richness
    pos_pool <- row$ne_pool
  }
  neg_pool <- row$eu_pool
  pos_se <- if ("pos_se" %in% names(row)) row$pos_se else NA_real_
  neg_se <- if ("neg_se" %in% names(row)) row$neg_se else NA_real_
  # An equal-n comparison is only as good as its smaller pool. A MISSING pool
  # must also return NA: the earlier guard was skipped when a pool was NA rather
  # than 0, and `pos` was then coerced to 0, so a lineage with no Near Eastern
  # samples anywhere up its spine scored (0 - neg)/(0 + neg) = -1, the strongest
  # possible European signal, instead of "not measurable".
  if (is.na(pos_pool) || is.na(neg_pool) || is.na(pos) || is.na(neg) ||
      min(pos_pool, neg_pool) < NEST_MIN_POOL) {
    return(NA_real_)
  }
  den <- pos + neg
  if (den == 0) return(0)
  frac <- (pos - neg) / den
  # Shrink toward zero by the channel's own uncertainty, so a near-tie does not
  # enter the synthesis with the same weight as a clear separation. Smooth,
  # bounded, and a no-op when |frac| is large relative to its SE.
  if (!is.na(pos_se) && !is.na(neg_se)) {
    se <- 2 / den^2 * sqrt(neg^2 * pos_se^2 + pos^2 * neg_se^2)
    if (is.finite(se) && se > 0) frac <- frac * frac^2 / (frac^2 + se^2)
  }
  frac
}

nesting_channels <- function(row) {
  if (!is.na(row$noneu_richness)) {
    return(list(pos = row$noneu_richness, neg = row$eu_richness,
                pos_pool = row$pos_pool, neg_pool = row$neg_pool,
                nest_frac = nesting_fraction(row)))
  }
  list(pos = row$ne_richness, neg = row$eu_richness,
       pos_pool = row$ne_pool, neg_pool = row$eu_pool,
       nest_frac = nesting_fraction(row))
}
