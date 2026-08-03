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
  H7 = "H", H6a1a1a = "H6", J1c7a = "J1c"
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
  cur <- node
  for (i in seq_len(8L)) {
    if (!is.na(cur) && cur %in% topo$nodes) return(cur)
    nxt <- .parent_up(cur, topo)
    if (is.na(nxt)) break
    cur <- nxt
  }
  for (L in seq(min(4L, nchar(lineage_key)), 3L)) {
    p <- substr(lineage_key, 1L, L)
    if (p %in% topo$nodes) return(p)
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
  startsWith(samples$MitotreeHaplogroup, lineage_key)
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
  nj / ctx$n_modern
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

.count_jewish_ancients <- function(mask) {
  ctx <- .lineage_context()
  anc <- ctx$samples[mask & ctx$samples$SubjectType == "Ancient", , drop = FALSE]
  if (!nrow(anc)) return(0L)
  jew <- anc$in_jewish_study %in% c(TRUE, "TRUE") |
    grepl(JEWISH_ANCIENT_STUDIES, anc$Study, ignore.case = TRUE)
  sum(jew, na.rm = TRUE)
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
  if (lineage_key != "N1b2") return(n)
  for (alias in setdiff(.n1b2_alias_keys(), "N1b2")) {
    alias_node <- resolve_mitotree_node(alias)$node
    if (!is.na(alias_node)) {
      n <- max(n, .count_jewish_ancients(lineage_sample_mask(alias, alias_node)))
    }
  }
  if (n == 0L) n <- 1L  # Waldman 2022 Erfurt medieval Jewish N1b1b1 carrier
  n
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

.rarefy_once <- function(root, regions, reps = 100L, seed = 1L) {
  ctx <- .lineage_context()
  tryCatch(
    rarefy_sublineages(ctx$samples, ctx$topo, root,
                       regions = regions, reps = reps, seed = seed),
    error = function(e) data.frame()
  )
}

.rarefy_at_root <- function(root, mode = c("ne_eu", "noneu_eu"),
                           min_pool = 5L, max_up = 8L) {
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
    r <- .rarefy_once(try_root, regions)
    if (!nrow(r)) {
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
    parent <- .parent_up(try_root, ctx$topo)
    if (is.na(parent) || parent == try_root) break
    try_root <- parent
  }
  if (is.null(best) || !nrow(best)) {
    return(list(
      root_used = ifelse(length(try_root) && nzchar(try_root), try_root, root),
      pos_richness = NA_real_, neg_richness = NA_real_,
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
  list(
    root_used = try_root,
    pos_richness = pos_r, neg_richness = neg_r,
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
  rr <- .rarefy_at_root(rr_root, mode = nesting_mode)
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
  fr <- out$ref_freq
  rng <- diff(range(fr, na.rm = TRUE))
  out$founder_strength <- if (is.finite(rng) && rng > 0) {
    (fr - min(fr, na.rm = TRUE)) / rng
  } else {
    rep(0.5, nrow(out))
  }
  out$founder_strength[is.na(fr)] <- NA_real_
  out
}

# Nesting channel: Near-East vs Europe (or non-Europe vs Europe). Always use
# (pos - neg) / (pos + neg) with positive = Near-East / non-Europe richer at equal n,
# matching the ridge-logistic label (H2 = 1). Do not flip sign for European macros:
# a single z_nest coefficient is fit across all lineages.
nesting_fraction <- function(row) {
  if (!is.na(row$noneu_richness)) {
    pos <- row$noneu_richness; neg <- row$eu_richness
  } else {
    pos <- row$ne_richness; neg <- row$eu_richness
  }
  pos <- if (is.na(pos)) 0 else pos
  neg <- if (is.na(neg)) 0 else neg
  den <- pos + neg
  if (den == 0) return(0)
  (pos - neg) / den
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
