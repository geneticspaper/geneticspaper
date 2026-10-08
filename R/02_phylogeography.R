# 02_phylogeography.R
# Module B -- phylogeography of the four major Ashkenazi founders using the
# Mitotree public sample set (S1) and tree topology (S10).
#
# Uses ape (tree construction / nesting plots) and pegas (haplotype network of
# reconstructed founder motifs). Addresses the Costa et al. (2013) claim that
# the founders nest within "solely European" sub-clades, and shows that the
# apparent European dominance of nesting sub-clades is a sample-size artifact
# (Livni-Skorecki section 3.3): rarefying to equal sample sizes removes it.

if (!exists("load_samples")) source(file.path("R", "utils.R"))

.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
suppressMessages({
  have_ape   <- requireNamespace("ape",   quietly = TRUE)
  have_pegas <- requireNamespace("pegas", quietly = TRUE)
})

# --------------------------------------------------------------------------
# B0. Topology helpers -- build parent/child maps and descendant sets.
# --------------------------------------------------------------------------
build_topology <- function(struct) {
  parent <- setNames(struct$ParentHaplogroup, struct$Haplogroup)
  tmrca  <- setNames(suppressWarnings(as.numeric(struct$TMRCA)), struct$Haplogroup)
  children <- split(struct$Haplogroup, struct$ParentHaplogroup)
  list(parent = parent, children = children, tmrca = tmrca,
       nodes = struct$Haplogroup)
}

descendants <- function(root, topo, include_self = TRUE) {
  out <- character(0)
  stack <- root
  while (length(stack)) {
    cur <- stack[[1]]; stack <- stack[-1]
    out <- c(out, cur)
    kids <- topo$children[[cur]]
    if (!is.null(kids)) stack <- c(stack, kids)
  }
  if (!include_self) out <- setdiff(out, root)
  unique(out)
}

# Cumulative mutation motif for a node (union of defining mutations up its spine).
motif_of <- function(node, struct_mut, parent) {
  muts <- character(0)
  cur <- node
  guard <- 0L
  while (!is.na(cur) && cur %in% names(parent) && guard < 200L) {
    m <- struct_mut[[cur]]
    if (!is.null(m) && !is.na(m) && nzchar(m) && m != ".")
      muts <- c(muts, strsplit(trimws(m), "\\s+")[[1]])
    nxt <- parent[[cur]]
    if (is.na(nxt) || nxt == cur || nxt == ".") break
    cur <- nxt; guard <- guard + 1L
  }
  unique(muts)
}

# --------------------------------------------------------------------------
# B1. Regional distribution of founders and their parent clades.
# --------------------------------------------------------------------------
regional_distribution <- function(samples, topo, clades) {
  modern <- samples[samples$SubjectType == "Modern", ]
  do.call(rbind, lapply(clades, function(cl) {
    desc <- descendants(cl, topo)
    sub <- modern[modern$MitotreeHaplogroup %in% desc, ]
    tab <- table(factor(sub$region, levels = REGION_LEVELS))
    data.frame(clade = cl, region = names(tab), n = as.integer(tab),
               row.names = NULL)
  }))
}

# --------------------------------------------------------------------------
# B2. Rarefaction: distinct sub-lineages vs sample size, per region.
#     Rebuts the Costa et al. sample-size confound (22,257 Eur vs 3,227 NE).
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# B0b. Resolution-robust rarefaction helpers.
# --------------------------------------------------------------------------
# Depth of `node` in edges below `root`; NA when `root` is not an ancestor.
depth_below <- function(node, root, topo, max_up = 200L) {
  cur <- node; d <- 0L
  while (!is.na(cur) && cur %in% names(topo$parent) && d <= max_up) {
    if (identical(cur, root)) return(d)
    p <- topo$parent[[cur]]
    if (is.na(p) || p == "." || p == cur) break
    cur <- p; d <- d + 1L
  }
  NA_integer_
}

# The ancestor of `node` exactly `k` edges below `root` (or `node` itself when
# it is shallower). NA when `root` is not an ancestor of `node`.
cap_to_level <- function(node, root, topo, k, max_up = 200L) {
  path <- character(0); cur <- node; g <- 0L
  repeat {
    if (is.na(cur) || !cur %in% names(topo$parent) || g > max_up)
      return(NA_character_)
    path <- c(path, cur)
    if (identical(cur, root)) break
    p <- topo$parent[[cur]]
    if (is.na(p) || p == "." || p == cur) return(NA_character_)
    cur <- p; g <- g + 1L
  }
  path <- rev(path)                       # root, ..., node
  path[[min(length(path), k + 1L)]]
}

# Exact Hurlbert (1971) rarefaction, E[S_n] = sum_i (1 - C(N-n_i, n)/C(N, n)).
# Replaces the Monte Carlo mean, which moved by up to 1.2 richness units across
# seeds at reps = 100 (sd 0.37 on the K1a European pool) -- larger than the
# between-region differences the channel is asked to adjudicate -- and whose
# set.seed() clobbered the caller's global RNG stream mid-pipeline.
rarefy_expected <- function(labels, n) {
  N <- length(labels)
  if (!N || is.na(n) || n < 1 || n > N) return(NA_real_)
  ni <- as.integer(table(labels))
  sum(1 - exp(lchoose(N - ni, n) - lchoose(N, n)))
}

# SE of `rarefy_expected`, treating the pool as ONE realisation of the regional
# population (Chao et al. 2014 / iNEXT bootstrap assemblage). The textbook Heck
# (1975) variance conditions on the pool and is identically 0 whenever
# n = pool size -- which, since the comparison runs at n = min(pool), is always
# the smaller Near Eastern pool, i.e. exactly the side the sign turns on.
rarefy_se <- function(labels, n, reps = 400L, seed = 1L) {
  N <- length(labels)
  if (!N || is.na(n) || n < 1 || n > N) return(NA_real_)
  ni <- as.integer(table(labels))
  f1 <- sum(ni == 1L); f2 <- sum(ni == 2L)
  chat <- if (f1 == 0L) 1 else
    1 - (f1 / N) * ((N - 1) * f1) / ((N - 1) * f1 + 2 * max(f2, 1))
  p <- ni / N
  den <- sum(p * (1 - p)^N)
  lam <- if (den > 0) (1 - chat) / den else 0
  padj <- p * (1 - lam * (1 - p)^N)
  f0 <- if (f2 > 0) ((N - 1) / N) * f1^2 / (2 * f2)
        else ((N - 1) / N) * f1 * (f1 - 1) / 2
  f0 <- max(0L, min(as.integer(ceiling(f0)), 10000L))
  if (f0 > 0L && chat < 1) padj <- c(padj, rep((1 - chat) / f0, f0))
  padj <- padj / sum(padj)
  old <- if (exists(".Random.seed", .GlobalEnv))
    get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  v <- vapply(seq_len(reps), function(i) {
    cnt <- stats::rmultinom(1L, N, padj)[, 1]; cnt <- cnt[cnt > 0]
    sum(1 - exp(lchoose(N - cnt, n) - lchoose(N, n)))
  }, numeric(1))
  stats::sd(v)
}

# Depth-matched rarefied richness: the resolution-robust replacement for
# "count distinct MitotreeHaplogroup".
#
# Mitotree does not name branches at a uniform depth by region -- it is built
# largely from European-ancestry consumer sequencing -- so a region whose
# samples sit deeper shows more distinct named nodes at the same sample size
# regardless of real diversity. Europe is the deeper-named side in 14 of the 16
# roots this channel uses (mean gap +0.28 levels, +1.19 under K1a), and the two
# exceptions are the two roots where the Near East is the better-sampled region.
#
# Pairing the two draws at random and capping BOTH members of each pair at
# min(depth_a, depth_b) edges below the root gives the two counted label sets
# identical depth distributions by construction. A single fixed cap is not
# portable: the roots sit at different depths and Mitotree has long single-child
# chains (U5 has 1 child, K has 3, K1a has 129), so "one level" is not a
# constant amount of evolution.
depth_matched_richness <- function(labels_a, labels_b, root, topo, n,
                                   reps = 400L, seed = 1L) {
  chain_of <- function(node) {
    path <- character(0); cur <- node; g <- 0L
    repeat {
      if (is.na(cur) || !cur %in% names(topo$parent) || g > 200L) return(NULL)
      path <- c(path, cur)
      if (identical(cur, root)) break
      p <- topo$parent[[cur]]
      if (is.na(p) || p == "." || p == cur) return(NULL)
      cur <- p; g <- g + 1L
    }
    rev(path)
  }
  uniq <- unique(c(labels_a, labels_b))
  chains <- lapply(uniq, chain_of)
  names(chains) <- uniq
  chains <- chains[!vapply(chains, is.null, logical(1))]
  if (!length(chains)) return(c(a = NA_real_, b = NA_real_, sd_a = NA_real_, sd_b = NA_real_))
  depth <- vapply(chains, length, integer(1)) - 1L
  maxd <- max(depth)
  CAP <- matrix(NA_character_, nrow = length(chains), ncol = maxd + 1L,
                dimnames = list(names(chains), NULL))
  for (i in seq_along(chains)) {
    ch <- chains[[i]]
    CAP[i, ] <- ch[pmin(seq_len(maxd + 1L), length(ch))]
  }
  ia_all <- match(labels_a, rownames(CAP)); ia_all <- ia_all[!is.na(ia_all)]
  ib_all <- match(labels_b, rownames(CAP)); ib_all <- ib_all[!is.na(ib_all)]
  if (!length(ia_all) || !length(ib_all) || is.na(n) ||
      n > min(length(ia_all), length(ib_all)))
    return(c(a = NA_real_, b = NA_real_, sd_a = NA_real_, sd_b = NA_real_))
  da <- depth[ia_all]; db <- depth[ib_all]
  old <- if (exists(".Random.seed", .GlobalEnv))
    get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  set.seed(seed)
  va <- numeric(reps); vb <- numeric(reps)
  for (r in seq_len(reps)) {
    pa <- sample.int(length(ia_all), n); pb <- sample.int(length(ib_all), n)
    # ONE cap level for the whole replicate, drawn from the paired-depth
    # distribution. Capping each PAIR at its own level is not a richness: two
    # records carrying the same label would be truncated to different depths,
    # so an ancestor and its own descendant get counted as two lineages, and
    # the result can exceed the number of distinct labels in the pool (the K1a
    # Near Eastern pool has 56 distinct labels and scored 56.20 at n = N, which
    # no rarefaction can do). A single level per replicate keeps both sides at
    # the same resolution AND keeps each replicate a genuine distinct count.
    lev <- min(da[pa[1L]], db[pb[1L]]) + 1L
    va[r] <- length(unique.default(CAP[cbind(ia_all[pa], lev)]))
    vb[r] <- length(unique.default(CAP[cbind(ib_all[pb], lev)]))
  }
  c(a = mean(va), b = mean(vb), sd_a = stats::sd(va), sd_b = stats::sd(vb))
}

# `exclude_clade` drops a focal clade's own carriers from both regional pools.
# Required whenever the question is "which region's variation does clade X nest
# inside?": X's own carriers are not evidence about where X nests, and because
# they concentrate on few named nodes they depress the richness of whichever
# region they are most common in. Left NULL for the standalone macro-clade
# rarefaction in Figure B6, which asks about the macro clade as a whole.
rarefy_sublineages <- function(samples, topo, root_clade,
                               regions = list(Europe = "Europe",
                                              `Near East` = NEAR_EAST_REGIONS),
                               n_grid = NULL, reps = 200L, seed = 1L,
                               exclude_clade = NULL,
                               statistic = c("depth_matched", "names", "haplotype"),
                               se_reps = 400L, match_reps = 400L) {
  statistic <- match.arg(statistic)
  modern <- samples[samples$SubjectType == "Modern", ]
  desc <- descendants(root_clade, topo)
  sub <- modern[modern$MitotreeHaplogroup %in% desc, ]
  if (!is.null(exclude_clade) && length(exclude_clade) &&
      !all(is.na(exclude_clade))) {
    drop <- unique(unlist(lapply(
      exclude_clade[!is.na(exclude_clade) & exclude_clade %in% topo$nodes],
      descendants, topo = topo)))
    if (length(drop)) sub <- sub[!(sub$MitotreeHaplogroup %in% drop), ]
  }
  value <- if (statistic == "haplotype") sub$Haplotype else sub$MitotreeHaplogroup
  idx <- lapply(regions, function(rg) which(as.character(sub$region) %in% rg))
  sizes <- vapply(idx, length, integer(1))
  # With every pool empty the old code did min(integer(0)) -> Inf and then
  # seq(1, Inf) -> error, which the caller's tryCatch swallowed into an empty
  # frame, hiding the failure. Return NULL explicitly.
  if (!any(sizes > 0)) return(NULL)
  nmax <- min(sizes[sizes > 0])
  if (is.null(n_grid)) n_grid <- unique(round(seq(1, nmax, length.out = 25)))
  n_grid <- n_grid[n_grid >= 1 & n_grid <= nmax]

  res <- do.call(rbind, lapply(names(idx), function(rg) {
    ii <- idx[[rg]]
    if (!length(ii)) return(NULL)
    pool <- value[ii]
    means <- vapply(n_grid, rarefy_expected, numeric(1), labels = pool)
    ses <- rep(NA_real_, length(n_grid))      # SE only at the reported n
    hit <- which(n_grid == nmax)
    if (length(hit)) ses[hit] <- rarefy_se(pool, nmax, se_reps, seed)
    data.frame(root = root_clade, region = rg, n = n_grid,
               distinct_lineages = means, se = ses,
               pool_size = length(pool), statistic = statistic,
               stringsAsFactors = FALSE)
  }))
  if (is.null(res) || !nrow(res)) return(NULL)

  # Depth matching needs the two pools jointly, so it is applied at n = nmax,
  # the only point the synthesis reads. The rest of the curve stays on raw node
  # counts for figures B6/B6b.
  # Only the n == nmax point is depth-matched, so label the rest honestly.
  res$statistic <- ifelse(res$n == nmax, statistic, "names")
  if (statistic == "depth_matched" && length(idx) == 2L && all(sizes > 0)) {
    dm <- depth_matched_richness(value[idx[[1]]], value[idx[[2]]],
                                 root_clade, topo, nmax,
                                 reps = match_reps, seed = seed)
    if (is.finite(dm[["a"]]) && is.finite(dm[["b"]])) {
      for (j in 1:2) {
        rw <- res$n == nmax & res$region == names(idx)[j]
        res$distinct_lineages[rw] <- dm[[c("a", "b")[j]]]
        # The bootstrap SE was computed for the RAW-name estimate; the reported
        # point estimate is the depth-matched one, so add its Monte Carlo SD in
        # quadrature rather than reporting a variance for a different estimator.
        res$se[rw] <- sqrt(res$se[rw]^2 + dm[[c("sd_a", "sd_b")[j]]]^2)
      }
    }
  }
  res
}

# --------------------------------------------------------------------------
# B3. ape: sibling-nesting tree for a founder within its parent clade.
# --------------------------------------------------------------------------
dominant_region <- function(clade, samples, topo) {
  desc <- descendants(clade, topo)
  sub <- samples[samples$SubjectType == "Modern" &
                   samples$MitotreeHaplogroup %in% desc, ]
  r <- as.character(sub$region)
  r <- r[r != "Unknown"]
  if (!length(r)) return("Unknown")
  names(sort(table(r), decreasing = TRUE))[1]
}

# Choose an informative ancestor clade to display siblings for: walk up from
# the founder (skipping intermediate 1-2 child "^"/combined nodes) until a clade
# with several children is reached.
choose_display_parent <- function(founder, topo, min_sibs = 3L, max_up = 5L) {
  p <- topo$parent[[founder]]; up <- 0L
  while (!is.na(p) && p != "." && up < max_up) {
    if (length(topo$children[[p]]) >= min_sibs) return(p)
    np <- topo$parent[[p]]
    if (is.na(np) || np == p || np == ".") break
    p <- np; up <- up + 1L
  }
  p
}

plot_sibling_nesting <- function(founder, samples, topo, file, parent = NULL) {
  if (!have_ape) return(invisible(NULL))
  if (is.null(parent)) parent <- choose_display_parent(founder, topo)
  sibs <- topo$children[[parent]]
  if (is.null(sibs) || length(sibs) < 2) return(invisible(NULL))
  page <- topo$tmrca[[parent]]
  if (is.na(page)) page <- max(topo$tmrca, na.rm = TRUE)
  info <- lapply(sibs, function(s) {
    age <- topo$tmrca[[s]]; if (is.na(age)) age <- page * 0.5
    bl  <- max(page - age, page * 0.02)
    dsc <- descendants(s, topo)
    reg <- dominant_region(s, samples, topo)
    n <- sum(samples$SubjectType == "Modern" & samples$MitotreeHaplogroup %in% dsc)
    # this sibling carries the founder if the founder nests within it
    is_line <- founder %in% dsc
    list(name = s, bl = bl, region = reg, n = n, is_line = is_line)
  })
  founder_line <- vapply(info, function(i) i$is_line, logical(1))
  labs <- vapply(info, function(i) sprintf("%s  [%s, n=%d]%s", i$name, i$region,
                 i$n, if (i$is_line) "  <- founder" else ""), "")
  bls  <- vapply(info, function(i) i$bl, numeric(1))
  regs <- vapply(info, function(i) i$region, "")
  clean <- gsub("'", "", labs)
  newick <- paste0("(", paste(sprintf("'%s':%f", clean, bls), collapse = ","),
                   ")", gsub("[^A-Za-z0-9]", "_", parent), ";")
  tr <- ape::read.tree(text = newick)
  ord <- match(gsub("'", "", tr$tip.label), clean)
  tip_reg <- regs[ord]
  tip_col <- PALETTE[tip_reg]; tip_col[is.na(tip_col)] <- PALETTE[["Unknown"]]
  is_f <- founder_line[ord]
  open_png(file, width = 1300, height = max(600, 24 * length(sibs) + 220))
  op <- par(mar = c(4, 1, 3, 1))
  ape::plot.phylo(tr, tip.color = tip_col, cex = 0.85, no.margin = FALSE,
                  font = ifelse(is_f, 2, 1),
                  main = sprintf("Nesting of %s among siblings of %s (Mitotree)",
                                 founder, parent))
  ape::tiplabels(pch = ifelse(is_f, 18, NA), col = "#c1121f", cex = 1.8)
  used <- intersect(names(PALETTE), unique(tip_reg))
  legend("bottomleft", legend = gsub("_", " ", used), col = PALETTE[used],
         pch = 15, bty = "n", cex = 0.8, title = "dominant region of clade")
  par(op); dev.off()
  invisible(data.frame(founder = founder, parent = parent,
                       sibling = vapply(info, `[[`, "", "name"), region = regs,
                       n = vapply(info, `[[`, 0L, "n"),
                       founder_line = founder_line, row.names = NULL))
}

# --------------------------------------------------------------------------
# B4. pegas: haplotype network of founder + parent motifs (reconstructed).
# --------------------------------------------------------------------------
founder_haplonet <- function(struct, topo, file) {
  if (!have_pegas || !have_ape) return(invisible(NULL))
  struct_mut <- setNames(struct$Mutations, struct$Haplogroup)
  nodes <- c("K1a", "K1a1b1", "K1a1b1a", "K1a9", "K2a", "K2a2a",
             "N1b", "N1b2", "K1a1b", "K1a1", "K2a2")
  nodes <- nodes[nodes %in% names(struct_mut)]
  motifs <- lapply(nodes, motif_of, struct_mut = struct_mut, parent = topo$parent)
  all_mut <- sort(unique(unlist(motifs)))
  if (length(all_mut) < 2) return(invisible(NULL))
  mat <- t(vapply(motifs, function(m) as.integer(all_mut %in% m), integer(length(all_mut))))
  rownames(mat) <- nodes
  # Encode as pseudo-DNA: ancestral 'a', derived 't'.
  seqmat <- ifelse(mat == 1L, "t", "a")
  db <- ape::as.DNAbin(seqmat)
  h <- pegas::haplotype(db)
  # Relabel haplotypes (roman numerals) with their haplogroup names.
  idx <- attr(h, "index")
  labs <- vapply(idx, function(i) paste(rownames(mat)[i], collapse = "/"), "")
  dimnames(h)[[1]] <- labs
  net <- pegas::haploNet(h)
  macro <- ifelse(grepl("^K1", labs), "K1 branch",
           ifelse(grepl("^K2", labs), "K2 branch", "N1b branch"))
  mcol <- c("K1 branch" = "#c1121f", "K2 branch" = "#e5793a",
            "N1b branch" = "#264653")
  open_png(file, width = 1100, height = 900)
  op <- par(mar = c(1, 1, 3, 1))
  plot(net, size = 2.4, show.mutation = 2, labels = TRUE, cex = 0.85,
       bg = mcol[macro], fast = FALSE,
       main = "Reconstructed haplotype network: K/N founders + parent clades")
  legend("bottomright", legend = names(mcol), pt.bg = mcol, pch = 21,
         pt.cex = 1.6, bty = "n")
  par(op); dev.off()
  invisible(data.frame(node = nodes, n_mutations = vapply(motifs, length, 0L)))
}

# --------------------------------------------------------------------------
# Runner
# --------------------------------------------------------------------------
run_phylogeography <- function() {
  message("== Module B: phylogeography (ape + pegas) ==")
  samples <- load_analysis_samples()
  struct  <- load_structure()
  topo    <- build_topology(struct)
  founders <- load_founders()

  clades <- unique(c(founders$founder, founders$parent_clade,
                     "K1a", "K1a1b1", "K2a", "N1b"))
  rd <- regional_distribution(samples, topo, clades)
  save_table(rd, "B_regional_distribution.csv")

  # ---- Figure B1: regional composition of parent clades (proportions) ----
  parent_clades <- c("K1a", "K1a1b1", "K2a", "N1b")
  keep_reg <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus",
                "North_Africa", "Europe")
  mat <- sapply(parent_clades, function(cl) {
    v <- rd$n[rd$clade == cl]; names(v) <- rd$region[rd$clade == cl]
    known <- v[keep_reg]; known[is.na(known)] <- 0
    known / sum(known)
  })
  open_png("B1_parentclade_region_props.png", width = 1350, height = 800)
  op <- par(mar = c(4.5, 5, 3, 12))
  barplot(mat, col = PALETTE[keep_reg], las = 1,
          ylab = "proportion of West-Eurasian / North-African samples",
          main = "Regional composition of founder parent-clades (Mitotree, modern)")
  legend(par("usr")[2] * 1.02, par("usr")[4], legend = gsub("_", " ", keep_reg),
         fill = PALETTE[keep_reg], xpd = NA, bty = "n", cex = 0.95)
  par(op); dev.off()

  # ---- Figure B2..B5: ape sibling-nesting trees per founder ----
  # Display parents chosen to give an informative sibling set (>=3 clades).
  nesting <- list()
  specs <- list(K1a1b1a = "K1a1b1", K1a9 = "K1a", K2a2a = NULL, N1b2 = "N1b")
  for (i in seq_along(specs)) {
    fnd <- names(specs)[i]
    f <- sprintf("B%d_nesting_%s.png", i + 1, fnd)
    nesting[[fnd]] <- plot_sibling_nesting(fnd, samples, topo, f,
                                           parent = specs[[i]])
  }
  nest_df <- do.call(rbind, nesting[!vapply(nesting, is.null, logical(1))])
  if (!is.null(nest_df)) save_table(nest_df, "B_nesting_siblings.csv")

  # ---- Figure B6: rarefaction (Costa sample-size rebuttal) ----
  # statistic = "names": these curves illustrate the raw sample-size argument
  # across the whole n grid, so every point must be the same statistic. Depth
  # matching is applied only at the single n the synthesis reads, so leaving it
  # on here would depth-match the last point and none of the others, putting a
  # spurious step at the right-hand end of the curve.
  rare <- do.call(rbind, lapply(c("K1a", "K2a", "N1b"),
                                rarefy_sublineages, samples = samples, topo = topo,
                                statistic = "names"))
  save_table(rare, "B_rarefaction.csv")
  open_png("B6_rarefaction.png", width = 1200, height = 800)
  op <- par(mar = c(4.5, 4.8, 3, 1))
  roots <- unique(rare$root)
  cols <- c(Europe = PALETTE[["Europe"]], `Near East` = PALETTE[["Levant"]])
  ltys <- setNames(seq_along(roots), roots)
  plot(NA, xlim = c(1, max(rare$n, na.rm = TRUE)),
       ylim = c(0, max(rare$distinct_lineages, na.rm = TRUE)),
       xlab = "samples drawn (rarefaction)", ylab = "distinct sub-lineages",
       main = "Sub-lineage richness at equal sample size: Europe vs Near East")
  for (rt in roots) for (rg in unique(rare$region)) {
    s <- rare[rare$root == rt & rare$region == rg, ]
    lines(s$n, s$distinct_lineages, col = cols[[rg]], lty = ltys[[rt]], lwd = 2)
  }
  legend("topleft", bty = "n", cex = 0.85,
         legend = c(paste("region:", names(cols)), paste("clade:", roots)),
         col = c(cols, rep("black", length(roots))),
         lty = c(1, 1, ltys), lwd = 2)
  par(op); dev.off()

  # ---- Figure B6b: rarefaction controls for the nesting/richness method ----
  # Keep these in a companion panel rather than adding them to B6: the founder
  # plot must remain readable. The selected roots have enough samples in both
  # Europe and the Near East to give interpretable equal-n curves. V and HV1b are
  # avoided here because their equal-sample caps are too small in this public
  # Mitotree subset, making rarefaction curves visually noisy.
  ctrl_roots <- data.frame(
    root = c("H", "U5", "J1c", "R0a", "U7", "U1"),
    expected = c(rep("European control", 3), rep("Near Eastern control", 3)),
    stringsAsFactors = FALSE
  )
  rare_ctrl <- do.call(rbind, lapply(ctrl_roots$root, rarefy_sublineages,
                                     samples = samples, topo = topo,
                                     statistic = "names"))
  rare_ctrl <- merge(rare_ctrl, ctrl_roots, by = "root", all.x = TRUE)
  save_table(rare_ctrl, "B_rarefaction_controls.csv")

  open_png("B6b_rarefaction_controls.png", width = 1400, height = 950)
  op <- par(mfrow = c(2, 3), mar = c(4, 4.4, 3, 1), oma = c(0, 0, 2, 0))
  ctrl_cols <- c(Europe = PALETTE[["Europe"]], `Near East` = PALETTE[["Levant"]])
  plot_ctrl_panel <- function(root) {
    d <- rare_ctrl[rare_ctrl$root == root, ]
    expectation <- unique(d$expected)
    plot(NA, xlim = c(1, max(d$n, na.rm = TRUE)),
         ylim = c(0, max(d$distinct_lineages, na.rm = TRUE)),
         xlab = "samples drawn (rarefaction)", ylab = "distinct sub-lineages",
         main = sprintf("%s (%s)", root, expectation))
    for (rg in unique(d$region)) {
      s <- d[d$region == rg, ]
      lines(s$n, s$distinct_lineages, col = ctrl_cols[[rg]], lwd = 2)
    }
    legend("topleft", bty = "n", cex = 0.78,
           legend = paste("region:", names(ctrl_cols)),
           col = ctrl_cols, lty = 1, lwd = 2)
  }
  for (rt in ctrl_roots$root) plot_ctrl_panel(rt)
  mtext("Rarefaction control panel: method moves in both directions", outer = TRUE,
        cex = 1.2, font = 2)
  par(op); dev.off()

  # ---- Figure B7: Table 7 non-Jewish rarity (founders; panel uses 1KG subclades) ----
  t7 <- load_table7()
  t7_livni <- t7[grepl("Livni", t7$source, fixed = TRUE), , drop = FALSE]
  kg1 <- tryCatch(load_1kg_clade_freq(), error = function(e) data.frame())
  open_png("B7_nonjewish_rarity.png", width = 1100, height = 780)
  op <- par(mar = c(4.5, 6, 3, 1))
  bp <- barplot(t7_livni$frequency * 100, names.arg = t7_livni$haplotype,
                horiz = TRUE, las = 1, col = "#264653",
                xlab = "frequency among non-Jews (%)",
                main = "Founder subclades: Livni-Skorecki Table 7 (n=27,651)")
  text(t7_livni$frequency * 100, bp,
       labels = sprintf("%.4f%%", t7_livni$frequency * 100),
       pos = 4, xpd = NA, cex = 0.85)
  mtext(sprintf("Panel/controls use 1000 Genomes EUR (CEU+GBR+FIN+IBS+TSI) subclade frequencies (n=%d); see D_lineage_computed_inputs.csv",
                if (nrow(kg1)) kg1$n_nonjews[1] else 503L),
        side = 1, line = 3.2, cex = 0.78)
  par(op); dev.off()

  # ---- pegas haplotype network ----
  hn <- founder_haplonet(struct, topo, "B8_founder_haplonet.png")
  if (!is.null(hn)) save_table(hn, "B_founder_motifs.csv")

  invisible(list(regional = rd, rarefaction = rare, nesting = nesting))
}

if (sys.nframe() == 0) run_phylogeography()
