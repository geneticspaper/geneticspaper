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
rarefy_sublineages <- function(samples, topo, root_clade,
                               regions = list(Europe = "Europe",
                                              `Near East` = NEAR_EAST_REGIONS),
                               n_grid = NULL, reps = 200L, seed = 1L) {
  set.seed(seed)
  modern <- samples[samples$SubjectType == "Modern", ]
  desc <- descendants(root_clade, topo)
  sub <- modern[modern$MitotreeHaplogroup %in% desc, ]
  pools <- lapply(regions, function(rg) {
    sub$MitotreeHaplogroup[as.character(sub$region) %in% rg]
  })
  sizes <- vapply(pools, length, integer(1))
  nmax <- min(sizes[sizes > 0])
  if (is.null(n_grid)) {
    n_grid <- unique(round(seq(1, nmax, length.out = 25)))
  }
  res <- do.call(rbind, lapply(names(pools), function(rg) {
    pool <- pools[[rg]]
    if (!length(pool)) return(NULL)
    means <- vapply(n_grid, function(n) {
      if (n > length(pool)) return(NA_real_)
      mean(replicate(reps, length(unique(sample(pool, n)))))
    }, numeric(1))
    data.frame(root = root_clade, region = rg, n = n_grid,
               distinct_lineages = means, pool_size = length(pool))
  }))
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
          ylab = "proportion of known-origin samples",
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
  rare <- do.call(rbind, lapply(c("K1a", "K2a", "N1b"),
                                rarefy_sublineages, samples = samples, topo = topo))
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
                                     samples = samples, topo = topo))
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
