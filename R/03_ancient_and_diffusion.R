# 03_ancient_and_diffusion.R
# Module C -- ancient DNA, coalescence ages, and the diffusion-law frequency
# argument for a Near Eastern / Levantine origin of the four major founders.
#
# Uses ape (TMRCA handling) and adegenet (glPca of reconstructed founder motifs).
#   * Ancient-DNA timeline of K1a / K2a / N1b by region and calendar age.
#   * Coalescence (TMRCA) ages of founders vs parent clades (Mitotree S10).
#   * Diffusion law (Fisher): parent clades peak in the Near East / N Africa,
#     so gene flow (high -> low concentration) points to a Near Eastern source.

if (!exists("load_samples")) source(file.path("R", "utils.R"))
if (!exists("build_topology")) source(file.path("R", "02_phylogeography.R"))

.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
suppressMessages({
  have_ape      <- requireNamespace("ape",      quietly = TRUE)
  have_adegenet <- requireNamespace("adegenet", quietly = TRUE)
})

# --------------------------------------------------------------------------
# C1. Ancient-DNA record for the founder super-clades.
# --------------------------------------------------------------------------
ancient_records <- function(samples, topo, roots = c("K1a", "K2a", "N1b")) {
  anc <- samples[samples$SubjectType == "Ancient", ]
  do.call(rbind, lapply(roots, function(rt) {
    desc <- descendants(rt, topo)
    sub <- anc[anc$MitotreeHaplogroup %in% desc & !is.na(anc$AgeEstimateMean), ]
    if (!nrow(sub)) return(NULL)
    data.frame(super_clade = rt,
               haplogroup = sub$MitotreeHaplogroup,
               country = sub$Country,
               region = as.character(sub$region),
               year = sub$AgeEstimateMean,
               row.names = NULL)
  }))
}

# --------------------------------------------------------------------------
# C2. Coalescence (TMRCA) ages of founders and parent clades.
# --------------------------------------------------------------------------
coalescence_ages <- function(struct) {
  keep <- c("K1a", "K1a1b", "K1a1b1", "K1a1b1a", "K1a9",
            "K2a", "K2a2a", "N1b", "N1b2")
  s <- struct[struct$Haplogroup %in% keep, ]
  data.frame(haplogroup = s$Haplogroup,
             tmrca = suppressWarnings(as.numeric(s$TMRCA)),
             tmrca_lo = suppressWarnings(as.numeric(s$TMRCA_95l)),
             tmrca_hi = suppressWarnings(as.numeric(s$TMRCA_95u)),
             ntips = suppressWarnings(as.integer(s$Ntips)))
}

# --------------------------------------------------------------------------
# C3. adegenet glPca of reconstructed founder / parent motifs.
# --------------------------------------------------------------------------
founder_glpca <- function(struct, topo, file) {
  if (!have_adegenet) return(invisible(NULL))
  struct_mut <- setNames(struct$Mutations, struct$Haplogroup)
  nodes <- c("K1a", "K1a1b", "K1a1b1", "K1a1b1a", "K1a9",
             "K2a", "K2a2", "K2a2a", "N1b", "N1b2")
  nodes <- nodes[nodes %in% names(struct_mut)]
  motifs <- lapply(nodes, motif_of, struct_mut = struct_mut, parent = topo$parent)
  all_mut <- sort(unique(unlist(motifs)))
  mat <- t(vapply(motifs, function(m) as.integer(all_mut %in% m),
                  integer(length(all_mut))))
  rownames(mat) <- nodes
  gl <- methods::new("genlight", mat, ploidy = 1)
  pca <- adegenet::glPca(gl, nf = 2, parallel = FALSE)
  scores <- as.data.frame(pca$scores)
  scores$node <- nodes
  macro <- ifelse(grepl("^K1", nodes), "K1 branch",
           ifelse(grepl("^K2", nodes), "K2 branch", "N1b branch"))
  mcol <- c("K1 branch" = "#c1121f", "K2 branch" = "#e5793a",
            "N1b branch" = "#264653")
  eig <- pca$eig / sum(pca$eig) * 100
  open_png(file, width = 1050, height = 850)
  op <- par(mar = c(4.5, 4.8, 3, 1))
  plot(scores$PC1, scores$PC2, pch = 21, bg = mcol[macro], cex = 2.2,
       xlab = sprintf("PC1 (%.1f%%)", eig[1]),
       ylab = sprintf("PC2 (%.1f%%)", eig[2]),
       main = "adegenet glPca of founder + parent-clade motifs")
  text(scores$PC1, scores$PC2, labels = nodes, pos = 3, cex = 0.85)
  legend("topright", legend = names(mcol), pt.bg = mcol, pch = 21,
         pt.cex = 1.6, bty = "n")
  par(op); dev.off()
  invisible(scores)
}

# Regional frequencies for control clades using the same broad buckets as the
# Livni-Skorecki Table 3 display. These are not part of Table 3 itself; they are
# a calibration companion showing that the frequency-gradient method can return
# European or Near Eastern profiles depending on the clade.
regional_control_frequencies <- function(samples, topo) {
  modern <- samples[samples$SubjectType == "Modern", ]
  region_sets <- list(
    Europe = "Europe",
    `North Africa` = "North_Africa",
    Caucasus = "Caucasus",
    `Middle East` = c("Levant", "Arabia_Mesopotamia"),
    Anatolia = "Anatolia"
  )
  controls <- data.frame(
    clade = c("H", "H6a1a", "U5", "J1c", "R0a", "U7", "U1"),
    expected = c("European control", "ambiguous H6 context",
                 "European control", "European control",
                 rep("Near Eastern control", 3)),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, lapply(seq_len(nrow(controls)), function(i) {
    cl <- controls$clade[i]
    desc <- descendants(cl, topo)
    do.call(rbind, lapply(names(region_sets), function(rg) {
      pool <- modern[modern$region %in% region_sets[[rg]], ]
      n_total <- nrow(pool)
      n_clade <- sum(pool$MitotreeHaplogroup %in% desc)
      data.frame(clade = cl, expected = controls$expected[i], region = rg,
                 n_clade = n_clade, n_total = n_total,
                 frequency = if (n_total) n_clade / n_total else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
  out
}

# --------------------------------------------------------------------------
# Runner
# --------------------------------------------------------------------------
run_ancient_and_diffusion <- function() {
  message("== Module C: ancient DNA + diffusion (ape + adegenet) ==")
  samples <- load_analysis_samples()
  struct  <- load_structure()
  topo    <- build_topology(struct)

  # ---- C1: ancient DNA timeline ----
  anc <- ancient_records(samples, topo)
  save_table(anc, "C_ancient_records.csv")
  focus_reg <- c("Levant", "Arabia_Mesopotamia", "Anatolia", "Caucasus",
                 "North_Africa", "Europe")
  anc$region[!(anc$region %in% focus_reg)] <- "Other/Unknown"
  reg_levels <- c(focus_reg, "Other/Unknown")
  anc$ry <- factor(anc$region, levels = rev(reg_levels))
  open_png("C1_ancient_timeline.png", width = 1300, height = 850)
  op <- par(mar = c(4.5, 8.5, 3, 1))
  pch_map <- c(K1a = 19, K2a = 17, N1b = 15)
  col_map <- c(K1a = "#c1121f", K2a = "#e5793a", N1b = "#264653")
  plot(anc$year, jitter(as.integer(anc$ry), 0.6), pch = pch_map[anc$super_clade],
       col = col_map[anc$super_clade], cex = 1,
       yaxt = "n", ylab = "", xlab = "calendar year (negative = BCE)",
       main = "Ancient DNA: K1a / K2a / N1b through time and space")
  axis(2, at = seq_along(reg_levels), labels = rev(gsub("_", " ", reg_levels)),
       las = 1)
  abline(v = 0, lty = 3, col = "grey60"); abline(h = seq_along(reg_levels),
         col = "grey92")
  legend("topleft", legend = names(pch_map), pch = pch_map, col = col_map,
         bty = "n", pt.cex = 1.2)
  par(op); dev.off()

  # Highlight table: oldest Near Eastern occurrences + founder footprint.
  ne <- anc[anc$region %in% c("Levant", "Arabia_Mesopotamia", "Anatolia",
                              "Caucasus"), ]
  # Keep only the reported columns; `anc` also carries a `ry` plotting factor
  # (a duplicate of region used for the timeline y-axis) that must not leak into
  # the saved table.
  ne_cols <- c("super_clade", "haplogroup", "country", "region", "year")
  oldest_ne <- ne[order(ne$year), ne_cols][seq_len(min(12, nrow(ne))), ]
  save_table(oldest_ne, "C_oldest_near_east.csv")
  k1a1b1a_anc <- anc[grepl("^K1a1b1a", anc$haplogroup), ]
  save_table(k1a1b1a_anc, "C_k1a1b1a_ancient_footprint.csv")

  # ---- C2: coalescence ages ----
  ca <- coalescence_ages(struct)
  save_table(ca, "C_coalescence_ages.csv")
  ord <- order(ca$tmrca)
  ca2 <- ca[ord, ]
  is_founder <- ca2$haplogroup %in% c("K1a1b1a", "K1a9", "K2a2a", "N1b2")
  open_png("C2_coalescence_ages.png", width = 1150, height = 800)
  op <- par(mar = c(4.5, 7, 3, 1))
  y <- seq_len(nrow(ca2))
  plot(ca2$tmrca, y, xlim = range(c(ca2$tmrca_lo, ca2$tmrca_hi), na.rm = TRUE),
       pch = 19, col = ifelse(is_founder, "#c1121f", "#264653"), cex = 1.3,
       yaxt = "n", ylab = "", xlab = "TMRCA (years before present)",
       main = "Coalescence ages: founders (red) vs parent clades")
  segments(ca2$tmrca_lo, y, ca2$tmrca_hi, y,
           col = ifelse(is_founder, "#c1121f", "#264653"))
  axis(2, at = y, labels = ca2$haplogroup, las = 1)
  par(op); dev.off()

  # ---- C3: diffusion-law regional frequencies (Table 3) ----
  t3 <- load_table3()
  save_table(t3, "C_table3_regional_freq.csv")
  m <- as.matrix(t3[, c("K1a", "K1a1b1", "K2a")]) * 100
  rownames(m) <- t3$region
  open_png("C3_diffusion_frequencies.png", width = 1200, height = 800)
  op <- par(mar = c(6.5, 4.8, 3, 8), xpd = NA)
  cols <- c(K1a = "#c1121f", K1a1b1 = "#f3a712", K2a = "#457b9d")
  bp <- barplot(t(m), beside = TRUE, col = cols, las = 2,
                ylab = "haplotype frequency (%)",
                main = "Parent-signature frequencies by region (Livni-Skorecki Table 3)")
  legend(par("usr")[2] * 1.01, par("usr")[4], legend = colnames(m),
         fill = cols, bty = "n")
  mtext("Diffusion law: gene flow runs high -> low concentration (NE/N.Africa -> Europe)",
        side = 1, line = 5, cex = 0.85)
  par(op); dev.off()

  # ---- C3b: Table-3-style frequency-gradient controls ----
  t3_ctrl <- regional_control_frequencies(samples, topo)
  save_table(t3_ctrl, "C_table3_control_freq.csv")
  open_png("C3b_table3_controls.png", width = 1500, height = 950)
  op <- par(mfrow = c(2, 4), mar = c(6.5, 4.4, 3, 1), oma = c(0, 0, 2, 0))
  reg_order <- c("Europe", "North Africa", "Caucasus", "Middle East", "Anatolia")
  reg_cols <- c("Europe" = PALETTE[["Europe"]],
                "North Africa" = PALETTE[["North_Africa"]],
                "Caucasus" = PALETTE[["Caucasus"]],
                "Middle East" = PALETTE[["Levant"]],
                "Anatolia" = PALETTE[["Anatolia"]])
  for (cl in unique(t3_ctrl$clade)) {
    d <- t3_ctrl[t3_ctrl$clade == cl, ]
    d <- d[match(reg_order, d$region), ]
    barplot(d$frequency * 100, names.arg = d$region, las = 2,
            col = reg_cols[d$region], ylab = "frequency among regional samples (%)",
            main = sprintf("%s (%s)", cl, unique(d$expected)))
  }
  mtext("Table-3-style control gradients from Mitotree/GenBank regional samples",
        outer = TRUE, cex = 1.15, font = 2)
  par(op); dev.off()

  # ---- adegenet glPca ----
  sc <- founder_glpca(struct, topo, "C4_glpca_founders.png")
  if (!is.null(sc)) save_table(sc, "C_glpca_scores.csv")

  message(sprintf("Ancient records: %d (K1a/K2a/N1b with dates); oldest NE year=%s",
                  nrow(anc), if (nrow(ne)) min(ne$year) else NA))
  invisible(list(ancient = anc, coalescence = ca, table3 = t3,
                 table3_controls = t3_ctrl))
}

if (sys.nframe() == 0) run_ancient_and_diffusion()
