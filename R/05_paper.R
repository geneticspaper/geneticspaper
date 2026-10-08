# 05_paper.R
# Second report format: a self-contained scientific-paper PDF.
#
# The Markdown report (R/04_report.R) is the quick, data-faithful dump. This
# module instead emits a LaTeX manuscript with proper prose, figures embedded
# beside the relevant text, and formatted tables, then compiles it to
# report/ashkenazi_mtdna_paper.pdf using tinytex::latexmk (no pandoc/quarto required).
#
# Design constraints:
#   * Base R only. No external R packages.
#   * ASCII-safe LaTeX (unicode is converted / escaped) so pdflatex succeeds
#     without fontspec or a system font.
#   * Reads the same generated tables/figures as the Markdown report so the two
#     formats never disagree.

if (!exists("ROOT")) source(file.path("R", "utils.R"))

# --------------------------------------------------------------------------
# LaTeX helpers.
# --------------------------------------------------------------------------

# Convert common unicode to LaTeX and escape special characters so the output
# compiles under pdflatex without fontspec.
.tex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  # Unicode -> LaTeX (do these before escaping backslash-producing commands).
  subs <- c(
    "\u2013" = "--", "\u2014" = "---", "\u2212" = "-",
    "\u2018" = "`", "\u2019" = "'", "\u201c" = "``", "\u201d" = "''",
    "\u2026" = "...", "\u00d7" = "$\\times$", "\u2248" = "$\\approx$",
    "\u2265" = "$\\geq$", "\u2264" = "$\\leq$", "\u00b1" = "$\\pm$",
    "\u00b5" = "$\\mu$", "\u03bc" = "$\\mu$", "\u03c3" = "$\\sigma$",
    "\u03bb" = "$\\lambda$", "\u03c1" = "$\\rho$", "\u00b0" = "$^{\\circ}$",
    "\u00e9" = "\\'e", "\u00e8" = "\\`e", "\u00ea" = "\\^e", "\u00eb" = '\\"e',
    "\u00e1" = "\\'a", "\u00e0" = "\\`a", "\u00e2" = "\\^a",
    "\u00ed" = "\\'i", "\u00ee" = "\\^i", "\u00ef" = '\\"i',
    "\u00f3" = "\\'o", "\u00f4" = "\\^o", "\u00f6" = '\\"o',
    "\u00fa" = "\\'u", "\u00fc" = '\\"u', "\u00f1" = "\\~n",
    "\u00e7" = "\\c{c}", "\u00c0" = "\\`A", "\u00c9" = "\\'E"
  )
  for (k in names(subs)) x <- gsub(k, subs[[k]], x, fixed = TRUE)
  # Escape LaTeX specials. Backslash first, but our unicode subs already added
  # intentional backslashes, so protect them with a placeholder.
  x <- gsub("\\", "\uE000", x, fixed = TRUE)               # temp-hide backslash
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)                    # escape specials
  # Use OT1 \char glyphs (not \textasciitilde/\textasciicircum, which resolve to
  # TS1 fonts that cannot be built on the fly in a sandbox without lmodern).
  x <- gsub("~", "\\\\char126{}", x)
  x <- gsub("\\^", "\\\\char94{}", x)
  x <- gsub("\uE000", "\\", x, fixed = TRUE)               # restore backslash
  # Drop any remaining non-ASCII (e.g. mojibake in imported study titles); all
  # intended output is already ASCII LaTeX, so nothing legitimate is lost.
  x <- gsub("[^\u0001-\u007F]", "", x)
  x
}

# Prose and captions are authored as valid LaTeX (ASCII, with explicit commands
# and $...$ math) and passed through verbatim. Only data-derived content (table
# cells and headers) is escaped via .tex_escape.

# Render a data.frame as a booktabs longtable-free tabular, sized to text width.
.tex_table <- function(df, caption, label, digits = 4, align = NULL,
                       colnames = NULL, fontsize = "\\small") {
  df <- as.data.frame(df)
  fmt <- function(v) {
    if (is.numeric(v)) formatC(v, format = "fg", digits = digits, big.mark = "")
    else as.character(v)
  }
  df[] <- lapply(df, fmt)
  header <- if (is.null(colnames)) names(df) else colnames
  # Guard against a colnames/ncol mismatch silently producing a misaligned table
  # with headerless trailing columns: pad the tail with the data-frame names and
  # truncate any excess so the header always has exactly ncol(df) cells.
  if (length(header) < ncol(df)) {
    header <- c(header, names(df)[(length(header) + 1L):ncol(df)])
  }
  header <- header[seq_len(ncol(df))]
  header <- vapply(header, .tex_escape, character(1))
  ncol <- ncol(df)
  if (is.null(align)) align <- paste0("l", strrep("l", ncol - 1))
  body <- apply(df, 1, function(r) paste(vapply(r, .tex_escape, character(1)), collapse = " & "))
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    fontsize,
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste0(paste(header, collapse = " & "), " \\\\"),
    "\\midrule",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

# Embed a figure with caption, placed near the calling text.
.tex_figure <- function(file, caption, label, width = 0.82) {
  path <- file.path("..", "outputs", "figures", file)
  c(
    "\\begin{figure}[htbp]",
    "\\centering",
    sprintf("\\includegraphics[width=%.2f\\linewidth]{%s}", width, path),
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\end{figure}"
  )
}

# Prose paragraph: authored LaTeX passed through verbatim, followed by a blank line.
.p <- function(...) c(paste(c(...), collapse = " "), "")

# --------------------------------------------------------------------------
# Manuscript body.
# --------------------------------------------------------------------------
build_paper <- function(compile = TRUE) {
  message("== Building LaTeX scientific paper ==")
  .read_tab <- function(name) read.csv(file.path(TAB_DIR, name), stringsAsFactors = FALSE)
  .read_ext <- function(name) {
    p <- file.path(ROOT, "data", "external", name)
    if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  }

  D0    <- .read_tab("00_enriched_dedup_summary.csv")
  A_t1  <- .read_tab("A_table1_single_descendant.csv")
  A_sens<- .read_tab("A_offspring_law_sensitivity.csv")
  A_mt  <- .read_tab("A_mitotree_founder_summary.csv")
  A_ma  <- .read_tab("A_mitotree_ancient_founders.csv")
  A_el  <- .read_tab("A_euro_levantine_feasibility.csv")
  A_els <- .read_tab("A_euro_levantine_sensitivity.csv")
  A_abs <- .read_tab("A_absorption_estimate.csv")
  A_absid <- .read_tab("A_absorption_identifiability.csv")
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
  D_fl  <- .read_tab("D_fitted_loo.csv")
  C_ca  <- .read_tab("C_coalescence_ages.csv")
  C_ne  <- .read_tab("C_oldest_near_east.csv")
  T7    <- load_table7()
  T7_founders <- T7[grepl("Livni", T7$source, fixed = TRUE), , drop = FALSE]
  KG1   <- tryCatch(load_1kg_clade_freq(), error = function(e) data.frame())
  D_li  <- .read_tab("D_lineage_computed_inputs.csv")
  PPNB  <- load_ppnb()
  NONEUR <- load_noneuropean()
  PUBLIC <- if (file.exists(file.path(ROOT, "data", "public_haplogroup_pages", "public_haplogroup_page_summary.csv"))) {
    load_public_haplogroup_pages()
  } else {
    data.frame()
  }
  G_mt  <- .read_ext("genes2026_mtdna.csv")

  elv <- setNames(A_el$value, A_el$quantity)
  p_ge3 <- as.numeric(elv[["prob_ge3_major"]]) * 100
  els_lo <- if (nrow(A_els)) 100 * min(A_els$prob_ge3_major) else p_ge3
  els_hi <- if (nrow(A_els)) 100 * max(A_els$prob_ge3_major) else p_ge3
  fs <- frequency_sources()
  ashk_src_label <- switch(
    fs$ashkenazi,
    "brook" = "Brook (2022) / Penninx (2019)",
    "mitotree" = "Mitotree",
    "livni_skorecki" = "Livni--Skorecki",
    fs$ashkenazi
  )
  nonjew_src_label <- switch(
    fs$nonjew,
    "1kg_eur" = "1000 Genomes EUR (CEU+GBR+FIN+IBS+TSI, n=503) subclades with Mitotree fallback",
    "livni_skorecki" = "Livni--Skorecki (+ Mitchell macro)",
    "mitotree" = "Mitotree",
    fs$nonjew
  )

  ppnb_assign <- PPNB[PPNB$haplogroup != "", ]
  ppnb_tab <- as.data.frame(table(ppnb_assign$haplogroup), stringsAsFactors = FALSE)
  names(ppnb_tab) <- c("Haplogroup", "Count")
  ppnb_tab <- ppnb_tab[order(-ppnb_tab$Count), ]

  noneur_tab <- data.frame(
    Lineage = sub(" \\(non-Eur\\)$", "", D_nec$founder),
    `Deep origin` = D_nec$deep_origin,
    `Ashk pct` = D_nec$ashkenazi_pct,
    `Erfurt` = D_nec$medieval_erfurt,
    `Post P(H2-side)` = sprintf("%.2f", D_nec$post_mean_H2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  noneur_afr_pct <- NONEUR$ashkenazi_pct[NONEUR$mtree_name == "L2a1l2a"]
  noneur_total_pct <- round(sum(NONEUR$ashkenazi_pct), 1)
  noneur_min_post <- min(D_nec$post_mean_H2)
  noneur_max_post <- max(D_nec$post_mean_H2)

  major_tab <- data.frame(
    Lineage = D_ml$founder,
    `Ashk pct` = D_ml$ashkenazi_pct,
    `Published` = D_ml$expected_origin,
    `Post P(H2)` = sprintf("%.2f", D_ml$post_mean_H2),
    `Model call` = ifelse(D_ml$post_mean_H2 >= 0.5, "Near Eastern", "European"),
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

  # data-calibrated robustness model (standardized channels + fitted logistic)
  fc_chan <- c("(Intercept)" = "Intercept", z_freq = "Frequency (non-Jewish rarity)",
               z_time = "Time (antiquity + carriers + public age)",
               # Data cell: goes through .tex_escape, so no math mode here.
               z_nest = "Nesting (equal-n richness)")
  fc_coef_tab <- data.frame(
    Channel = unname(fc_chan[D_fc$term]),
    `Posterior weight` = sprintf("%.2f [%.2f, %.2f]", D_fc$mean, D_fc$lo, D_fc$hi),
    check.names = FALSE, stringsAsFactors = FALSE)
  fc_loo_acc <- sum(D_fl$correct); fc_loo_n <- nrow(D_fl)
  v7_li <- D_li[D_li$lineage_key == "V7a2c1b", , drop = FALSE]
  u5_li <- D_li[D_li$lineage_key == "U5a1f1a3", , drop = FALSE]
  nc_post <- function(l, col = "post_mean_H2") .posterior_lineage(D_nc, l, col = col)
  v7_match <- if (nrow(v7_li)) v7_li$nonjew_freq_1kg_match[1] else "V"
  v7_pct <- if (nrow(v7_li)) 100 * v7_li$nonjew_freq[1] else 0
  u5_match <- if (nrow(u5_li)) u5_li$nonjew_freq_1kg_match[1] else "U5"
  u5_pct <- if (nrow(u5_li)) 100 * u5_li$nonjew_freq[1] else 0

  first_n <- function(x, n = 1) {
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
      Lineage = PUBLIC$lineage_key,
      Group = PUBLIC$group,
      `YFull TMRCA` = PUBLIC$yfull_tmrca,
      `YFull n` = PUBLIC$yfull_public_ids_n,
      `FTDNA TMRCA` = vapply(PUBLIC$ftdna_tmrca_mean, cal_year, character(1)),
      `First ancient` = vapply(PUBLIC$ftdna_ancient_codes, first_n, character(1)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    public_tab$`FTDNA TMRCA`[public_tab$`FTDNA TMRCA` == ""] <- "not exposed"
    public_tab$`First ancient`[public_tab$`First ancient` == ""] <- "not exposed"
  } else {
    public_tab <- data.frame()
  }

  # Public-enrichment coverage counts (dynamic, so the prose tracks whatever was
  # fetched rather than a hard-coded lineage count).
  n_public <- nrow(PUBLIC)
  n_yfull_ok <- if (n_public && "yfull_status" %in% names(PUBLIC)) sum(PUBLIC$yfull_status == "ok", na.rm = TRUE) else n_public
  n_ftdna_ok <- if (n_public && "ftdna_status" %in% names(PUBLIC)) sum(PUBLIC$ftdna_status == "ok", na.rm = TRUE) else n_public

  n_modern <- A_mt$n_modern_analysis[1]
  k_ancient <- A_mt$ancient_count[A_mt$founder == "K1a1b1a"]

  # Compact negative-binomial table (three representative growth ratios).
  D_nb_c <- D_nb[D_nb$growth_rate %in% c(1.15, 1.20, 1.28),
                 c("offspring_law", "growth_rate", "var_mean_ratio",
                   "prob_extinct", "prob_one_descendant")]
  # Collapse the 15-generation growth schedule to its three phases.
  D_sch_ph <- do.call(rbind, lapply(split(D_sch, D_sch$phase), function(s) {
    data.frame(phase = s$phase[1],
               generations = sprintf("%d-%d", min(s$generation), max(s$generation)),
               growth_ratio = s$growth_ratio[1], stringsAsFactors = FALSE)
  }))
  D_sch_ph <- D_sch_ph[match(c("medieval_stasis", "early_modern",
                               "great_expansion", "modern"), D_sch_ph$phase), ]
  D_sch_ph <- D_sch_ph[!is.na(D_sch_ph$phase), , drop = FALSE]

  # Posterior convenience accessors (Bayesian origin synthesis, Module D).
  # Fail loudly on an unknown key: a zero-length lookup here collapses the whole
  # enclosing sprintf() to character(0) and silently deletes the paragraph.
  pfmt  <- function(f, col = "post_mean_H2") {
    v <- D_bp[[col]][D_bp$founder == f]
    if (!length(v)) stop(sprintf("pfmt(): '%s' is not a founder in D_bp", f))
    .fmt_post(v)
  }
  # Control posteriors live in the negative-control table, not D_bp (which holds
  # the four founders only). Looking a control up via pfmt() silently yields
  # character(0), which collapses the whole enclosing sprintf() to character(0)
  # and drops the paragraph from the paper without any error -- which is exactly
  # what happened to the graded-synthesis paragraph in the Discussion.
  cfmt  <- function(f, col = "post_mean_H2") {
    v <- D_nc[[col]][D_nc$lineage == f]
    if (!length(v)) stop(sprintf("cfmt(): no control posterior for '%s'", f))
    .fmt_post(v)
  }
  a_eps_mt <- A_abs$absorbed_fraction_a_eps[A_abs$sample == "enriched_mitotree_modern"][1]

  # ---- preamble ----
  pre <- c(
    "\\documentclass[11pt]{article}",
    # Make the PDF text layer extractable (copy-paste / pdftotext / screen
    # readers). Without a ToUnicode map the subsetted Type1 fonts render
    # correctly on screen but export garbage, which mangles the verbatim code
    # excerpts in Appendix A when copied out. cmap adds the CMap tables; the
    # glyphtounicode fallback covers glyphs cmap misses. Both are guarded so a
    # missing file cannot break the build.
    "\\IfFileExists{glyphtounicode.tex}{\\input{glyphtounicode}\\pdfgentounicode=1}{}",
    "\\IfFileExists{cmap.sty}{\\usepackage{cmap}}{}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{lmodern}",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{graphicx}",
    "\\usepackage{booktabs}",
    "\\usepackage{amsmath}",
    "\\usepackage[hidelinks]{hyperref}",
    "\\IfFileExists{xurl.sty}{\\usepackage{xurl}}{\\usepackage{url}}",
    "\\hypersetup{breaklinks=true}",
    "\\Urlmuskip=0mu plus 2mu",
    "\\usepackage{caption}",
    "\\captionsetup{font=small,labelfont=bf}",
    "\\usepackage{parskip}",
    "\\graphicspath{{../outputs/figures/}}",
    "\\emergencystretch=3em",
    "\\sloppy",
    "\\title{\\bfseries Evidence for a Near Eastern origin of the four major Ashkenazi mitochondrial founder lineages}",
    "\\author{}",
    "\\date{\\today}",
    "\\begin{document}",
    "\\maketitle"
  )

  abstract <- c(
    "\\begin{abstract}",
    sprintf(paste(
      "Four maternal founder lineages -- K1a1b1a, K1a9, K2a2a, and N1b2 -- account for roughly",
      "%.0f\\%% of Ashkenazi Jewish mitochondrial DNA (mtDNA), but their geographic origin has been",
      "disputed between a prehistoric European assimilation model (H1; Costa et al., 2013) and a",
      "Near Eastern / Levantine origin model (H2; Behar et al., 2006; Livni and Skorecki, 2025).",
      "We test these hypotheses by re-implementing the Livni--Skorecki branching-process and",
      "founder-versus-host framework and integrating it with a phylogeographic, ancient-DNA, and",
      "Bayesian analysis of a deduplicated, enriched reference dataset (%s modern records from the",
      "Mitotree universal mtDNA phylogeny (Maier et al., 2026) and GenBank, with medieval Jewish",
      "ancient DNA). Five analyses favor a Near Eastern origin. (i) The",
      "estimated absorbed host-lineage fraction ($a_\\epsilon \\approx %.2f$) places all four founders",
      "far out in the multi-copy founder tail rather than among absorbed singletons. (ii) The four",
      "founders are rare among 27,651 non-Jews, which is difficult to reconcile with simple recent",
      "European-host assimilation. (iii) The Costa 'solely European nesting' argument is",
      "substantially weakened by rarefaction to equal sample size. (iv) Haplogroup K is the most frequent lineage (42.8\\%%, 6 of the 14 assigned skeletons) in Pre-Pottery",
      "Neolithic Near Eastern farmers (Fernandez et al., 2014), and the founder clades occur in",
      "medieval Jewish individuals from Erfurt (Waldman et al., 2022) and Tarrega",
      "(Pallares-Vina et al., 2026); Norwich",
      "(Brace et al., 2022) provides an earlier independent anchor for an Ashkenazi-like medieval",
      "maternal pool, though not direct K/N founder carriers. (v) A ridge-logistic origin synthesis",
      "on standardized data channels assigns every founder a posterior probability of Near",
      "Eastern origin above 0.5: the K1a founders most strongly (K1a1b1a %s, K1a9 %s), then",
      "N1b2 (%s), whose Europe-leaning equal-$n$ nesting tempers its rarity and antiquity",
      "signals, with K2a2a the most moderate (%s) and the only founder whose 90\\%% credible",
      "interval still includes parity -- it has no usable Near Eastern sister pool within K2a and",
      "no Jewish carrier dated before 1800, so its score rests on non-Jewish rarity almost alone.",
      "Run across a literature-coded broader-lineage panel, the same synthesis reproduces the",
      "evaluable published-origin calls in both directions -- scoring high the frequent lineages",
      "Costa et al. (2013) classify as Near Eastern or ultimately Near Eastern (HV1b2, R0a, U7,",
      "U1) while placing the clean H7 and J1c controls on the European side and flagging",
      "H6a1a1a as ambiguous. We",
      "conclude that all four major founders are best explained by a Near Eastern origin, with",
      "the strongest support for the K1a founders and the most tentative signal for K2a2a;",
      "mtDNA cannot fix an origin with autosomal-level certainty, but under these data a recent",
      "European-host origin is disfavored for every founder."),
      sum(load_founders()$ashkenazi_pct),
      format(n_modern, big.mark = ","), a_eps_mt,
      pfmt("K1a1b1a"), pfmt("K1a9"), pfmt("N1b2"), pfmt("K2a2a")),
    "\\end{abstract}"
  )

  intro <- c(
    "\\section{Introduction}",
    .p(paste(
      "Ashkenazi Jews descend from a small late-medieval founder population, a bottleneck",
      "that left strong signatures in both autosomal and uniparental markers. On the maternal",
      "side, a striking observation is that roughly 38\\% of Ashkenazi individuals carry one of",
      "four mitochondrial DNA (mtDNA) founder lineages: K1a1b1a, K1a9, K2a2a, and N1b2. The",
      "origin of these founders is contested.")),
    .p(paste(
      "Under hypothesis H1, articulated by Costa et al. (2013), the maternal founders are",
      "assimilated prehistoric European women, inferred from the observation that the founders",
      "appear to nest within predominantly European mtDNA sub-clades. Under hypothesis H2,",
      "supported by Behar et al. (2006) and formalized by Livni \\& Skorecki (2025), the",
      "founders are Levantine / Near Eastern, consistent with a Middle Eastern origin of the",
      "Ashkenazi maternal gene pool followed by demographic expansion.")),
    .p(paste(
      "This manuscript has two goals. First, it reproduces the Livni--Skorecki branching-process",
      "and founder-versus-host detection model in a transparent, auditable form. Second, it",
      "tests the competing hypotheses against a large, deduplicated reference dataset built from",
      "the Mitotree universal mtDNA phylogeny (Maier et al. 2026) and GenBank, enriched with",
      "medieval Jewish ancient DNA and framed by recent autosomal and ancient-genome studies.",
      "Throughout, we separate model-based inference from direct observation, and we",
      "avoid overclaiming from a reference phylogeny that is not a random population sample."))
  )

  methods <- c(
    "\\section{Materials and Methods}",
    "\\subsection{Data set}",
    .p(paste(
      "The analysis set is built from Mitotree S1 (the public release of the universal human",
      "mtDNA reference phylogeny), enriched with GenBank founder-query metadata and then",
      "deduplicated by accession and sample key. Behar (2006) and Costa (2013) records are not",
      "added as separate frequency tables; they are already present within Mitotree and are",
      "carried as study/provenance flags. Medieval Jewish ancient mtDNA (Erfurt; the Tarrega",
      "Roquetes pogrom, Pallares-Vina et al., 2026) is incorporated so that founder clades can",
      "be anchored in time. Table~\\ref{tab:data} summarizes the composition.")),
    .tex_table(D0, "Composition of the enriched, deduplicated reference data set.",
               "tab:data", colnames = c("Metric", "Value")),
    "\\subsection{Study design and data provenance}",
    .p(paste(
      "The analysis constructs an enriched, deduplicated reference sample set and applies five",
      "methodologically independent tests to it: a branching-process founder-versus-host model, an",
      "equal-sample-size phylogeographic rarefaction analysis, an ancient-DNA and diffusion-law",
      "summary, a non-Jewish rarity benchmark, and a Bayesian evidence synthesis. Every reported",
      "quantity is computed directly from the sample tables, tree topology, and cited study data,",
      "so each figure and posterior is traceable to a specific input rather than asserted.")),
    .p(paste(
      "Each module corresponds to a falsifiable claim in the literature. The branching module",
      "tests Livni--Skorecki's founder-versus-host logic: after a bottleneck, true founders should",
      "fall in the multi-copy descendant tail, whereas recently absorbed host lineages should",
      "appear mainly as singletons. The Behar sample is used as an Ashkenazi-specific benchmark",
      "for that founder spectrum. The phylogeography module tests Costa et al.'s central nesting",
      "claim by counting descendant sub-lineages in Europe and the Near East only after rarefying",
      "both regions to equal sample size. The ancient-DNA module tests whether the relevant parent",
      "macro-clades and medieval Jewish carriers exist in the right temporal and geographic",
      "context. Livni--Skorecki Table 7 supplies the non-Jewish rarity benchmark, and Table 3 is",
      "used as an external frequency-gradient benchmark. The Bayesian module then combines these",
      "pre-computed channels into a transparent evidence synthesis; it does not create new data.")),
    .p(paste(
      "The main inputs are data-based. Modern and ancient mtDNA rows, haplogroup names, country",
      "labels, TMRCA values and tree parent-child relationships are read from Mitotree and GenBank",
      "exports. Medieval Jewish records are keyed to their published studies. Behar, Costa,",
      "Livni--Skorecki, Fernandez, Waldman, Brace, Pallares-Vina, Brook/Penninx and Shamoon-Pour",
      "contribute explicit reference tables, frequencies or literature-coded classifications.",
      "Where the analysis uses manually parameterized controls (for example V7a2c1b, U5a1f1a3",
      "and the deep non-European lineages), they are labelled as calibration controls, not as",
      "independent discoveries from the Mitotree pipeline. Their role is to check whether the same",
      "scoring machinery can recover known European, Near Eastern and non-European cases.")),
    .p(paste(
      "This design makes the test meaningful and hard to tune post hoc. The same code must",
      "simultaneously reproduce Livni--Skorecki's branching predictions, keep Behar/Costa",
      "benchmarks separate from the Mitotree reference set, weaken or preserve Costa's nesting",
      "claim under equal-sample rarefaction, place European controls on the European side, place",
      "Near Eastern controls on the Near Eastern side, and still score the four founders from",
      "their observed rarity, ancient context and nesting. A contrived model could be made to",
      "favor one founder or one diagram; it is much harder for a single auditable pipeline to",
      "match this pattern across branching, rarity, rarefaction, diffusion, ancient DNA, negative",
      "controls and positive controls unless the underlying evidence is coherent.")),
    "\\subsection{Data structures and algorithms}",
    .p(paste(
      "All core analysis uses base R, with \\texttt{ape}, \\texttt{pegas} and \\texttt{adegenet}",
      "used only for the phylogenetic/network and ordination displays. The analysis draws on four",
      "data objects: (i) a sample table containing subject type, haplogroup, country, age and",
      "provenance flags; (ii) a Mitotree structure table containing each node, parent node, TMRCA",
      "and tip count; (iii) external reference tables for Livni--Skorecki, Behar, Costa,",
      "Brook/Penninx and ancient-DNA records; and (iv) a public YFull/FTDNA enrichment table. All",
      "joins are key-based: samples are assigned to clades by exact membership in the descendant set",
      "of a Mitotree node, not by free-text matching. Appendix~\\ref{app:code} gives short excerpts",
      "of the key algorithmic steps.")),
    .p(paste(
      "Topology is represented as a parent map and a child adjacency list. Descendant sets are",
      "computed by depth-first traversal: starting from a root haplogroup, the algorithm pushes",
      "children from the adjacency list onto a stack until all reachable descendants are visited.",
      "This operation is used identically for founder counts, parent-clade regional counts, ancient",
      "record extraction, rarefaction, and control panels. Deduplication is performed before",
      "analysis by accession/sample keys so that the same sequence is not double-counted when it",
      "appears through both Mitotree and GenBank/provenance enrichment.")),
    .p(paste(
      "The Galton--Watson distribution is computed from the probability-generating function",
      "$f_t(s)=\\sum_j p_{tj}s^j$. For a scalar growth ratio $m$ the same offspring law is",
      "composed $k$ times; for piecewise history we compose generation-specific functions",
      "$g_k(s)=f_1(f_2(\\cdots f_k(s)))$. The coefficients of $g_k$ are recovered on a finite",
      "support by evaluating the composed PGF on roots of unity and applying an inverse FFT. An",
      "independent direct polynomial-convolution calculation confirms this transform: the largest",
      "discrepancy across the reported growth grid is below $10^{-4}$, so numerical error",
      "is negligible relative to model uncertainty.")),
    .p(paste(
      "Rarefaction addresses the largest computational confound in Costa-style nesting arguments.",
      "For each clade and region, descendant sub-lineages are sampled without replacement to a",
      "common sample size $n$; the reported curve is the expected number of distinct descendant",
      "haplogroups after equalizing $n$ between Europe and the Near East. Thus a larger European",
      "reference panel cannot mechanically create a larger European nesting count. Regional",
      "frequency-gradient plots are separate descriptive summaries: they report",
      "$n_{clade,region}/n_{region}$ and are interpreted as directional context, not a stand-alone",
      "origin likelihood.")),
    .p(paste(
      "The origin synthesis combines three data-computed channels into one posterior",
      "probability of Near Eastern / Levantine origin (H2). Non-Jewish rarity uses",
      "Livni--Skorecki Table~7 subclade frequencies for the four founders ($n=27{,}651$) and",
      "the full 1000 Genomes Phase~3 European super-population (CEU+GBR+FIN+IBS+TSI) subclade frequencies (Haplogrep3; $n=503$) for all other",
      "scored lineages. Where the 1000 Genomes panel is too small to contain a scored lineage's",
      "macro background and clade-matching would otherwise climb to an unrelated broader clade",
      "(e.g.\\ HV1b2), the rarity is instead measured as that macro background's frequency among",
      "the $\\sim$14.7k European non-Jews in the enriched Mitotree pool, rather than a spurious",
      "climbed value. Rarity is combined by rank-based normal scores so these heterogeneous",
      "reference frames enter on a common scale; Brook/Penninx Ashkenazi percentages are",
      "recorded separately for audit only. Time combines Near Eastern macro-clade depth,",
      "medieval Jewish carrier counts, and a bounded YFull/FTDNA public-age adjustment.",
      "Nesting uses equal-sample-size rarefaction richness (positive = Near-East richer",
      "at equal $n$). Richness is counted on DEPTH-MATCHED labels: Mitotree does not name",
      "branches at a uniform depth by region --- it is built largely from European-ancestry",
      "consumer sequencing, and Europe is the more deeply named side in 14 of the 16 root",
      "clades used here --- so a plain count of distinct named nodes rewards whichever region",
      "the tree resolves more finely. Pairing the two draws at random and truncating both",
      "members of each pair to the shallower one's depth gives the two counted sets identical",
      "depth distributions by construction. Rarefied richness itself uses the exact Hurlbert",
      "(1971) expectation rather than a Monte Carlo mean, whose seed-to-seed spread exceeded",
      "the between-region differences being compared, and each value carries a Chao et al.",
      "(2014) bootstrap standard error; the channel is shrunk toward zero by its own",
      "uncertainty so a near-tie cannot enter the synthesis with the weight of a clear",
      "separation. The measurement is made on the lineage's SISTER variation: its own descendants,",
      "and those of any sibling founder sharing the same root clade, are excluded from both",
      "regional pools, since its own carriers are not",
      "evidence about which regional variation it nests inside, and because they concentrate",
      "on few named nodes they would otherwise depress the measured richness of whichever",
      "region they are most common in. The comparison root climbs to a parent clade until",
      "both regional sister pools reach at least 10 samples, and the root actually used is",
      "reported per lineage. Each channel is z-scored across scored lineages; a",
      "ridge-logistic model is fit to labeled controls (European $=0$; Near Eastern /",
      "non-European diaspora $=1$) with founders held out. Coefficient uncertainty is",
      "propagated by Laplace-approximation Monte Carlo; leave-one-out validation on the",
      "labeled panel guards against overfitting.")),
    "\\subsection{Non-Jewish European reference frequencies (1000 Genomes)}",
    .p(paste(
      "For non-founder lineages, non-Jewish rarity is estimated from the finest matching",
      "haplogroup prefix observed among the European reference individuals in the",
      "1000 Genomes Project Phase~3 callMom full-mtDNA call set --- all five EUR populations",
      "(CEU, $n=99$; GBR, $n=91$; FIN, $n=99$; IBS, $n=107$; TSI, $n=107$; $n=503$ in total),",
      "spanning northern/western and southern Europe (1000 Genomes Project Consortium, 2015). Sequences were classified with Haplogrep3",
      "(Phylotree~17 / rCRS; \\texttt{scripts/compute\\_1kg\\_eu\\_freq.py}). Each modeled",
      "lineage is matched by longest observed prefix on its YFull/FTDNA aliases, then on its",
      "Mitotree macro root; single-letter macro buckets are used only when no finer prefix is",
      "present in the panel (e.g.\\ V for V7 when no V7 carriers occur). This replaces Mitchell",
      "et al.\\ (2014) NHANES macro-haplogroup frequencies, which had inflated U and V counts for",
      "European controls. Diroma et al.\\ (2014) reconstructed mtDNA from 1000 Genomes whole-exome",
      "off-target reads; our benchmark uses the dedicated Phase~3 mitochondrial FASTA resource",
      "rather than that WES extraction pipeline.")),
    "\\subsection{Public YFull and FamilyTreeDNA enrichment}",
    .p(sprintf(paste(
      "We also added a reproducible public-page enrichment layer for every modeled haplogroup.",
      "The script \\texttt{scripts/03\\_fetch\\_public\\_haplogroup\\_pages.py} reads",
      "\\texttt{data/reference/public\\_haplogroup\\_page\\_targets.csv}, fetches and caches",
      "YFull MTree pages and FamilyTreeDNA Discover public mtDNA JSON, and extracts only fields",
      "visible without login: aliases, public sample counts, branch age estimates, country tags,",
      "ancient connections, variants and public project counts. Public pages were compiled for",
      "%d lineages (K1a4a has none);",
      "YFull pages were available for %d and FTDNA JSON for %d",
      "(Table~\\ref{tab:public-enrichment}); U1b1a1 resolves through the public parent U1b1 and",
      "A-a1b3a1 sits below macro-A on FTDNA. The public branch age adjusts the time channel only",
      "for lineages that also carry a Jewish or ancient anchor, so the European calibration",
      "lineages enter with no public-age bonus."),
      n_public, n_yfull_ok, n_ftdna_ok)),
    .p(paste(
      "YFull and FTDNA are genealogy and testing databases, so their modern country and project",
      "counts reflect participation, self-reporting and project enrollment. The model therefore",
      "does not use tester country proportions as population frequencies. It does use the least",
      "sample-sensitive fields to make the time channel more realistic: supported aliases, public",
      "branch/TMRCA estimates and public ancient-sample links. Very shallow public TMRCA estimates modestly",
      "penalize the time channel; mature public branch ages modestly strengthen it only when paired",
      "with an existing Jewish/ancient anchor. That anchor is restricted to pre-modern Jewish sites",
      "and studies (Erfurt/Waldman, Tarrega/Roquetes/Pallares-Vina, Chapelfield): the Sobibor series",
      "(Diepenbroek, 2021) is excluded, because those individuals were born 1893--1923 and so cannot",
      "evidence a lineage's presence in the Jewish maternal pool before the modern era -- the same",
      "date restriction applied to the medieval-carrier count. Table~\\ref{tab:public-enrichment} lists the",
      "public fields extracted for each modeled haplogroup.")),
    if (nrow(public_tab)) .tex_table(public_tab,
      "Public YFull/FTDNA enrichment for the modeled haplogroups. Calendar years are FTDNA public TMRCA means where exposed; YFull counts are public identifiers parsed from the MTree page.",
      "tab:public-enrichment", fontsize = "\\scriptsize") else character(0),
    "\\subsection{Branching-process model}",
    .p(paste(
      "We model each matriline as a Galton--Watson branching process. Each mother has a random",
      "number of daughters with mean equal to the per-generation growth ratio $m$, with the",
      "standard Galton--Watson assumption that daughter counts are independent conditional on",
      "the offspring law for that generation. The analysis uses the \\emph{finite-time}",
      "descendant-count distribution after the founder bottleneck, not only the ultimate",
      "extinction fixed point. The probability that a founder leaves exactly one matrilineal",
      "descendant after $k$ generations is obtained by $k$-fold composition of the offspring",
      "probability generating function, evaluated on the roots of unity and inverted by FFT. We",
      "verify this fast implementation against a slow, direct polynomial-convolution",
      "implementation; the two agree to better than $10^{-4}$.")),
    .p(paste(
      "The founder-versus-host detection argument is Poisson: a lineage of population frequency",
      "$\\nu$ is expected to appear $\\lambda = N\\nu$ times in a sample of size $N$, and the",
      "probability of observing it in at least $j$ copies follows the upper-tail Poisson",
      "distribution. Founder lineages (high $\\nu$) are seen in many copies; host lineages",
      "absorbed by conversion (low $\\nu$) appear as singletons.")),
    .p(paste(
      "We use this machinery to \\emph{estimate} the host/convert absorption rate. The",
      "individual-level singleton fraction $a_\\epsilon$ -- the share of sampled individuals who",
      "are the sole carrier of their fully resolved lineage -- estimates the absorbed fraction,",
      "and the per-generation rate follows Livni--Skorecki as $\\rho = 1 - (1-a_\\epsilon)^{1/K}$",
      "over founder-event depth $K$. We estimate it from the enriched Mitotree+GenBank set,",
      "treated as the largest available maternal sample, and cross-check against the dedicated",
      "Behar 2006 Ashkenazi sample; the derived surviving-absorbed-lineage count uses the",
      "founder-event growth, size, and extinction parameters so it is internally consistent.")),
    "\\subsection{Phylogeography and ancient DNA}",
    .p(paste(
      "We reconstruct the founder topology from the Mitotree structure table and track all",
      "descendant nodes of each founder clade. To address the Costa 'solely European nesting'",
      "argument, which is confounded by far larger European than Near Eastern sampling, we",
      "rarefy both regions to a common sample size before comparing sub-lineage richness.",
      "Ancient records with dates and country metadata provide time anchors, and coalescence",
      "ages are read from Mitotree. All tables and figures are produced by \\texttt{run\\_all.R}."))
  )

  # Results -----------------------------------------------------------------
  results <- c(
    "\\section{Results}",
    "\\subsection{The branching model makes small surviving lineages rare and is robust to the offspring law}",
    .p(sprintf(paste(
      "Reproducing Livni--Skorecki Table 1, the probability that a founder leaves exactly one",
      "matrilineal descendant after %d generations is far below 1\\%% across plausible growth",
      "ratios (Table~\\ref{tab:t1}, Figure~\\ref{fig:a1}). Because the original model assumes a",
      "Poisson offspring law, we add a sensitivity analysis under a geometric offspring law (a",
      "strongly overdispersed negative-binomial special case). The qualitative conclusion is",
      "unchanged: exact single-descendant survival remains a low-probability tail event",
      "(Table~\\ref{tab:sens}, Figure~\\ref{fig:a2b}). The same Poisson copy-number logic",
      "separates true founders from recently absorbed host lineages in a realistically sized",
      "sample, in which genuine founders appear in many copies while absorbed matrilines remain",
      "singletons (Figure~\\ref{fig:a2})."), GW_K)),
    .tex_table(A_t1, sprintf("Probability of exactly one matrilineal descendant after %d generations, by growth ratio.", GW_K),
               "tab:t1", colnames = c("Growth ratio", "P(one descendant)")),
    .tex_figure("A1_single_descendant.png",
                sprintf("Probability of a single surviving descendant at generation %d versus growth ratio.", GW_K),
                "fig:a1"),
    .tex_table(A_sens, "Sensitivity of the branching model to the offspring law (Poisson vs. geometric).",
               "tab:sens",
               colnames = c("Offspring model", "Growth", "P(extinct)", "P(one)", "P(>=3 desc)")),
    .tex_figure("A2b_offspring_law_sensitivity.png",
                "Single-descendant probability under Poisson and geometric (overdispersed) offspring laws.",
                "fig:a2b"),
    .tex_figure("A2_detection_curves.png",
                "Detectability of founder versus absorbed lineages in a sample of $N=600$ (cumulative Poisson).",
                "fig:a2"),
    "\\subsection{A simple extension: intermediate reproductive overdispersion}",
    .p(paste(
      "The Poisson-versus-geometric contrast spans only the extremes of a single axis,",
      "reproductive overdispersion. A negative-binomial offspring law with dispersion",
      "$\\phi$ has variance $m + m^2/\\phi$, so $\\phi=1$ is geometric and $\\phi\\to\\infty$ is",
      "Poisson. Sweeping intermediate dispersions (Table~\\ref{tab:nb},",
      "Figure~\\ref{fig:d1}) shows the single-descendant and extinction probabilities vary",
      "smoothly between the extremes: no member of this family removes the conclusion that",
      "exact single-descendant survival is a rare tail event.")),
    .tex_table(D_nb_c,
               "Single-descendant and extinction probabilities across the negative-binomial size ladder (geometric $=1$, Poisson $=\\infty$), at three growth ratios.",
               "tab:nb",
               colnames = c("Offspring law", "Growth", "Var/mean", "P(extinct)", "P(one desc)"),
               fontsize = "\\footnotesize"),
    .tex_figure("D1_nbinom_size.png",
                "Single-descendant probability across the negative-binomial dispersion ladder.",
                "fig:d1"),
    "\\subsection{A complex extension: piecewise historical growth}",
    .p(sprintf(paste(
      "A single constant growth ratio is demographically unrealistic. We replace it with a",
      "piecewise per-generation schedule (Table~\\ref{tab:sch}): a near-stationary medieval",
      "founder phase followed by rapid modern expansion, and compare it against a",
      "constant-growth process matched on \\emph{net} expansion. Holding total growth fixed at",
      "about %.1f-fold over %d generations, the realistic schedule raises extinction to %.1f\\%%",
      "(from %.1f\\%%) and lowers the singleton probability to %.2f\\%% (from %.2f\\%%;",
      "Figure~\\ref{fig:d2}). Early stasis prunes more lineages while rapid late expansion",
      "leaves survivors in many copies, so realistic history sharpens rather than weakens the",
      "founder-versus-host distinction."),
      D_pw$net_expansion[1], GW_K, 100 * D_pw$prob_extinct[2], 100 * D_pw$prob_extinct[1],
      100 * D_pw$prob_one_descendant[2], 100 * D_pw$prob_one_descendant[1])),
    .tex_table(D_sch_ph, sprintf("Calendar-anchored piecewise growth schedule over %d generations.", GW_K),
               "tab:sch", colnames = c("Phase", "Generations", "Growth ratio")),
    .tex_figure("D2_piecewise_growth.png",
                "Piecewise historical growth versus a constant-growth process matched on net expansion.",
                "fig:d2"),
    "\\subsection{A complex extension: origin synthesis (single fitted model)}",
    .p(sprintf(paste(
      "We combine the separate lines of evidence into a single posterior probability that each",
      "lineage is Near Eastern / Levantine (H2). Active inputs: Ashkenazi frequencies from %s;",
      "non-Jewish frequencies from %s (Figure~\\ref{fig:d3}). Three raw channels --- non-Jewish",
      "rarity, antiquity/time (depth, medieval carriers, public age), and equal-$n$ nesting ---",
      "are z-scored and combined by a ridge-logistic model fit to labeled controls with founders",
      "held out; the fitted channel weights are reported in Table~\\ref{tab:fcoef} and the",
      "standardized founder channels in Table~\\ref{tab:bscore}."),
      ashk_src_label, nonjew_src_label)),
    .p(sprintf(paste(
      "Leave-one-out validation on the %d labeled lineages classifies \\textbf{%d of %d correctly",
      "(%.0f\\%%)}, and the resulting per-founder posteriors are collected in",
      "Table~\\ref{tab:bayes}."),
      fc_loo_n, fc_loo_acc, fc_loo_n, 100 * fc_loo_acc / fc_loo_n)),
    .tex_table(fc_coef_tab,
               "Fitted channel weights (posterior mean and 90 percent credible interval).",
               "tab:fcoef", colnames = c("Channel", "Posterior weight"),
               fontsize = "\\small"),
    .tex_table(D_bs[, c("founder", "z_freq", "z_time", "z_nest", "rarity", "nest_frac")],
               "Standardized channel values for founders.",
               "tab:bscore",
               colnames = c("Founder", "z(freq)", "z(time)", "z(nest)", "rarity", "nest frac"),
               fontsize = "\\footnotesize"),
    .tex_table(D_bp, "Posterior probability of a Near Eastern / Levantine origin (H2) per founder, with 90 percent credible intervals.",
               "tab:bayes",
               colnames = c("Founder", "Post mean", "q05", "Median", "q95")),
    .tex_figure("D3_bayes_origin.png",
                "Origin synthesis: founders (red), European controls (blue), deep non-European controls (green); 90\\% credible intervals from the fitted logistic on standardized data channels.",
                "fig:d3"),
    .p(sprintf(paste(
      "All four founders sit above 0.5: K1a1b1a %s, K1a9 %s, K2a2a %s and N1b2 %s. K2a2a's",
      "interval is by far the widest and is the only one of the four that still includes parity.",
      "Its sister pool in the Near East is a single sample once its own carriers are excluded, so",
      "its comparison climbs to K, where the sister richness is near parity and the nesting",
      "channel contributes nothing; it also has no Jewish carrier dated before 1800, its only",
      "Jewish ancient records being 20th-century Sobibor victims, so its antiquity channel is",
      "negative. Its score therefore rests on non-Jewish rarity almost alone, and is best read as",
      "a single-channel result rather than as convergent evidence. The nesting channel is a",
      "negative contributor for N1b2 as well -- at the N1b macro background its sister lineages",
      "remain Europe-richer at equal $n$ -- but its rarity and time channels keep its interval",
      "clear of parity. (We report three decimals because rounding would overstate certainty.)"),
      pfmt("K1a1b1a"), pfmt("K1a9"), pfmt("K2a2a"), pfmt("N1b2"))),
    "\\subsection{Negative controls: European (absorbed) Ashkenazi lineages}",
    .p(paste(
      "A synthesis that returned a Levantine origin for \\emph{every} Ashkenazi lineage would be",
      "uninformative, so we applied the identical model to a small manually parameterized",
      "calibration panel of lineages that are Ashkenazi-carried but European in phylogeographic",
      "origin, and required the model to place them well below 0.5. The inputs are literature-coded",
      "rather than estimated de novo from the Mitotree pipeline. The first, V7a2c1b, is not a rare oddity: it is a common, established",
      "Ashkenazi maternal cluster at about 2.4\\% (Penninx, 2019; compiled in Brook, 2022;",
      "defining mutation A12753G), comparable to HV5a and not far below K2a2a. Its carriers",
      "include Yiddish-speaking individuals from Lithuania (YFull YF083232, YF108999), Galicia",
      "(GenBank OR803751, Dynow) and Ukraine (GenBank PQ337267, Berezna). Yet on every diagnostic",
      "that distinguishes origin it is the mirror image of the four founders: it nests inside an",
      "entirely European branch of haplogroup V (a post-glacial European expansion) rather than",
      "K1a/K2a/N1b; within the same subclade it is shared with non-Jewish Europeans (Danish",
      "KF162774, Russian JQ703830, and German and Belarusian samples), so the",
      "near-absence-among-non-Jews signal that flags the founders is absent; and it has no deep",
      "Near Eastern ancient record and no medieval Jewish (Erfurt/Tarrega) carrier.")),
    .p(paste(
      "We deliberately avoid a \\emph{too-young} argument, because the age of this branch is",
      "uncertain: YFull dates the tight subclade to about 100 years, whereas FamilyTreeDNA's",
      "Discover tool (which names the same branch V7a12a1) suggests a TMRCA near 371 BCE",
      "(95\\% CI 696--66 BCE; about 2,400 years). Coalescence age alone does not indicate a",
      "Near Eastern origin, and",
      "granting V7a2c1b even a Roman-era age does not move it toward H2, because its European",
      "signal comes entirely from nesting and non-Jewish sharing. V7a2c1b",
      "has only ten Near Eastern samples in its comparison pool even after the climb to HV0,",
      "sitting exactly on the minimum-pool floor and in the root with the largest",
      "naming-depth asymmetry in the panel, so we do not read its nesting value as evidence",
      "in either direction. The second control, U5a1f1a3,",
      "belongs to haplogroup U5 --- the signature European Mesolithic hunter-gatherer lineage,",
      "which dominates pre-Neolithic Europe and is essentially absent from the Pre-Pottery",
      "Neolithic Near East (Fernandez et al., 2014, recover K, not U5, in PPNB farmers). Any",
      "Ashkenazi U5 lineage is thus European by deep phylogeography, and because U5 is common in",
      "non-Jewish Europeans its non-Jewish-rarity signal is strongly negative. It provides an",
      "independent European control with a different origin story (Paleolithic hunter-gatherer",
      "versus post-glacial V), guarding against the objection that the V7a2c1b result is a",
      "one-off.")),
    .p(sprintf(paste(
      "Under the same fitted synthesis, the two controls return $P(\\text{H2}) = %.2f$",
      "(90\\%% CI %.2f--%.2f) for V7a2c1b and $%.2f$ (%.2f--%.2f) for U5a1f1a3 --- far below",
      "all four founders (%.2f--%.2f; Table~\\ref{tab:nc}, Figure~\\ref{fig:d3}). 1000 Genomes EUR subclade",
      "frequencies (%s %.1f\\%%, %s %.1f\\%%) supply the non-Jewish-rarity contrast that",
      "Table~7 alone cannot provide for non-founder lineages."),
      nc_post("V7a2c1b"), nc_post("V7a2c1b", "post_lo"), nc_post("V7a2c1b", "post_hi"),
      nc_post("U5a1f1a3"), nc_post("U5a1f1a3", "post_lo"), nc_post("U5a1f1a3", "post_hi"),
      min(D_bp$post_mean_H2), max(D_bp$post_mean_H2),
      v7_match, v7_pct, u5_match, u5_pct)),
    .tex_table(D_nc[, c("founder", "nonjew_freq", "z_freq", "z_time", "z_nest",
                         "post_mean_H2", "post_lo", "post_hi")],
               "Negative-control lineages: non-Jewish frequency, standardized channels, and posterior origin probability.",
               "tab:nc",
               colnames = c("Lineage", "non-Jew freq", "z(freq)", "z(time)", "z(nest)",
                            "Post mean", "q05", "q95"),
               fontsize = "\\footnotesize"),
    .p(paste(
      "The two negative controls are intentionally sparse but mechanistically different. V7a2c1b",
      "tests a common Ashkenazi cluster that is nevertheless shared with non-Jewish Europeans;",
      "U5a1f1a3 tests a deep European hunter-gatherer macro-clade. This guards against two",
      "failure modes: a classifier that simply rewards Ashkenazi frequency, and a classifier that",
      "cannot recognize older European lineages.")),
    "\\subsection{Positive controls: deep non-European lineages through the same model}",
    .p(sprintf(paste(
      "The founders and the European host lineages are not the whole story. The Ashkenazi",
      "maternal pool also contains lineages whose \\emph{deep} origin is neither European nor",
      "Levantine but African or Asian, and their presence is hard to reconcile with the strong",
      "form of H1 (recent assimilation from European hosts), because a European source cannot",
      "contribute African or East Asian mtDNA. The largest is L2a1l2a (about %s\\%% of Ashkenazi",
      "maternal lines; Penninx, 2019; Brook, 2022), a lineage of deep Sub-Saharan/East African",
      "origin already present in the 14th-century Erfurt Jewish cemetery (individual I13865).",
      "Several smaller lineages of African or Asian deep origin tell the same story."),
      noneur_afr_pct)),
    .p(paste(
      "We run these through the \\emph{identical} synthesis as a manually parameterized",
      "positive-control panel, which is instructive once the axis is read correctly. The model",
      "contrasts a recent \\emph{European-host} origin (H1) with a",
      "Near Eastern / Jewish-diaspora origin (H2); for the four founders the H2 pole is",
      "specifically Levantine, whereas for these lineages it means entry through the non-European",
      "diaspora rather than a literal Levantine birthplace. A high posterior is therefore a",
      "rejection of European-host assimilation. Because these clades are essentially absent among",
      "non-Jewish Europeans and nest entirely outside European diversity (with medieval Jewish",
      "carriers where documented), the model places them firmly on the non-European side",
      "(Table~\\ref{tab:noneur}).")),
    .tex_table(noneur_tab,
               "Ashkenazi maternal lineages of deep African or Asian origin, scored under the identical synthesis (frequencies Penninx 2019 / Brook 2022; Erfurt from Waldman et al. 2022).",
               "tab:noneur",
               colnames = c("Lineage", "Deep origin", "Ashk pct", "Erfurt", "Post P(H2-side)"),
               fontsize = "\\footnotesize"),
    .p(sprintf(paste(
      "All five land above 0.5 (posterior means %.2f--%.2f), \\emph{with} the founders and",
      "opposite the European controls (Figure~\\ref{fig:d3}); the two lineages with a documented",
      "Erfurt carrier (L2a1l2a and N9a3a1b1) are scored with that medieval Jewish anchor",
      "contributing to their time channel. This is",
      "exactly the intended discriminant behaviour: the same machinery scores genuinely European",
      "lineages low and genuinely non-European lineages high. Together these deep African/Asian",
      "lineages account for roughly %s\\%% of Ashkenazi maternal lines, and their complementary",
      "message is that the Ashkenazi maternal pool was assembled by a population that acquired",
      "lineages from Africa and Asia through the cosmopolitan diaspora history of a Near Eastern /",
      "Mediterranean Jewish community, not by recent absorption from northern-European hosts. For",
      "these lineages a high posterior means non-European, not specifically Levantine; their deep",
      "origin is African or Asian."), noneur_min_post, noneur_max_post, noneur_total_pct)),
    "\\subsection{Calibration on the broader frequent Ashkenazi pool}",
    .p(paste(
      "The four founders account for a large minority, not the entirety, of the Ashkenazi",
      "maternal pool. To test whether the synthesis tracks the published literature in both",
      "directions, we scored a small literature-coded panel of frequent non-founder lineages",
      "spanning the origin spectrum with the same transformations and posterior",
      "(Table~\\ref{tab:major}). This is a calibration analysis, not a new population-frequency",
      "estimate: a useful classifier should place lineages with published European assignments",
      "on the European side while retaining lineages with published Near Eastern assignments on",
      "the Near Eastern side, and should expose disputed cases rather than counting them as",
      "validation successes.")),
    .p(paste(
      "The Near Eastern side of this panel is deliberately conservative: Costa et al. (2013),",
      "despite arguing for substantial European maternal ancestry, classify HV1b2, R0a-associated",
      "lineages and U7 as Near Eastern sources and describe U1 as ultimately Near Eastern.",
      "Shamoon-Pour et al. (2019) further place HV1b2 near Assyrian HV1b branches in northern",
      "Mesopotamia and the South Caucasus, a plausible context for ancient Jewish communities.",
      "HV1b2 and U1b1a1 are also documented in the 14th-century Erfurt Jewish cemetery. Against",
      "these we set frequent lineages with published European assignments (H7 and J1c7a) and",
      "two deliberately ambiguous test cases, H6a1a1a and K1a4a. Brook/Wexler appear to treat the Ashkenazi",
      "H6a1a1a/H6a1a1a1 line as Middle Eastern, whereas the public Mitotree/FTDNA context and the",
      "supported parent H6a1a look Europe-rich. It is therefore not counted as a clean European",
      "validation control.")),
    .tex_table(major_tab,
               "Broader frequent Ashkenazi lineages scored under the identical synthesis, with published origin (Costa et al. 2013; Shamoon-Pour et al. 2019) and the model's call.",
               "tab:major",
               colnames = c("Lineage", "Ashk pct", "Published", "Post P(H2)", "Model call"),
               fontsize = "\\footnotesize"),
    .p(sprintf(paste(
      "The model agrees with the published origin in \\textbf{%d of %d evaluable} panel cases",
      "(Figure~\\ref{fig:d4}); %s are retained separately as ambiguous/disputed. The Near Eastern",
      "members all score high (%.2f--%.2f), whereas the clean European members (%s) score low",
      "(%.2f--%.2f). The ambiguous H6a1a1a case also scores low (%.2f), which is useful as a",
      "data-driven disagreement with a possible Brook Middle Eastern coding rather than as a",
      "validation point. This matters for interpretation in two ways. First, the four founders are",
      "not an isolated signal; several additional frequent Ashkenazi lineages, totaling roughly",
      "%s\\%% in the Brook/Penninx frequency compilation, are also Near Eastern under published",
      "phylogeographic assignments. Second, the low posteriors for the clean European controls",
      "show that the classifier is not biased toward H2: it reproduces European calls when the",
      "input evidence points that way, while exposing ambiguous cases."),
      major_conc, major_eval_n, major_amb_names,
      min(major_ne$post_mean_H2), max(major_ne$post_mean_H2),
      major_eu_names, min(major_eu$post_mean_H2), max(major_eu$post_mean_H2),
      if (nrow(major_amb)) major_amb$post_mean_H2[1] else NA_real_,
      major_ne_pct)),
    .tex_figure("D4_broader_pool.png",
                "Broader frequent Ashkenazi pool versus published origin: four founders (grey), lineages with published Near Eastern assignments (red diamonds), lineages with published European assignments (blue diamonds), and ambiguous/disputed lineages (grey squares).",
                "fig:d4"),
    "\\subsection{Founder clades are present and rare in the reference set}",
    .p(sprintf(paste(
      "All four founder clades are present in the enriched dataset and remain rare among modern",
      "reference samples (Table~\\ref{tab:mt}, Figure~\\ref{fig:a3}). K1a1b1a is the most",
      "frequent of the four and is anchored by %d ancient/historical records, including medieval",
      "Jewish individuals. We interpret modern frequencies as rarity/context only, because",
      "modern country reflects diaspora residence rather than lineage origin."), k_ancient)),
    .tex_table(A_mt[, c("founder", "modern_count", "modern_frequency", "ancient_count",
                        "distinct_modern_haplogroups", "known_country_count")],
               "Founder-clade counts in the enriched, deduplicated reference set.",
               "tab:mt",
               colnames = c("Founder", "Modern n", "Modern freq", "Ancient n", "Distinct hg", "Known country")),
    .tex_figure("A3_mitotree_founder_frequencies.png",
                "Frequency of the four major founder clades among modern Mitotree reference samples.",
                "fig:a3"),
    "\\subsection{Estimated host-lineage absorption rate}",
    .p(sprintf(paste(
      "Treating the enriched Mitotree+GenBank set as the largest available maternal sample, the",
      "individual-level singleton fraction is $a_\\epsilon \\approx %.2f$ (%s individuals), giving",
      "a per-generation absorption rate of roughly %.2f\\%%--%.2f\\%% across founder-event depths",
      "$K = 25$--$32$ (Table~\\ref{tab:abs}, Figure~\\ref{fig:a4}). The dedicated Behar 2006",
      "Ashkenazi sample yields a lower absorbed fraction ($\\approx %.2f$) and rate, as expected:",
      "the enriched set is global, only a small fraction is explicitly Ashkenazi, and it over-counts",
      "rare non-founder lineages. Its singleton fraction is, however, strongly a function of how",
      "many samples were drawn and how finely the tree labels them rather than of absorption:",
      "rarefying the same pool gives $a_\\epsilon$ from 0.22 at $n = 57{,}531$ to 0.96 at $n = 100$,",
      "and truncating the clade label from full depth to one character takes it from 0.22 to 0.00",
      "(Table~\\ref{tab:absid}). The enriched set's LOW value therefore reflects its large $n$, not",
      "less absorption, and the two estimates are comparable only at equal $n$ and equal clade",
      "resolution. We report both with Behar as the",
      "population-appropriate anchor. Under either estimate the four major founders lie far out in",
      "the multi-copy tail rather than among the absorbed singletons."),
      A_abs$absorbed_fraction_a_eps[A_abs$sample == "enriched_mitotree_modern"][1],
      format(A_abs$n_individuals[A_abs$sample == "enriched_mitotree_modern"][1], big.mark = ","),
      100 * min(A_abs$absorption_rate_per_gen[A_abs$sample == "enriched_mitotree_modern"]),
      100 * max(A_abs$absorption_rate_per_gen[A_abs$sample == "enriched_mitotree_modern"]),
      A_abs$absorbed_fraction_a_eps[A_abs$sample == "behar2006_ashkenazi"][1])),
    .tex_table(A_abs[, c("sample", "n_individuals", "absorbed_fraction_a_eps",
                         "K_generations", "absorption_rate_per_gen", "absorbed_lineages_surviving")],
               "Estimated host-lineage absorption from the largest available sample and the Behar 2006 benchmark.",
               "tab:abs",
               colnames = c("Sample", "n", "a_eps", "K", "Rate/gen", "Absorbed surviving"),
               fontsize = "\\footnotesize"),
    if (nrow(A_absid)) .tex_table(
      A_absid[, c("varied", "setting", "n", "clade_resolution", "a_eps",
                  "absorption_rate_per_gen")],
      paste("Identifiability of the absorbed-fraction estimate: the singleton fraction of one",
            "unchanging pool, varied only by sample size and by clade-label depth."),
      "tab:absid",
      colnames = c("Varied", "Setting", "n", "Resolution", "a_eps", "Rate/gen"),
      fontsize = "\\footnotesize") else character(0),
    .tex_figure("A4_absorption_rate.png",
                "Estimated per-generation absorption rate versus founder-event depth.",
                "fig:a4"),
    .p(sprintf(paste(
      "Under the Livni--Skorecki 'Euro-Levantine' reconciliation scenario, a large Roman-era",
      "Jewish population produces three or more of the four major founders with probability",
      "%.2f\\%% at our base assumptions (Figure~\\ref{fig:a6}). The per-lineage major-emergence",
      "probability, the private-mutation rate and the number of founder families are explicit",
      "assumptions rather than fitted quantities; varying each across a plausible range keeps",
      "this probability between %.2f\\%% and %.2f\\%% (Table~\\ref{tab:elsens}), so the",
      "model-disfavored conclusion does not rest on any single value. We report this as a model",
      "sensitivity result, not as an empirical measurement."), p_ge3, els_lo, els_hi)),
    .tex_figure("A6_euro_levantine.png",
                "Probability of $\\geq 3$ major founders arising in a large Euro-Levantine pool.",
                "fig:a6"),
    if (nrow(A_els)) .tex_table(
      transform(A_els, prob_ge3_major = round(100 * prob_ge3_major, 3),
                prob_ge2_major = round(100 * prob_ge2_major, 3)),
      "Sensitivity of the Euro-Levantine feasibility to its three explicit assumptions (per-lineage major-emergence probability, private-mutation rate, founder families). P values are percentages.",
      "tab:elsens",
      colnames = c("Assumption varied", "Value", "lambda", "P(>=2) %", "P(>=3) %"),
      fontsize = "\\footnotesize") else character(0),
    "\\subsection{The 'solely European nesting' argument weakens after sampling control}",
    .p(paste(
      "The Costa argument rests on counting more European than Near Eastern sub-lineages, but that",
      "count is confounded: the European reference sample is several-fold larger. When Europe and",
      "the Near East are rarefied to a common sample size the argument's premise weakens sharply",
      "(Figure~\\ref{fig:b6}). For K1a, which contains K1a1b1a, sub-lineage richness is in fact",
      "higher in the Near East at equal sampling; for K2a and N1b it remains higher in Europe, so",
      "the nesting evidence ranges from positively Near Eastern (K1a) to inconclusive (K2a, N1b)",
      "but in no case supports the strong European-nesting claim once sampling is controlled.",
      "As a calibration check, the same rarefaction line-chart method was applied to literature-coded",
      "control clades in a separate companion panel: European-assigned controls (H, U5, J1c) are",
      "Europe-equal or Europe-richer, whereas Near Eastern-assigned controls (R0a, U7, U1) are",
      "Near-East-richer (Figure~\\ref{fig:b6b}). The overall regional composition of the founder",
      "parent clades among reference samples with known country is shown in Figure~\\ref{fig:b1},",
      "and the immediate sibling clades of the founders are likewise not uniformly European",
      "(Figure~\\ref{fig:b2}).")),
    .tex_figure("B6_rarefaction.png",
                "Sub-lineage richness at equal sample size, Europe versus Near East.",
                "fig:b6"),
    .tex_figure("B6b_rarefaction_controls.png",
                "Rarefaction control panel: European-assigned controls (H, U5, J1c) and Near Eastern-assigned controls (R0a, U7, U1), each plotted as Europe versus Near East at equal sample size. The companion panel shows the method can move in both directions and is not mechanically biased toward a Near Eastern result.",
                "fig:b6b"),
    .tex_figure("B1_parentclade_region_props.png",
                "Regional composition of the parent clades (modern samples with known country).",
                "fig:b1"),
    .tex_figure("B2_nesting_K1a1b1a.png",
                "K1a1b1a among the sibling clades of K1a1b1.",
                "fig:b2"),
    .p(paste(
      "Decisively, the four founders are essentially absent among 27,651 non-Jews in the",
      "Livni--Skorecki / FTDNA benchmark (Table~\\ref{tab:t7}, Figure~\\ref{fig:b7}), with",
      "frequencies of order $10^{-4}$. If these lineages had entered the Ashkenazi gene pool by",
      "recent assimilation of local European women, they should be detectable at appreciable",
      "frequency in the surrounding non-Jewish populations; their near-absence is the single",
      "hardest observation for H1 to accommodate and is naturally explained if the founders",
      "entered with a Near Eastern source population. A reconstructed haplotype network of the",
      "founders together with their parent clades is shown in Figure~\\ref{fig:b8}.")),
    .tex_table(T7_founders[, c("haplotype", "count", "frequency", "n_nonjews")],
               "Frequency of the four major founders among 27,651 non-Jews (Livni--Skorecki Table 7 reference benchmark).",
               "tab:t7", colnames = c("Haplotype", "Count", "Frequency", "n non-Jews")),
    .tex_figure("B7_nonjewish_rarity.png",
                "Frequency of the four major founders among non-Jews.",
                "fig:b7"),
    .tex_figure("B8_founder_haplonet.png",
                "Reconstructed haplotype network of the founders and their parent clades.",
                "fig:b8"),
    "\\subsection{Ancient DNA places the parent macro-clades deep in the Near East}",
    .p(paste(
      "The oldest ancient K1a lineages in the enriched dataset are Near Eastern (Anatolia and",
      "the Levant), millennia before any Ashkenazi founder event (Table~\\ref{tab:ne},",
      "Figure~\\ref{fig:c1}). Independently, Fernandez et al. (2014) report that haplogroup K is",
      "the most frequent lineage among Pre-Pottery Neolithic B farmers of Syria (K in 6 of 15",
      "HVR1-typed individuals, about 40\\%, in the profiles compiled here; Table~\\ref{tab:ppnb}),",
      "with rare N* also present. These data are HVR1-level and therefore support the",
      "Near Eastern context of the K/N macro-clades rather than the derived Ashkenazi founder",
      "motifs themselves; nonetheless, they are difficult to reconcile with a purely European",
      "macro-clade narrative. The founders coalesce recently relative to their parent clades",
      "(Table~\\ref{tab:coal}, Figure~\\ref{fig:c2}), the expected signature of founder lineages.")),
    .tex_table(head(C_ne, 8), "Oldest Near Eastern K1a records in the enriched dataset.",
               "tab:ne", colnames = names(C_ne)),
    .tex_table(ppnb_tab, "PPNB Near Eastern Neolithic haplogroup composition (Fernandez et al. 2014).",
               "tab:ppnb"),
    .tex_figure("C1_ancient_timeline.png",
                "Ancient K1a / K2a / N1b records through time and space.",
                "fig:c1"),
    .tex_table(C_ca, "Coalescence ages (TMRCA) of founders versus parent clades.",
               "tab:coal", colnames = names(C_ca)),
    .tex_figure("C2_coalescence_ages.png",
                "TMRCA of founders (recent) versus their much older parent clades.",
                "fig:c2"),
    .tex_figure("C3_diffusion_frequencies.png",
                "Parent-signature frequencies by region (diffusion-law context).",
                "fig:c3"),
    .p(paste(
      "Regional parent-signature frequencies for the founder macro-clades are summarized in",
      "Figure~\\ref{fig:c3}. We also apply the same regional-frequency-gradient display to control",
      "clades (Figure~\\ref{fig:c3b}). European-assigned controls (H, U5, J1c) peak in Europe or are",
      "Europe-rich; Near Eastern-assigned controls (R0a, U7, U1) peak in Middle East, Caucasus",
      "or Anatolia in this Mitotree/GenBank reference set. H6a1a is shown separately as the",
      "Mitotree-supported parent/context for the disputed Ashkenazi H6a1a1a lineage, not as a clean",
      "published-European control. This supports using Table 3 as directional context while keeping",
      "it a supporting benchmark rather than a stand-alone origin proof.")),
    .tex_figure("C3b_table3_controls.png",
                "Table-3-style frequency-gradient controls from Mitotree/GenBank regional samples: European controls (H, U5, J1c), Near Eastern controls (R0a, U7, U1), and H6a1a as ambiguous H6 context.",
                "fig:c3b"),
    "\\subsection{Medieval Jewish anchors}",
    .p(paste(
      "Founder clades appear directly in medieval Jewish individuals. In Erfurt (14th century),",
      "11 of 31 unrelated individuals carried K1a1b1a, with additional K1a9 and N1b2/N1b1b1 carriers.",
      "In the Tarrega Roquetes pogrom sample (1348; Pallares-Vina et al., 2026), individual ROQ12 is K1a1b1a. These records",
      "move the argument from modern-frequency inference to demonstrated medieval Jewish",
      "presence. (The Erfurt K1a9 and N1b1b1 individuals come from the AADR compilation used in",
      "the stress tests, not from the Mitotree ancient subset tabulated in",
      "Table~\\ref{tab:anc}.)")),
    .p(paste(
      "An independent and earlier anchor comes from Norwich, England: Brace et al. (2022)",
      "sequenced six individuals from a medieval well (radiocarbon 1161--1216 calCE, consistent",
      "with the historically attested antisemitic massacre of 1190 CE) and found strong genetic",
      "affinity to modern Ashkenazi Jews (a qpAdm model of 100\\% Chapelfield fits present-day",
      "Ashkenazim), together with Ashkenazi-associated disease alleles already near modern",
      "frequency. Three of the Norwich individuals carried mitochondrial haplogroup H5c2, of",
      "which Ashkenazi Jews are the majority of modern carriers. Norwich does not sample the four",
      "K/N founders directly, but it independently demonstrates that a distinctively Ashkenazi",
      "maternal pool -- and the founder event that shaped it -- pre-dates the 12th century,",
      "corroborating the deep continuity that the founder hypothesis requires.")),
    .p(paste(
      "A nomenclature note: the Ashkenazi founder here labelled N1b2 (following Behar and Costa)",
      "is resolved as N1b1b1 on Family Tree DNA's February 2025 mitochondrial tree and in Brook",
      "(2022); 23andMe still reports it as N1b2, and the true PhyloTree N1b2 is essentially",
      "absent from Ashkenazim. The medieval Erfurt carrier belongs to this N1b1b1 clade. We",
      "retain the widely used N1b2 label for continuity with the prior literature.")),
    .tex_table(A_ma[, c("founder", "subject", "haplogroup", "country", "year", "study")],
               "Ancient / historical founder-clade records in the enriched dataset.",
               "tab:anc",
               colnames = c("Founder", "Subject", "Haplogroup", "Country", "Year", "Study"),
               fontsize = "\\footnotesize")
  )

  discussion <- c(
    "\\section{Discussion}",
    .p(paste(
      "The five analyses assembled here point consistently in one direction. Each attacks the",
      "problem differently -- a branching-process treatment of lineage survival, a",
      "sampling-controlled test of the nesting argument, an ancient-DNA record, a non-Jewish",
      "rarity benchmark, and a Bayesian synthesis -- yet none supports a European origin for any",
      "of the four founders, and every one either favors or is consistent with a Near Eastern",
      "origin; N1b2's Europe-leaning equal-$n$ nesting tempers an otherwise Near Eastern signal,",
      "but it stays clear of parity. K2a2a is the most tentative and rests on the narrowest",
      "evidential base: once its own carriers are excluded there is a single Near Eastern sister",
      "sample within K2a, so its nesting is measured at K and near parity, and it has no pre-1800",
      "Jewish carrier, leaving non-Jewish rarity as effectively its only positive channel -- and",
      "it is the one founder whose credible interval still admits parity.",
      "The convergence of",
      "methodologically independent lines of evidence is itself the central result: it is",
      "difficult to construct a European-assimilation scenario that simultaneously reproduces the",
      "founder-tail position of all four lineages, their near-total absence among non-Jews, the",
      "equal-sample-size nesting, the deep Neolithic Near Eastern K frequency, and the medieval",
      "Jewish carriers. A further, orthogonal observation points the",
      "same way: the Ashkenazi maternal pool also carries lineages of deep African or Asian origin",
      "(notably L2a1l2a, about 2.2\\%, already present at Erfurt), which a purely European source",
      "could not have supplied and which instead reflect the cosmopolitan diaspora history of a",
      "Near Eastern / Mediterranean Jewish population.")),
    .p(sprintf(paste(
      "The origin synthesis makes the conclusion quantitative and graded. The K1a founders are",
      "supported strongly: K1a1b1a (posterior mean %s; 90\\%% credible interval %s--%s) and",
      "K1a9 (%s; %s--%s) both have credible intervals that exclude parity, as does",
      "N1b2 -- scored as a single lineage together with its FTDNA-tree synonym N1b1b1 --",
      "(%s; %s--%s), whose posterior mean holds up on the strength of its rarity and time",
      "channels even though its equal-sample-size nesting in the N1b macro-clade leans European",
      "in this reference set on a heavily weighted channel. K2a2a (%s; %s--%s) is the most",
      "tentative of the four and the only one whose interval still admits parity, on a markedly",
      "narrower evidential base: excluding its own",
      "carriers leaves a single Near Eastern sister sample within K2a, so its nesting is measured",
      "at K, where sister richness is near parity, and it has no Jewish carrier dated before 1800",
      "(its Sobibor records are 20th-century and are not treated as pre-modern anchors). Its",
      "score therefore rests on non-Jewish rarity almost alone, and its point estimate should be",
      "read against that single-channel support rather than as convergent evidence.",
      "The honest reading is",
      "therefore a clear gradient: firm support for the K1a founders, a positive signal for N1b2,",
      "and a directionally Near Eastern but parity-admitting result for K2a2a. Crucially, no",
      "founder favors H1. That this gradient is not an artifact of a model tuned to return H2 is",
      "shown by the negative controls: the same synthesis assigns the European Ashkenazi lineages",
      "V7a2c1b and U5a1f1a3 posteriors of only %s and %s --- below all four founders --- and by",
      "leave-one-out validation on the labeled panel (%d/%d correct)."),
      pfmt("K1a1b1a"), pfmt("K1a1b1a", "post_lo"), pfmt("K1a1b1a", "post_hi"),
      pfmt("K1a9"), pfmt("K1a9", "post_lo"), pfmt("K1a9", "post_hi"),
      pfmt("N1b2"), pfmt("N1b2", "post_lo"), pfmt("N1b2", "post_hi"),
      pfmt("K2a2a"), pfmt("K2a2a", "post_lo"), pfmt("K2a2a", "post_hi"),
      cfmt("V7a2c1b"), cfmt("U5a1f1a3"),
      fc_loo_acc, fc_loo_n)),
    .p(sprintf(paste(
      "The branching model underpins this interpretation. Under realistic reproduction -- Poisson,",
      "geometric, or any intermediate negative-binomial law, and under piecewise historical growth",
      "as well as constant growth -- a lineage that survived the founder bottleneck and rose to",
      "the observed frequencies is overwhelmingly likely to appear in many copies, whereas a",
      "recently absorbed host matriline appears as a singleton. The measured absorbed fraction",
      "($a_\\epsilon \\approx %.2f$ in the largest available sample) confirms that the four major",
      "founders sit in the founder tail, not the host tail. Consistent with this, the",
      "Euro-Levantine reconciliation -- in which a large Roman-era population generates three or",
      "more of the founders by private mutation -- carries only a %.2f\\%% probability."),
      a_eps_mt, p_ge3)),
    .p(paste(
      "Genome-wide studies frame and, on balance, reinforce a Near Eastern origin for the Jewish",
      "autosomal core. Atzmon et al. (2010) showed that the major Jewish Diaspora groups form",
      "distinct clusters with shared Middle Eastern ancestry, and framed the four mtDNA founders",
      "(about 40\\% of the Ashkenazi maternal pool) as Middle Eastern in origin. Agranat-Tamir et",
      "al. (2020) supply the deep ancient-DNA anchor: Bronze Age 'Canaanite' genomes of the",
      "Southern Levant to which present-day Jewish groups, including Ashkenazi Jews, trace 50\\% or",
      "more of their ancestry. Carmi et al. (2014) document the severe late-medieval bottleneck and",
      "roughly even European / Middle Eastern autosomal ancestry.")),
    .p(paste(
      "The complementary European component is real and well characterised, but it is",
      "predominantly Southern European and does not bear directly on the maternal founders. Xue et",
      "al. (2017) date this gene flow to at least two events around the founder bottleneck, and",
      "recent local-ancestry inference (Lerga-Jaso et al., 2025) assigns Ashkenazi Jews a large",
      "Italian / Southern European component alongside a Levantine one. Because Southern European",
      "populations themselves carry substantial Near Eastern ancestry, and because autosomal",
      "ancestry is not uniparental, a large Southern European autosomal fraction is fully",
      "compatible with --- and indeed expected under --- a Levantine maternal core (H2) followed by",
      "later European admixture; it does not imply that the maternal founders were European.")),
    .p(paste(
      "Two caveats bound but do not overturn the conclusion. First, the Pre-Pottery Neolithic K",
      "samples (Fernandez et al., 2014) were typed only on the control region and cannot be",
      "resolved to the specific K1a1b1a or K2a2a sub-motifs; they establish deep Near Eastern K at",
      "the macro-clade level rather than the founder sub-clade level. Second, lineage loss through",
      "drift means the absence of a derived sub-clade in a small ancient panel is not proof of its",
      "past absence. Both caveats limit the resolution of the ancient-DNA channel. They weaken",
      "over-specific claims at the founder-subclade level, but they do not supply positive",
      "evidence for the recent European-host origin model."))
  )

  conclusions <- c(
    "\\section{Conclusions}",
    .p(sprintf(paste(
      "Integrating a branching-process model, a phylogeographic rarefaction analysis, an",
      "ancient-DNA record, a non-Jewish rarity benchmark, and a ridge-logistic origin synthesis over",
      "the largest available maternal reference dataset, we find that all four major Ashkenazi mtDNA",
      "founders -- K1a1b1a, K1a9, K2a2a and N1b2 -- are best explained by a Near Eastern /",
      "Levantine origin. The support is strongest for the two K1a founders (K1a9 %s,",
      "K1a1b1a %s, credible intervals excluding parity), then N1b2 (%s), whose Europe-leaning",
      "nesting tempers its rarity and antiquity signals but leaves its interval clear of parity,",
      "and is positive though most tentative for K2a2a (%s), the one founder whose interval still",
      "includes parity and whose score rests on non-Jewish rarity almost alone.",
      "No founder is better explained by a simple prehistoric",
      "European-host origin under the analyses performed here."),
      pfmt("K1a9"), pfmt("K1a1b1a"), pfmt("N1b2"), pfmt("K2a2a"))),
    .p(paste(
      "These results support and extend the Livni--Skorecki (2025) and Behar et al. (2006)",
      "Near Eastern model relative to the Costa et al. (2013) European-assimilation model, and they",
      "do so with an explicit, reproducible quantification of the uncertainty. Mitochondrial DNA",
      "alone cannot fix a geographic origin with the precision of genome-wide data, and finer",
      "resolution -- coding-region typing of ancient Near Eastern K and N carriers, and dense",
      "Levantine population sampling -- would sharpen the N1b2 estimate in particular.",
      "On present evidence, however, the Near Eastern origin of the Ashkenazi maternal founder",
      "core is the parsimonious and best-supported conclusion.")),
    "\\vspace{0.5em}",
    "\\noindent\\textit{Reproducibility.} All results, tables, and figures are reproduced end-to-end from the public Mitotree release, GenBank queries, and the cited study data by a single analysis script.")

  limitations <- c(
    "\\section{Limitations}",
    "\\begin{itemize}",
    "\\item \\textbf{Reference-set sampling.} Mitotree/GenBank is a cross-study reference set; we treat it as the largest available maternal sample and estimate absorption from its singleton fraction, but that fraction is sensitive to lineage-calling resolution and global (non-Ashkenazi) diversity, so the enriched-set estimate is an upper bound and Behar 2006 is the Ashkenazi-specific anchor.",
    "\\item \\textbf{Diaspora bias.} Modern country reflects residence, not origin, and is not used as direct evidence of origin.",
    "\\item \\textbf{Public subset.} The public Mitotree release is a subset of the full phylogeny; ancient Near Eastern coverage is sparse.",
    "\\item \\textbf{YFull/FTDNA public enrichment.} Public genealogy pages improve alias resolution, branch-age context and ancient-anchor visibility, but tester country/project counts are participation-biased and are not treated as population frequencies.",
    "\\item \\textbf{Literature-coded controls.} Some non-founder controls use published classifications and manually curated frequencies rather than being discovered de novo from Mitotree. They are used for calibration and falsification, not as independent population estimates.",
    "\\item \\textbf{Regularized logistic synthesis.} Channel weights are fit to a small labeled panel with ridge priors and leave-one-out validation; the posterior is an evidence synthesis, not a fully generative demographic likelihood. European controls (V7, U5) use 1000 Genomes EUR (n=503) subclade frequencies and score far below all four founders.",
    "\\item \\textbf{Branching independence.} The Galton--Watson layer assumes independent offspring counts conditional on the chosen law and generation. Negative-binomial overdispersion and piecewise growth test sensitivity to realistic variance and timing, but they do not model full family-level correlation or overlapping generations.",
    "\\item \\textbf{Reconstructed motifs.} The haplotype network and glPCA use motifs reconstructed from Mitotree defining-mutation strings, not per-sample alignments.",
    "\\item \\textbf{Model framework.} The branching equations are used both as a sensitivity framework and, via the singleton fraction, as an absorption-rate estimator; the estimate is a moment-style estimate, not a full likelihood fit with per-individual sampling weights.",
    "\\item \\textbf{Not a primary study.} This pipeline integrates public data to test H2 against H1; it is not a substitute for peer-reviewed primary analysis.",
    "\\end{itemize}"
  )

  code_appendix <- c(
    "\\appendix",
    "\\section{Computational Appendix: Core Algorithm Excerpts}",
    "\\label{app:code}",
    .p(paste(
      "The full source code is the reproducible record. The excerpts below are included because",
      "they define the computational genetics operations most likely to affect interpretation:",
      "tree traversal, branching-process coefficient recovery, rarefaction, public-page time",
      "adjustment, and Bayesian score construction. Boilerplate plotting, file I/O and table",
      "formatting are intentionally omitted.")),
    "\\subsection{Mitotree descendant traversal}",
    .p(paste(
      "All clade counts use graph traversal over the Mitotree parent-child table. A sample is",
      "inside a clade if its resolved haplogroup is a member of the descendant set.")),
    "\\begingroup\\footnotesize",
    "\\begin{verbatim}",
    "build_topology <- function(struct) {",
    "  parent <- setNames(struct$ParentHaplogroup, struct$Haplogroup)",
    "  children <- split(struct$Haplogroup, struct$ParentHaplogroup)",
    "  list(parent = parent, children = children, nodes = struct$Haplogroup)",
    "}",
    "",
    "descendants <- function(root, topo, include_self = TRUE) {",
    "  out <- character(0); stack <- root",
    "  while (length(stack)) {",
    "    node <- stack[[1]]; stack <- stack[-1]",
    "    if (node %in% out) next",
    "    out <- c(out, node)",
    "    stack <- c(stack, topo$children[[node]])",
    "  }",
    "  if (!include_self) setdiff(out, root) else out",
    "}",
    "\\end{verbatim}",
    "\\endgroup",
    "\\subsection{Branching-process PGF composition}",
    .p(paste(
      "The descendant-count distribution after $k$ generations is recovered by composing",
      "offspring probability-generating functions and inverting the coefficients by FFT. This",
      "avoids slow repeated convolutions while preserving an audit path: the paper also writes a",
      "direct-convolution check to \\texttt{A\\_branching\\_math\\_audit.csv}.")),
    "\\begingroup\\footnotesize",
    "\\begin{verbatim}",
    "gw_descendant_pmf <- function(m, k = GW_K, jmax = NULL, model, size = 1) {",
    "  if (length(m) == 1L) m <- rep(m, k)",
    "  # The inverse FFT ALIASES mass above the grid back onto low j rather than",
    "  # truncating it, so the grid is sized from E[Z_k] = prod(m).",
    "  N <- if (is.null(jmax)) 2^ceiling(log2(max(512, 64 + 40 * prod(m)))) else jmax + 1L",
    "  omega <- exp(-2i * pi * (0:(N - 1)) / N)",
    "  eval_pgf <- function(z, off) {",
    "    val <- complex(length(z))",
    "    for (jj in seq(length(off) - 1, 0)) val <- val * z + off[jj + 1]",
    "    val",
    "  }",
    "  gz <- omega",
    "  for (t in rev(seq_len(k))) {",
    "    off <- gw_offspring_pmf(m[t], jmax = jmax, model = model, size = size)",
    "    gz <- eval_pgf(gz, off)",
    "  }",
    "  p <- Re(fft(gz, inverse = TRUE) / N)",
    "  p[p < 0] <- 0; p / sum(p)",
    "}",
    "\\end{verbatim}",
    "\\endgroup",
    "\\subsection{Rarefied nesting score}",
    .p(paste(
      "The Costa nesting claim is evaluated after equalizing sample size. The model uses",
      "expected sub-lineage richness at matched $n$, not raw European and Near Eastern counts.")),
    "\\begingroup\\footnotesize",
    "\\begin{verbatim}",
    "# .rarefy_at_root() in R/lineage_data.R; sketch:",
    "rarefy_at_root <- function(root) {",
    "  s <- B_rar[B_rar$root == root, ]",
    "  nmax <- max(s$n[!is.na(s$distinct_lineages)])",
    "  w <- s[s$n == nmax, ]",
    "  list(ne = w$distinct_lineages[w$region == 'Near East'],",
    "       eu = w$distinct_lineages[w$region == 'Europe'],",
    "       ne_pool = s$pool_size[s$region == 'Near East'][1],",
    "       eu_pool = s$pool_size[s$region == 'Europe'][1])",
    "}",
    "",
    "# positive = Near East (or non-Europe) richer at equal n; a single z_nest",
    "# weight is fit across all lineages, so the sign is NOT flipped by macro-clade.",
    "s_nest <- (ne_richness - eu_richness) / (ne_richness + eu_richness)",
    "# richness is DEPTH-MATCHED, and s_nest is then shrunk toward 0 by its own",
    "# bootstrap SE:  s_nest * s_nest^2 / (s_nest^2 + se^2)",
    "\\end{verbatim}",
    "\\endgroup",
    "\\subsection{Origin synthesis (standardized channels + ridge logistic)}",
    .p(paste(
      "Raw channels are z-scored across scored lineages; a ridge-logistic model is fit to labeled",
      "controls with founders held out. Coefficient uncertainty is propagated by Laplace-approximation",
      "Monte Carlo.")),
    "\\begingroup\\footnotesize",
    "\\begin{verbatim}",
    "z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)",
    "X <- cbind(1, z_freq, z_time, z_nest)",
    "fit <- fit_logistic_bayes(train)   # ridge priors tau = 2, tau0 = 5",
    "draws <- fit$draws                 # M = 40000 Laplace-approximation draws",
    "p <- plogis(X %*% t(draws))",
    "\\end{verbatim}",
    "\\endgroup",
    "\\subsection{Public YFull/FTDNA time adjustment}",
    .p(paste(
      "Public genealogy data enter only through bounded branch-age and ancient-anchor context.",
      "Tester country/project counts remain descriptive and do not enter the frequency channel.")),
    "\\begingroup\\footnotesize",
    "\\begin{verbatim}",
    "public_time_adjustment <- function(tmrca_kyr, confidence, ancient_anchor,",
    "                                   medieval_carriers) {",
    "  if (is.na(tmrca_kyr) || confidence <= 0) return(0)",
    "  recent_penalty <- if (tmrca_kyr < 1) -0.25 * confidence else 0",
    "  # carriers are their own sub-channel of z_time, so they must NOT also",
    "  # set `anchored` here -- that would count one record twice in one channel.",
    "  anchored <- public_jewish_ancient_anchor > 0",
    "  mature_bonus <- if (anchored && tmrca_kyr >= 1.5) 0.12 * confidence else 0",
    "  recent_penalty + mature_bonus",
    "}",
    "\\end{verbatim}",
    "\\endgroup"
  )

  refs <- c(
    "\\section*{References}",
    "\\begin{itemize}\\small",
    "\\item Livni J. \\& Skorecki K. (2025) Distinguishing between founder and host population mtDNA lineages in the Ashkenazi population. Human Gene; SSRN 5035272.",
    "\\item Costa M.D. et al. (2013) A substantial prehistoric European ancestry amongst Ashkenazi maternal lineages. Nat. Commun. 4:2543.",
    "\\item Behar D.M. et al. (2006) The matrilineal ancestry of Ashkenazi Jewry. Am. J. Hum. Genet. 78:487--497.",
    "\\item Atzmon G. et al. (2010) Abraham's children in the genome era: major Jewish Diaspora populations comprise distinct genetic clusters with shared Middle Eastern ancestry. Am. J. Hum. Genet. 86:850--859.",
    "\\item Carmi S. et al. (2014) Sequencing an Ashkenazi reference panel supports population-targeted personal genomics and illuminates Jewish and European origins. Nat. Commun. 5:4835.",
    "\\item Agranat-Tamir L. et al. (2020) The genomic history of the Bronze Age Southern Levant. Cell 181:1146--1157.",
    "\\item Lerga-Jaso J. et al. (2025) Tracing human genetic histories and natural selection with precise local ancestry inference. Nat. Commun. 16:4576.",
    "\\item Xue J. et al. (2017) The time and place of European admixture in Ashkenazi Jewish history. PLoS Genet. 13(4):e1006644.",
    "\\item Waldman S. et al. (2022) Genome-wide data from medieval German Jews show that the Ashkenazi founder event pre-dated the 14th century. Cell 185:4703--4716.",
    "\\item Brace S. et al. (2022) Genomes from a medieval mass burial show Ashkenazi-associated hereditary diseases pre-date the 12th century. Curr. Biol. 32:4350--4359.",
    "\\item Fernandez E. et al. (2014) Ancient DNA analysis of 8000 B.C. Near Eastern farmers. PLoS Genet. 10(6):e1004401.",
    "\\item Shamoon-Pour M., Li M. \\& Merriwether D.A. (2019) Rare human mitochondrial HV lineages spread from the Near East and Caucasus during post-LGM and Neolithic expansions. Sci. Rep. 9:14751.",
    "\\item Pallares-Vina L. et al. (2026) Uncovering a medieval pogrom: genetic history of a Jewish community in Catalonia (Spain). Genes 17(3):358.",
    "\\item Maier P.A. et al. (2026) Mitotree: the universal human mitochondrial reference phylogeny. bioRxiv.",
    "\\item 1000 Genomes Project Consortium (2015) A global reference for human genetic variation. Nature 526:68--74.",
    "\\item Diroma M.A. et al. (2014) Extraction and annotation of human mitochondrial genomes from 1000 Genomes Whole Exome Sequencing data. BMC Genomics 15(Suppl 3):S2.",
    "\\item Mitchell S.L. et al. (2014) Characterization of mitochondrial haplogroups in a large population-based sample from the United States. Hum. Genet. 133:861--868.",
    "\\item Brook K.A. (2022) The Maternal Genetic Lineages of Ashkenazic Jews. Academic Studies Press.",
    "\\item Penninx W. (2019) Analysis of mtDNA haplogroup frequencies in the Ashkenazi Jewish population (compiled in Brook 2022).",
    "\\item Wexler J.D. (2019--2026) Ashkenazi Y-DNA and mtDNA (compilation). \\url{https://sites.google.com/view/ashkenazi-y-dna-and-mtdna}",
    "\\item YFull MTree (accessed by public haplogroup pages, 2026). \\url{https://www.yfull.com/mtree/}",
    "\\item FamilyTreeDNA Discover public mtDNA haplogroup pages and JSON resources (accessed 2026). \\url{https://discover.familytreedna.com/mtdna/}",
    "\\item FamilyTreeDNA Discover, mtDNA haplogroup V7a12a1 age estimate. \\url{https://discover.familytreedna.com/mtdna/V7a12a1/classic}",
    "\\end{itemize}",
    "\\end{document}"
  )

  tex <- c(pre, abstract, intro, methods, results, discussion, conclusions,
           limitations, code_appendix, refs)
  out_tex <- file.path(ROOT, "report", "ashkenazi_mtdna_paper.tex")
  writeLines(tex, out_tex)
  message(sprintf("Wrote %s (%d lines)", out_tex, length(tex)))

  if (compile) compile_paper(out_tex)
  invisible(out_tex)
}

compile_paper <- function(tex_path) {
  wd <- dirname(tex_path)
  base <- basename(tex_path)
  pdf <- sub("\\.tex$", ".pdf", base)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)

  # Make a project-local library visible so tinytex is found without a global
  # install, and put the TinyTeX binaries on PATH.
  rlib <- file.path(ROOT, ".Rlib")
  if (dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))
  for (bindir in c("~/.local/bin", "~/bin", "~/.TinyTeX/bin/x86_64-linux")) {
    bindir <- path.expand(bindir)
    if (dir.exists(bindir) && !grepl(bindir, Sys.getenv("PATH"), fixed = TRUE)) {
      Sys.setenv(PATH = paste(bindir, Sys.getenv("PATH"), sep = .Platform$path.sep))
    }
  }

  # Preferred, idiomatic R path: tinytex::latexmk() auto-installs any missing
  # LaTeX packages/fonts (e.g. lmodern) on demand and runs the needed passes.
  if (requireNamespace("tinytex", quietly = TRUE)) {
    ok <- tryCatch({
      tinytex::latexmk(base, engine = "pdflatex")
      TRUE
    }, error = function(e) { message("tinytex::latexmk failed: ", conditionMessage(e)); FALSE })
    if (ok && file.exists(pdf)) {
      message(sprintf("Built PDF (tinytex): %s", file.path(wd, pdf)))
      return(invisible(file.path(wd, pdf)))
    }
  }

  # Fallback: system LaTeX engine, two passes for cross-references.
  engine <- Sys.which(c("pdflatex", "lualatex", "xelatex"))
  engine <- engine[nzchar(engine)]
  if (!length(engine)) {
    message("No LaTeX engine available; wrote .tex only.")
    return(invisible(NULL))
  }
  for (i in 1:2) {
    suppressWarnings(system2(engine[[1]],
      c("-interaction=nonstopmode", "-halt-on-error", shQuote(base)),
      stdout = FALSE, stderr = FALSE))
  }
  if (file.exists(pdf)) {
    for (ext in c(".aux", ".log", ".out", ".fls", ".fdb_latexmk")) {
      f <- sub("\\.tex$", ext, base); if (file.exists(f)) unlink(f)
    }
    message(sprintf("Built PDF: %s", file.path(wd, pdf)))
  } else {
    message("LaTeX compile did not produce a PDF; see .log in report/.")
  }
  invisible(file.path(wd, pdf))
}

if (sys.nframe() == 0) { source(file.path("R", "utils.R")); build_paper() }
