# =============================================================================
# make_figures.R
#
# Draws Figures 1 and 2 (network diagrams) from the pipeline's outputs.
#
# Usage:
#   Rscript make_figures.R <run_dir> [output_folder] [figure=1|2|both] [threshold=0.01|perm|<p>]
#                          [rule=mean|pop] [prob=0.95] [z=z_score|z_child|z_block]
#                          [show=all|strong] [perm=auto|0] [alpha=0.05]
#
#   run_dir        the pipeline's RUN_DIR (holds GRANGER/ and MODELS-STEP1/)
#   output_folder  default: <run_dir>/FIGURES
#   figure=        which figure to draw: 1, 2 or both (default both)
#   threshold=     which Figure 2 blocks are drawn:
#                    0.01 (default) — every block with a final inclusion
#                         probability of at least 0.01, i.e. all blocks with
#                         non-negligible weight; line width shows the weight.
#                         This is a display floor, not an evidential threshold.
#                    perm — only blocks above the permutation-calibrated
#                         threshold: the 95th percentile, over the permutation
#                         runs, of the LARGEST final inclusion probability of
#                         any block (needs permutation_null.R output). A block
#                         above it has a familywise permutation p < .05.
#                    <number> — any other cut-off (the original paper used 1/3)
#   perm=          Figure 2 line style from the permutation null. Default 0:
#                  no significance marking, all lines solid. perm=1 (or auto:
#                  when permutation results exist): solid = significant,
#                  dashed = not. Uses the
#                  single-block test (PERMUTATIONS/single_block_null.csv,
#                  p_single_max < alpha) if present, otherwise the two-level
#                  results (p_final_max and p_inner both < alpha)
#   alpha=         significance level for perm (default 0.05)
#   weight=        what Figure 2's links show (default inclusion):
#                    inclusion — final inclusion probability (two-level model
#                                averaging); blocks chosen by threshold=
#                    gain — strength of evidence from the single-block test
#                           (single_block_test.R): every block that on its
#                           own improves leave-one-child-out prediction of
#                           vocabulary (ELPD gain > min_gain, default 0), line
#                           width = ELPD gain over the intercept-only model
#                    perm — as gain, but line width = -log10(permutation p of
#                           the block, p_single)
#                  With gain/perm, colour = sign of the block's effect in its
#                  single-block model, a node ring = the mean rate (ring width
#                  = its evidence); with perm=1 also solid lines = p_single <
#                  alpha (uncorrected), dashed = not.
#   min_gain=      minimum ELPD gain for weight=gain|perm (default 0)
#   rule=          Figure 1 "positive relation" rule:
#                    mean (default) — pooled mean lag coefficient > 0
#                    pop  — population-level effect (A_pop summed over the 5 lags)
#                           with P(effect > 0) >= prob; needs granger_direction.R output
#   prob=          probability cut-off for rule=pop (default 0.95)
#   z=             which z-score sets "strong support" (z > 2), arrow width
#                  and colour: z_score (default; observation-level SE),
#                  z_child (child-clustered SE) or z_block (block SE); the
#                  last two need the current granger_analysis.R output
#   show=          all (default): every edge passing the rule, strong edges
#                  (z > 2) solid and the others dashed; strong: strong edges only
#
# Figure 1 (Step 1b, granger_analysis.R)
#   Input: GRANGER/granger_loocv_pairwise_elpd.csv
#   An arrow x -> y for every directed pair where adding x's lags improved
#   predictions of y (ELPD difference > 0) AND the average lag coefficient
#   is positive. Arrow width and colour show the z-score (ELPD difference /
#   its SE); z > 2 is "strong support".
#
# Figure 2 (Step 2b, final_term_weighting.R)
#   Input: MODELS-STEP1/dsem_selection_from_single_behaviour_sets_stacking_coefs.csv
#   An arrow x -> y for every lag block A[y, x, lags 1-5] whose final
#   inclusion weight is >= threshold (self-lags as loops); red = higher
#   vocabulary, blue = lower. Arrow width shows the inclusion weight. A
#   coloured ring round a node means that behaviour's mean (intercept) is
#   included, coloured by direction. (The paper's log-residual-variance
#   text has no counterpart: those terms are no longer in the model.)
#
# Outputs (in output_folder): figure1_granger_network.{pdf,png},
#   figure2_vocabulary_network.{pdf,png}, and CSVs of exactly what was drawn.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript make_figures.R <run_dir> [output_folder] [figure=1|2|both] [threshold=0.333]")
keyed     <- args[grepl("=", args)]
pos_args  <- args[!grepl("=", args)]
run_dir   <- pos_args[1]
out_dir   <- if (length(pos_args) >= 2) pos_args[2] else file.path(run_dir, "FIGURES")
thr_arg   <- sub("^threshold=", "", keyed[grepl("^threshold=", keyed)])
thr_mode  <- if (length(thr_arg)) thr_arg[length(thr_arg)] else "0.01"
threshold <- if (thr_mode == "perm") NA_real_ else suppressWarnings(as.numeric(thr_mode))
if (thr_mode != "perm" && (is.na(threshold) || threshold < 0 || threshold > 1))
  stop("threshold= must be perm or a number between 0 and 1")
perm_arg  <- sub("^perm=", "", keyed[grepl("^perm=", keyed)])
perm_mode <- if (length(perm_arg)) tolower(perm_arg[length(perm_arg)]) else "0"
if (!perm_mode %in% c("auto", "0", "1")) stop("perm= must be auto, 0 or 1")
w2_arg    <- sub("^weight=", "", keyed[grepl("^weight=", keyed)])
fig2_weight <- if (length(w2_arg)) tolower(w2_arg[length(w2_arg)]) else "inclusion"
if (!fig2_weight %in% c("inclusion", "gain", "perm")) stop("weight= must be inclusion, gain or perm")
mg_arg    <- sub("^min_gain=", "", keyed[grepl("^min_gain=", keyed)])
fig2_min_gain <- if (length(mg_arg)) as.numeric(mg_arg) else 0
alpha_arg <- sub("^alpha=", "", keyed[grepl("^alpha=", keyed)])
perm_alpha <- if (length(alpha_arg)) as.numeric(alpha_arg) else 0.05
fig_arg   <- sub("^figure=", "", keyed[grepl("^figure=", keyed)])
which_fig <- if (length(fig_arg)) tolower(fig_arg[length(fig_arg)]) else "both"
if (!which_fig %in% c("1", "2", "both")) stop("figure= must be 1, 2 or both")
do_fig1 <- which_fig %in% c("1", "both")
do_fig2 <- which_fig %in% c("2", "both")
rule_arg  <- sub("^rule=", "", keyed[grepl("^rule=", keyed)])
fig1_rule <- if (length(rule_arg)) tolower(rule_arg[length(rule_arg)]) else "mean"
if (!fig1_rule %in% c("mean", "pop")) stop("rule= must be mean or pop")
prob_arg  <- sub("^prob=", "", keyed[grepl("^prob=", keyed)])
fig1_prob <- if (length(prob_arg)) as.numeric(prob_arg) else 0.95
z_arg     <- sub("^z=", "", keyed[grepl("^z=", keyed)])
fig1_z    <- if (length(z_arg)) z_arg[length(z_arg)] else "z_score"
if (!fig1_z %in% c("z_score", "z_child", "z_block")) stop("z= must be z_score, z_child or z_block")
show_arg  <- sub("^show=", "", keyed[grepl("^show=", keyed)])
fig1_show <- if (length(show_arg)) tolower(show_arg[length(show_arg)]) else "all"
if (!fig1_show %in% c("all", "strong")) stop("show= must be all or strong")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(ggplot2))

granger_csv <- file.path(run_dir, "GRANGER",
                         if (fig1_rule == "pop") "granger_direction.csv" else "granger_loocv_pairwise_elpd.csv")
coefs_csv   <- file.path(run_dir, "MODELS-STEP1",
                         "dsem_selection_from_single_behaviour_sets_stacking_coefs.csv")
perm_blocks_csv <- file.path(run_dir, "PERMUTATIONS", "permutation_null_blocks.csv")
perm_global_csv <- file.path(run_dir, "PERMUTATIONS", "permutation_null_global.csv")

# -----------------------------------------------------------------------------
# Layout: same node positions as the paper, clockwise from the top.
# -----------------------------------------------------------------------------
node_order <- c("gaze", "CV", "contingent", "show", "other_gesture",
                "open_point", "non_CV", "index_point", "give")
node_label <- c(contingent = "contingent", gaze = "gaze", CV = "CV", non_CV = "non CV",
                index_point = "index point", open_point = "open point", give = "give",
                show = "show", other_gesture = "other gesture")
ang   <- pi / 2 - 2 * pi * (seq_along(node_order) - 1) / length(node_order)
nodes <- data.frame(name = node_order, x = cos(ang), y = sin(ang),
                    label = node_label[node_order], stringsAsFactors = FALSE)
NODE_R <- 0.065           # node radius (data units; layout circle has radius 1)
LABEL_GAP <- 0.11         # labels sit just outside each node, away from the centre

RED  <- "#C8102E"
BLUE <- "#2B6CB0"

# Edge geometry: shorten each edge so it starts and ends at the node rims,
# and bend edges that have a partner in the opposite direction so the two
# don't overlap.
edge_geometry <- function(e) {
  if (nrow(e) == 0) return(e)
  fx <- nodes$x[match(e$from, nodes$name)]; fy <- nodes$y[match(e$from, nodes$name)]
  tx <- nodes$x[match(e$to,   nodes$name)]; ty <- nodes$y[match(e$to,   nodes$name)]
  d  <- sqrt((tx - fx)^2 + (ty - fy)^2)
  ux <- (tx - fx) / d; uy <- (ty - fy) / d
  e$x    <- fx + ux * NODE_R;        e$y    <- fy + uy * NODE_R
  e$xend <- tx - ux * NODE_R * 1.15; e$yend <- ty - uy * NODE_R * 1.15
  key  <- paste(e$from, e$to); rkey <- paste(e$to, e$from)
  e$bidirectional <- rkey %in% key
  e
}

# Self-loop: a circle sitting just outside the node, pointing away from the
# centre of the layout, drawn as a path with an arrow at its end.
loop_paths <- function(loops) {
  if (nrow(loops) == 0) return(NULL)
  out <- lapply(seq_len(nrow(loops)), function(i) {
    n  <- nodes[nodes$name == loops$from[i], ]
    r  <- 0.085
    th <- atan2(n$y, n$x)                       # outward direction
    cx <- n$x + cos(th) * (NODE_R + r * 0.55)
    cy <- n$y + sin(th) * (NODE_R + r * 0.55)
    a0 <- th + pi + 1.0; a1 <- th + pi - 1.0 + 2 * pi   # skip the arc inside the node
    a  <- seq(a0, a1, length.out = 60)
    data.frame(loops[rep(i, 60), , drop = FALSE], px = cx + r * cos(a), py = cy + r * sin(a),
               row.names = NULL)
  })
  do.call(rbind, out)
}

base_plot <- function() {
  ggplot() +
    coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.3, 1.35), clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.position = "right", plot.margin = margin(10, 10, 10, 10),
          legend.box.spacing = unit(1.6, "cm"),
          plot.background = element_rect(fill = "white", colour = NA))
}

# Nodes are drawn as circles in data units, so arrows end exactly at the rim.
node_circles <- function(nd) {
  a <- seq(0, 2 * pi, length.out = 60)
  do.call(rbind, lapply(seq_len(nrow(nd)), function(i) data.frame(
    name = nd$name[i], ring = nd$ring[i],
    cx = nd$x[i] + NODE_R * cos(a), cy = nd$y[i] + NODE_R * sin(a))))
}

add_nodes <- function(p, ring = NULL, loop_nodes = character(0)) {
  nd <- nodes
  nd$ring <- if (is.null(ring)) "none" else ifelse(is.na(ring[nd$name]), "none", ring[nd$name])
  circ <- node_circles(nd)
  th <- atan2(nd$y, nd$x)
  gap <- ifelse(nd$name %in% loop_nodes, LABEL_GAP + 0.17, LABEL_GAP)  # clear self-loops
  lab <- data.frame(label = nd$label,
                    lx = nd$x + cos(th) * (NODE_R + gap),
                    ly = nd$y + sin(th) * (NODE_R + gap),
                    hjust = ifelse(abs(cos(th)) < 0.3, 0.5, ifelse(cos(th) > 0, 0, 1)))
  p <- p + geom_polygon(data = circ[circ$ring == "none", ], aes(cx, cy, group = name),
                        fill = "#F4F6FB", colour = "#1F2A44", linewidth = 0.5)
  if (any(circ$ring != "none"))
    p <- p + geom_polygon(data = circ[circ$ring != "none", ],
                          aes(cx, cy, group = name, colour = ring),
                          fill = "#F4F6FB", linewidth = 1.8, show.legend = FALSE)
  p + geom_text(data = lab, aes(lx, ly, label = label, hjust = hjust), size = 3.6,
                fontface = "bold", colour = "#1F2A44")
}

# PNG output on a cluster node: the default png() device needs an X11
# display, which batch jobs don't have. Use ragg if installed, otherwise
# the cairo-based png(); if neither is available, write the PDF only.
png_device <- function() {
  if (requireNamespace("ragg", quietly = TRUE)) return(ragg::agg_png)
  if (isTRUE(capabilities("cairo")))
    return(function(filename, width, height, res, units, ...)
      grDevices::png(filename = filename, width = width, height = height,
                     res = res, units = units, type = "cairo", ...))
  NULL
}

save_fig <- function(p, stem, w = 9, h = 7) {
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  ggsave(pdf_file, p, width = w, height = h, device = grDevices::pdf)
  message("Saved ", pdf_file)
  dev_png <- png_device()
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  if (is.null(dev_png)) {
    message("PNG skipped: no display-free PNG device (install the 'ragg' package to enable it)")
  } else {
    ok <- tryCatch({ ggsave(png_file, p, width = w, height = h, dpi = 300, device = dev_png); TRUE },
                   error = function(e) { message("PNG skipped: ", conditionMessage(e)); FALSE })
    if (ok) message("Saved ", png_file)
  }
}

arrow_spec <- arrow(length = unit(0.16, "inches"), type = "closed", angle = 22)

# Legend key: a plain line (the default key would also draw the arrowhead)
key_line <- function(data, params, size) {
  params$arrow <- NULL
  draw_key_path(data, params, size)
}

# =============================================================================
# Figure 1
# =============================================================================
if (!do_fig1) {
  message("Figure 1 not requested (figure=", which_fig, ")")
} else if (file.exists(granger_csv)) {
  g <- read.csv(granger_csv, stringsAsFactors = FALSE)
  if (!fig1_z %in% names(g)) stop("Column ", fig1_z, " not in ", granger_csv,
                                  " (re-run granger_analysis.R with the current code)")
  g$z_use <- g[[fig1_z]]
  positive <- if (fig1_rule == "pop") g$pop_prob_positive >= fig1_prob else g$mean_coef_x_to_y > 0
  rule_txt <- if (fig1_rule == "pop") sprintf("population effect P(> 0) >= %.2f", fig1_prob) else "mean coefficient > 0"
  shown <- g[g$elpd_diff_full_minus_reduced > 0 & positive, ]
  shown$strong <- !is.na(shown$z_use) & shown$z_use > 2
  if (fig1_show == "strong") shown <- shown[shown$strong, ]
  shown <- shown[order(shown$z_use), ]
  message(sprintf("Figure 1: %d of %d directed pairs shown (ELPD improved and %s%s); %d with %s > 2",
                  nrow(shown), nrow(g), rule_txt, if (fig1_show == "strong") ", strong only" else "",
                  sum(shown$strong), fig1_z))
  write.csv(shown, file.path(out_dir, "figure1_edges.csv"), row.names = FALSE)

  e <- edge_geometry(data.frame(from = shown$predictor_x, to = shown$outcome_y,
                                z = pmax(shown$z_use, 0),
                                support = factor(ifelse(shown$strong, "strong", "weak"),
                                                 levels = c("strong", "weak")),
                                stringsAsFactors = FALSE))
  p1 <- base_plot()
  if (nrow(e) > 0) {
    zmax <- max(e$z, 2)
    zname <- c(z_score = "z", z_child = "z\n(child-\nclustered)", z_block = "z\n(block SE)")[[fig1_z]]
    sc_col <- scale_colour_gradient(low = "#F9C9B6", high = "#B2182B", limits = c(0, zmax), name = zname)
    sc_lw  <- scale_linewidth(range = c(0.4, 2.6), limits = c(0, zmax), guide = "none")
    sc_lt  <- scale_linetype_manual(values = c(strong = "solid", weak = "22"), drop = TRUE,
                                    name = "Support",
                                    labels = c(strong = "strong (z > 2)", weak = "weak (0 < z <= 2)"))
    if (any(e$bidirectional))
      p1 <- p1 + geom_curve(data = e[e$bidirectional, ],
                            aes(x, y, xend = xend, yend = yend, colour = z, linewidth = z, linetype = support),
                            curvature = 0.25, arrow = arrow_spec, lineend = "butt", key_glyph = key_line)
    if (any(!e$bidirectional))
      p1 <- p1 + geom_segment(data = e[!e$bidirectional, ],
                              aes(x, y, xend = xend, yend = yend, colour = z, linewidth = z, linetype = support),
                              arrow = arrow_spec, lineend = "butt", key_glyph = key_line)
    p1 <- p1 + sc_col + sc_lw + sc_lt +
      guides(linetype = guide_legend(override.aes = list(colour = "grey35", linewidth = 0.9), order = 2),
             colour = guide_colourbar(order = 1))
  }
  p1 <- add_nodes(p1)
  save_fig(p1, "figure1_granger_network")
} else {
  message("Figure 1 skipped: ", granger_csv, " not found (run granger + granger_combine first",
          if (fig1_rule == "pop") ", then granger_direction" else "", ")")
}

# =============================================================================
# Figure 2
# =============================================================================
if (!do_fig2) {
  message("Figure 2 not requested (figure=", which_fig, ")")
} else if (file.exists(coefs_csv)) {
  sb_csv <- file.path(run_dir, "PERMUTATIONS", "single_block_null.csv")
  if (fig2_weight != "inclusion") {
    # ---------------- Figure 2 weighted by single-block evidence ----------------
    if (!file.exists(sb_csv)) stop("weight=", fig2_weight, " needs ", sb_csv, " (run the single_block stage)")
    sb <- read.csv(sb_csv, stringsAsFactors = FALSE)
    sb <- sb[!is.na(sb$elpd_gain) & sb$elpd_gain > fig2_min_gain, ]
    sb$w <- if (fig2_weight == "gain") sb$elpd_gain else -log10(sb$p_single)
    sb$direction <- ifelse(sb$effect_mean > 0, "positive", "negative")
    sb$sig <- factor(ifelse(!is.na(sb$p_single) & sb$p_single < perm_alpha, "sig", "ns"), levels = c("sig", "ns"))
    is_mean <- grepl("^mean ", sb$block)
    edges <- data.frame(from = sub(" -> .*$", "", sb$block[!is_mean]), to = sub("^.* -> ", "", sb$block[!is_mean]),
                        weight = sb$w[!is_mean], direction = sb$direction[!is_mean], sig = sb$sig[!is_mean],
                        elpd_gain = sb$elpd_gain[!is_mean], p_single = sb$p_single[!is_mean],
                        stringsAsFactors = FALSE)
    edges <- edges[order(edges$weight), ]
    means <- data.frame(behaviour = sub("^mean ", "", sb$block[is_mean]), weight = sb$w[is_mean],
                        direction = sb$direction[is_mean], elpd_gain = sb$elpd_gain[is_mean],
                        p_single = sb$p_single[is_mean], stringsAsFactors = FALSE)
    message(sprintf("Figure 2 (weight=%s): %d lag blocks and %d means with ELPD gain > %g",
                    fig2_weight, nrow(edges), nrow(means), fig2_min_gain))
    write.csv(rbind(
      data.frame(type = ifelse(edges$from == edges$to, "self-lag", "cross-lag"), from = edges$from, to = edges$to,
                 weight = edges$weight, elpd_gain = edges$elpd_gain, p_single = edges$p_single, direction = edges$direction),
      data.frame(type = rep("mean", nrow(means)), from = means$behaviour, to = means$behaviour, weight = means$weight,
                 elpd_gain = means$elpd_gain, p_single = means$p_single, direction = means$direction)),
      file.path(out_dir, "figure2_terms.csv"), row.names = FALSE)

    wmax <- max(c(edges$weight, means$weight, 1e-6))
    wlim <- c(0, wmax)
    lw_of <- function(w) 0.3 + 2.7 * w / wmax           # same mapping as the line-width scale
    cross <- edge_geometry(edges[edges$from != edges$to, ])
    loops <- loop_paths(edges[edges$from == edges$to, ])
    mark_sig <- perm_mode %in% c("1", "auto")
    edge_aes <- if (mark_sig) aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight, linetype = sig) else
      aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight)
    loop_aes <- if (mark_sig) aes(px, py, group = from, colour = direction, linewidth = weight, linetype = sig) else
      aes(px, py, group = from, colour = direction, linewidth = weight)
    lend <- if (mark_sig) "butt" else "round"
    p2 <- base_plot()
    if (!is.null(cross) && nrow(cross) > 0) {
      if (any(cross$bidirectional))
        p2 <- p2 + geom_curve(data = cross[cross$bidirectional, ], mapping = edge_aes,
                              curvature = 0.25, arrow = arrow_spec, lineend = lend, key_glyph = key_line)
      if (any(!cross$bidirectional))
        p2 <- p2 + geom_segment(data = cross[!cross$bidirectional, ], mapping = edge_aes,
                                arrow = arrow_spec, lineend = lend, key_glyph = key_line)
    }
    if (!is.null(loops))
      p2 <- p2 + geom_path(data = loops, mapping = loop_aes,
                           arrow = arrow(length = unit(0.13, "inches"), type = "closed", angle = 22),
                           lineend = lend, key_glyph = key_line)
    p2 <- add_nodes(p2, loop_nodes = edges$from[edges$from == edges$to])
    # Mean rates: a ring round the node, coloured by direction, width = evidence.
    dir_col <- c(positive = RED, negative = BLUE)
    a <- seq(0, 2 * pi, length.out = 60)
    for (i in seq_len(nrow(means))) {
      n <- nodes[nodes$name == means$behaviour[i], ]
      ring_df <- data.frame(cx = n$x + NODE_R * cos(a), cy = n$y + NODE_R * sin(a))
      p2 <- p2 + geom_polygon(data = ring_df, aes(cx, cy), fill = NA, colour = dir_col[[means$direction[i]]],
                              linewidth = lw_of(means$weight[i]),
                              linetype = if (!mark_sig || means$p_single[i] < perm_alpha) "solid" else "22")
    }
    w_name <- if (fig2_weight == "gain") "ELPD gain\nover intercept-\nonly model" else "Evidence\n(-log10 p)"
    brks <- pretty(wlim, 3); brks <- brks[brks > 0 & brks <= wmax]
    p2 <- p2 +
      scale_colour_manual(values = dir_col, name = "Vocabulary",
                          labels = c(positive = "higher (positive)", negative = "lower (negative)"), drop = FALSE) +
      scale_linewidth(range = c(0.3, 3), limits = wlim, name = w_name, breaks = brks) +
      guides(linewidth = guide_legend(override.aes = list(colour = "grey35")))
    if (mark_sig)
      p2 <- p2 + scale_linetype_manual(values = c(sig = "solid", ns = "22"), drop = FALSE,
                                       name = "Permutation test\n(single block)",
                                       labels = c(sig = sprintf("p < %.2f (uncorrected)", perm_alpha),
                                                  ns = "not significant")) +
        guides(linetype = guide_legend(override.aes = list(colour = "grey35", linewidth = 0.9)))
    save_fig(p2, "figure2_vocabulary_network")
  } else {
  cf <- read.csv(coefs_csv, stringsAsFactors = FALSE)
  is_A  <- grepl("^A\\[", cf$term_label)
  is_mu <- grepl("^mu\\[", cf$term_label)

  # Lag blocks: all five lags of one x -> y block share one weight and
  # direction, so collapse to one edge per block.
  A <- cf[is_A, ]
  A$y <- sub("^A\\[y=([^,]+),x=([^,]+),lag=([0-9]+)\\]$", "\\1", A$term_label)
  A$x <- sub("^A\\[y=([^,]+),x=([^,]+),lag=([0-9]+)\\]$", "\\2", A$term_label)
  blocks <- do.call(rbind, lapply(split(A, paste(A$x, A$y)), function(b) data.frame(
    from = b$x[1], to = b$y[1], weight = max(b$inclusion_weight),
    direction = b$effect_direction[which.max(b$inclusion_weight)], stringsAsFactors = FALSE)))
  blocks <- if (is.null(blocks)) data.frame(from = character(0), to = character(0),
                                            weight = numeric(0), direction = character(0)) else blocks

  # Permutation-calibrated threshold: 95th percentile of the largest final
  # inclusion probability in the permutation runs.
  if (thr_mode == "perm") {
    if (!file.exists(perm_global_csv)) stop("threshold=perm needs ", perm_global_csv, " (run perm_summary)")
    pg <- read.csv(perm_global_csv, stringsAsFactors = FALSE)
    threshold <- pg$null_q95[pg$statistic == "max_final_inclusion"]
    if (length(threshold) != 1 || !is.finite(threshold))
      stop("No max_final_inclusion row in ", perm_global_csv, " (re-run permutation_null.R)")
    message(sprintf("Permutation-calibrated threshold: %.3f", threshold))
  }

  # Permutation p-values for the line style. Preferred: the single-block test
  # (each block alone vs intercept-only; familywise p_single_max). Otherwise
  # the two-level results, where a block counts only if BOTH its familywise
  # final p and its within-group p are below alpha.
  have_any <- file.exists(sb_csv) || file.exists(perm_blocks_csv)
  use_perm <- perm_mode == "1" || (perm_mode == "auto" && have_any)
  if (use_perm && !have_any) stop("perm=1 needs ", sb_csv, " or ", perm_blocks_csv)
  p_lookup <- NULL; perm_label <- ""
  if (use_perm && file.exists(sb_csv)) {
    sb <- read.csv(sb_csv, stringsAsFactors = FALSE)
    p_lookup <- setNames(sb$p_single_max, sb$block)
    perm_label <- "single-block, familywise"
  } else if (use_perm) {
    pb <- read.csv(perm_blocks_csv, stringsAsFactors = FALSE)
    p_lookup <- setNames(pmax(pb$p_final_max, pb$p_inner), pb$block)
    perm_label <- "final familywise and within-group"
  }
  if (use_perm) message("Figure 2 line style from permutation test: ", perm_label)
  block_key <- function(from, to) paste0(from, " -> ", to)

  edges <- blocks[blocks$weight >= threshold, ]
  edges <- edges[order(edges$weight), ]
  edges$p_perm <- if (use_perm && nrow(edges) > 0) unname(p_lookup[block_key(edges$from, edges$to)]) else rep(NA_real_, nrow(edges))
  edges$sig <- factor(ifelse(!is.na(edges$p_perm) & edges$p_perm < perm_alpha, "sig", "ns"),
                      levels = c("sig", "ns"))

  mu <- cf[is_mu, ]
  mu$behaviour <- sub("^mu\\[(.+)\\]$", "\\1", mu$term_label)
  mu_in <- mu[mu$inclusion_weight >= threshold, ]
  mu_in$p_perm <- if (use_perm && nrow(mu_in) > 0) unname(p_lookup[paste0("mean ", mu_in$behaviour)]) else rep(NA_real_, nrow(mu_in))

  message(sprintf("Figure 2: %d lag blocks (of %d) and %d means with final inclusion >= %.3f%s",
                  nrow(edges), nrow(blocks), nrow(mu_in), threshold,
                  if (use_perm) sprintf("; %d lag blocks and %d means with permutation p < %.2f (familywise)",
                                        sum(edges$sig == "sig"), sum(!is.na(mu_in$p_perm) & mu_in$p_perm < perm_alpha),
                                        perm_alpha) else ""))
  write.csv(rbind(
    data.frame(type = ifelse(edges$from == edges$to, "self-lag", "cross-lag"),
               from = edges$from, to = edges$to, weight = edges$weight, direction = edges$direction,
               p_perm_familywise = edges$p_perm),
    data.frame(type = rep("mean", nrow(mu_in)), from = mu_in$behaviour, to = mu_in$behaviour,
               weight = mu_in$inclusion_weight, direction = mu_in$effect_direction,
               p_perm_familywise = mu_in$p_perm)),
    file.path(out_dir, "figure2_terms.csv"), row.names = FALSE)

  dir_col <- c(positive = RED, negative = BLUE)
  cross   <- edge_geometry(edges[edges$from != edges$to, ])
  loops   <- loop_paths(edges[edges$from == edges$to, ])
  wlim    <- c(0, 1)
  # With the permutation null, line type marks familywise significance.
  edge_aes <- if (use_perm)
    aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight, linetype = sig) else
    aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight)
  loop_aes <- if (use_perm)
    aes(px, py, group = from, colour = direction, linewidth = weight, linetype = sig) else
    aes(px, py, group = from, colour = direction, linewidth = weight)
  lend <- if (use_perm) "butt" else "round"
  p2 <- base_plot()
  if (!is.null(cross) && nrow(cross) > 0) {
    if (any(cross$bidirectional))
      p2 <- p2 + geom_curve(data = cross[cross$bidirectional, ], mapping = edge_aes,
                            curvature = 0.25, arrow = arrow_spec, lineend = lend, key_glyph = key_line)
    if (any(!cross$bidirectional))
      p2 <- p2 + geom_segment(data = cross[!cross$bidirectional, ], mapping = edge_aes,
                              arrow = arrow_spec, lineend = lend, key_glyph = key_line)
  }
  if (!is.null(loops))
    p2 <- p2 + geom_path(data = loops, mapping = loop_aes,
                         arrow = arrow(length = unit(0.13, "inches"), type = "closed", angle = 22),
                         lineend = lend, key_glyph = key_line)
  ring <- setNames(mu_in$effect_direction, mu_in$behaviour)
  p2 <- add_nodes(p2, ring = ring, loop_nodes = edges$from[edges$from == edges$to]) +
    scale_colour_manual(values = dir_col, name = "Vocabulary",
                        labels = c(positive = "higher (positive)", negative = "lower (negative)"),
                        drop = FALSE) +
    scale_linewidth(range = c(0.3, 3), limits = wlim, name = "Final\ninclusion\nprobability",
                    breaks = c(0.1, 0.5, 1)) +
    guides(linewidth = guide_legend(override.aes = list(colour = "grey35")))
  if (use_perm)
    p2 <- p2 + scale_linetype_manual(values = c(sig = "solid", ns = "22"), drop = FALSE,
                                     name = paste0("Permutation test\n(", perm_label, ")"),
                                     labels = c(sig = sprintf("p < %.2f", perm_alpha),
                                                ns = "not significant")) +
      guides(linetype = guide_legend(override.aes = list(colour = "grey35", linewidth = 0.9)))
  save_fig(p2, "figure2_vocabulary_network")
  }
} else {
  message("Figure 2 skipped: ", coefs_csv, " not found (run step2a + step2b first)")
}
