# =============================================================================
# make_figures.R
#
# Draws Figures 1 and 2 (network diagrams) from the pipeline's outputs.
#
# Usage:
#   Rscript make_figures.R <run_dir> [output_folder] [threshold=0.333]
#
#   run_dir        the pipeline's RUN_DIR (holds GRANGER/ and MODELS-STEP1/)
#   output_folder  default: <run_dir>/FIGURES
#   threshold=     Figure 2 inclusion-weight cut-off (default 1/3, as in the paper)
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
if (length(args) < 1) stop("Usage: Rscript make_figures.R <run_dir> [output_folder] [threshold=0.333]")
keyed     <- args[grepl("=", args)]
pos_args  <- args[!grepl("=", args)]
run_dir   <- pos_args[1]
out_dir   <- if (length(pos_args) >= 2) pos_args[2] else file.path(run_dir, "FIGURES")
thr_arg   <- sub("^threshold=", "", keyed[grepl("^threshold=", keyed)])
threshold <- if (length(thr_arg)) as.numeric(thr_arg) else 1/3
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(ggplot2))

granger_csv <- file.path(run_dir, "GRANGER", "granger_loocv_pairwise_elpd.csv")
coefs_csv   <- file.path(run_dir, "MODELS-STEP1",
                         "dsem_selection_from_single_behaviour_sets_stacking_coefs.csv")

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

save_fig <- function(p, stem, w = 9, h = 7) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  message("Saved ", file.path(out_dir, paste0(stem, ".pdf/.png")))
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
if (file.exists(granger_csv)) {
  g <- read.csv(granger_csv, stringsAsFactors = FALSE)
  shown <- g[g$elpd_diff_full_minus_reduced > 0 & g$mean_coef_x_to_y > 0, ]
  shown <- shown[order(shown$z_score), ]
  message(sprintf("Figure 1: %d of %d directed pairs shown (ELPD improved and positive coefficient); %d with z > 2",
                  nrow(shown), nrow(g), sum(shown$z_score > 2, na.rm = TRUE)))
  write.csv(shown, file.path(out_dir, "figure1_edges.csv"), row.names = FALSE)

  e <- edge_geometry(data.frame(from = shown$predictor_x, to = shown$outcome_y,
                                z = pmax(shown$z_score, 0), stringsAsFactors = FALSE))
  p1 <- base_plot()
  if (nrow(e) > 0) {
    zmax <- max(e$z, 2)
    sc_col <- scale_colour_gradient(low = "#F9C9B6", high = "#B2182B", limits = c(0, zmax), name = "z")
    sc_lw  <- scale_linewidth(range = c(0.4, 2.6), limits = c(0, zmax), guide = "none")
    if (any(e$bidirectional))
      p1 <- p1 + geom_curve(data = e[e$bidirectional, ],
                            aes(x, y, xend = xend, yend = yend, colour = z, linewidth = z),
                            curvature = 0.25, arrow = arrow_spec, lineend = "round", key_glyph = key_line)
    if (any(!e$bidirectional))
      p1 <- p1 + geom_segment(data = e[!e$bidirectional, ],
                              aes(x, y, xend = xend, yend = yend, colour = z, linewidth = z),
                              arrow = arrow_spec, lineend = "round", key_glyph = key_line)
    p1 <- p1 + sc_col + sc_lw
  }
  p1 <- add_nodes(p1)
  save_fig(p1, "figure1_granger_network")
} else {
  message("Figure 1 skipped: ", granger_csv, " not found (run granger + granger_combine first)")
}

# =============================================================================
# Figure 2
# =============================================================================
if (file.exists(coefs_csv)) {
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
  edges <- blocks[blocks$weight >= threshold, ]
  edges <- edges[order(edges$weight), ]

  mu <- cf[is_mu, ]
  mu$behaviour <- sub("^mu\\[(.+)\\]$", "\\1", mu$term_label)
  mu_in <- mu[mu$inclusion_weight >= threshold, ]

  message(sprintf("Figure 2: %d lag blocks (of %d) and %d means with inclusion weight >= %.3f",
                  nrow(edges), nrow(blocks), nrow(mu_in), threshold))
  write.csv(rbind(
    data.frame(type = ifelse(edges$from == edges$to, "self-lag", "cross-lag"),
               from = edges$from, to = edges$to, weight = edges$weight, direction = edges$direction),
    data.frame(type = rep("mean", nrow(mu_in)), from = mu_in$behaviour, to = mu_in$behaviour,
               weight = mu_in$inclusion_weight, direction = mu_in$effect_direction)),
    file.path(out_dir, "figure2_terms.csv"), row.names = FALSE)

  dir_col <- c(positive = RED, negative = BLUE)
  cross   <- edge_geometry(edges[edges$from != edges$to, ])
  loops   <- loop_paths(edges[edges$from == edges$to, ])
  wlim    <- c(threshold, 1)
  p2 <- base_plot()
  if (!is.null(cross) && nrow(cross) > 0) {
    if (any(cross$bidirectional))
      p2 <- p2 + geom_curve(data = cross[cross$bidirectional, ],
                            aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight),
                            curvature = 0.25, arrow = arrow_spec, lineend = "round", key_glyph = key_line)
    if (any(!cross$bidirectional))
      p2 <- p2 + geom_segment(data = cross[!cross$bidirectional, ],
                              aes(x, y, xend = xend, yend = yend, colour = direction, linewidth = weight),
                              arrow = arrow_spec, lineend = "round", key_glyph = key_line)
  }
  if (!is.null(loops))
    p2 <- p2 + geom_path(data = loops, aes(px, py, group = from, colour = direction, linewidth = weight),
                         arrow = arrow(length = unit(0.13, "inches"), type = "closed", angle = 22),
                         lineend = "round", key_glyph = key_line)
  ring <- setNames(mu_in$effect_direction, mu_in$behaviour)
  p2 <- add_nodes(p2, ring = ring, loop_nodes = edges$from[edges$from == edges$to]) +
    scale_colour_manual(values = dir_col, name = "Vocabulary",
                        labels = c(positive = "higher (positive)", negative = "lower (negative)"),
                        drop = FALSE) +
    scale_linewidth(range = c(0.5, 3), limits = wlim, name = "Inclusion\nweight",
                    breaks = pretty(wlim, 3)) +
    guides(linewidth = guide_legend(override.aes = list(colour = "grey35")))
  save_fig(p2, "figure2_vocabulary_network")
} else {
  message("Figure 2 skipped: ", coefs_csv, " not found (run step2a + step2b first)")
}
