# =============================================================================
# make_tables.R
#
# Supplementary tables for the statistical support behind Figures 1 and 2.
#
# Usage:
#   Rscript make_tables.R <run_dir> [output_folder] [z=z_score|z_child|z_block]
#     output_folder default: <run_dir>/TABLES
#     z=  which z-score defines "strong" in Table S1 (default: z_child if
#         present, else z_score)
#
# Table S1 (table_S1_granger_pairs.csv) — all 72 directed pairs x -> y:
#   ELPD difference (full minus reduced model) with three SEs and z-scores
#   (observation-level, child-clustered, block; see granger_analysis.R),
#   pooled mean lag coefficient, population-level effect (A_pop summed over
#   the lags: mean, 90% interval, P(> 0); from granger_direction.R if run),
#   and whether the pair is drawn in Figure 1 and whether it is strong.
#
# Table S2 (table_S2_step2_blocks.csv) — every Step 2 predictor block
#   (lag block x -> y, or the mean of y):
#   inner (within-group) inclusion, group (outer) inclusion and final
#   inclusion probability (from final_term_weighting.R); the block's effect
#   on log vocabulary (sum of the coefficients of its terms) averaged over
#   the Step 2a models containing it — mean, 90% interval and P(> 0) — and
#   the 90% interval of the fully model-averaged effect; and, when available,
#   permutation results: the group's p (p_group, p_group_max; from
#   permutation_null.R), the single-block test (single_elpd_gain, its SE,
#   p_single, p_single_max; from single_block_test.R) and the two-level block
#   p-values (p_inner, p_final, p_final_max). Effects need Step 2a output from the current
#   select_terms_for_individual_behaviours.R ($block_effects).
#   The weighting (pseudo-BMA+ or stacking) follows the Step 2b run.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript make_tables.R <run_dir> [output_folder] [z=...]")
keyed    <- args[grepl("=", args)]
pos_args <- args[!grepl("=", args)]
run_dir  <- pos_args[1]
out_dir  <- if (length(pos_args) >= 2) pos_args[2] else file.path(run_dir, "TABLES")
z_arg    <- sub("^z=", "", keyed[grepl("^z=", keyed)])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
prefix   <- "dsem_selection_from_single_behaviour_sets"

# ---------------------------------------------------------------- Table S1
g_dir <- file.path(run_dir, "GRANGER")
g_file <- if (file.exists(file.path(g_dir, "granger_direction.csv"))) file.path(g_dir, "granger_direction.csv") else
  file.path(g_dir, "granger_loocv_pairwise_elpd.csv")
if (file.exists(g_file)) {
  g <- read.csv(g_file, stringsAsFactors = FALSE)
  zc <- if (length(z_arg)) z_arg[length(z_arg)] else if ("z_child" %in% names(g)) "z_child" else "z_score"
  if (!zc %in% names(g)) stop("Column ", zc, " not in ", g_file)
  keep <- intersect(c("predictor_x", "outcome_y", "n_obs", "elpd_diff_full_minus_reduced",
                      "se_diff", "z_score", "se_diff_child", "z_child", "se_diff_block", "z_block",
                      "mean_coef_x_to_y", "pop_effect_mean", "pop_effect_q05", "pop_effect_q95",
                      "pop_prob_positive"), names(g))
  s1 <- g[, keep]
  s1$in_figure1 <- g$elpd_diff_full_minus_reduced > 0 & g$mean_coef_x_to_y > 0
  s1$strong     <- s1$in_figure1 & !is.na(g[[zc]]) & g[[zc]] > 2
  s1 <- s1[order(-s1$in_figure1, -g[[zc]]), ]
  names(s1)[names(s1) == "elpd_diff_full_minus_reduced"] <- "elpd_diff"
  write.csv(s1, file.path(out_dir, "table_S1_granger_pairs.csv"), row.names = FALSE)
  message(sprintf("Table S1: %d pairs; %d in Figure 1, %d strong (%s > 2)", nrow(s1),
                  sum(s1$in_figure1), sum(s1$strong), zc))
  for (z in intersect(c("z_score", "z_child", "z_block"), names(g)))
    message(sprintf("  strong by %-8s: %d", z, sum(s1$in_figure1 & g[[z]][match(paste(s1$predictor_x, s1$outcome_y),
                                                                                 paste(g$predictor_x, g$outcome_y))] > 2, na.rm = TRUE)))
} else message("Table S1 skipped: no Granger output in ", g_dir)

# ---------------------------------------------------------------- Table S2
step1_dir <- file.path(run_dir, "MODELS-STEP1")
s2b_file  <- file.path(step1_dir, paste0(prefix, ".rds"))
if (file.exists(s2b_file)) {
  r  <- readRDS(s2b_file)
  wm <- if (is.null(r$weight_mode)) "stacking" else r$weight_mode
  cf <- r$stacking_coefs
  is_A <- grepl("^A\\[", cf$term_label)
  cf$group <- ifelse(is_A, sub("^A\\[y=([^,]+),.*$", "\\1", cf$term_label), sub("^mu\\[(.+)\\]$", "\\1", cf$term_label))
  cf$predictor <- ifelse(is_A, sub("^A\\[y=[^,]+,x=([^,]+),.*$", "\\1", cf$term_label), NA_character_)
  cf$block_label <- ifelse(is_A, paste0("Lags: ", cf$predictor), "Mean")
  cf$block <- ifelse(is_A, paste0(cf$predictor, " -> ", cf$group), paste0("mean ", cf$group))
  s2 <- do.call(rbind, lapply(split(cf, cf$block), function(b) data.frame(
    block = b$block[1], group = b$group[1], block_label = b$block_label[1],
    inner_inclusion = max(b$inner_inclusion), group_inclusion = max(b$outer_inclusion),
    final_inclusion = max(b$inclusion_weight), stringsAsFactors = FALSE)))

  # Model-averaged effects from Step 2a.
  pre <- if (wm == "pseudobma") "pbma_" else "stack_"
  eff_cols <- paste0(pre, c("cond_mean", "cond_q05", "cond_q95", "cond_p_positive", "avg_q05", "avg_q95"))
  eff <- do.call(rbind, lapply(unique(s2$group), function(b) {
    f <- file.path(step1_dir, paste0("dsem_selection_", b, "_single_behaviour.rds"))
    if (!file.exists(f)) return(NULL)
    be <- readRDS(f)$block_effects
    if (is.null(be)) return(NULL)
    data.frame(group = b, be[, c("block_label", eff_cols)], stringsAsFactors = FALSE)
  }))
  if (is.null(eff)) {
    message("Table S2: Step 2a files have no block effects (re-run step2a); ",
            "effect columns left out")
  } else {
    names(eff) <- sub(paste0("^", pre), "", names(eff))
    names(eff)[names(eff) == "cond_mean"] <- "effect_mean_if_included"
    names(eff)[names(eff) == "cond_q05"]  <- "effect_q05_if_included"
    names(eff)[names(eff) == "cond_q95"]  <- "effect_q95_if_included"
    names(eff)[names(eff) == "cond_p_positive"] <- "p_effect_positive_if_included"
    names(eff)[names(eff) == "avg_q05"]   <- "model_averaged_q05"
    names(eff)[names(eff) == "avg_q95"]   <- "model_averaged_q95"
    s2 <- merge(s2, eff, by = c("group", "block_label"), all.x = TRUE, sort = FALSE)
  }
  pg <- file.path(run_dir, "PERMUTATIONS", "permutation_null_groups.csv")
  if (file.exists(pg)) {
    g <- read.csv(pg, stringsAsFactors = FALSE)
    names(g)[names(g) == "p"] <- "p_group"; names(g)[names(g) == "p_max"] <- "p_group_max"
    s2 <- merge(s2, g[, c("group", "p_group", "p_group_max")], by = "group", all.x = TRUE, sort = FALSE)
    message("Table S2: group permutation p-values added from ", pg)
  }
  sb <- file.path(run_dir, "PERMUTATIONS", "single_block_null.csv")
  if (file.exists(sb)) {
    b <- read.csv(sb, stringsAsFactors = FALSE)
    b <- b[, intersect(c("block", "elpd_gain", "elpd_gain_se", "effect_mean", "effect_q05", "effect_q95",
                         "p_single", "p_single_max", "q_single"), names(b))]
    names(b)[names(b) %in% c("effect_mean", "effect_q05", "effect_q95")] <-
      paste0("single_", names(b)[names(b) %in% c("effect_mean", "effect_q05", "effect_q95")])
    names(b)[2:3] <- c("single_elpd_gain", "single_elpd_gain_se")
    s2 <- merge(s2, b, by = "block", all.x = TRUE, sort = FALSE)
    message("Table S2: single-block test added from ", sb)
  }
  pn <- file.path(run_dir, "PERMUTATIONS", "permutation_null_blocks.csv")
  if (file.exists(pn)) {
    p <- read.csv(pn, stringsAsFactors = FALSE)
    s2 <- merge(s2, p[, c("block", "p_inner", "p_final", "p_final_max")], by = "block", all.x = TRUE, sort = FALSE)
    message("Table S2: two-level permutation p-values added from ", pn)
  }
  s2$block_label <- NULL
  s2 <- s2[order(-s2$final_inclusion), ]
  write.csv(s2, file.path(out_dir, "table_S2_step2_blocks.csv"), row.names = FALSE)
  message(sprintf("Table S2: %d blocks (weights: %s)", nrow(s2), wm))
} else message("Table S2 skipped: ", s2b_file, " not found")
message("Tables written to ", out_dir)
