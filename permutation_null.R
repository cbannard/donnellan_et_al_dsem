# =============================================================================
# permutation_null.R
#
# Empirical (permutation) null for the Step 2 results.
#
# Each permutation run (stages perm2a + perm2b of run_dsem_pipeline.slurm)
# repeats Step 2a and Step 2b exactly as in the observed run, but with the
# children's vocabulary scores shuffled across children. Shuffling breaks any
# link between a child's interaction dynamics and their vocabulary while
# keeping everything else (the Step 1 estimates and their correlations, the
# model search, the weighting). The distribution of each statistic over the
# permutations is what that statistic looks like when there is nothing to
# find, given this search over many models and 58 children.
#
# This script compares the observed run with those permutation runs:
#   blocks  — for every predictor block (a lag block x -> y, or the mean of y):
#             within-group (inner) inclusion, group (outer) inclusion and
#             final inclusion; p = share of permutations with a value at least
#             as large (with the usual +1 correction); p_final_max compares the
#             observed final inclusion with the LARGEST final inclusion in each
#             permutation run, which corrects for looking at all blocks.
#   groups  — outer inclusion per behaviour group.
#   global  — weight of the intercept-only model (small = signal; p = share
#             of permutations with weight at most as large), the ELPD
#             difference of the best model over the intercept-only model, and
#             the LARGEST final / inner / group inclusion probability in the
#             run. The 95th percentile of the null for max_final_inclusion is
#             the permutation-calibrated inclusion threshold: in only 5% of
#             runs with no real association does any block reach it
#             (make_figures.R threshold=perm uses it).
#
# Usage:
#   Rscript permutation_null.R <run_dir> [prefix]
#     <run_dir>/MODELS-STEP1/   observed Step 2a/2b outputs
#     <run_dir>/PERMUTATIONS/perm_<k>/   permutation outputs
#   prefix: Step 2b output prefix (default dsem_selection_from_single_behaviour_sets)
#
# Outputs in <run_dir>/PERMUTATIONS/: permutation_null_blocks.csv,
#   permutation_null_groups.csv, permutation_null_global.csv,
#   permutation_null_report.txt
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript permutation_null.R <run_dir> [prefix]")
run_dir <- args[1]
prefix  <- if (length(args) >= 2) args[2] else "dsem_selection_from_single_behaviour_sets"
obs_dir  <- file.path(run_dir, "MODELS-STEP1")
# Variant analyses (other Step 2b prefixes) get their own output files.
tag <- sub("^dsem_selection_from_single_behaviour_sets_?", "", prefix)
out_name <- function(stem, ext) file.path(file.path(run_dir, "PERMUTATIONS"),
  paste0("permutation_null", if (nzchar(tag)) paste0("_", tag) else "", "_", stem, ".", ext))
perm_dir <- file.path(run_dir, "PERMUTATIONS")

# One Step 2b result -> list(blocks, groups, global)
read_run <- function(dir) {
  f <- file.path(dir, paste0(prefix, ".rds"))
  if (!file.exists(f)) return(NULL)
  r  <- readRDS(f)
  wm <- if (is.null(r$weight_mode)) "stacking" else r$weight_mode
  wcol <- if (wm == "pseudobma") "pseudobma_weight" else "stacking_weight"
  cf <- r$stacking_coefs
  is_A <- grepl("^A\\[", cf$term_label)
  cf$group <- ifelse(is_A, sub("^A\\[y=([^,]+),.*$", "\\1", cf$term_label),
                     sub("^mu\\[(.+)\\]$", "\\1", cf$term_label))
  cf$block <- ifelse(is_A,
                     sub("^A\\[y=([^,]+),x=([^,]+),lag=[0-9]+\\]$", "\\2 -> \\1", cf$term_label),
                     paste0("mean ", cf$group))
  blocks <- do.call(rbind, lapply(split(cf, cf$block), function(b) data.frame(
    block = b$block[1], group = b$group[1], inner = max(b$inner_inclusion),
    outer = max(b$outer_inclusion), final = max(b$inclusion_weight),
    direction = b$effect_direction[which.max(b$inclusion_weight)], stringsAsFactors = FALSE)))
  # Step 2b-only runs (no_selection=1) have no within-group weights: their
  # block-level inner/final values are unavailable and only the group and
  # global statistics are used from them.
  if (isTRUE(r$no_selection)) { blocks$inner <- NA_real_; blocks$final <- NA_real_ }
  ss <- r$stacking_summary
  list(blocks = blocks,
       groups = data.frame(group = names(r$outer_inclusion), outer = as.numeric(r$outer_inclusion),
                           stringsAsFactors = FALSE),
       global = c(intercept_weight = ss[[wcol]][ss$behaviours == "intercept_only"][1],
                  best_elpd_diff   = max(ss$elpd_diff[ss$behaviours != "intercept_only"], na.rm = TRUE)),
       weight_mode = wm, permute_seed = r$permute_seed, block_mode = r$block_mode)
}

obs <- read_run(obs_dir)
if (is.null(obs)) stop("Observed Step 2b output not found in ", obs_dir)
perm_dirs <- list.files(perm_dir, pattern = "^perm_[0-9]+$", full.names = TRUE)
perms <- lapply(perm_dirs, read_run)
ok <- !vapply(perms, is.null, logical(1))
message(sum(ok), " of ", length(perm_dirs), " permutation runs complete")
perms <- perms[ok]
if (length(perms) == 0) stop("No completed permutation runs under ", perm_dir)
modes <- unique(vapply(perms, `[[`, "", "weight_mode"))
if (!all(modes == obs$weight_mode))
  warning("Weight mode differs between observed (", obs$weight_mode, ") and permutation runs (",
          paste(modes, collapse = ", "), ")")
n <- length(perms)
p_ge <- function(obs_val, null_vals) (1 + sum(null_vals >= obs_val - 1e-12, na.rm = TRUE)) / (1 + sum(!is.na(null_vals)))
p_le <- function(obs_val, null_vals) (1 + sum(null_vals <= obs_val + 1e-12, na.rm = TRUE)) / (1 + sum(!is.na(null_vals)))

# Null values for each block (a block missing from a run counts as 0).
null_block <- function(blk, what) vapply(perms, function(p) {
  v <- p$blocks[[what]][p$blocks$block == blk]
  if (length(v)) v[1] else 0
}, numeric(1))
null_max_final <- vapply(perms, function(p)
  if (all(is.na(p$blocks$final))) NA_real_ else max(c(0, p$blocks$final), na.rm = TRUE), numeric(1))
n_block_level <- sum(!is.na(null_max_final))
q_na <- function(x, pr) if (all(is.na(x))) NA_real_ else unname(quantile(x, pr, na.rm = TRUE))

ob <- obs$blocks
blocks_out <- do.call(rbind, lapply(seq_len(nrow(ob)), function(j) {
  nb_in <- null_block(ob$block[j], "inner"); nb_fi <- null_block(ob$block[j], "final")
  data.frame(ob[j, c("block", "group", "direction")],
             inner = ob$inner[j], inner_null_q95 = q_na(nb_in, 0.95), p_inner = p_ge(ob$inner[j], nb_in),
             outer = ob$outer[j],
             final = ob$final[j], final_null_q95 = q_na(nb_fi, 0.95), p_final = p_ge(ob$final[j], nb_fi),
             p_final_max = p_ge(ob$final[j], null_max_final),
             stringsAsFactors = FALSE, row.names = NULL)
}))
blocks_out <- blocks_out[order(-blocks_out$final), ]

groups_out <- do.call(rbind, lapply(seq_len(nrow(obs$groups)), function(j) {
  g <- obs$groups$group[j]
  nv <- vapply(perms, function(p) { v <- p$groups$outer[p$groups$group == g]; if (length(v)) v[1] else 0 }, numeric(1))
  data.frame(group = g, outer = obs$groups$outer[j], null_q95 = unname(quantile(nv, 0.95)),
             p = p_ge(obs$groups$outer[j], nv), stringsAsFactors = FALSE)
}))
null_max_outer <- vapply(perms, function(p) max(c(0, p$groups$outer)), numeric(1))
groups_out$p_max <- vapply(groups_out$outer, p_ge, numeric(1), null_vals = null_max_outer)
groups_out <- groups_out[order(-groups_out$outer), ]

null_int  <- vapply(perms, function(p) p$global[["intercept_weight"]], numeric(1))
null_best <- vapply(perms, function(p) p$global[["best_elpd_diff"]], numeric(1))
null_max_inner <- vapply(perms, function(p)
  if (all(is.na(p$blocks$inner))) NA_real_ else max(c(0, p$blocks$inner), na.rm = TRUE), numeric(1))
g_row <- function(stat, obs_val, nv, lower_is_signal, note) data.frame(
  statistic = stat, observed = obs_val,
  null_q05 = q_na(nv, 0.05), null_median = q_na(nv, 0.5), null_q95 = q_na(nv, 0.95),
  n_perm = sum(!is.na(nv)),
  p = if (lower_is_signal) p_le(obs_val, nv) else p_ge(obs_val, nv), note = note,
  stringsAsFactors = FALSE, row.names = NULL)
global_out <- rbind(
  g_row("intercept_only_weight", obs$global[["intercept_weight"]], null_int, TRUE,
        "p = share of permutations with intercept-only weight <= observed (small weight = signal)"),
  g_row("best_model_elpd_diff", obs$global[["best_elpd_diff"]], null_best, FALSE,
        "p = share of permutations with best-model ELPD gain >= observed"),
  g_row("max_final_inclusion", max(c(0, ob$final)), null_max_final, FALSE,
        "largest final inclusion of any block; null_q95 = permutation-calibrated threshold"),
  g_row("max_inner_inclusion", max(c(0, ob$inner)), null_max_inner, FALSE,
        "largest within-group inclusion of any block"),
  g_row("max_group_inclusion", max(c(0, obs$groups$outer)), null_max_outer, FALSE,
        "largest group (outer) inclusion"))

write.csv(blocks_out, out_name("blocks", "csv"), row.names = FALSE)
write.csv(groups_out, out_name("groups", "csv"), row.names = FALSE)
write.csv(global_out, out_name("global", "csv"), row.names = FALSE)

thr <- global_out$null_q95[global_out$statistic == "max_final_inclusion"]
block_level <- is.finite(thr) && any(is.finite(blocks_out$final))
top <- blocks_out[1, ]
calib_txt <- if (!block_level) "Block-level statistics are not available for this run (Step 2b only, no Step 2a).\n" else sprintf(paste0(
  "Permutation-calibrated threshold for final inclusion (95th percentile of the largest final\n",
  "inclusion in each permutation run): %.3f. Blocks above it: %d.\n",
  "Top observed block: %s (final inclusion %.3f); the largest block in a permutation run\n",
  "reached at least this in %.0f%% of runs (familywise p = %.3f).\n"),
  thr, sum(blocks_out$final > thr), top$block, top$final,
  100 * mean(null_max_final >= top$final - 1e-12, na.rm = TRUE), top$p_final_max)

fmt <- function(df) paste(capture.output(print(format(df, digits = 3), row.names = FALSE)), collapse = "\n")
report <- paste0(
  "Permutation null for Step 2", if (nzchar(tag)) paste0(" [", tag, "]") else "", " (", n, " permutations; weights = ", obs$weight_mode,
  ", blocks = ", obs$block_mode, ")\n",
  "Smallest attainable p-value with ", n, " permutations: ", signif(1 / (n + 1), 3), "\n",
  "Block-level statistics (inner, final) use the ", n_block_level, " runs with Step 2a;\n",
  "group and global statistics use all ", n, " runs.\n\n",
  calib_txt, "\n",
  "GLOBAL\n", fmt(global_out[, c("statistic", "observed", "null_q05", "null_median", "null_q95", "p", "n_perm")]), "\n\n",
  "BEHAVIOUR GROUPS (outer inclusion)\n", fmt(groups_out), "\n\n",
  "BLOCKS (top 20 by final inclusion)\n", fmt(head(blocks_out, 20)), "\n\n",
  "p_inner / p_final: vs the same block in the permutation runs;\n",
  "p_final_max: vs the largest final inclusion of any block in each run (corrects for all blocks);\n",
  "p_max (groups): vs the largest group inclusion in each run.\n")
writeLines(report, out_name("report", "txt"))
cat(report)
