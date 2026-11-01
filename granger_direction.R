# =============================================================================
# granger_direction.R
#
# Adds a population-level direction test to the Granger results.
#
# For each directed pair x -> y it takes the population lag coefficients
# A_pop[y, x, lag 1..5] from the Step 1 fits and, draw by draw, sums them over
# the five lags. That sum is the overall effect of x on y (change in log-odds
# of y when x was present in all five previous samples, for the average
# child). It reports the posterior mean and 90% interval of that sum and the
# posterior probability that it is above 0.
#
# This is a better-founded test of "x has a positive effect on y" than the
# share of positive draws pooled over lags, children and draws
# (prop_positive_x_to_y), which is pulled towards 50% by between-child
# variation and by lags with opposite signs.
#
# Usage:
#   Rscript granger_direction.R <run_dir>
#
# Input:   <run_dir>/MODELS-STEP1/dsem_step1_*.rds
#          <run_dir>/GRANGER/granger_loocv_pairwise_elpd.csv   (from granger_combine)
# Output:  <run_dir>/GRANGER/granger_direction.csv  — every pair, with the
#          Granger columns plus:
#            pop_effect_mean, pop_effect_q05, pop_effect_q95,
#            pop_prob_positive
#          and a table (printed) of how many relations pass at different
#          thresholds under each rule.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript granger_direction.R <run_dir>")
run_dir   <- args[1]
step1_dir <- file.path(run_dir, "MODELS-STEP1")
g_csv     <- file.path(run_dir, "GRANGER", "granger_loocv_pairwise_elpd.csv")
out_csv   <- file.path(run_dir, "GRANGER", "granger_direction.csv")
if (!file.exists(g_csv)) stop("Not found: ", g_csv, " (run granger and granger_combine first)")

suppressPackageStartupMessages(library(posterior))

files <- list.files(step1_dir, pattern = "^dsem_step1_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No dsem_step1_<n>.rds files in ", step1_dir)

# -----------------------------------------------------------------------------
# Read A_pop draws from every chain (only A_pop is kept from each file).
# -----------------------------------------------------------------------------
prep <- NULL
apop_list <- vector("list", length(files))
for (j in seq_along(files)) {
  message("Reading ", basename(files[j]))
  obj <- readRDS(files[j])
  if (is.null(prep)) prep <- obj$prep
  v <- variables(obj$draws)
  apop_vars <- v[grepl("^A_pop\\[", v)]
  apop_list[[j]] <- as_draws_matrix(subset_draws(obj$draws, variable = apop_vars))
  rm(obj); invisible(gc())
}
apop <- do.call(rbind, lapply(apop_list, unclass))     # draws x variables
rm(apop_list)
message("A_pop draws: ", nrow(apop), " from ", length(files), " chains")

behaviours <- prep$behaviours
lag_names  <- prep$lag_col_names                         # e.g. "gaze.l3", length K
V <- length(behaviours); K <- length(lag_names)
lag_beh <- sub("\\.l[0-9]+$", "", lag_names)             # predictor behaviour of each column
idx <- regmatches(colnames(apop), regexpr("\\[(\\d+),(\\d+)\\]", colnames(apop)))
row_i <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx))   # outcome y
col_k <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx))   # lag column

pop <- do.call(rbind, lapply(seq_len(V), function(y) do.call(rbind, lapply(seq_len(V), function(x) {
  if (x == y) return(NULL)
  cols <- which(row_i == y & lag_beh[col_k] == behaviours[x])
  if (length(cols) == 0) stop("No A_pop columns for ", behaviours[x], " -> ", behaviours[y])
  s <- rowSums(apop[, cols, drop = FALSE])               # sum over the lags, per draw
  data.frame(predictor_x = behaviours[x], outcome_y = behaviours[y],
             n_lags = length(cols),
             pop_effect_mean   = mean(s),
             pop_effect_q05    = unname(quantile(s, 0.05)),
             pop_effect_q95    = unname(quantile(s, 0.95)),
             pop_prob_positive = mean(s > 0),
             stringsAsFactors = FALSE)
}))))

g <- read.csv(g_csv, stringsAsFactors = FALSE)
out <- merge(g, pop, by = c("predictor_x", "outcome_y"), all.x = TRUE, sort = FALSE)
out <- out[order(-out$z_score), ]
if (any(is.na(out$pop_prob_positive))) warning("Some Granger pairs had no matching A_pop columns")
write.csv(out, out_csv, row.names = FALSE)

# -----------------------------------------------------------------------------
# How many relations pass under each rule?
#   All rules also require the Granger criterion: ELPD difference > 0.
# -----------------------------------------------------------------------------
imp <- out$elpd_diff_full_minus_reduced > 0
strong <- out$z_score > 2
count <- function(keep) c(shown = sum(imp & keep), strong = sum(imp & keep & strong))

tab <- rbind(
  "Mean pooled coefficient > 0 (current rule)"  = count(out$mean_coef_x_to_y > 0),
  "Pooled share of draws positive >= 60%"        = count(out$prop_positive_x_to_y >= 0.60),
  "Population effect: P(sum of lags > 0) >= 80%" = count(out$pop_prob_positive >= 0.80),
  "Population effect: P(sum of lags > 0) >= 90%" = count(out$pop_prob_positive >= 0.90),
  "Population effect: P(sum of lags > 0) >= 95%" = count(out$pop_prob_positive >= 0.95),
  "Population effect: P(sum of lags > 0) >= 99%" = count(out$pop_prob_positive >= 0.99))

cat("\nPairs improving prediction (ELPD difference > 0): ", sum(imp),
    "   of which z > 2: ", sum(imp & strong), "\n\n", sep = "")
cat("Relations included under each 'positive' rule (all also need ELPD difference > 0):\n")
print(tab)

cat("\nRelations that improve prediction, by population-level P(effect > 0):\n")
show <- out[imp, c("predictor_x", "outcome_y", "z_score", "pop_effect_mean",
                   "pop_effect_q05", "pop_effect_q95", "pop_prob_positive",
                   "mean_coef_x_to_y", "prop_positive_x_to_y")]
show <- show[order(-show$pop_prob_positive, -show$z_score), ]
print(format(show, digits = 3), row.names = FALSE)
cat("\nSaved: ", out_csv, "\n", sep = "")
