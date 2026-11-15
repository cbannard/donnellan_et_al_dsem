# =============================================================================
# single_block_test.R
#
# Per-predictor test of usefulness for vocabulary, with a permutation null.
#
# Question: does this relation (one predictor block — the five lags of
# x -> y, or the mean rate of y), on its own, predict later vocabulary better
# than an intercept-only model, beyond what chance would give?
#
# For every one of the 90 blocks (9 behaviour groups x [9 lag blocks + mean])
# this fits the Step 2 error-in-variables regression with that block alone and
# compares it with the intercept-only model by PSIS-LOO (leave-one-child-out):
#   elpd_gain = ELPD(block alone) - ELPD(intercept only)
# It does this for the observed vocabulary (run 0) and for the same shuffles
# of vocabulary as the Step 2 permutation runs (run k uses permute_seed =
# perm_seed_base + k, so run k matches PERMUTATIONS/perm_<k>).
#
# Unlike the two-level inclusion probabilities, the result for a block does
# not depend on its behaviour group's support or on other blocks sharing
# its information. Two correlated blocks can therefore both be useful while
# carrying the same information.
#
# Modes:
#   Rscript single_block_test.R fit <run_dir> <k> [perm_seed_base=20261112] [chains=4] [iter=500]
#       k = 0: observed vocabulary; k >= 1: shuffle k.
#       Writes <run_dir>/PERMUTATIONS/single_block/run_<kkk>.csv
#   Rscript single_block_test.R summary <run_dir>
#       Writes <run_dir>/PERMUTATIONS/single_block_null.csv and
#       single_block_report.txt:
#         elpd_gain, its SE and z (observed);
#         effect_mean / effect_q05 / effect_q95 / p_effect_positive — the
#           block's effect on log vocabulary (sum of its coefficients) in the
#           single-block model;
#         null_q95 — 95th percentile of the same block's gain under shuffling;
#         p_single — share of shuffles where this block's gain >= observed;
#         p_single_max — share of shuffles where the LARGEST gain of any of the
#           90 blocks >= observed (familywise: corrects for testing 90 blocks);
#         q_single — permutation false discovery rate: if this block and every
#           block with a larger gain were declared real, the expected share of
#           them that would be chance findings (average number of shuffled
#           blocks with a gain at least this large / number of observed ones),
#           made monotone (Storey & Tibshirani, 2003).
#       The report also compares the number of blocks with a positive gain and
#       the sum of positive gains with their permutation distributions.
#
# Needs: <run_dir>/PERMUTATIONS/bridge.rds (make_bridge.R) and ./CDI_data.csv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript single_block_test.R fit|summary <run_dir> [k] [...]")
mode    <- args[1]
run_dir <- args[2]
out_dir <- file.path(run_dir, "PERMUTATIONS", "single_block")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
keyed_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", key, "="), "", hit[length(hit)]) else default
}

# =============================================================================
# FIT one run (observed or one shuffle)
# =============================================================================
if (mode == "fit") {
  k <- as.integer(args[3])
  if (is.na(k) || k < 0) stop("k must be 0 (observed) or a positive permutation number")
  perm_seed_base <- as.integer(keyed_arg("perm_seed_base", 20261112))
  n_chains <- as.integer(keyed_arg("chains", 4))
  n_iter   <- as.integer(keyed_arg("iter", 500))
  out_csv  <- file.path(out_dir, sprintf("run_%03d.csv", k))
  if (file.exists(out_csv)) { message("Exists, skipping: ", out_csv); quit(save = "no") }

  suppressPackageStartupMessages({ library(cmdstanr); library(posterior); library(loo); library(dplyr) })

  bridge <- readRDS(file.path(run_dir, "PERMUTATIONS", "bridge.rds"))
  prep <- bridge$prep
  behaviours <- prep$behaviours; V <- prep$V; K <- prep$K
  K_per_pair <- as.integer(K / V)

  cdi <- read.csv("CDI_data.csv", stringsAsFactors = FALSE) |>
    mutate(participant = as.character(participant)) |>
    rename(participant_id = participant)
  merged <- left_join(data.frame(participant_id = as.character(prep$participants), stringsAsFactors = FALSE),
                      cdi[, c("participant_id", "us18_imp")], by = "participant_id")
  has_outcome <- !is.na(merged$us18_imp)
  y <- log1p(as.numeric(merged$us18_imp[has_outcome]))
  Zm <- bridge$z_hat[has_outcome, , drop = FALSE]
  Zs <- bridge$z_sd[has_outcome, , drop = FALSE]

  # Same shuffle as select_terms_for_individual_behaviours.R / final_term_weighting.R
  if (k > 0) {
    set.seed(perm_seed_base + k)
    y <- y[sample.int(length(y))]
    message("Shuffle ", k, " (permute_seed=", perm_seed_base + k, ")")
  } else message("Observed vocabulary")

  # The 90 blocks: q-columns of the z_person vector.
  blocks <- do.call(rbind, lapply(seq_len(V), function(v) {
    lag_rows <- lapply(seq_len(V), function(x) {
      ks <- (x - 1L) * K_per_pair + seq_len(K_per_pair)
      data.frame(group = behaviours[v], block = paste0(behaviours[x], " -> ", behaviours[v]),
                 cols = I(list((ks - 1L) * V + v)), stringsAsFactors = FALSE)
    })
    rbind(do.call(rbind, lag_rows),
          data.frame(group = behaviours[v], block = paste0("mean ", behaviours[v]),
                     cols = I(list(V * K + v)), stringsAsFactors = FALSE))
  }))

  stan_intercept <- "
data { int<lower=1> N; vector[N] y; real<lower=0> prior_sd_sigma; }
parameters { real alpha; real<lower=0> sigma; }
model { alpha ~ normal(0, 10); sigma ~ normal(0, prior_sd_sigma); y ~ normal(alpha, sigma); }
generated quantities { vector[N] log_lik; for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | alpha, sigma); }
"
  # Same error-in-variables model as Step 2a/2b.
  stan_normal <- "
data { int<lower=1> N; int<lower=1> P; vector[N] y; matrix[N, P] X_mean; matrix[N, P] X_sd;
       real<lower=0> prior_sd_sigma; }
parameters { real alpha; vector[P] beta; real<lower=0> sigma; }
model {
  beta ~ normal(0, 0.5); alpha ~ normal(0, 10); sigma ~ normal(0, prior_sd_sigma);
  for (n in 1:N) {
    real mu_n = alpha + dot_product(X_mean[n], beta);
    real sd_n = sqrt(square(sigma) + dot_product(square(X_sd[n]), square(beta)));
    y[n] ~ normal(mu_n, sd_n);
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n = alpha + dot_product(X_mean[n], beta);
    real sd_n = sqrt(square(sigma) + dot_product(square(X_sd[n]), square(beta)));
    log_lik[n] = normal_lpdf(y[n] | mu_n, sd_n);
  }
}
"
  mod_int <- cmdstan_model(write_stan_file(stan_intercept))
  mod_nrm <- cmdstan_model(write_stan_file(stan_normal))
  fit_seed <- (perm_seed_base %% 100000L) * 1000L + k * 100L
  fit_loo <- function(mod, data, seed) {
    fit <- mod$sample(data = data, seed = seed, chains = n_chains, iter_warmup = n_iter,
                      iter_sampling = n_iter, refresh = 0, show_messages = FALSE)
    ll <- fit$draws("log_lik", format = "matrix")
    list(fit = fit, loo = suppressWarnings(loo::loo(ll, cores = 1)))
  }
  N <- length(y)
  base <- fit_loo(mod_int, list(N = N, y = y, prior_sd_sigma = 2.5), fit_seed)
  base_pw <- base$loo$pointwise[, "elpd_loo"]

  res <- vector("list", nrow(blocks))
  for (j in seq_len(nrow(blocks))) {
    cols <- blocks$cols[[j]]
    r <- tryCatch(fit_loo(mod_nrm, list(N = N, P = length(cols), y = y,
                                        X_mean = Zm[, cols, drop = FALSE], X_sd = Zs[, cols, drop = FALSE],
                                        prior_sd_sigma = 2.5), fit_seed + j),
                  error = function(e) NULL)
    if (is.null(r)) {
      res[[j]] <- data.frame(elpd_gain = NA_real_, elpd_gain_se = NA_real_, effect_mean = NA_real_,
                             effect_q05 = NA_real_, effect_q95 = NA_real_, p_effect_positive = NA_real_)
      next
    }
    d <- r$loo$pointwise[, "elpd_loo"] - base_pw
    eff <- rowSums(r$fit$draws("beta", format = "matrix"))
    res[[j]] <- data.frame(elpd_gain = sum(d), elpd_gain_se = sqrt(length(d) * var(d)),
                           effect_mean = mean(eff), effect_q05 = unname(quantile(eff, 0.05)),
                           effect_q95 = unname(quantile(eff, 0.95)), p_effect_positive = mean(eff > 0))
    if (j %% 10 == 0) message("Progress ", j, "/", nrow(blocks))
  }
  out <- cbind(run = k, blocks[, c("group", "block")], do.call(rbind, res))
  write.csv(out, out_csv, row.names = FALSE)
  message("Saved ", out_csv)
  quit(save = "no")
}

# =============================================================================
# SUMMARY: observed vs shuffles
# =============================================================================
if (mode == "summary") {
  files <- list.files(out_dir, pattern = "^run_[0-9]+\\.csv$", full.names = TRUE)
  runs <- lapply(files, read.csv, stringsAsFactors = FALSE)
  ks <- vapply(runs, function(r) r$run[1], numeric(1))
  if (!0 %in% ks) stop("Observed run (run_000.csv) not found in ", out_dir)
  obs <- runs[[which(ks == 0)]]
  perms <- runs[ks > 0]
  n <- length(perms)
  if (n == 0) stop("No shuffled runs found in ", out_dir)
  message(n, " shuffled runs")
  null_mat <- sapply(perms, function(r) r$elpd_gain[match(obs$block, r$block)])   # blocks x runs
  null_max <- apply(null_mat, 2, max, na.rm = TRUE)
  p_ge <- function(o, nv) (1 + sum(nv >= o - 1e-12, na.rm = TRUE)) / (1 + sum(!is.na(nv)))
  out <- obs[, c("group", "block", "elpd_gain", "elpd_gain_se", "effect_mean", "effect_q05",
                 "effect_q95", "p_effect_positive")]
  out$z <- out$elpd_gain / out$elpd_gain_se
  out$null_q95 <- apply(null_mat, 1, quantile, 0.95, na.rm = TRUE)
  out$p_single <- vapply(seq_len(nrow(out)), function(j) p_ge(out$elpd_gain[j], null_mat[j, ]), numeric(1))
  out$p_single_max <- vapply(out$elpd_gain, p_ge, numeric(1), nv = null_max)
  out <- out[order(-out$elpd_gain), ]
  # Permutation FDR / q-values over the ranking by gain.
  null_exceed <- vapply(out$elpd_gain, function(g) mean(colSums(null_mat >= g - 1e-12, na.rm = TRUE)), numeric(1))
  fdr <- pmin(1, null_exceed / seq_len(nrow(out)))
  out$q_single <- rev(cummin(rev(fdr)))
  write.csv(out, file.path(run_dir, "PERMUTATIONS", "single_block_null.csv"), row.names = FALSE)
  n_pos <- sum(out$elpd_gain > 0, na.rm = TRUE)
  null_pos <- colSums(null_mat > 0, na.rm = TRUE)
  sum_pos <- sum(pmax(out$elpd_gain, 0), na.rm = TRUE)
  null_sum <- colSums(pmax(null_mat, 0), na.rm = TRUE)

  thr <- unname(quantile(null_max, 0.95))
  fmt <- function(df) paste(capture.output(print(format(df, digits = 3), row.names = FALSE)), collapse = "\n")
  rep <- paste0(
    "Single-block test: each block alone vs intercept-only (ELPD gain, leave-one-child-out)\n",
    n, " shuffles of vocabulary; smallest attainable p = ", signif(1 / (n + 1), 3), "\n\n",
    sprintf("Familywise threshold (95th percentile of the largest gain of any block per shuffle): %.2f\n", thr),
    sprintf("Blocks above it (p_single_max < .05): %d\n", sum(out$p_single_max < 0.05)),
    sprintf("Blocks with p_single < .05 (uncorrected): %d (about %.1f expected by chance)\n",
            sum(out$p_single < 0.05), 0.05 * nrow(out)),
    sprintf("Blocks with a positive gain: %d (shuffles: median %.0f, 95th percentile %.0f; p = %.3f)\n",
            n_pos, median(null_pos), quantile(null_pos, 0.95), p_ge(n_pos, null_pos)),
    sprintf("Sum of positive gains: %.2f (shuffles: median %.2f, 95th percentile %.2f; p = %.3f)\n",
            sum_pos, median(null_sum), quantile(null_sum, 0.95), p_ge(sum_pos, null_sum)),
    sprintf("Blocks with q_single < .25 / < .5: %d / %d\n\n", sum(out$q_single < 0.25), sum(out$q_single < 0.5)),
    "TOP 25 BLOCKS\n",
    fmt(head(out[, c("block", "elpd_gain", "elpd_gain_se", "effect_mean", "effect_q05", "effect_q95",
                     "null_q95", "p_single", "p_single_max", "q_single")], 25)), "\n\n",
    "elpd_gain: ELPD(block alone) - ELPD(intercept only); positive = better held-out prediction.\n",
    "effect_*: block effect on log vocabulary (sum of its coefficients) in the single-block model.\n",
    "p_single: vs the same block under shuffling; p_single_max: vs the best of all 90 blocks per shuffle;\n",
    "q_single: expected share of chance findings among this block and all blocks ranked above it.\n")
  writeLines(rep, file.path(run_dir, "PERMUTATIONS", "single_block_report.txt"))
  cat(rep)
  quit(save = "no")
}
stop("mode must be fit or summary")
