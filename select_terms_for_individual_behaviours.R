# =============================================================================
# select_terms_for_individual_behaviours.R
#
# Finds best sets of predictors from VAR for a given behaviours
#
# Usage:
#   Rscript select_terms_for_individual_behaviours.R <behaviour> [input_folder]
#
# Example:
#   Rscript select_terms_for_individual_behaviours.R gaze
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript select_terms_for_individual_behaviours.R <behaviour> [input_folder]")
}

target_behaviour <- args[1]
input_folder <- if (length(args) >= 2) args[2] else "MODELS-AGGREGATE-MODELSONLY"

is_dir <- function(path) {
  info <- suppressWarnings(file.info(path)$isdir)
  if (is.na(info)) FALSE else info
}

combine_step1_results <- function(rds_list) {
  if (length(rds_list) == 1) return(rds_list[[1]])
  draws_list <- lapply(rds_list, function(x) posterior::as_draws_matrix(x$draws))
  prep <- rds_list[[1]]$prep
  draws_combined <- posterior::bind_draws(draws_list, along = "draw")
  list(draws = draws_combined, prep = prep)
}

if (!is_dir(input_folder)) stop("Provided path is not a directory: ", input_folder)
rds_files <- list.files(input_folder, pattern = "^dsem_step1_.*\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
  rds_files <- list.files(input_folder, pattern = "^dsem_result_.*\\.rds$", full.names = TRUE)
}
if (length(rds_files) == 0) stop("No step-1 RDS files (dsem_step1_*.rds) found in folder: ", input_folder)
message("Loading and combining RDS files from folder: ", input_folder)
step1_results <- lapply(rds_files, readRDS)
step1_result  <- combine_step1_results(step1_results)

library(cmdstanr)
library(posterior)
library(loo)
library(dplyr)

stan_intercept <- "
data {
  int<lower=1> N;
  vector[N] y;
  real<lower=0> prior_sd_sigma;
}
parameters {
  real alpha;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 10);
  sigma ~ normal(0, prior_sd_sigma);
  y ~ normal(alpha, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | alpha, sigma);
}
"

stan_normal <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] y;
  matrix[N, P] X_mean;
  matrix[N, P] X_sd;
  real<lower=0> prior_sd_sigma;
}
parameters {
  real           alpha;
  vector[P]      beta;
  real<lower=0>  sigma;
}
model {
  beta  ~ normal(0, 0.5);
  alpha ~ normal(0, 10);
  sigma ~ normal(0, prior_sd_sigma);
  for (n in 1:N) {
    real mu_n  = alpha + dot_product(X_mean[n], beta);
    real sd_n  = sqrt(square(sigma) +
                      dot_product(square(X_sd[n]), square(beta)));
    y[n] ~ normal(mu_n, sd_n);
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    real mu_n  = alpha + dot_product(X_mean[n], beta);
    real sd_n  = sqrt(square(sigma) +
                      dot_product(square(X_sd[n]), square(beta)));
    log_lik[n] = normal_lpdf(y[n] | mu_n, sd_n);
  }
}
"

step1_prep <- step1_result$prep
behaviours <- step1_prep$behaviours
lag_names  <- step1_prep$lag_col_names
V          <- step1_prep$V
K          <- step1_prep$K
Q          <- step1_prep$Q

if (!(target_behaviour %in% behaviours)) {
  stop(
    "Unknown behaviour '", target_behaviour, "'. Available behaviours: ",
    paste(behaviours, collapse = ", ")
  )
}

target_b_idx <- match(target_behaviour, behaviours)

stopifnot("K must be divisible by V" = K %% V == 0L)
K_per_pair <- as.integer(K / V)

lag_blocks <- lapply(seq_len(V), function(pred_v) {
  seq_len(K_per_pair) + (pred_v - 1L) * K_per_pair
})
names(lag_blocks) <- behaviours
block_names <- c(paste0("Lags: ", behaviours[seq_len(V)]), "Mean", "Log-var")

message(
  "Target behaviour: ", target_behaviour,
  "  N=", step1_prep$stan_data$N,
  "  V=", V,
  "  K_total=", K,
  "  K_per_pair=", K_per_pair,
  "  Q=", Q
)

bridge_step1 <- function(step1_result) {
  draws <- step1_result$draws
  prep  <- step1_result$prep
  N     <- prep$stan_data$N
  Q     <- prep$Q

  all_vars <- posterior::variables(draws)
  z_vars   <- all_vars[grepl("^z_person\\[", all_vars)]
  stopifnot(length(z_vars) == N * Q)

  idx_mat <- regmatches(z_vars, regexpr("\\[(\\d+),(\\d+)\\]", z_vars))
  i_idx   <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx_mat))
  q_idx   <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx_mat))

  z_mat <- posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = z_vars)
  )
  message("  Bridge: S=", nrow(z_mat), " draws  N=", N, "  Q=", Q)

  z_hat <- matrix(0, nrow = N, ncol = Q)
  z_sd  <- matrix(0, nrow = N, ncol = Q)
  for (col_j in seq_len(ncol(z_mat))) {
    ii <- i_idx[col_j]
    qq <- q_idx[col_j]
    z_hat[ii, qq] <- mean(z_mat[, col_j])
    z_sd[ii, qq]  <- sd(z_mat[, col_j])
  }

  z_sd <- pmax(z_sd, 1e-4)
  list(z_hat = z_hat, z_sd = z_sd)
}

bridge <- bridge_step1(step1_result)

cdi_imp <- read.csv(
  "CDI_data.csv",
  header = TRUE,
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(participant = as.character(participant)) |>
  dplyr::rename(participant_id = participant)

participants <- step1_prep$participants

pid_df <- data.frame(
  participant_id = as.character(participants),
  step1_row      = seq_along(participants),
  stringsAsFactors = FALSE
)
merged <- dplyr::left_join(
  pid_df,
  cdi_imp[, c("participant_id", "us18_imp")],
  by = "participant_id"
)

has_outcome <- !is.na(merged$us18_imp)
out_idx     <- which(has_outcome)        # indexes into merged / step1 rows
N_out       <- length(out_idx)
y_out       <- log1p(as.numeric(merged$us18_imp[has_outcome]))

Z_mean <- bridge$z_hat[out_idx, , drop = FALSE]
Z_sd   <- bridge$z_sd[out_idx,  , drop = FALSE]

message(
  "Outcome sample: N_out=", N_out,
  "  y range=[", round(min(y_out), 2), ", ", round(max(y_out), 2), "]"
)

group_indices <- lapply(seq_len(V), function(v) {
  A_cols <- as.integer(unlist(lapply(seq_len(V), function(pred_v) {
    sapply(seq_len(K_per_pair), function(l) {
      k <- (pred_v - 1L) * K_per_pair + l
      (k - 1L) * V + v
    })
  })))
  list(
    A   = A_cols,
    mu  = V * K + v,
    lv  = V * K + V + v,
    all = c(A_cols, V * K + v, V * K + V + v)
  )
})
names(group_indices) <- behaviours

message("Compiling Stan models...")
mod_intercept <- cmdstan_model(write_stan_file(stan_intercept))
mod_normal    <- cmdstan_model(write_stan_file(stan_normal))
message("  Done.")

fit_and_loo <- function(mod, data_list,
                        chains        = 4,
                        iter_warmup   = 500,
                        iter_sampling = 500,
                        refresh       = 0) {
  fit <- mod$sample(
    data          = data_list,
    chains        = chains,
    iter_warmup   = iter_warmup,
    iter_sampling = iter_sampling,
    refresh       = refresh,
    show_messages = FALSE
  )
  ll <- fit$draws("log_lik", format = "matrix")
  loo_obj <- loo::loo(ll, cores = 1)
  list(fit = fit, loo = loo_obj)
}

prior_sd_sigma <- 2.5
min_n_frac <- 0.80

group_selection <- list()

b <- target_behaviour
q_cols <- group_indices[[b]]$all
n_terms <- length(q_cols)

message("\n", strrep("=", 60))
message(
  "Group: ", b, "  (", n_terms, " candidate terms: ",
  K_per_pair, " lag-coef per predictor + 1 mean + 1 log-var)"
)
message(strrep("=", 60))

block_defs <- c(
  lapply(seq_len(V), function(i) lag_blocks[[i]]),
  list(K + 1L),
  list(K + 2L)
)
n_blocks <- length(block_defs)

all_combos <- unlist(lapply(1:n_blocks, function(k) {
  combn(seq_len(n_blocks), k, simplify = FALSE)
}), recursive = FALSE)

# ---------------------------------------------------------------------------
# Shared complete-case mask across all combinations.
# stacking_weights() requires all models to be evaluated on the same
# observation set so that pointwise log-lik matrices are conformable.
# ---------------------------------------------------------------------------
all_q_cols_for_behaviour <- sort(unique(unlist(block_defs)))
all_trial_cols <- q_cols[all_q_cols_for_behaviour]
shared_mask <- complete.cases(Z_mean[, all_trial_cols, drop = FALSE]) &
               complete.cases(Z_sd[,   all_trial_cols, drop = FALSE])
N_shared <- sum(shared_mask)
message(sprintf("  Shared complete-case N: %d / %d", N_shared, N_out))

if (N_shared < ceiling(min_n_frac * N_out)) {
  stop(sprintf(
    "Shared complete-case N (%d) is below min_n_frac threshold (%.0f%% of %d).",
    N_shared, 100 * min_n_frac, N_out
  ))
}

y_shared   <- y_out[shared_mask]
Zm_shared  <- Z_mean[shared_mask, , drop = FALSE]
Zsd_shared <- Z_sd[shared_mask,  , drop = FALSE]

# Fit intercept-only baseline on shared sample.
message("  Fitting intercept-only baseline on shared sample...")
baseline_shared_res <- fit_and_loo(
  mod_intercept,
  list(N = N_shared, y = y_shared, prior_sd_sigma = prior_sd_sigma)
)
baseline_shared_elpd <- as.numeric(
  baseline_shared_res$loo$estimates["elpd_loo", "Estimate"]
)

combo_results <- list()

for (combo in all_combos) {
  trial_local <- sort(unique(unlist(block_defs[combo])))
  trial_cols  <- q_cols[trial_local]
  P_trial     <- length(trial_local)

  trial_data <- list(
    N              = N_shared,
    P              = P_trial,
    y              = y_shared,
    X_mean         = Zm_shared[,  trial_cols, drop = FALSE],
    X_sd           = Zsd_shared[, trial_cols, drop = FALSE],
    prior_sd_sigma = prior_sd_sigma
  )

  res <- tryCatch(
    fit_and_loo(mod_normal, trial_data),
    error = function(e) NULL
  )
  if (is.null(res)) next

  elpd      <- as.numeric(res$loo$estimates["elpd_loo", "Estimate"])
  elpd_se   <- as.numeric(res$loo$estimates["elpd_loo", "SE"])
  elpd_diff <- elpd - baseline_shared_elpd

  combo_results[[length(combo_results) + 1L]] <- list(
    combo      = combo,
    loo        = res$loo,
    elpd       = elpd,
    elpd_se    = elpd_se,
    elpd_diff  = elpd_diff,
    trial_cols = trial_cols
  )
}

# ---------------------------------------------------------------------------
# Stacking weights and per-block inclusion weights.
# ---------------------------------------------------------------------------
if (length(combo_results) == 0) {
  message("  No combinations fitted successfully for '", b, "'.")
  group_selection[[b]] <- list(
    behaviour         = b,
    block_inclusion   = NULL,
    stacking_weights  = NULL,
    baseline_elpd     = baseline_shared_elpd,
    n_cc              = N_shared
  )
} else {
  message(sprintf("  Computing stacking weights over %d fitted combinations...",
                  length(combo_results)))

  # Build named list of loo objects: intercept-only first, then combos.
  combo_labels <- vapply(combo_results, function(cr) {
    paste(block_names[cr$combo], collapse = "+")
  }, character(1))

  loo_list <- c(
    list(intercept_only = baseline_shared_res$loo),
    setNames(lapply(combo_results, `[[`, "loo"), combo_labels)
  )

  lpd_point   <- do.call(cbind, lapply(loo_list, function(l) l$pointwise[, "elpd_loo"]))
  stacking_w  <- as.numeric(loo::stacking_weights(lpd_point))
  pseudobma_w <- as.numeric(loo::pseudobma_weights(lpd_point, BB = TRUE))
  names(stacking_w)  <- names(loo_list)
  names(pseudobma_w) <- names(loo_list)

  # Per-block inclusion weight: sum stacking weights over models containing
  # each block (intercept-only contains no blocks so never contributes).
  # block_defs indices 1..n_blocks map to block_names.
  block_inclusion <- data.frame(
    block_idx        = seq_len(n_blocks),
    block_label      = block_names,
    inclusion_weight = vapply(seq_len(n_blocks), function(blk) {
      sum(vapply(seq_along(combo_results), function(ci) {
        if (blk %in% combo_results[[ci]]$combo) stacking_w[ci + 1L] else 0
      }, numeric(1)))
    }, numeric(1)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(inclusion_weight))

  # Effect direction from correlation — valid for zero-centred symmetric prior.
  block_inclusion$effect_direction <- vapply(block_inclusion$block_idx, function(blk) {
    blk_local <- sort(unique(unlist(block_defs[blk])))
    blk_cols  <- q_cols[blk_local]
    # Use mean correlation across all terms in the block.
    cors <- vapply(blk_cols, function(col) {
      stats::cor(Zm_shared[, col], y_shared)
    }, numeric(1))
    r <- mean(cors, na.rm = TRUE)
    if (is.na(r))       NA_character_
    else if (r > 0)     "positive"
    else                "negative"
  }, character(1))

  # Summary of all combinations for diagnostics.
  combo_summary <- data.frame(
    combo_label      = combo_labels,
    n_blocks         = vapply(combo_results, function(cr) length(cr$combo), integer(1)),
    elpd             = vapply(combo_results, `[[`, numeric(1), "elpd"),
    elpd_se          = vapply(combo_results, `[[`, numeric(1), "elpd_se"),
    elpd_diff        = vapply(combo_results, `[[`, numeric(1), "elpd_diff"),
    stacking_weight  = stacking_w[-1],   # drop intercept-only
    pseudobma_weight = pseudobma_w[-1],
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(stacking_weight))

  message(sprintf(
    "  Intercept-only stacking weight: %.4f", stacking_w["intercept_only"]
  ))
  message(sprintf(
    "  Top block by inclusion weight: %s  (w=%.4f  dir=%s)",
    block_inclusion$block_label[1],
    block_inclusion$inclusion_weight[1],
    block_inclusion$effect_direction[1]
  ))

  group_selection[[b]] <- list(
    behaviour         = b,
    block_inclusion   = block_inclusion,
    combo_summary     = combo_summary,
    stacking_weights  = stacking_w,
    pseudobma_weights = pseudobma_w,
    baseline_elpd     = baseline_shared_elpd,
    n_cc              = N_shared
  )
}

result_file <- file.path(
  input_folder,
  paste0("dsem_selection_", target_behaviour, "_single_behaviour.rds")
)
saveRDS(group_selection[[target_behaviour]], file = result_file)
message("\nResults saved to ", result_file)

