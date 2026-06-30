# =============================================================================
# final_term_weighting.R
#
# Bayesian model averaging via two-level stacking over behaviour-set
# combinations.
#
# Step 2a (select_terms_for_individual_behaviours.R) computes within-behaviour
# stacking weights across all block combinations for each behaviour, producing
# per-block inner inclusion weights.
#
# This script (Step 2b) forms combinations of the 9 behaviour blocks, fits
# every combination on the shared complete-case CDI sample, and computes
# outer stacking weights via loo::stacking_weights().
#
# The final per-term inclusion weight is the product of the outer inclusion
# weight for the behaviour and the inner inclusion weight for the block
# within that behaviour:
#
#   w_final(q) = outer_inclusion(behaviour(q)) × inner_inclusion(q | behaviour)
#
# This gives fine-grained inclusion probabilities at the level of individual
# x→y predictor pairs, means, and log-variances, without requiring an
# exhaustive search over the full 2^99 predictor space.
#
# Usage:
#   Rscript dsem_step2_selection_from_single_behaviour_sets.R [selection_folder] [step1_folder] [output_prefix]
# =============================================================================#

args <- commandArgs(trailingOnly = TRUE)
selection_folder <- if (length(args) >= 1) args[1] else "DSEM-MODELS-STEP1-PERCULATE"
step1_folder <- if (length(args) >= 2) args[2] else selection_folder
output_prefix <- if (length(args) >= 3) args[3] else "dsem_selection_from_single_behaviour_sets"

is_dir <- function(path) {
  info <- suppressWarnings(file.info(path)$isdir)
  if (is.na(info)) FALSE else info
}
if (!is_dir(selection_folder)) stop("Selection folder not found: ", selection_folder)
if (!is_dir(step1_folder)) stop("Step-1 folder not found: ", step1_folder)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(dplyr)
})

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

combine_step1_results <- function(rds_list) {
  if (length(rds_list) == 1) return(rds_list[[1]])
  draws_list <- lapply(rds_list, function(x) posterior::as_draws_matrix(x$draws))
  prep <- rds_list[[1]]$prep
  draws_combined <- posterior::bind_draws(draws_list, along = "draw")
  list(draws = draws_combined, prep = prep)
}

extract_single_behaviour_record <- function(obj) {
  records <- list()

  collect_record <- function(gs) {
    if (!is.list(gs)) return(NULL)
    if (is.null(gs$behaviour)) return(NULL)
    # New format: block_inclusion data frame with inclusion_weight per block.
    if (!is.null(gs$block_inclusion) && is.data.frame(gs$block_inclusion)) {
      return(list(
        behaviour       = as.character(gs$behaviour),
        block_inclusion = gs$block_inclusion,  # data frame: block_idx, block_label, inclusion_weight, effect_direction
        n_cc            = as.integer(gs$n_cc)
      ))
    }
    NULL
  }

  one <- collect_record(obj)
  if (!is.null(one)) records[[length(records) + 1L]] <- one

  if (is.list(obj) && length(obj) > 0) {
    for (el in obj) {
      rec <- collect_record(el)
      if (!is.null(rec)) records[[length(records) + 1L]] <- rec
    }
  }

  if (length(records) == 0) return(NULL)

  rec_map <- list()
  for (r in records) rec_map[[r$behaviour]] <- r
  unname(rec_map)
}

selection_files <- list.files(
  selection_folder,
  pattern = "dsem_selection_.*single_behaviour\\.rds$",
  full.names = TRUE
)
if (length(selection_files) == 0) {
  stop("No single-behaviour selection files found in: ", selection_folder)
}
message("Loading ", length(selection_files), " single-behaviour selection RDS files")

single_records <- list()
for (f in selection_files) {
  obj <- readRDS(f)
  recs <- extract_single_behaviour_record(obj)
  if (is.null(recs)) next
  for (r in recs) single_records[[r$behaviour]] <- r
}
if (length(single_records) == 0) {
  stop("No valid {behaviour, selected_blocks} records found in selection RDS files")
}

step1_files <- list.files(step1_folder, pattern = "^dsem_step1_.*\\.rds$", full.names = TRUE)
if (length(step1_files) == 0) {
  step1_files <- list.files(step1_folder, pattern = "^dsem_result_.*\\.rds$", full.names = TRUE)
}
if (length(step1_files) == 0) {
  stop("No step-1 files found (expected dsem_step1_*.rds or dsem_result_*.rds) in: ", step1_folder)
}
message("Loading and combining ", length(step1_files), " step-1 files")
step1_result <- combine_step1_results(lapply(step1_files, readRDS))

step1_prep <- step1_result$prep
behaviours <- step1_prep$behaviours
V <- step1_prep$V
K <- step1_prep$K
Q <- step1_prep$Q
N <- step1_prep$stan_data$N

stopifnot("K must be divisible by V" = K %% V == 0L)
K_per_pair <- as.integer(K / V)

lag_blocks <- lapply(seq_len(V), function(pred_v) {
  seq_len(K_per_pair) + (pred_v - 1L) * K_per_pair
})
names(lag_blocks) <- behaviours

group_indices <- lapply(seq_len(V), function(v) {
  A_cols <- as.integer(unlist(lapply(seq_len(V), function(pred_v) {
    sapply(seq_len(K_per_pair), function(l) {
      k <- (pred_v - 1L) * K_per_pair + l
      (k - 1L) * V + v
    })
  })))
  list(
    A = A_cols,
    mu = V * K + v,
    lv = V * K + V + v,
    all = c(A_cols, V * K + v, V * K + V + v)
  )
})
names(group_indices) <- behaviours

label_q_term <- function(q_col, behaviours, V, K, K_per_pair) {
  if (q_col <= V * K) {
    v <- ((q_col - 1L) %% V) + 1L
    k <- ((q_col - 1L) %/% V) + 1L
    pred_v <- ((k - 1L) %/% K_per_pair) + 1L
    lag <- ((k - 1L) %% K_per_pair) + 1L
    return(sprintf("A[y=%s,x=%s,lag=%d]", behaviours[v], behaviours[pred_v], lag))
  }

  if (q_col <= V * K + V) {
    v <- q_col - V * K
    return(sprintf("mu[%s]", behaviours[v]))
  }

  v <- q_col - (V * K + V)
  sprintf("lv[%s]", behaviours[v])
}

decode_q_term <- function(q_col, behaviours, V, K, K_per_pair) {
  if (q_col <= V * K) {
    y_idx <- ((q_col - 1L) %% V) + 1L
    k <- ((q_col - 1L) %/% V) + 1L
    x_idx <- ((k - 1L) %/% K_per_pair) + 1L
    lag <- ((k - 1L) %% K_per_pair) + 1L
    return(list(
      family = "A",
      y = behaviours[y_idx],
      x = behaviours[x_idx],
      lag = lag
    ))
  }

  if (q_col <= V * K + V) {
    y_idx <- q_col - V * K
    return(list(family = "mu", y = behaviours[y_idx], x = NA_character_, lag = NA_integer_))
  }

  y_idx <- q_col - (V * K + V)
  list(family = "lv", y = behaviours[y_idx], x = NA_character_, lag = NA_integer_)
}

bridge_step1 <- function(step1_obj) {
  draws <- step1_obj$draws
  prep <- step1_obj$prep
  N_local <- prep$stan_data$N
  Q_local <- prep$Q

  all_vars <- posterior::variables(draws)
  z_vars <- all_vars[grepl("^z_person\\[", all_vars)]
  stopifnot(length(z_vars) == N_local * Q_local)

  idx_mat <- regmatches(z_vars, regexpr("\\[(\\d+),(\\d+)\\]", z_vars))
  i_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx_mat))
  q_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx_mat))

  z_mat <- posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = z_vars)
  )

  z_hat <- matrix(0, nrow = N_local, ncol = Q_local)
  z_sd <- matrix(0, nrow = N_local, ncol = Q_local)
  for (col_j in seq_len(ncol(z_mat))) {
    ii <- i_idx[col_j]
    qq <- q_idx[col_j]
    z_hat[ii, qq] <- mean(z_mat[, col_j])
    z_sd[ii, qq] <- sd(z_mat[, col_j])
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
pid_df <- data.frame(participant_id = as.character(participants), stringsAsFactors = FALSE)
merged <- dplyr::left_join(pid_df, cdi_imp[, c("participant_id", "us18_imp")], by = "participant_id")

has_outcome <- !is.na(merged$us18_imp)
out_idx <- which(has_outcome)
N_out <- length(out_idx)
y_out <- log1p(as.numeric(merged$us18_imp[has_outcome]))

Z_mean <- bridge$z_hat[out_idx, , drop = FALSE]
Z_sd <- bridge$z_sd[out_idx, , drop = FALSE]

message("Outcome sample: N_out=", N_out)

# Build per-behaviour q-column sets and inner block inclusion weights.
# Every block with any nonzero inner inclusion weight enters the outer search.
# selected_q_by_behaviour: all q_cols for blocks with inclusion_weight > 0.
# inner_block_inclusion: keyed by behaviour, then by q_col, giving the
#   within-behaviour stacking inclusion weight for that block.
selected_q_by_behaviour  <- list()
inner_block_inclusion    <- list()  # behaviour -> named numeric vector (q_col -> inner weight)
inner_block_direction    <- list()  # behaviour -> named character vector (q_col -> direction)

for (b in names(single_records)) {
  if (!(b %in% behaviours)) {
    warning("Skipping behaviour not found in step-1 behaviours: ", b)
    next
  }

  rec        <- single_records[[b]]
  bi         <- rec$block_inclusion   # data frame: block_idx, block_label, inclusion_weight, effect_direction
  q_cols_b   <- group_indices[[b]]$all
  block_defs_b <- c(
    lapply(seq_len(V), function(i) lag_blocks[[i]]),
    list(K + 1L),
    list(K + 2L)
  )

  # Only consider blocks that exist and have nonzero inclusion weight.
  bi <- bi[is.finite(bi$block_idx) &
           bi$block_idx >= 1L &
           bi$block_idx <= length(block_defs_b) &
           bi$inclusion_weight > 0, , drop = FALSE]

  if (nrow(bi) == 0) {
    selected_q_by_behaviour[[b]] <- integer(0)
    next
  }

  # Union of all q_cols across nonzero-weight blocks.
  local_idx <- sort(unique(unlist(block_defs_b[bi$block_idx])))
  b_q_cols  <- as.integer(q_cols_b[local_idx])
  selected_q_by_behaviour[[b]] <- b_q_cols

  # Map each q_col to its block's inner inclusion weight and direction.
  # q_cols within the same block share the block's weight.
  w_vec   <- setNames(numeric(length(b_q_cols)),   as.character(b_q_cols))
  dir_vec <- setNames(character(length(b_q_cols)), as.character(b_q_cols))

  for (row_i in seq_len(nrow(bi))) {
    blk_local <- sort(unique(unlist(block_defs_b[bi$block_idx[row_i]])))
    blk_q     <- as.integer(q_cols_b[blk_local])
    blk_q_chr <- as.character(blk_q)
    blk_q_chr <- blk_q_chr[blk_q_chr %in% names(w_vec)]
    w_vec[blk_q_chr]   <- bi$inclusion_weight[row_i]
    dir_vec[blk_q_chr] <- as.character(bi$effect_direction[row_i])
  }

  inner_block_inclusion[[b]] <- w_vec
  inner_block_direction[[b]] <- dir_vec
}

selected_q_by_behaviour <- selected_q_by_behaviour[lengths(selected_q_by_behaviour) > 0]
if (length(selected_q_by_behaviour) == 0) {
  stop("No behaviours had blocks with nonzero inner inclusion weight")
}

candidate_behaviours <- sort(names(selected_q_by_behaviour))
message("Candidate behaviours with nonzero inner inclusion weight: ",
        paste(candidate_behaviours, collapse = ", "))

message("Compiling Stan models...")
mod_intercept <- cmdstan_model(write_stan_file(stan_intercept))
mod_normal <- cmdstan_model(write_stan_file(stan_normal))

fit_and_loo <- function(mod, data_list, chains = 4, iter_warmup = 500, iter_sampling = 500, refresh = 0) {
  fit <- mod$sample(
    data = data_list,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    show_messages = FALSE
  )
  ll <- fit$draws("log_lik", format = "matrix")
  loo_obj <- loo::loo(ll, cores = 1)
  list(fit = fit, loo = loo_obj)
}

prior_sd_sigma <- 2.5
min_n_frac <- 0.80

combos <- unlist(lapply(seq_along(candidate_behaviours), function(k) {
  combn(candidate_behaviours, k, simplify = FALSE)
}), recursive = FALSE)

# ---------------------------------------------------------------------------
# Determine the shared complete-case mask across ALL combinations.
# Stacking requires every model to be evaluated on exactly the same set of
# observations so that the N x M pointwise log-lik matrix is conformable.
# ---------------------------------------------------------------------------
message("Computing shared complete-case mask across all ", length(combos), " combinations...")
all_cols_used <- sort(unique(unlist(lapply(combos, function(bset) {
  unlist(selected_q_by_behaviour[bset])
}))))
shared_mask <- complete.cases(Z_mean[, all_cols_used, drop = FALSE]) &
               complete.cases(Z_sd[,   all_cols_used, drop = FALSE])
N_shared <- sum(shared_mask)
message(sprintf("Shared complete-case N: %d / %d", N_shared, N_out))

if (N_shared < ceiling(min_n_frac * N_out)) {
  stop(sprintf(
    "Shared complete-case N (%d) is below the minimum fraction %.0f%% of N_out (%d). ",
    N_shared, 100 * min_n_frac, N_out,
    "Consider lowering min_n_frac or reviewing missingness in the predictor matrix."
  ))
}

y_shared   <- y_out[shared_mask]
Zm_shared  <- Z_mean[shared_mask, , drop = FALSE]
Zsd_shared <- Z_sd[shared_mask,  , drop = FALSE]

# ---------------------------------------------------------------------------
# Checkpoint setup: resume from a previous partial run if available.
# ---------------------------------------------------------------------------
checkpoint_file <- file.path(selection_folder, paste0(output_prefix, "_checkpoint.rds"))
checkpoint_every <- 25L

# ---------------------------------------------------------------------------
# Fit the intercept-only baseline on the shared sample (included in stacking).
# ---------------------------------------------------------------------------
message("Fitting intercept-only baseline on shared sample...")
baseline_res <- fit_and_loo(
  mod_intercept,
  list(N = N_shared, y = y_shared, prior_sd_sigma = prior_sd_sigma)
)
baseline_elpd_full <- as.numeric(baseline_res$loo$estimates["elpd_loo", "Estimate"])

# ---------------------------------------------------------------------------
# Fit every combination model and retain their LOO objects.
# ---------------------------------------------------------------------------
message("Fitting ", length(combos), " behaviour-set combination models on shared sample...")
results      <- vector("list", length(combos))
loo_list     <- vector("list", length(combos) + 1L)  # +1 for intercept-only
loo_list[[1]] <- baseline_res$loo  # intercept-only is model 1

resume_from <- 1L
if (file.exists(checkpoint_file)) {
  message("Checkpoint found — resuming from: ", checkpoint_file)
  chk <- readRDS(checkpoint_file)
  if (!is.null(chk$results))  results  <- chk$results
  if (!is.null(chk$loo_list)) loo_list <- chk$loo_list
  resume_from <- chk$last_completed + 1L
  message(sprintf("  Resuming from combo %d / %d", resume_from, length(combos)))
}

progress_every <- checkpoint_every
loop_start <- proc.time()[["elapsed"]]

for (ii in seq(resume_from, length(combos))) {
  if (ii == resume_from || ii %% progress_every == 0L || ii == length(combos)) {
    elapsed <- proc.time()[["elapsed"]] - loop_start
    n_done  <- ii - resume_from + 1L
    avg_per_combo <- elapsed / n_done
    remaining <- avg_per_combo * (length(combos) - ii)
    message(sprintf(
      "Progress %d/%d (%.1f%%) | elapsed %.1f min | avg %.2f sec/combo | ETA %.1f min",
      ii, length(combos), 100 * ii / length(combos),
      elapsed / 60, avg_per_combo, remaining / 60
    ))
    saveRDS(
      list(last_completed = ii, results = results, loo_list = loo_list),
      checkpoint_file
    )
  }

  bset       <- combos[[ii]]
  trial_cols <- sort(unique(unlist(selected_q_by_behaviour[bset])))

  trial_data <- list(
    N            = N_shared,
    P            = length(trial_cols),
    y            = y_shared,
    X_mean       = Zm_shared[,  trial_cols, drop = FALSE],
    X_sd         = Zsd_shared[, trial_cols, drop = FALSE],
    prior_sd_sigma = prior_sd_sigma
  )

  res <- tryCatch(fit_and_loo(mod_normal, trial_data), error = function(e) NULL)

  if (is.null(res)) {
    results[[ii]] <- data.frame(
      combo_id     = ii,
      behaviours   = paste(bset, collapse = "+"),
      n_behaviours = length(bset),
      n_terms      = length(trial_cols),
      n_cc         = N_shared,
      elpd         = NA_real_,
      elpd_se      = NA_real_,
      baseline_elpd = baseline_elpd_full,
      elpd_diff    = NA_real_,
      elpd_diff_se = NA_real_,
      stacking_weight   = NA_real_,
      pseudobma_weight  = NA_real_,
      stringsAsFactors = FALSE
    )
    loo_list[[ii + 1L]] <- NULL
    next
  }

  elpd     <- as.numeric(res$loo$estimates["elpd_loo", "Estimate"])
  elpd_se  <- as.numeric(res$loo$estimates["elpd_loo", "SE"])
  pt_diff  <- res$loo$pointwise[, "elpd_loo"] - baseline_res$loo$pointwise[, "elpd_loo"]
  elpd_diff    <- elpd - baseline_elpd_full
  elpd_diff_se <- sqrt(length(pt_diff) * stats::var(pt_diff))

  results[[ii]] <- data.frame(
    combo_id      = ii,
    behaviours    = paste(bset, collapse = "+"),
    n_behaviours  = length(bset),
    n_terms       = length(trial_cols),
    n_cc          = N_shared,
    elpd          = elpd,
    elpd_se       = elpd_se,
    baseline_elpd = baseline_elpd_full,
    elpd_diff     = elpd_diff,
    elpd_diff_se  = elpd_diff_se,
    stacking_weight  = NA_real_,   # filled in after stacking
    pseudobma_weight = NA_real_,
    stringsAsFactors = FALSE
  )
  loo_list[[ii + 1L]] <- res$loo
}

total_loop_elapsed <- proc.time()[["elapsed"]] - loop_start
message(sprintf("Finished fitting all models in %.1f minutes", total_loop_elapsed / 60))

# ---------------------------------------------------------------------------
# Stacking of predictive distributions (Yao et al., 2018).
# Only models that fitted successfully contribute; failed models get weight 0.
# ---------------------------------------------------------------------------
message("Computing stacking weights (Yao et al., 2018)...")

# Identify models with valid LOO objects (intercept-only is index 1 of loo_list).
valid_idx <- which(!vapply(loo_list, is.null, logical(1)))
if (length(valid_idx) < 2L) {
  stop("Fewer than 2 models fitted successfully; cannot compute stacking weights.")
}
loo_list_valid <- loo_list[valid_idx]

# loo::stacking_weights() accepts a named list of loo objects.
model_labels <- c(
  "intercept_only",
  vapply(combos, function(bset) paste(bset, collapse = "+"), character(1))
)
names(loo_list_valid) <- model_labels[valid_idx]

lpd_point   <- do.call(cbind, lapply(loo_list_valid, function(l) l$pointwise[, "elpd_loo"]))
stacking_w  <- loo::stacking_weights(lpd_point)
pseudobma_w <- loo::pseudobma_weights(lpd_point, BB = TRUE)

# Map weights back to the results data frame rows.
# Row ii in results corresponds to loo_list index ii+1 (index 1 = intercept).
stacking_weights_named  <- as.numeric(stacking_w)
pseudobma_weights_named <- as.numeric(pseudobma_w)
names(stacking_weights_named)  <- names(loo_list_valid)
names(pseudobma_weights_named) <- names(loo_list_valid)

for (ii in seq_along(combos)) {
  list_idx <- ii + 1L
  if (!list_idx %in% valid_idx) next
  lbl <- model_labels[list_idx]
  results[[ii]]$stacking_weight  <- stacking_weights_named[lbl]
  results[[ii]]$pseudobma_weight <- pseudobma_weights_named[lbl]
}

results_df <- dplyr::bind_rows(results) |>
  dplyr::arrange(dplyr::desc(stacking_weight))

intercept_stacking_w  <- stacking_weights_named["intercept_only"]
intercept_pseudobma_w <- pseudobma_weights_named["intercept_only"]
message(sprintf(
  "Intercept-only stacking weight: %.4f  pseudo-BMA+ weight: %.4f",
  intercept_stacking_w, intercept_pseudobma_w
))
message(sprintf(
  "Top combination by stacking weight: %s  w=%.4f  elpd=%.2f",
  results_df$behaviours[1],
  results_df$stacking_weight[1],
  results_df$elpd[1]
))

# ---------------------------------------------------------------------------
# Two-level term inclusion weights.
#
# The final inclusion weight for a q_col is:
#
#   w_final(q) = outer_inclusion(behaviour) × inner_inclusion(q | behaviour)
#
# where:
#   outer_inclusion(b) = sum of outer stacking weights over combinations
#                        that include behaviour b.
#   inner_inclusion(q | b) = within-behaviour stacking weight for the block
#                            containing q, as computed by Step 2a.
#
# Effect direction comes from Step 2a (sign of mean correlation within the
# block on the shared CDI sample), which is the same outcome and sample
# as the outer search.
# ---------------------------------------------------------------------------
message("Computing two-level term inclusion weights...")

# Step 1: outer inclusion weight per behaviour — sum outer stacking weights
# over all combos containing that behaviour.
nonzero_combo_idx <- which(
  !vapply(results, function(r) is.null(r) || is.na(r$stacking_weight) || r$stacking_weight == 0,
          logical(1))
)

outer_inclusion <- setNames(numeric(length(candidate_behaviours)), candidate_behaviours)
for (ii in nonzero_combo_idx) {
  bset  <- combos[[ii]]
  w_raw <- results[[ii]]$stacking_weight
  for (b in bset) outer_inclusion[b] <- outer_inclusion[b] + w_raw
}

# Step 2: multiply outer × inner for each q_col.
# Collect all q_cols across all candidate behaviours.
all_q_chr <- sort(unique(unlist(lapply(candidate_behaviours, function(b) {
  as.character(selected_q_by_behaviour[[b]])
}))))

coef_rows <- lapply(all_q_chr, function(q_chr) {
  q <- as.integer(q_chr)

  # Find which behaviour this q_col belongs to (may appear in multiple if
  # behaviours share q_cols — unlikely but handled by taking the max).
  w_final <- 0
  direction <- NA_character_

  for (b in candidate_behaviours) {
    if (!q_chr %in% names(inner_block_inclusion[[b]])) next
    w_outer <- outer_inclusion[b]
    w_inner <- inner_block_inclusion[[b]][q_chr]
    candidate_w <- w_outer * w_inner
    if (candidate_w > w_final) {
      w_final   <- candidate_w
      direction <- inner_block_direction[[b]][q_chr]
    }
  }

  data.frame(
    q_col            = q,
    term_label       = label_q_term(q, behaviours, V, K, K_per_pair),
    outer_inclusion  = max(vapply(candidate_behaviours, function(b) {
                         if (q_chr %in% names(inner_block_inclusion[[b]]))
                           outer_inclusion[b] else 0
                       }, numeric(1))),
    inner_inclusion  = max(vapply(candidate_behaviours, function(b) {
                         if (q_chr %in% names(inner_block_inclusion[[b]]))
                           inner_block_inclusion[[b]][q_chr] else 0
                       }, numeric(1))),
    inclusion_weight = w_final,
    effect_direction = direction,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
})

coef_df <- dplyr::bind_rows(coef_rows) |>
  dplyr::arrange(dplyr::desc(inclusion_weight))

message(sprintf(
  "Two-level inclusion weights computed for %d unique terms across %d behaviours.",
  nrow(coef_df), length(candidate_behaviours)
))
message(sprintf(
  "Top term: %s  outer=%.4f  inner=%.4f  combined=%.4f  dir=%s",
  coef_df$term_label[1],
  coef_df$outer_inclusion[1],
  coef_df$inner_inclusion[1],
  coef_df$inclusion_weight[1],
  coef_df$effect_direction[1]
))

# Output file paths
coef_csv     <- file.path(selection_folder, paste0(output_prefix, "_stacking_coefs.csv"))
coef_rds     <- file.path(selection_folder, paste0(output_prefix, "_stacking_coefs.rds"))
stacking_csv <- file.path(selection_folder, paste0(output_prefix, "_stacking_weights.csv"))
stacking_rds <- file.path(selection_folder, paste0(output_prefix, "_stacking_weights.rds"))
out_rds      <- file.path(selection_folder, paste0(output_prefix, ".rds"))
out_csv      <- file.path(selection_folder, paste0(output_prefix, ".csv"))

# Stacking summary table: one row per model (including intercept-only).
intercept_row <- data.frame(
  combo_id         = 0L,
  behaviours       = "intercept_only",
  n_behaviours     = 0L,
  n_terms          = 0L,
  n_cc             = N_shared,
  elpd             = baseline_elpd_full,
  elpd_se          = as.numeric(baseline_res$loo$estimates["elpd_loo", "SE"]),
  baseline_elpd    = baseline_elpd_full,
  elpd_diff        = 0,
  elpd_diff_se     = 0,
  stacking_weight  = intercept_stacking_w,
  pseudobma_weight = intercept_pseudobma_w,
  stringsAsFactors = FALSE
)

stacking_df <- dplyr::bind_rows(intercept_row, results_df) |>
  dplyr::arrange(dplyr::desc(stacking_weight))

write.csv(stacking_df, stacking_csv, row.names = FALSE)
saveRDS(stacking_df, stacking_rds)

write.csv(coef_df, coef_csv, row.names = FALSE)
saveRDS(coef_df, coef_rds)

saveRDS(
  list(
    summary                 = results_df,
    stacking_summary        = stacking_df,
    stacking_weights        = stacking_weights_named,
    pseudobma_weights       = pseudobma_weights_named,
    stacking_coefs          = coef_df,
    outer_inclusion         = outer_inclusion,
    inner_block_inclusion   = inner_block_inclusion,
    inner_block_direction   = inner_block_direction,
    intercept_only_elpd     = baseline_elpd_full,
    n_shared                = N_shared,
    candidate_behaviours    = candidate_behaviours,
    selected_q_by_behaviour = selected_q_by_behaviour,
    prior_sd_sigma          = prior_sd_sigma,
    min_n_frac              = min_n_frac
  ),
  out_rds
)
write.csv(results_df, out_csv, row.names = FALSE)

message("Saved: ", stacking_csv)
message("Saved: ", stacking_rds)
message("Saved: ", coef_csv)
message("Saved: ", coef_rds)
message("Saved: ", out_rds)
message("Saved: ", out_csv)

if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
  message("Checkpoint deleted: ", checkpoint_file)
}

