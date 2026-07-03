# =============================================================================
# granger_analysis.R
#
# Usage:
#   Rscript granger_analysis.R <step1_rds_or_folder> [output_prefix]
#   Rscript granger_analysis.R <step1_rds_or_folder> [output_prefix] array_task_id=1 array_task_count=9
#
# Examples:
#   Rscript granger_analysis.R dsem_step1_1.rds
#   Rscript granger_analysis.R MODELS-STEP1 granger_loocv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript granger_analysis.R <step1_rds_or_folder> [output_prefix]")
}

input_path    <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "granger_loocv"

# Parse optional keyed arguments for array job sharding only.
opt_args <- if (length(args) >= 3) args[3:length(args)] else character(0)
keyed    <- opt_args[grepl("=", opt_args)]

parse_keyed_int <- function(keyed_args, key) {
  pat <- paste0("^(--)?", key, "=")
  hit <- keyed_args[grepl(pat, keyed_args)]
  if (length(hit) == 0) return(NA_integer_)
  as.integer(sub(pat, "", hit[length(hit)]))
}

array_task_id_keyed    <- parse_keyed_int(keyed, "array_task_id")
array_task_count_keyed <- parse_keyed_int(keyed, "array_task_count")

array_task_id_env    <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",    unset = "")))
array_task_count_env <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "")))

array_task_id <- if (!is.na(array_task_id_keyed)) {
  array_task_id_keyed
} else if (!is.na(array_task_id_env)) {
  array_task_id_env
} else {
  1L
}

array_task_count <- if (!is.na(array_task_count_keyed)) {
  array_task_count_keyed
} else if (!is.na(array_task_count_env)) {
  array_task_count_env
} else {
  1L
}

if (!is.finite(array_task_id)    || array_task_id    < 1L) stop("array_task_id must be a positive integer")
if (!is.finite(array_task_count) || array_task_count < 1L) stop("array_task_count must be a positive integer")
if (array_task_id > array_task_count)                      stop("array_task_id cannot exceed array_task_count")

message("Array shard settings: task_id=", array_task_id, ", task_count=", array_task_count)

suppressPackageStartupMessages({
  library(posterior)
  library(loo)
})


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

is_dir <- function(path) {
  info <- suppressWarnings(file.info(path)$isdir)
  if (is.na(info)) FALSE else info
}

combine_step1_results <- function(rds_list) {
  if (length(rds_list) == 1) return(rds_list[[1]])
  draws_list     <- lapply(rds_list, function(x) posterior::as_draws_matrix(x$draws))
  prep           <- rds_list[[1]]$prep
  draws_combined <- posterior::bind_draws(draws_list, along = "draw")
  list(draws = draws_combined, prep = prep)
}

read_step1_input <- function(path) {
  if (is_dir(path)) {
    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
    if (length(files) == 0) stop("No .rds files found in folder: ", path)
    message("Reading and combining ", length(files), " RDS files from folder: ", path)
    step1   <- combine_step1_results(lapply(files, readRDS))
    out_dir <- path
  } else {
    if (!file.exists(path)) stop("Input file not found: ", path)
    message("Reading single RDS file: ", path)
    step1   <- readRDS(path)
    out_dir <- dirname(path)
  }
  if (is.null(step1$draws) || is.null(step1$prep))
    stop("Input object must contain $draws and $prep")
  list(step1 = step1, out_dir = out_dir)
}

softplus <- function(x) {
  log1p(exp(-abs(x))) + pmax(x, 0)
}

extract_z_person_matrix <- function(draws, N, Q) {
  var_names <- posterior::variables(draws)
  z_vars    <- var_names[grepl("^z_person\\[", var_names)]
  if (length(z_vars) != N * Q)
    stop("Expected ", N * Q, " z_person variables, found ", length(z_vars))

  idx_txt <- regmatches(z_vars, regexpr("\\[(\\d+),(\\d+)\\]", z_vars))
  i_idx   <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx_txt))
  q_idx   <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx_txt))

  z_mat      <- posterior::as_draws_matrix(posterior::subset_draws(draws, variable = z_vars))
  col_lookup <- matrix(NA_integer_, nrow = N, ncol = Q)
  for (j in seq_along(z_vars)) col_lookup[i_idx[j], q_idx[j]] <- j

  list(z_mat = z_mat, col_lookup = col_lookup)
}

build_obs_index <- function(stan_data) {
  idx <- which(stan_data$obs_mask == 1L, arr.ind = TRUE)
  colnames(idx) <- c("i", "t")
  idx
}

build_y_block_data <- function(stan_data, obs_idx, y_idx, K) {
  n_obs <- nrow(obs_idx)
  i_idx <- obs_idx[, "i"]
  t_idx <- obs_idx[, "t"]
  y_vec <- as.numeric(stan_data$Y[cbind(i_idx, t_idx, rep(y_idx, n_obs))])
  X_mat <- matrix(0, nrow = n_obs, ncol = K)
  for (k in seq_len(K))
    X_mat[, k] <- as.numeric(stan_data$X_lag[cbind(i_idx, t_idx, rep(k, n_obs))])
  list(person_idx = i_idx, y = y_vec, X = X_mat)
}

loglik_binary <- function(eta, y_vec) {
  y_mat <- matrix(rep(y_vec, each = nrow(eta)), nrow = nrow(eta))
  ll1   <- -softplus(-eta)
  ll0   <- -softplus(eta)
  y_mat * ll1 + (1 - y_mat) * ll0
}

loglik_gaussian <- function(eta, y_vec, sigma_vec) {
  y_mat  <- matrix(rep(y_vec, each = nrow(eta)), nrow = nrow(eta))
  sigma2 <- sigma_vec^2
  ll     <- sweep(-(y_mat - eta)^2, 1, 2 * sigma2, "/")
  sweep(ll, 1, -0.5 * log(2 * pi * sigma2), "+")
}

# Compute PSIS-LOO directly on the full matrix — no chunking.
run_loo <- function(ll) {
  loo::loo(ll, cores = 1L)
}

# Precompute full-model log-likelihood and LOO for outcome y.
prepare_y_cache <- function(y_dat, y_idx, z_mat, col_lookup, V, K, is_binary_y) {
  S     <- nrow(z_mat)
  n_obs <- nrow(y_dat$X)

  ll_full      <- matrix(NA_real_, nrow = S, ncol = n_obs)
  person_cache <- vector("list", length(unique(y_dat$person_idx)))
  pc_id        <- 1L

  A_q_idx <- as.integer(((seq_len(K) - 1L) * V) + y_idx)
  mu_q    <- V * K + y_idx
  lv_q    <- V * K + V + y_idx

  for (i in unique(y_dat$person_idx)) {
    rows <- which(y_dat$person_idx == i)
    Xi   <- y_dat$X[rows, , drop = FALSE]
    yi   <- y_dat$y[rows]

    Ai  <- z_mat[, col_lookup[i, A_q_idx], drop = FALSE]
    mui <- z_mat[, col_lookup[i, mu_q]]
    lvi <- z_mat[, col_lookup[i, lv_q]]

    eta_full <- Ai %*% t(Xi)
    eta_full <- sweep(eta_full, 1, mui, "+")

    if (is_binary_y) {
      ll_full_i <- loglik_binary(eta_full, yi)
      sigma_i   <- NULL
    } else {
      sigma_i   <- sqrt(exp(lvi))
      ll_full_i <- loglik_gaussian(eta_full, yi, sigma_i)
    }

    ll_full[, rows] <- ll_full_i

    person_cache[[pc_id]] <- list(
      rows     = rows,
      Xi       = Xi,
      yi       = yi,
      Ai       = Ai,
      eta_full = eta_full,
      sigma_i  = sigma_i
    )
    pc_id <- pc_id + 1L
  }

  loo_full <- run_loo(ll_full)

  list(
    person_cache   = person_cache,
    ll_full        = ll_full,
    loo_full       = loo_full,
    pointwise_full = loo_full$pointwise[, "elpd_loo"],
    n_obs          = n_obs,
    is_binary_y    = is_binary_y
  )
}

# Score a single x -> y pair against the cached full model.
score_pair_from_cache <- function(y_cache, x_idx, K_per_pair) {
  S      <- nrow(y_cache$ll_full)
  n_obs  <- y_cache$n_obs
  ll_red <- matrix(NA_real_, nrow = S, ncol = n_obs)
  block_k    <- as.integer(seq_len(K_per_pair) + (x_idx - 1L) * K_per_pair)
  coef_pool  <- numeric(0)

  for (pc in y_cache$person_cache) {
    contrib_block <- pc$Ai[, block_k, drop = FALSE] %*% t(pc$Xi[, block_k, drop = FALSE])
    eta_red       <- pc$eta_full - contrib_block

    if (y_cache$is_binary_y) {
      ll_red_i <- loglik_binary(eta_red, pc$yi)
    } else {
      ll_red_i <- loglik_gaussian(eta_red, pc$yi, pc$sigma_i)
    }
    ll_red[, pc$rows] <- ll_red_i
    coef_pool <- c(coef_pool, as.vector(pc$Ai[, block_k, drop = FALSE]))
  }

  mean_coef_x_to_y     <- mean(coef_pool)
  prop_positive_x_to_y <- mean(coef_pool > 0)
  coef_direction_x_to_y <- if (prop_positive_x_to_y >= 0.95) {
    "positive"
  } else if (prop_positive_x_to_y <= 0.05) {
    "negative"
  } else {
    "mixed"
  }

  loo_red   <- run_loo(ll_red)
  elpd_full <- y_cache$loo_full$estimates["elpd_loo", "Estimate"]
  elpd_red  <- loo_red$estimates["elpd_loo", "Estimate"]
  elpd_diff <- elpd_full - elpd_red
  pt_diff   <- y_cache$pointwise_full - loo_red$pointwise[, "elpd_loo"]
  se_diff   <- sqrt(length(pt_diff) * stats::var(pt_diff))

  list(
    elpd_full             = elpd_full,
    elpd_reduced          = elpd_red,
    elpd_diff             = elpd_diff,
    se_diff               = se_diff,
    n_obs                 = n_obs,
    mean_coef_x_to_y      = mean_coef_x_to_y,
    prop_positive_x_to_y  = prop_positive_x_to_y,
    coef_direction_x_to_y = coef_direction_x_to_y
  )
}


# =============================================================================
# LOAD DATA
# =============================================================================

loaded    <- read_step1_input(input_path)
step1     <- loaded$step1
out_dir   <- loaded$out_dir
prep      <- step1$prep
stan_data <- prep$stan_data
behaviours <- prep$behaviours
V          <- if (!is.null(prep$V)) prep$V else stan_data$V
K          <- if (!is.null(prep$K)) prep$K else stan_data$K
Q          <- prep$Q
N          <- stan_data$N
K_per_pair <- as.integer(K / V)
if (K %% V != 0L) stop("K must be divisible by V")
is_binary <- as.logical(stan_data$is_binary)


# =============================================================================
# SHARD ASSIGNMENT
# =============================================================================

all_y_idx      <- seq_len(V)
y_idx_assigned <- all_y_idx[((all_y_idx - 1L) %% array_task_count) + 1L == array_task_id]
message("Assigned outcomes: ", paste(behaviours[y_idx_assigned], collapse = ", "))

output_prefix_effective <- if (array_task_count > 1L) {
  paste0(output_prefix,
         "_task", sprintf("%03d", array_task_id),
         "of",   sprintf("%03d", array_task_count))
} else {
  output_prefix
}

run_start_time <- proc.time()[["elapsed"]]


# =============================================================================
# EXTRACT POSTERIOR DRAWS
# =============================================================================

message("Extracting z_person draws...")
z_info     <- extract_z_person_matrix(step1$draws, N = N, Q = Q)
z_mat      <- z_info$z_mat
col_lookup <- z_info$col_lookup


# =============================================================================
# BUILD OBSERVATION INDEX AND PER-Y DATA BLOCKS
# =============================================================================

message("Building observation index...")
obs_idx <- build_obs_index(stan_data)

message("Precomputing y-specific design blocks...")
y_blocks <- vector("list", V)
for (y_idx in seq_len(V))
  y_blocks[[y_idx]] <- build_y_block_data(stan_data, obs_idx, y_idx, K)


# =============================================================================
# MAIN SCORING LOOP
# =============================================================================

results <- list()
message("Scoring directed pairs x -> y with ELPD...")

for (y_idx in y_idx_assigned) {
  outcome_start_time <- proc.time()[["elapsed"]]
  y_name      <- behaviours[y_idx]
  y_dat       <- y_blocks[[y_idx]]
  is_binary_y <- is_binary[y_idx]

  message("Outcome ", y_name, ": computing full-model LOO...")
  y_cache <- prepare_y_cache(
    y_dat       = y_dat,
    y_idx       = y_idx,
    z_mat       = z_mat,
    col_lookup  = col_lookup,
    V           = V,
    K           = K,
    is_binary_y = is_binary_y
  )

  x_candidates <- setdiff(seq_len(V), y_idx)

  for (x_idx in x_candidates) {
    x_name <- behaviours[x_idx]
    message("  Pair ", x_name, " -> ", y_name)

    pair_score <- score_pair_from_cache(y_cache, x_idx, K_per_pair)

    results[[length(results) + 1L]] <- data.frame(
      predictor_x                  = x_name,
      outcome_y                    = y_name,
      outcome_is_binary            = is_binary_y,
      n_obs                        = pair_score$n_obs,
      elpd_full                    = pair_score$elpd_full,
      elpd_reduced_drop_xlags      = pair_score$elpd_reduced,
      elpd_diff_full_minus_reduced = pair_score$elpd_diff,
      se_diff                      = pair_score$se_diff,
      z_score = if (is.finite(pair_score$se_diff) && pair_score$se_diff > 0)
        pair_score$elpd_diff / pair_score$se_diff else NA_real_,
      supports_x_granger_causing_y = pair_score$elpd_diff > 0,
      strong_support_diff_gt_2se   = if (is.finite(pair_score$se_diff) && pair_score$se_diff > 0)
        pair_score$elpd_diff > (2 * pair_score$se_diff) else FALSE,
      mean_coef_x_to_y      = pair_score$mean_coef_x_to_y,
      prop_positive_x_to_y  = pair_score$prop_positive_x_to_y,
      coef_direction_x_to_y = pair_score$coef_direction_x_to_y,
      stringsAsFactors = FALSE
    )
  }

  message(
    "Outcome ", y_name, " done in ",
    format(round(proc.time()[["elapsed"]] - outcome_start_time, 1), nsmall = 1), " sec"
  )
}


# =============================================================================
# SAVE RESULTS
# =============================================================================

if (length(results) > 0) {
  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$elpd_diff_full_minus_reduced, decreasing = TRUE), , drop = FALSE]
} else {
  results_df <- data.frame(
    predictor_x                  = character(0),
    outcome_y                    = character(0),
    outcome_is_binary            = logical(0),
    n_obs                        = integer(0),
    elpd_full                    = numeric(0),
    elpd_reduced_drop_xlags      = numeric(0),
    elpd_diff_full_minus_reduced = numeric(0),
    se_diff                      = numeric(0),
    z_score                      = numeric(0),
    supports_x_granger_causing_y = logical(0),
    strong_support_diff_gt_2se   = logical(0),
    mean_coef_x_to_y             = numeric(0),
    prop_positive_x_to_y         = numeric(0),
    coef_direction_x_to_y        = character(0),
    stringsAsFactors = FALSE
  )
}

csv_file <- file.path(out_dir, paste0(output_prefix_effective, "_pairwise_elpd.csv"))
rds_file <- file.path(out_dir, paste0(output_prefix_effective, "_pairwise_elpd.rds"))
write.csv(results_df, csv_file, row.names = FALSE)
saveRDS(results_df, rds_file)

message("Total runtime (sec): ", format(round(proc.time()[["elapsed"]] - run_start_time, 1), nsmall = 1))
message("Saved: ", csv_file)
message("Saved: ", rds_file)

cat("\nTop pairs by ELPD advantage:\n")
print(utils::head(results_df, 15))
