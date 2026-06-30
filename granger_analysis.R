# =============================================================================
# granger_analysis.R
#
# Granger-style pairwise relationship scoring from fitted DSEM step-1 VAR output,
# using PSIS-LOO ELPD as the comparison metric.
#
# Usage:
#   Rscript dsem_var_granger_loocv_array.R <step1_rds_or_folder> [output_prefix] [loo_cores] [pair_cores] [chunk_size]
#   Rscript dsem_var_granger_loocv_array.R <step1_rds_or_folder> [output_prefix] <pair_cores>
#   Rscript dsem_var_granger_loocv_array.R <step1_rds_or_folder> [output_prefix] loo_cores=1 pair_cores=6 chunk_size=2000
#   Rscript dsem_var_granger_loocv_array.R <step1_rds_or_folder> [output_prefix] loo_cores=1 pair_cores=6 chunk_size=2000 array_task_id=1 array_task_count=10
#
# Examples:
#   Rscript dsem_var_granger_loocv_array.R dsem_step1_unknownjob.rds
#   Rscript dsem_var_granger_loocv_array.R DSEM-JOINT-1500subsample granger_pairs
#
# Notes:
#   - This script reads the already-fitted step-1 VAR posterior and does not
#     refit separate reduced models for each pair.
#   - For each directed pair x -> y, it compares:
#       full predictive log-likelihood for y,
#       versus a reduced version where all lag coefficients from x to y are
#       set to 0 in the linear predictor.
#   - The result is an approximate Granger-style ELPD contribution score for the
#     x-lag block in predicting y under the fitted DSEM VAR posterior.
# =============================================================================

# Capture command-line arguments supplied to this script.
args <- commandArgs(trailingOnly = TRUE)
# Require at least one positional argument (input RDS or folder path).
if (length(args) < 1) {
  stop("Usage: Rscript dsem_var_granger_loocv_array.R <step1_rds_or_folder> [output_prefix] [loo_cores] [pair_cores] [chunk_size]")
}

# First argument is the input source (single file or folder of files).
input_path <- args[1]
# Optional output prefix; default keeps output names predictable.
output_prefix <- if (length(args) >= 2) args[2] else "granger_loocv"

# Parse optional core arguments robustly.
# Accepted forms:
#   positional: [loo_cores] [pair_cores]
#   keyed: loo_cores=1 pair_cores=6 chunk_size=2000 (or --loo_cores=1 --pair_cores=6)
#   keyed for sharding: array_task_id=1 array_task_count=10
opt_args <- if (length(args) >= 3) args[3:length(args)] else character(0)
keyed <- opt_args[grepl("=", opt_args)]
positional <- opt_args[!grepl("=", opt_args)]

parse_keyed_int <- function(keyed_args, key) {
  pat <- paste0("^(--)?", key, "=")
  hit <- keyed_args[grepl(pat, keyed_args)]
  if (length(hit) == 0) return(NA_integer_)
  as.integer(sub(pat, "", hit[length(hit)]))
}

loo_cores_keyed <- parse_keyed_int(keyed, "loo_cores")
pair_cores_keyed <- parse_keyed_int(keyed, "pair_cores")
chunk_size_keyed <- parse_keyed_int(keyed, "chunk_size")
array_task_id_keyed <- parse_keyed_int(keyed, "array_task_id")
array_task_count_keyed <- parse_keyed_int(keyed, "array_task_count")

# Optional scheduler fallback (primarily SLURM array jobs).
array_task_id_env <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")))
array_task_count_env <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "")))

array_task_id <- if (!is.na(array_task_id_keyed)) {
  array_task_id_keyed
} else if (!is.na(array_task_id_env)) {
  array_task_id_env
} else if (!is.na(array_task_id_pos)) {
  array_task_id_pos
} else {
  1L
}

array_task_count <- if (!is.na(array_task_count_keyed)) {
  array_task_count_keyed
} else if (!is.na(array_task_count_env)) {
  array_task_count_env
} else if (!is.na(array_task_count_pos)) {
  array_task_count_pos
} else {
  1L
}

if (!is.finite(array_task_id) || array_task_id < 1L) {
  stop("array_task_id must be a positive integer")
}
if (!is.finite(array_task_count) || array_task_count < 1L) {
  stop("array_task_count must be a positive integer")
}
if (array_task_id > array_task_count) {
  stop("array_task_id cannot exceed array_task_count")
}

# Backward/ergonomic behavior:
#   - one positional numeric argument => treat as pair_cores (common intent)
#   - two positional numerics => [loo_cores, pair_cores]
#   - three positional numerics => [loo_cores, pair_cores, chunk_size]
if (length(positional) == 1) {
  loo_cores_pos <- NA_integer_
  pair_cores_pos <- as.integer(positional[1])
} else {
  loo_cores_pos <- if (length(positional) >= 1) as.integer(positional[1]) else NA_integer_
  pair_cores_pos <- if (length(positional) >= 2) as.integer(positional[2]) else NA_integer_
}

chunk_size_pos <- if (length(positional) >= 3) as.integer(positional[3]) else NA_integer_
array_task_id_pos <- if (length(positional) >= 4) as.integer(positional[4]) else NA_integer_
array_task_count_pos <- if (length(positional) >= 5) as.integer(positional[5]) else NA_integer_

# Optional core count for loo::loo calls; defaults to 1 for reproducibility.
loo_cores <- if (!is.na(loo_cores_keyed)) loo_cores_keyed else if (!is.na(loo_cores_pos)) loo_cores_pos else 1L
if (!is.finite(loo_cores) || loo_cores < 1L) {
  stop("loo_cores must be a positive integer")
}
# Detect available physical cores for automatic defaults.
detected_physical_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_physical_cores) || detected_physical_cores < 1L) {
  detected_physical_cores <- 1L
}
# Optional core count for parallel x->y pair evaluation per outcome y.
# If not provided, default to physical cores minus one (at least 1).
pair_cores <- if (!is.na(pair_cores_keyed)) {
  pair_cores_keyed
} else if (!is.na(pair_cores_pos)) {
  pair_cores_pos
} else {
  max(1L, as.integer(detected_physical_cores) - 1L)
}
if (!is.finite(pair_cores) || pair_cores < 1L) {
  stop("pair_cores must be a positive integer")
}
if (pair_cores > 1L && .Platform$OS.type == "windows") {
  warning("pair_cores > 1 is not supported on Windows in this script; using pair_cores=1")
  pair_cores <- 1L
}
if (pair_cores > 1L && loo_cores > 1L) {
  warning("Both pair_cores and loo_cores are >1; consider setting one of them to 1 to avoid CPU oversubscription")
}

# Chunk size for PSIS-LOO evaluation over observation columns.
# Smaller values reduce peak memory at the cost of more loo() calls.
loo_chunk_size <- if (!is.na(chunk_size_keyed)) {
  chunk_size_keyed
} else if (!is.na(chunk_size_pos)) {
  chunk_size_pos
} else {
  2000L
}
if (!is.finite(loo_chunk_size) || loo_chunk_size < 1L) {
  stop("chunk_size must be a positive integer")
}
message(
  "Core settings: detected_physical_cores=", detected_physical_cores,
  ", loo_cores=", loo_cores,
  ", pair_cores=", pair_cores,
  " (set via args 3 and 4 to override)"
)
message(
  "Array shard settings: task_id=", array_task_id,
  ", task_count=", array_task_count,
  " (keyed args or SLURM env)"
)

# If loo_cores is left at 1, use pair_cores for the per-outcome full-model
# LOO cache so the early scoring phase does not bottleneck on one core.
cache_loo_cores <- if (loo_cores > 1L) loo_cores else pair_cores
if (!is.finite(cache_loo_cores) || cache_loo_cores < 1L) {
  cache_loo_cores <- 1L
}

message(
  "Execution cores: cache_loo_cores=", cache_loo_cores,
  ", pair_eval_cores=", pair_cores,
  ", loo_cores_within_pair=", loo_cores,
  ", loo_chunk_size=", loo_chunk_size
)

if (pair_cores > 1L) {
  message("Parallel backend: parallel::mclapply (forking) with mc.cores=", pair_cores)
} else {
  message("Parallel backend: sequential (pair_cores=1)")
}

# Load only the libraries needed for posterior extraction and PSIS-LOO.
suppressPackageStartupMessages({
  library(posterior)
  library(loo)
})

# Helper: robust directory check that treats invalid paths as FALSE.
is_dir <- function(path) {
  # file.info may emit warnings/NA for missing paths.
  info <- suppressWarnings(file.info(path)$isdir)
  # Convert NA to FALSE to simplify downstream checks.
  if (is.na(info)) FALSE else info
}

# Combine multiple step-1 results when runs were sharded across files.
combine_step1_results <- function(rds_list) {
  # Fast path: no combination needed when exactly one object is supplied.
  if (length(rds_list) == 1) return(rds_list[[1]])
  # Convert each draw object to matrix form for row-wise draw binding.
  draws_list <- lapply(rds_list, function(x) posterior::as_draws_matrix(x$draws))
  # Keep prep metadata from the first object; should be identical across shards.
  prep <- rds_list[[1]]$prep
  # Stack posterior draws along the draw dimension.
  draws_combined <- posterior::bind_draws(draws_list, along = "draw")
  # Return object in the same structure expected by downstream code.
  list(draws = draws_combined, prep = prep)
}

# Read either a single step-1 RDS or all RDS files in a directory.
read_step1_input <- function(path) {
  # Branch: folder input means multiple shard files that need combination.
  if (is_dir(path)) {
    # Find all serialized result files in the provided directory.
    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
    # Stop early if directory does not contain usable files.
    if (length(files) == 0) stop("No .rds files found in folder: ", path)
    # Informative progress message for larger folders.
    message("Reading and combining ", length(files), " RDS files from folder: ", path)
    # Read each RDS file.
    rds_list <- lapply(files, readRDS)
    # Merge shard objects to one analysis object.
    step1 <- combine_step1_results(rds_list)
    # Write outputs beside input folder contents.
    out_dir <- path
  } else {
    # Branch: path should be a single file.
    if (!file.exists(path)) stop("Input file not found: ", path)
    message("Reading single RDS file: ", path)
    # Read single serialized step-1 object.
    step1 <- readRDS(path)
    # Write outputs next to this file.
    out_dir <- dirname(path)
  }

  # Validate required fields used by this script.
  if (is.null(step1$draws) || is.null(step1$prep)) {
    stop("Input object must contain $draws and $prep")
  }
  # Return both loaded object and output destination directory.
  list(step1 = step1, out_dir = out_dir)
}

# Numerically stable softplus for Bernoulli log-likelihood calculations.
softplus <- function(x) {
  # Equivalent to log(1 + exp(x)) without overflow/underflow issues.
  log1p(exp(-abs(x))) + pmax(x, 0)
}

# Extract z_person draws and build a lookup from (person, q-index) to matrix col.
extract_z_person_matrix <- function(draws, N, Q) {
  # Enumerate variable names present in posterior draws.
  var_names <- posterior::variables(draws)
  # Keep only person-level latent vectors z_person[i,q].
  z_vars <- var_names[grepl("^z_person\\[", var_names)]
  # Guard against malformed inputs.
  if (length(z_vars) != N * Q) {
    stop("Expected ", N * Q, " z_person variables, found ", length(z_vars))
  }

  # Parse the i,q indices out of variable names.
  idx_txt <- regmatches(z_vars, regexpr("\\[(\\d+),(\\d+)\\]", z_vars))
  i_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx_txt))
  q_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx_txt))

  # Materialize selected variables into a draws matrix [S x (N*Q)].
  z_mat <- posterior::as_draws_matrix(
    posterior::subset_draws(draws, variable = z_vars)
  )

  # Preallocate lookup matrix where lookup[i,q] = column index in z_mat.
  col_lookup <- matrix(NA_integer_, nrow = N, ncol = Q)
  # Fill lookup so later code can pull any person/parameter quickly.
  for (j in seq_along(z_vars)) {
    col_lookup[i_idx[j], q_idx[j]] <- j
  }

  # Return posterior draw matrix and integer lookup table.
  list(z_mat = z_mat, col_lookup = col_lookup)
}

# Build global observation index of observed time points across all persons.
build_obs_index <- function(stan_data) {
  # obs_mask==1 marks valid observations in the original Stan data.
  idx <- which(stan_data$obs_mask == 1L, arr.ind = TRUE)
  # Name index columns for readability in downstream code.
  colnames(idx) <- c("i", "t")
  # Return matrix with one row per observed (person,time) pair.
  idx
}

# Build y-specific response vector and lag design matrix for a target y.
build_y_block_data <- function(stan_data, obs_idx, y_idx, K) {
  # Number of observed rows for this target variable.
  n_obs <- nrow(obs_idx)
  # Person indices for each observed row.
  i_idx <- obs_idx[, "i"]
  # Time indices for each observed row.
  t_idx <- obs_idx[, "t"]

  # Extract response y at each observed row.
  y_vec <- as.numeric(stan_data$Y[cbind(i_idx, t_idx, rep(y_idx, n_obs))])
  # Allocate lag design matrix with K lag predictors.
  X_mat <- matrix(0, nrow = n_obs, ncol = K)
  # Fill each lag column from the 3D Stan array X_lag.
  for (k in seq_len(K)) {
    X_mat[, k] <- as.numeric(stan_data$X_lag[cbind(i_idx, t_idx, rep(k, n_obs))])
  }

  # Return compact data object for this target y.
  list(person_idx = i_idx, y = y_vec, X = X_mat)
}

# Bernoulli-logit pointwise log-likelihood for matrix eta [S x Nobs].
loglik_binary <- function(eta, y_vec) {
  # Broadcast y over draws dimension for vectorized operations.
  y_mat <- matrix(rep(y_vec, each = nrow(eta)), nrow = nrow(eta))
  # log(sigmoid(eta)).
  ll1 <- -softplus(-eta)
  # log(1 - sigmoid(eta)).
  ll0 <- -softplus(eta)
  # Select ll1 for y=1 and ll0 for y=0.
  y_mat * ll1 + (1 - y_mat) * ll0
}

# Gaussian pointwise log-likelihood for matrix eta [S x Nobs].
loglik_gaussian <- function(eta, y_vec, sigma_vec) {
  # Broadcast y over draws dimension for vectorized operations.
  y_mat <- matrix(rep(y_vec, each = nrow(eta)), nrow = nrow(eta))
  # Per-draw variance.
  sigma2 <- sigma_vec^2
  # Quadratic term in Gaussian log-density.
  ll <- sweep(-(y_mat - eta)^2, 1, 2 * sigma2, "/")
  # Add normalization constant term.
  ll <- sweep(ll, 1, -0.5 * log(2 * pi * sigma2), "+")
  # Return pointwise log-likelihood matrix [S x Nobs].
  ll
}

# Compute loo::loo() on a large log-likelihood matrix in column chunks so we
# never materialize an excessively wide intermediate object inside loo.
loo_matrix_chunked <- function(ll, cores = 1L, chunk_size = 2000L) {
  if (!is.matrix(ll)) {
    ll <- as.matrix(ll)
  }

  n_obs <- ncol(ll)
  if (n_obs == 0L) {
    stop("log-likelihood matrix has zero columns")
  }

  chunk_size <- as.integer(chunk_size)
  if (!is.finite(chunk_size) || chunk_size < 1L) {
    chunk_size <- n_obs
  }
  chunk_size <- min(chunk_size, n_obs)

  chunk_starts <- seq(1L, n_obs, by = chunk_size)
  pointwise_chunks <- vector("list", length(chunk_starts))

  for (chunk_i in seq_along(chunk_starts)) {
    start_col <- chunk_starts[chunk_i]
    end_col <- min(n_obs, start_col + chunk_size - 1L)
    ll_chunk <- ll[, start_col:end_col, drop = FALSE]

    loo_chunk <- loo::loo(ll_chunk, cores = cores)
    pointwise_chunks[[chunk_i]] <- loo_chunk$pointwise
  }

  pointwise <- do.call(rbind, pointwise_chunks)
  rownames(pointwise) <- NULL

  # Rebuild a loo object from the combined pointwise summaries.
  loo:::importance_sampling_loo_object(
    pointwise = pointwise,
    diagnostics = list(),
    dims = dim(ll),
    is_method = "psis",
    is_object = NULL
  )
}

# Precompute y-specific quantities that are shared across all x -> y scores.
prepare_y_cache <- function(y_dat, y_idx, z_mat, col_lookup, V, K, is_binary_y, loo_cores) {
  # Number of posterior draws.
  S <- nrow(z_mat)
  # Number of pointwise observations for this target y.
  n_obs <- nrow(y_dat$X)

  # Container for full-model pointwise log-likelihood draws.
  ll_full <- matrix(NA_real_, nrow = S, ncol = n_obs)
  # Cache of person-level objects reused for each predictor x.
  person_cache <- vector("list", length(unique(y_dat$person_idx)))
  pc_id <- 1L

  # q-indices for all lag coefficients feeding target y.
  A_q_idx <- as.integer(((seq_len(K) - 1L) * V) + y_idx)
  # q-index for y-specific intercept/mean term.
  mu_q <- V * K + y_idx
  # q-index for y-specific log-variance term.
  lv_q <- V * K + V + y_idx

  # Process one person at a time so person-specific parameters are matched correctly.
  for (i in unique(y_dat$person_idx)) {
    # Observation rows belonging to person i.
    rows <- which(y_dat$person_idx == i)
    # Person-specific lag predictors and outcomes.
    Xi <- y_dat$X[rows, , drop = FALSE]
    yi <- y_dat$y[rows]

    # Extract all A coefficients for target y for person i across draws.
    Ai <- z_mat[, col_lookup[i, A_q_idx], drop = FALSE]
    # Extract person i mean term for target y across draws.
    mui <- z_mat[, col_lookup[i, mu_q]]
    # Extract person i log-variance term for target y across draws.
    lvi <- z_mat[, col_lookup[i, lv_q]]

    # Full linear predictor: mu + A * X.
    eta_full <- Ai %*% t(Xi)
    eta_full <- sweep(eta_full, 1, mui, "+")

    # Compute full-model log-likelihood under the correct family for outcome y.
    if (is_binary_y) {
      ll_full_i <- loglik_binary(eta_full, yi)
      sigma_i <- NULL
    } else {
      # Convert log-variance to standard deviation.
      sigma_i <- sqrt(exp(lvi))
      ll_full_i <- loglik_gaussian(eta_full, yi, sigma_i)
    }

    # Store full-model blocks back into global matrix.
    ll_full[, rows] <- ll_full_i

    # Persist person-level pieces needed for fast reduced-model reconstruction.
    person_cache[[pc_id]] <- list(
      rows = rows,
      Xi = Xi,
      yi = yi,
      Ai = Ai,
      eta_full = eta_full,
      sigma_i = sigma_i
    )
    pc_id <- pc_id + 1L
  }

  # PSIS-LOO for the full predictive matrix (computed once per y).
  loo_full <- loo_matrix_chunked(ll_full, cores = loo_cores, chunk_size = loo_chunk_size)
  # Return cache used by all x -> y pair evaluations.
  list(
    person_cache = person_cache,
    ll_full = ll_full,
    loo_full = loo_full,
    pointwise_full = loo_full$pointwise[, "elpd_loo"],
    n_obs = n_obs,
    is_binary_y = is_binary_y
  )
}

# Score directed pair x -> y using cached full-model pieces for y.
score_pair_from_cache <- function(y_cache, x_idx, K_per_pair, loo_cores) {
  S <- nrow(y_cache$ll_full)
  n_obs <- y_cache$n_obs
  ll_red <- matrix(NA_real_, nrow = S, ncol = n_obs)

  # k-indices of lag coefficients belonging to predictor block x.
  block_k <- as.integer(seq_len(K_per_pair) + (x_idx - 1L) * K_per_pair)

  # Pool the x -> y lag coefficients (across all posterior draws, all
  # people, and all K_per_pair lags) so we can summarize their overall
  # sign/direction alongside the ELPD comparison.
  coef_pool <- numeric(0)

  for (pc in y_cache$person_cache) {
    # Reconstruct reduced predictor by subtracting x-block contribution.
    contrib_block <- pc$Ai[, block_k, drop = FALSE] %*% t(pc$Xi[, block_k, drop = FALSE])
    eta_red <- pc$eta_full - contrib_block

    if (y_cache$is_binary_y) {
      ll_red_i <- loglik_binary(eta_red, pc$yi)
    } else {
      ll_red_i <- loglik_gaussian(eta_red, pc$yi, pc$sigma_i)
    }
    ll_red[, pc$rows] <- ll_red_i

    coef_pool <- c(coef_pool, as.vector(pc$Ai[, block_k, drop = FALSE]))
  }

  mean_coef_x_to_y <- mean(coef_pool)
  prop_positive_x_to_y <- mean(coef_pool > 0)
  coef_direction_x_to_y <- if (prop_positive_x_to_y >= 0.95) {
    "positive"
  } else if (prop_positive_x_to_y <= 0.05) {
    "negative"
  } else {
    "mixed"
  }

  # Only reduced-model loo is recomputed per x -> y pair.
  loo_red <- loo_matrix_chunked(ll_red, cores = loo_cores, chunk_size = loo_chunk_size)

  elpd_full <- y_cache$loo_full$estimates["elpd_loo", "Estimate"]
  elpd_red <- loo_red$estimates["elpd_loo", "Estimate"]
  elpd_diff <- elpd_full - elpd_red

  pt_diff <- y_cache$pointwise_full - loo_red$pointwise[, "elpd_loo"]
  se_diff <- sqrt(length(pt_diff) * stats::var(pt_diff))

  list(
    elpd_full = elpd_full,
    elpd_reduced = elpd_red,
    elpd_diff = elpd_diff,
    se_diff = se_diff,
    n_obs = n_obs,
    loo_red = loo_red,
    mean_coef_x_to_y = mean_coef_x_to_y,
    prop_positive_x_to_y = prop_positive_x_to_y,
    coef_direction_x_to_y = coef_direction_x_to_y
  )
}

# Load and validate input object(s).
loaded <- read_step1_input(input_path)
# Extract the step-1 object.
step1 <- loaded$step1
# Directory where output files will be written.
out_dir <- loaded$out_dir

# Extract prep metadata and stan data arrays.
prep <- step1$prep
stan_data <- prep$stan_data
behaviours <- prep$behaviours
V <- if (!is.null(prep$V)) prep$V else stan_data$V
K <- if (!is.null(prep$K)) prep$K else stan_data$K
Q <- prep$Q
N <- stan_data$N
# Number of lags per predictor block.
K_per_pair <- as.integer(K / V)

# Sanity check on lag-block structure.
if (K %% V != 0L) stop("K must be divisible by V")

# Outcome-type indicator from step-1 data prep.
is_binary <- as.logical(stan_data$is_binary)

# Assign this task a subset of outcomes in round-robin fashion for array jobs.
all_y_idx <- seq_len(V)
y_idx_assigned <- all_y_idx[((all_y_idx - 1L) %% array_task_count) + 1L == array_task_id]
message(
  "Assigned outcomes for this shard: ", length(y_idx_assigned),
  " of ", V, " total"
)
if (length(y_idx_assigned) > 0) {
  message("Assigned outcome names: ", paste(behaviours[y_idx_assigned], collapse = ", "))
}

# Use a shard-specific output prefix when running as an array.
output_prefix_effective <- if (array_task_count > 1L) {
  paste0(
    output_prefix,
    "_task", sprintf("%03d", array_task_id),
    "of", sprintf("%03d", array_task_count)
  )
} else {
  output_prefix
}

# Wall-clock timer for lightweight runtime reporting.
run_start_time <- proc.time()[["elapsed"]]

message("Extracting z_person draws...")
# Extract z-person posterior draws plus lookup table.
z_info <- extract_z_person_matrix(step1$draws, N = N, Q = Q)
z_mat <- z_info$z_mat
col_lookup <- z_info$col_lookup

message("Building observation index...")
# Build common observation index from mask.
obs_idx <- build_obs_index(stan_data)

message("Precomputing y-specific design blocks...")
# Precompute y-specific data so each pair score can reuse these efficiently.
y_blocks <- vector("list", V)
for (y_idx in seq_len(V)) {
  y_blocks[[y_idx]] <- build_y_block_data(stan_data, obs_idx, y_idx, K)
}

# Collect pairwise result rows in a list and bind once at the end.
results <- list()

message("Scoring directed pairs x -> y with ELPD...")
# Loop over each possible outcome y.
for (y_idx in y_idx_assigned) {
  outcome_start_time <- proc.time()[["elapsed"]]
  y_name <- behaviours[y_idx]
  y_dat <- y_blocks[[y_idx]]
  is_binary_y <- is_binary[y_idx]

  message(
    "Outcome ", y_name,
    ": precomputing full-model LOO with cache_loo_cores=", cache_loo_cores,
    "; then evaluating ", V - 1L, " pair(s) with pair_cores=", pair_cores,
    "."
  )

  # Build y-specific full-model cache once, reused for all x -> y comparisons.
  cache_start_time <- proc.time()[["elapsed"]]
  y_cache <- prepare_y_cache(
    y_dat = y_dat,
    y_idx = y_idx,
    z_mat = z_mat,
    col_lookup = col_lookup,
    V = V,
    K = K,
    is_binary_y = is_binary_y,
    loo_cores = cache_loo_cores
  )
  cache_elapsed <- proc.time()[["elapsed"]] - cache_start_time

  # Predictor candidates for this fixed outcome y.
  x_candidates <- setdiff(seq_len(V), y_idx)

  # Build one result row for a single x -> y pair.
  eval_one_pair <- function(x_idx) {
    x_name <- behaviours[x_idx]
    message("  Pair ", x_name, " -> ", y_name)

    pair_score <- score_pair_from_cache(
      y_cache = y_cache,
      x_idx = x_idx,
      K_per_pair = K_per_pair,
      loo_cores = loo_cores
    )

    data.frame(
      predictor_x = x_name,
      outcome_y = y_name,
      outcome_is_binary = is_binary_y,
      n_obs = pair_score$n_obs,
      elpd_full = pair_score$elpd_full,
      elpd_reduced_drop_xlags = pair_score$elpd_reduced,
      elpd_diff_full_minus_reduced = pair_score$elpd_diff,
      se_diff = pair_score$se_diff,
      # z-score style standardized effect size for ranking/interpretation.
      z_score = if (is.finite(pair_score$se_diff) && pair_score$se_diff > 0) {
        pair_score$elpd_diff / pair_score$se_diff
      } else {
        NA_real_
      },
      # Basic support criterion: any positive ELPD difference.
      supports_x_granger_causing_y = pair_score$elpd_diff > 0,
      # Stronger criterion: ELPD difference exceeds 2*SE.
      strong_support_diff_gt_2se = if (is.finite(pair_score$se_diff) && pair_score$se_diff > 0) {
        pair_score$elpd_diff > (2 * pair_score$se_diff)
      } else {
        FALSE
      },
      # Summary of the sign of the x -> y lag coefficients (A[y,x,lag])
      # pooled across all posterior draws, all people, and all lags.
      mean_coef_x_to_y = pair_score$mean_coef_x_to_y,
      prop_positive_x_to_y = pair_score$prop_positive_x_to_y,
      coef_direction_x_to_y = pair_score$coef_direction_x_to_y,
      stringsAsFactors = FALSE
    )
  }

  # Evaluate all x -> y pairs for this y, optionally in parallel.
  y_rows <- if (pair_cores > 1L) {
    parallel::mclapply(
      x_candidates,
      eval_one_pair,
      mc.cores = pair_cores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(x_candidates, eval_one_pair)
  }
  pair_elapsed <- proc.time()[["elapsed"]] - (cache_start_time + cache_elapsed)

  message(
    "Outcome ", y_name,
    " timing (sec): cache=", format(round(cache_elapsed, 2), nsmall = 2),
    ", pair_eval=", format(round(pair_elapsed, 2), nsmall = 2),
    ", total=", format(round(proc.time()[["elapsed"]] - outcome_start_time, 2), nsmall = 2)
  )

  results <- c(results, y_rows)
}

# Bind list of row data frames into one result table.
if (length(results) > 0) {
  results_df <- do.call(rbind, results)
} else {
  results_df <- data.frame(
    predictor_x = character(0),
    outcome_y = character(0),
    outcome_is_binary = logical(0),
    n_obs = integer(0),
    elpd_full = numeric(0),
    elpd_reduced_drop_xlags = numeric(0),
    elpd_diff_full_minus_reduced = numeric(0),
    se_diff = numeric(0),
    z_score = numeric(0),
    supports_x_granger_causing_y = logical(0),
    strong_support_diff_gt_2se = logical(0),
    mean_coef_x_to_y = numeric(0),
    prop_positive_x_to_y = numeric(0),
    coef_direction_x_to_y = character(0),
    stringsAsFactors = FALSE
  )
}
# Rank by largest predictive gain when including x-lag block.
if (nrow(results_df) > 0) {
  results_df <- results_df[order(results_df$elpd_diff_full_minus_reduced, decreasing = TRUE), , drop = FALSE]
}

# Output file paths.
csv_file <- file.path(out_dir, paste0(output_prefix_effective, "_pairwise_elpd.csv"))
rds_file <- file.path(out_dir, paste0(output_prefix_effective, "_pairwise_elpd.rds"))

# Save results as CSV for easy inspection and RDS for exact reload.
write.csv(results_df, csv_file, row.names = FALSE)
saveRDS(results_df, rds_file)

total_elapsed <- proc.time()[["elapsed"]] - run_start_time
message("Total runtime (sec): ", format(round(total_elapsed, 2), nsmall = 2))

# Report output locations.
message("\nSaved pairwise ELPD results:")
message("  ", csv_file)
message("  ", rds_file)

# Print top ranked directed pairs for quick terminal review.
cat("\nTop pairs by ELPD advantage (full - drop x lags):\n")
print(utils::head(results_df, 15))
