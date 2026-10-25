# =============================================================================
# dsem_step1.R
#
# Step 1 of the two-step DSEM pipeline.
#
# Fits a multilevel VAR to the infant behaviour time-series data and saves
# the full posterior draws to an RDS file for use by dsem_step2.R.
#
# Usage (interactive):
#   source("dsem_step1.R")
#
# Usage (cluster, with job ID passed as argument):
#   Rscript dsem_step1.R <jobid>
#
# Output:
#   dsem_step1_<jobid>.rds   — list with elements:
#     $fit      CmdStanFit object (CSV files still on disk)
#     $draws    posterior::draws_array  ← eagerly extracted before job ends
#     $prep     Stan data + metadata (participants, behaviours, lag names, etc.)
#
# NOTE: $draws is extracted immediately after sampling so that the RDS is
# self-contained even if the CmdStan CSV files are later cleaned up.
# dsem_step2.R reads only $draws and $prep — it does not need $fit.
# =============================================================================

args  <- commandArgs(trailingOnly = TRUE)
jobid <- if (length(args) >= 1) args[1] else "unknownjob"


# =============================================================================
# STAN MODEL
# =============================================================================

step1_stan <- "
// dsem_step1_var.stan
// Multilevel VAR model — no outcome regression.
//
// All behaviours are binary (0/1), so the likelihood is Bernoulli-logit
// only. There is no residual-variance parameter: for a binary outcome the
// variance is fixed by the probability, p(1 - p), so a per-person
// log-variance would never enter the likelihood and would just reproduce
// its prior.
//
// generated quantities block exposes z_person[i], the Q-vector
//   z_i = [ vec(A[i]) | mu_person[i] ]
// which the bridge step in dsem_step2.R uses to construct the Step 2
// measurement model.

functions {
  // Per-person contribution to the VAR log-likelihood, summed over all
  // observed timepoints. Vectorized over behaviours (all binary, logit
  // link) so each (person, timepoint) costs one vectorized statement
  // instead of V scalar ones. Designed for use with reduce_sum to
  // parallelize over people.
  real partial_sum_lpmf(array[] int i_idx,
                         int start, int end,
                         data array[,,] real Y,
                         data array[,,] real X_lag,
                         data array[] int T_i,
                         data array[,] int obs_mask,
                         array[] matrix A,
                         array[] vector mu_person) {
    real lp = 0;

    for (i in i_idx) {
      for (t in 1:T_i[i]) {
        if (obs_mask[i, t] == 0) continue;

        vector[cols(A[i])] x_t = to_vector(X_lag[i, t]);
        vector[rows(A[i])] eta = mu_person[i] + A[i] * x_t;

        vector[rows(A[i])] y_t = to_vector(Y[i, t]);
        // Manual Bernoulli-logit log-density (avoids bernoulli_logit_lpmf,
        // whose templated gradient code triggers a GCC 14 ICE via
        // Boost's unchecked_factorial table). Bernoulli has no
        // combinatorial term, so this is exact, not an approximation:
        //   log p(y|eta) = y*log_inv_logit(eta) + (1-y)*log1m_inv_logit(eta)
        lp += sum(y_t       .* log_inv_logit(eta)
                  + (1 - y_t) .* log1m_inv_logit(eta));
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> V;
  int<lower=1> K;
  int<lower=1> T_max;

  array[N] int<lower=1>            T_i;
  array[N, T_max, V]    real        Y;
  array[N, T_max, K]    real        X_lag;
  array[N, T_max]       int<lower=0,upper=1> obs_mask;
  array[V]              int<lower=0,upper=1> is_binary;

  real<lower=0> prior_sd_A;
  real<lower=0> prior_sd_intercept;
  real<lower=0> prior_sd_sigma;

  int<lower=1> grainsize;
}

transformed data {
  // The model has no residual-variance parameters, so it is only valid for
  // binary outcomes.
  for (v in 1:V)
    if (is_binary[v] != 1)
      reject(\"Behaviour \", v, \" is not binary; this model is Bernoulli-logit only.\");

  array[N] int seq_N;
  for (i in 1:N) seq_N[i] = i;
}

parameters {
  matrix[V, K]           A_pop;
  matrix<lower=0>[V, K]  sigma_A;
  vector[V]              mu_pop;
  vector<lower=0>[V]     sigma_mu;

  // Non-centred person-level deviations
  array[N] matrix[V, K]  A_z;
  array[N] vector[V]     mu_z;
}

transformed parameters {
  array[N] matrix[V, K]        A;
  array[N] vector[V]           mu_person;

  for (i in 1:N) {
    A[i]             = A_pop + sigma_A .* A_z[i];
    mu_person[i]     = mu_pop + sigma_mu .* mu_z[i];
  }
}

model {
  // Population priors
  to_vector(A_pop)   ~ normal(0, prior_sd_A);
  to_vector(sigma_A) ~ normal(0, prior_sd_A);
  mu_pop             ~ normal(0, prior_sd_intercept);
  sigma_mu           ~ normal(0, prior_sd_intercept);
  // Non-centred person-level priors
  for (i in 1:N) {
    to_vector(A_z[i]) ~ std_normal();
    mu_z[i]           ~ std_normal();
  }

  // VAR likelihood (vectorized over behaviours; called directly without
  // reduce_sum/threading to avoid a GCC 14 + TBB compiler ICE on this
  // cluster). start/end are unused by partial_sum_lpmf's body (it just
  // loops over i_idx), so passing the full 1:N range here is equivalent
  // to evaluating all N people's contributions.
  target += partial_sum_lpmf(seq_N | 1, N,
                              Y, X_lag, T_i, obs_mask,
                              A, mu_person);
}

generated quantities {
  // z_person[i] is the full Q-vector of person i's VAR parameters.
  // The bridge step in dsem_step2.R reads these draws to form the
  // Gaussian measurement model that propagates uncertainty into Step 2.
  array[N] vector[V * K + V] z_person;
  for (i in 1:N) {
    z_person[i] = append_row(to_vector(A[i]), mu_person[i]);
  }
  matrix[V, K] A_pop_out = A_pop;
}
"


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

import_data <- function(file_name) {
  message(file_name)
  dat <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
  dat <- t(dat)
  data.frame(
    PID  = basename(file_name),
    TIME = seq_len(nrow(dat)),
    dat,
    stringsAsFactors = FALSE
  )
}

extract_child_id <- function(pid) sub("_(11m|12m)\\.csv$", "", pid)

make_lagged_comp <- function(df_person, vars, p = 1, min_T = 30) {
  df_person[vars] <- lapply(df_person[vars], as.numeric)
  n_total <- nrow(df_person)
  if (n_total <= p) return(NULL)
  lagged_list <- list()
  lag_names   <- character(0)
  for (lag in seq_len(p)) {
    for (v in vars) {
      nm <- paste0(v, ".l", lag)
      lagged_list[[nm]] <- c(rep(NA, lag),
                             if (lag < n_total) df_person[[v]][1:(n_total - lag)]
                             else rep(NA, n_total))
      lag_names <- c(lag_names, nm)
    }
  }
  combined    <- cbind(df_person[vars],
                       as.data.frame(lagged_list, stringsAsFactors = FALSE))
  complete_ix <- complete.cases(combined)
  if (sum(complete_ix) < min_T) return(NULL)
  combined[complete_ix, , drop = FALSE]
}


# =============================================================================
# prepare_step1_data
# =============================================================================

prepare_step1_data <- function(
  combined_data,
  behaviours,
  p,
  id_col             = "participant_id",
  prior_sd_A         = 1,
  prior_sd_intercept = 5,
  prior_sd_sigma     = 2.5,
  ts_method          = c("subsample", "aggregate", "none"),
  T_subsample        = 500,
  window_size        = 10,
  grainsize          = 1
) {
  ts_method    <- match.arg(ts_method)
  participants <- unique(combined_data[[id_col]])
  V            <- length(behaviours)
  K            <- p * V

  all_comps <- lapply(participants, function(pid) {
    dfp <- combined_data[combined_data[[id_col]] == pid, , drop = FALSE]
    make_lagged_comp(dfp, behaviours, p = p, min_T = p + 1)
  })
  names(all_comps) <- as.character(participants)

  valid        <- !vapply(all_comps, is.null, logical(1))
  all_comps    <- all_comps[valid]
  participants <- participants[valid]
  N            <- length(participants)

  lag_col_names <- as.vector(
    outer(seq_len(p), behaviours, function(j, v) paste0(v, ".l", j))
  )

  reduce_ts <- function(comp, pid) {
    Ti <- nrow(comp)
    if (ts_method == "none") return(comp)
    if (ts_method == "subsample") {
      if (Ti <= T_subsample) return(comp)
      idx <- round(seq(1, Ti, length.out = T_subsample))
      message("  [", pid, "] subsampling: ", Ti, " -> ", length(idx))
      return(comp[idx, , drop = FALSE])
    }
    if (ts_method == "aggregate") {
      n_windows <- floor(Ti / window_size)
      if (n_windows == 0L) return(comp)
      message("  [", pid, "] aggregating: ", Ti, " -> ", n_windows, " windows")
      rows <- lapply(seq_len(n_windows), function(w) {
        idx  <- ((w - 1) * window_size + 1):(w * window_size)
        resp <- colMeans(comp[idx, behaviours,    drop = FALSE], na.rm = TRUE)
        lags <- colMeans(comp[idx, lag_col_names, drop = FALSE], na.rm = TRUE)
        as.data.frame(as.list(c(resp, lags)))
      })
      return(do.call(rbind, rows))
    }
  }

  message("Applying ts_method='", ts_method, "'...")
  all_comps <- mapply(reduce_ts, all_comps, names(all_comps), SIMPLIFY = FALSE)

  T_i   <- vapply(all_comps, nrow, integer(1))
  T_max <- max(T_i)
  message("T_max=", T_max, "  N=", N, "  total obs=", sum(T_i))

  is_binary <- vapply(behaviours, function(b) {
    vals <- na.omit(unique(combined_data[[b]]))
    as.integer(length(vals) == 2L && all(vals %in% c(0, 1)))
  }, integer(1))
  message("Binary behaviours: ",
          paste(behaviours[is_binary == 1], collapse = ", "))
  if (any(is_binary == 0L))
    stop("Non-binary behaviour(s): ",
         paste(behaviours[is_binary == 0L], collapse = ", "),
         ". The Step 1 model is Bernoulli-logit only (no residual variance).")

  X_all_flat <- matrix(NA_real_, nrow = sum(T_i), ncol = K)
  row_ptr    <- 1L
  for (i in seq_along(participants)) {
    Ti <- T_i[i]
    for (k in seq_len(K))
      X_all_flat[row_ptr:(row_ptr + Ti - 1L), k] <-
        as.numeric(all_comps[[i]][[lag_col_names[k]]])
    row_ptr <- row_ptr + Ti
  }
  col_mn <- colMeans(X_all_flat, na.rm = TRUE)
  col_sd <- apply(X_all_flat, 2, sd, na.rm = TRUE)
  col_sd[col_sd == 0 | !is.finite(col_sd)] <- 1

  Y_arr    <- array(0,  dim = c(N, T_max, V))
  X_arr    <- array(0,  dim = c(N, T_max, K))
  obs_mask <- array(0L, dim = c(N, T_max))

  for (i in seq_along(participants)) {
    comp <- all_comps[[i]]
    Ti   <- T_i[i]
    for (v in seq_len(V)) {
      y_raw <- as.numeric(comp[[behaviours[v]]])
      y_raw[!is.finite(y_raw)] <- 0
      Y_arr[i, seq_len(Ti), v] <- y_raw
    }
    for (k in seq_len(K)) {
      x_raw <- (as.numeric(comp[[lag_col_names[k]]]) - col_mn[k]) / col_sd[k]
      x_raw[!is.finite(x_raw)] <- 0
      X_arr[i, seq_len(Ti), k] <- x_raw
    }
    obs_mask[i, seq_len(Ti)] <- 1L
  }

  list(
    stan_data = list(
      N        = N,
      V        = V,
      K        = K,
      T_max    = T_max,
      T_i      = as.integer(T_i),
      Y        = Y_arr,
      X_lag    = X_arr,
      obs_mask = obs_mask,
      is_binary          = as.integer(is_binary),
      prior_sd_A         = prior_sd_A,
      prior_sd_intercept = prior_sd_intercept,
      prior_sd_sigma     = prior_sd_sigma,
      grainsize          = as.integer(grainsize)
    ),
    participants  = as.character(participants),
    behaviours    = behaviours,
    lag_col_names = lag_col_names,
    X_lag_center  = col_mn,
    X_lag_scale   = col_sd,
    V             = V,
    K             = K,
    Q             = V * K + V
  )
}


# =============================================================================
# fit_step1
# =============================================================================

fit_step1 <- function(
  combined_data,
  behaviours,
  p                  = 5,
  id_col             = "participant_id",
  prior_sd_A         = 1,
  prior_sd_intercept = 5,
  prior_sd_sigma     = 2.5,
  ts_method          = "subsample",
  T_subsample        = 500,
  window_size        = 10,
  chains             = 2,
  iter_warmup        = 500,
  iter_sampling      = 500,
  cores              = 4,
  grainsize          = 1,
  adapt_delta        = 0.95,
  max_treedepth      = 10,
  seed               = NULL,
  output_dir         = tempdir(),
  use_pathfinder     = TRUE,
  pathfinder_paths   = 4,
  pathfinder_draws   = 200,
  ...
) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Install 'cmdstanr'")

  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  message("=== STEP 1: Writing and compiling Stan model ===")
  stan_file <- cmdstanr::write_stan_file(step1_stan)
  mod       <- cmdstanr::cmdstan_model(stan_file, compile = TRUE,
                                        force_recompile = TRUE)

  message("=== STEP 1: Preparing data ===")
  prep <- prepare_step1_data(
    combined_data      = combined_data,
    behaviours         = behaviours,
    p                  = p,
    id_col             = id_col,
    prior_sd_A         = prior_sd_A,
    prior_sd_intercept = prior_sd_intercept,
    prior_sd_sigma     = prior_sd_sigma,
    ts_method          = ts_method,
    T_subsample        = T_subsample,
    window_size        = window_size,
    grainsize          = grainsize
  )

  V <- prep$stan_data$V
  K <- prep$stan_data$K
  N <- prep$stan_data$N

  init_fn <- function() list(
    A_pop        = matrix(0,   nrow = V, ncol = K),
    sigma_A      = matrix(0.1, nrow = V, ncol = K),
    A_z          = array(0,    dim  = c(N, V, K)),
    mu_pop       = rep(0,   V),
    sigma_mu     = rep(0.1, V),
    mu_z         = matrix(0, nrow = N, ncol = V)
  )

  # Run Pathfinder to get a fast variational approximation, then use its
  # draws to initialize HMC. This gives HMC a starting point much closer
  # to the typical set (and, via init-from-draws, an implicit starting
  # mass-matrix scale), which can substantially shorten the iter_warmup
  # needed for step-size/metric adaptation.
  init_arg <- init_fn
  if (use_pathfinder) {
    message("=== STEP 1: Running Pathfinder for initialization ===")
    pf <- tryCatch(
      mod$pathfinder(
        data       = prep$stan_data,
        num_paths  = pathfinder_paths,
        draws      = pathfinder_draws,
        seed       = seed
      ),
      error = function(e) {
        warning("Pathfinder failed (", conditionMessage(e),
                "); falling back to default init.")
        NULL
      }
    )
    if (!is.null(pf)) {
      init_arg <- pf  # cmdstanr accepts a CmdStanPathfinder object as init=
    }
  }

  message("=== STEP 1: Sampling (chains=", chains,
          ", warmup=", iter_warmup, ", sampling=", iter_sampling,
          ", seed=", seed, ") ===")
  fit <- mod$sample(
    data              = prep$stan_data,
    chains            = chains,
    parallel_chains   = cores,
    iter_warmup       = iter_warmup,
    iter_sampling     = iter_sampling,
    adapt_delta       = adapt_delta,
    max_treedepth     = max_treedepth,
    seed              = seed,
    init              = init_arg,
    output_dir        = output_dir,
    metric      = "diag_e",
    save_warmup = TRUE,
    ...
  )

  list(fit = fit, prep = prep, draws = fit$draws())
}


# =============================================================================
# DATA IMPORT
# =============================================================================
library(dplyr)

all_vars <- c("contingent", "gaze", "CV", "non_CV", "index_point",
              "open_point", "give", "show", "other_gesture")

file_names <- dir("data_full", pattern = "_11m\\.csv$",
                  full.names = TRUE, recursive = TRUE)
data_list  <- lapply(file_names, import_data)
df         <- data.frame(do.call(rbind, data_list), row.names = NULL)
names(df)  <- c("PID", "TIME", all_vars)
df$month   <- 11

file_names_12m <- dir("data_full", pattern = "_12m\\.csv$",
                      full.names = TRUE, recursive = TRUE)
data_list_12m  <- lapply(file_names_12m, import_data)
df_12m         <- data.frame(do.call(rbind, data_list_12m), row.names = NULL)
names(df_12m)  <- c("PID", "TIME", all_vars)
df_12m$month   <- 12

df$child_id     <- extract_child_id(df$PID)
df_12m$child_id <- extract_child_id(df_12m$PID)
df_all          <- bind_rows(df, df_12m)
parts           <- unique(df_all$child_id)

behaviour_cols <- all_vars
behaviours     <- all_vars
combined_data  <- data.frame()

for (i in seq_along(parts)) {
  thispart_data <- df_all |>
    filter(child_id == parts[i]) |>
    arrange(month, TIME) |>
    select(all_of(behaviour_cols))
  if (nrow(thispart_data) > 25) {
    combined_data <- rbind(
      combined_data,
      cbind(thispart_data,
            participant_id = gsub("p([0-9]+)", "\\1", parts[i]))
    )
  }
}
combined_data <- as.data.frame(lapply(combined_data, as.numeric))
message("combined_data: ", nrow(combined_data), " rows, ",
        length(unique(combined_data$participant_id)), " participants")


# =============================================================================
# RUN
# =============================================================================

step1_result <- fit_step1(
  combined_data      = combined_data,
  behaviours         = behaviours,
  p                  = 5,
  # Priors as reported in the paper: normal(0, 5) for the population
  # intercepts and for children's variation around them (sigma_mu), and
  # normal(0, 1) for the lag coefficients and their between-child SDs
  # (sigma_A). prior_sd_sigma is not used by the model.
  prior_sd_A         = 1,
  prior_sd_intercept = 5,
  prior_sd_sigma     = 2.5,
  ts_method          = "subsample",
  T_subsample 	     = 3165,
  chains             = 1,
  iter_warmup        = 200,
  iter_sampling      = 400,
  cores              = 1,
  adapt_delta        = 0.95,
  refresh            = 1, 
  seed=as.integer(jobid),
)

# $draws is already extracted inside fit_step1(); saving here makes the RDS
# self-contained so dsem_step2.R does not need the CmdStan CSV files.
saveRDS(step1_result, file = paste0("dsem_step1_", jobid, ".rds"))
message("Step 1 complete. Saved to dsem_step1_", jobid, ".rds")
