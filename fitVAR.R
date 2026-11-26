# =============================================================================
# fitVAR.R
#
# Step 1 of the two-step DSEM pipeline.
#
# Fits a multilevel VAR to the infant behaviour time-series data and saves
# the full posterior draws to an RDS file for use by the later stages
# (check_convergence.R, granger_analysis.R, granger_direction.R,
# select_terms_for_individual_behaviours.R, final_term_weighting.R and
# make_bridge.R). Run it from a folder that contains (or links to) data_full/.
#
# Usage (interactive):
#   source("fitVAR.R")
#
# Usage (one chain; the argument is both the seed and the file suffix):
#   Rscript fitVAR.R <chain>
#
# Output:
#   dsem_step1_<jobid>.rds   — list with elements:
#     $fit      CmdStanFit object (CSV files still on disk)
#     $draws    posterior::draws_array  ← eagerly extracted before job ends
#     $prep     Stan data + metadata (participants, behaviours, lag names, etc.)
#
# NOTE: $draws is extracted immediately after sampling so that the RDS is
# self-contained even if the CmdStan CSV files are later cleaned up.
# The later stages read only $draws and $prep — they do not need $fit.
#
# OPTIONS (environment variables; the defaults are the settings of the paper)
#   STEP1_SESSION=1    (default) Treat the 11- and 12-month recordings as separate
#                      series: lags are built within each session (the first
#                      samples of a session never take lags from the end of
#                      the other session), and each behaviour gets a session
#                      effect on its intercept (population mean + child-level
#                      deviation, both on the log-odds scale). STEP1_SESSION=0:
#                      the two sessions are concatenated into one series per child.
#   STEP1_ONSET=1      Onset (transition) model: each behaviour is predicted
#                      only at samples where it was absent in the previous
#                      200 ms sample, so the model predicts the START of a
#                      behaviour rather than its continuation. All behaviours'
#                      lags are still predictors. (The lag-1 self-lag is then
#                      uninformed by the data and stays at its prior.)
#   STEP1_LAG_STEP=s   Spacing of the lags in 200 ms samples (default 1). With
#                      s = 2 the five lags are t-2, t-4, ..., t-10, i.e. a 2 s
#                      window instead of 1 s.
#   STEP1_T_SUBSAMPLE, STEP1_ITER_WARMUP, STEP1_ITER_SAMPLING, STEP1_N_CHILDREN
#                      run-size settings (defaults 3165, 1000, 1000, all);
#                      only change these for quick tests.
# The options used are saved in $prep$options and are read by the
# downstream scripts (granger_analysis.R needs them).
# =============================================================================

args  <- commandArgs(trailingOnly = TRUE)
jobid <- if (length(args) >= 1) args[1] else "unknownjob"

env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(as.integer(default))
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) stop(name, " must be an integer (got '", v, "')")
  out
}
opt_session   <- env_int("STEP1_SESSION", 1L) == 1L
opt_onset     <- env_int("STEP1_ONSET", 0L) == 1L
opt_lag_step  <- env_int("STEP1_LAG_STEP", 1L)
opt_T_sub     <- env_int("STEP1_T_SUBSAMPLE", 3165L)
opt_warmup    <- env_int("STEP1_ITER_WARMUP", 1000L)
opt_sampling  <- env_int("STEP1_ITER_SAMPLING", 1000L)
opt_n_child   <- env_int("STEP1_N_CHILDREN", 0L)     # 0 = all children
if (opt_lag_step < 1L) stop("STEP1_LAG_STEP must be >= 1")
message("Step 1 options: session=", opt_session, "  onset=", opt_onset,
        "  lag_step=", opt_lag_step, "  T_subsample=", opt_T_sub,
        "  warmup=", opt_warmup, "  sampling=", opt_sampling,
        if (opt_n_child > 0) paste0("  n_children=", opt_n_child) else "")


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
// Likelihood is computed with one matrix product per child (all of that
// child's observed time points at once) rather than one small
// matrix-vector product per time point. The two are equivalent; the matrix
// form needs far less automatic differentiation and runs much faster.
//
// Optional data flags (paper: use_session = 1, use_onset = 0):
//   use_session = 1: each behaviour's intercept differs between the two
//     recording sessions by delta_person[i][v] (population mean delta_pop,
//     child-level SD sigma_delta), multiplying a centred session indicator.
//   use_onset = 1: y[i, t, v] only enters the likelihood where
//     y_include[i, t, v] = 1 (behaviour v absent at the previous sample), so
//     the model is for the onset of each behaviour.
//
// generated quantities block exposes z_person[i], the Q-vector
//   z_i = [ vec(A[i]) | mu_person[i] ]
// which the Step 2 scripts (and make_bridge.R) use to construct the Step 2
// measurement model.

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
  real<lower=0> prior_sd_sigma_mu;   // SD of half-normal prior on sigma_mu

  int<lower=1> grainsize;

  int<lower=0, upper=1> use_session;
  array[N, T_max]       real session_x;           // centred session indicator
  int<lower=0, upper=1> use_onset;
  array[N, T_max, V]    int<lower=0, upper=1> y_include;
}

transformed data {
  // The model has no residual-variance parameters, so it is only valid for
  // binary outcomes.
  for (v in 1:V)
    if (is_binary[v] != 1)
      reject(\"Behaviour \", v, \" is not binary; this model is Bernoulli-logit only.\");

  // Stack every child's observed time points into one design matrix,
  // with a leading column of 1s for the intercept, so each child's linear
  // predictor is a single product: X1[rows of i] * [mu_i | A_i]'.
  array[N] int n_obs;
  for (i in 1:N) {
    n_obs[i] = 0;
    for (t in 1:T_i[i]) n_obs[i] += obs_mask[i, t];
  }
  int n_total = sum(n_obs);
  array[N] int row_start;
  array[N] int row_end;
  matrix[n_total, K + 1] X1;
  matrix[n_total, V]     Y_obs;
  vector[n_total]        S_obs;    // session indicator (used if use_session)
  matrix[n_total, V]     W_obs;    // 1 = enters the likelihood (used if use_onset)
  {
    int r = 0;
    for (i in 1:N) {
      row_start[i] = r + 1;
      for (t in 1:T_i[i]) {
        if (obs_mask[i, t] == 0) continue;
        r += 1;
        X1[r, 1] = 1;
        for (k in 1:K) X1[r, k + 1] = X_lag[i, t, k];
        for (v in 1:V) Y_obs[r, v] = Y[i, t, v];
        S_obs[r] = session_x[i, t];
        for (v in 1:V) W_obs[r, v] = y_include[i, t, v];
      }
      row_end[i] = r;
    }
  }
}

parameters {
  matrix[V, K]           A_pop;
  matrix<lower=0>[V, K]  sigma_A;
  vector[V]              mu_pop;
  vector<lower=0>[V]     sigma_mu;

  // Non-centred person-level deviations
  array[N] matrix[V, K]  A_z;
  array[N] vector[V]     mu_z;

  // Session effects (zero-length unless use_session = 1)
  vector[use_session * V]           delta_pop;
  vector<lower=0>[use_session * V]  sigma_delta;
  array[N] vector[use_session * V]  delta_z;
}

transformed parameters {
  array[N] matrix[V, K]        A;
  array[N] vector[V]           mu_person;
  array[N] vector[use_session * V] delta_person;

  for (i in 1:N) {
    A[i]             = A_pop + sigma_A .* A_z[i];
    mu_person[i]     = mu_pop + sigma_mu .* mu_z[i];
    if (use_session == 1)
      delta_person[i] = delta_pop + sigma_delta .* delta_z[i];
  }
}

model {
  // Population priors
  to_vector(A_pop)   ~ normal(0, prior_sd_A);
  to_vector(sigma_A) ~ normal(0, prior_sd_A);
  mu_pop             ~ normal(0, prior_sd_intercept);
  sigma_mu           ~ normal(0, prior_sd_sigma_mu);
  // Non-centred person-level priors
  for (i in 1:N) {
    to_vector(A_z[i]) ~ std_normal();
    mu_z[i]           ~ std_normal();
  }
  if (use_session == 1) {
    delta_pop   ~ normal(0, 1);
    sigma_delta ~ normal(0, 1);
    for (i in 1:N) delta_z[i] ~ std_normal();
  }

  // VAR likelihood, one matrix product per child. Bernoulli-logit
  // log-density written as y*eta - log(1 + exp(eta)), which equals
  // y*log_inv_logit(eta) + (1-y)*log1m_inv_logit(eta) exactly. (Written by
  // hand rather than with bernoulli_logit_lpmf, whose gradient code
  // triggers a GCC 14 internal compiler error; no threading/reduce_sum for
  // the same reason.)
  for (i in 1:N) {
    if (n_obs[i] == 0) continue;
    matrix[n_obs[i], V] eta =
      X1[row_start[i]:row_end[i]] * append_col(mu_person[i], A[i])';
    if (use_session == 1)
      eta += S_obs[row_start[i]:row_end[i]] * delta_person[i]';
    if (use_onset == 1)
      target += sum(W_obs[row_start[i]:row_end[i]] .*
                    (Y_obs[row_start[i]:row_end[i]] .* eta - log1p_exp(eta)));
    else
      target += sum(Y_obs[row_start[i]:row_end[i]] .* eta - log1p_exp(eta));
  }
}

generated quantities {
  // z_person[i] is the full Q-vector of person i's VAR parameters.
  // The Step 2 scripts read these draws to form the
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

# Lagged design for one child. Lag j is the value lag_step * j samples
# earlier; column <v>.prev is the value one sample earlier (used to define
# onsets). With session_col given, lags are built within each session, so no
# lag crosses from one recording into the next; the session is kept as a
# column.
shift_back <- function(x, d) {
  n <- length(x)
  if (d >= n) return(rep(NA_real_, n))
  c(rep(NA_real_, d), x[1:(n - d)])
}
make_lagged_comp <- function(df_person, vars, p = 1, min_T = 30,
                             lag_step = 1L, session_col = NULL) {
  df_person[vars] <- lapply(df_person[vars], as.numeric)
  if (!is.null(session_col)) {
    pieces <- lapply(split(df_person, df_person[[session_col]]), function(d)
      make_lagged_comp(d, vars, p = p, min_T = 1, lag_step = lag_step))
    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    if (length(pieces) == 0) return(NULL)
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    if (nrow(out) < min_T) return(NULL)
    return(out)
  }
  n_total <- nrow(df_person)
  if (n_total <= p * lag_step) return(NULL)
  lagged_list <- list()
  for (lag in seq_len(p))
    for (v in vars)
      lagged_list[[paste0(v, ".l", lag)]] <- shift_back(df_person[[v]], lag * lag_step)
  for (v in vars)
    lagged_list[[paste0(v, ".prev")]] <- shift_back(df_person[[v]], 1L)
  extra <- setdiff(names(df_person), vars)
  combined    <- cbind(df_person[vars],
                       as.data.frame(lagged_list, stringsAsFactors = FALSE),
                       df_person[extra])
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
  prior_sd_sigma_mu  = 2,
  ts_method          = c("subsample", "aggregate", "none"),
  T_subsample        = 500,
  window_size        = 10,
  grainsize          = 1,
  lag_step           = 1L,
  use_session        = FALSE,
  use_onset          = FALSE,
  session_col        = "session"
) {
  ts_method    <- match.arg(ts_method)
  if (ts_method == "aggregate" && (use_session || use_onset))
    stop("ts_method='aggregate' cannot be combined with the session or onset options")
  if (use_session && !session_col %in% names(combined_data))
    stop("use_session needs a '", session_col, "' column in combined_data")
  participants <- unique(combined_data[[id_col]])
  V            <- length(behaviours)
  K            <- p * V

  all_comps <- lapply(participants, function(pid) {
    dfp <- combined_data[combined_data[[id_col]] == pid, , drop = FALSE]
    keep <- c(behaviours, if (use_session) session_col)
    make_lagged_comp(dfp[, keep, drop = FALSE], behaviours, p = p, min_T = p + 1,
                     lag_step = lag_step, session_col = if (use_session) session_col)
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
  # Lagged predictors are 0/1 indicators. They are centred (so mu is the
  # log-odds at average lagged values) but NOT divided by their SD: for rare
  # behaviours (e.g. ~0.1% of samples) scaling would turn one occurrence into
  # a value of ~+30, making the normal(0, prior_sd_A) prior on A effectively
  # flat. Unscaled, A is the change in log-odds when the lagged behaviour is
  # present, and the prior regularises it on that scale.
  col_mn <- colMeans(X_all_flat, na.rm = TRUE)
  col_sd <- rep(1, K)

  Y_arr    <- array(0,  dim = c(N, T_max, V))
  X_arr    <- array(0,  dim = c(N, T_max, K))
  obs_mask <- array(0L, dim = c(N, T_max))
  S_arr    <- array(0,  dim = c(N, T_max))
  W_arr    <- array(1L, dim = c(N, T_max, V))

  # Session indicator: 0/1 (second session = 1), centred over all rows so
  # mu stays the log-odds at the average session.
  if (use_session) {
    sess_levels <- sort(unique(unlist(lapply(all_comps, function(cmp) cmp[[session_col]]))))
    if (length(sess_levels) != 2)
      stop("use_session expects exactly 2 sessions, found: ", paste(sess_levels, collapse = ", "))
    sess_all  <- unlist(lapply(all_comps, function(cmp) as.numeric(cmp[[session_col]] == sess_levels[2])))
    sess_mean <- mean(sess_all)
  }

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
    if (use_session)
      S_arr[i, seq_len(Ti)] <- as.numeric(comp[[session_col]] == sess_levels[2]) - sess_mean
    if (use_onset)
      for (v in seq_len(V))
        W_arr[i, seq_len(Ti), v] <- as.integer(as.numeric(comp[[paste0(behaviours[v], ".prev")]]) == 0)
  }
  if (use_onset)
    message("Onset model: outcome samples entering the likelihood per behaviour: ",
            paste(sprintf("%s %d", behaviours,
                          vapply(seq_len(V), function(v) sum(W_arr[, , v] * obs_mask), numeric(1))),
                  collapse = ", "))

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
      prior_sd_sigma_mu  = prior_sd_sigma_mu,
      grainsize          = as.integer(grainsize),
      use_session        = as.integer(use_session),
      session_x          = S_arr,
      use_onset          = as.integer(use_onset),
      y_include          = W_arr
    ),
    participants  = as.character(participants),
    behaviours    = behaviours,
    lag_col_names = lag_col_names,
    X_lag_center  = col_mn,
    X_lag_scale   = col_sd,
    V             = V,
    K             = K,
    Q             = V * K + V,
    options       = list(session = use_session, onset = use_onset,
                         lag_step = as.integer(lag_step), p = as.integer(p),
                         ts_method = ts_method, T_subsample = T_subsample,
                         session_levels = if (use_session) sess_levels else NULL)
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
  prior_sd_sigma_mu  = 2,
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
  lag_step           = 1L,
  use_session        = FALSE,
  use_onset          = FALSE,
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
    prior_sd_sigma_mu  = prior_sd_sigma_mu,
    ts_method          = ts_method,
    T_subsample        = T_subsample,
    window_size        = window_size,
    grainsize          = grainsize,
    lag_step           = lag_step,
    use_session        = use_session,
    use_onset          = use_onset
  )

  V <- prep$stan_data$V
  K <- prep$stan_data$K
  N <- prep$stan_data$N

  init_fn <- function() {
    init <- list(
      A_pop        = matrix(0,   nrow = V, ncol = K),
      sigma_A      = matrix(0.1, nrow = V, ncol = K),
      A_z          = array(0,    dim  = c(N, V, K)),
      mu_pop       = rep(0,   V),
      sigma_mu     = rep(0.1, V),
      mu_z         = matrix(0, nrow = N, ncol = V)
    )
    if (use_session) {
      init$delta_pop   <- rep(0,   V)
      init$sigma_delta <- rep(0.1, V)
      init$delta_z     <- matrix(0, nrow = N, ncol = V)
    }
    init
  }

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

if (opt_n_child > 0) parts <- head(parts, opt_n_child)   # quick tests only

for (i in seq_along(parts)) {
  thispart_data <- df_all |>
    filter(child_id == parts[i]) |>
    arrange(month, TIME) |>
    mutate(session = month) |>
    select(all_of(c(behaviour_cols, "session")))
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
  # Priors: normal(0, 5) for the population intercepts (mu_pop);
  # half-normal(0, 2) for children's variation around them (sigma_mu), which
  # stops children who never produce a rare behaviour from dragging the
  # intercept hierarchy towards -infinity; normal(0, 1) for the lag
  # coefficients (A_pop, on the 0/1 log-odds scale) and their between-child
  # SDs (sigma_A).
  prior_sd_A         = 1,
  prior_sd_intercept = 5,
  prior_sd_sigma_mu  = 2,
  ts_method          = "subsample",
  T_subsample 	     = opt_T_sub,
  chains             = 1,        # one chain per job; run 16 jobs
  iter_warmup        = opt_warmup,
  iter_sampling      = opt_sampling,
  lag_step           = opt_lag_step,
  use_session        = opt_session,
  use_onset          = opt_onset,
  cores              = 1,
  adapt_delta        = 0.95,
  refresh            = 1, 
  seed=as.integer(jobid),
)

# $draws is already extracted inside fit_step1(); saving here makes the RDS
# self-contained so the later stages do not need the CmdStan CSV files.
saveRDS(step1_result, file = paste0("dsem_step1_", jobid, ".rds"))
message("Step 1 complete. Saved to dsem_step1_", jobid, ".rds")
