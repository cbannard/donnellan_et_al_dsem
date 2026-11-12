# =============================================================================
# make_bridge.R
#
# Summarises the Step 1 posterior once, for the Step 2 scripts: for every
# child and every Step 1 parameter used in Step 2 (the z_person vector: lag
# coefficients and intercepts) it stores the posterior mean and SD. These are
# exactly the quantities select_terms_for_individual_behaviours.R and
# final_term_weighting.R compute from the draws themselves, so passing
# bridge=<file> to those scripts gives the same results while needing a few MB
# of memory instead of loading every Step 1 chain. Used by the permutation
# runs (permutation_null.R), which run Step 2 hundreds of times.
#
# Usage:
#   Rscript make_bridge.R <step1_folder> <output_file.rds>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript make_bridge.R <step1_folder> <output_file.rds>")
step1_folder <- args[1]
out_file     <- args[2]

suppressPackageStartupMessages(library(posterior))

files <- list.files(step1_folder, pattern = "^dsem_step1_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No dsem_step1_<n>.rds files in ", step1_folder)

prep <- NULL
z_list <- vector("list", length(files))
for (j in seq_along(files)) {
  message("Reading ", basename(files[j]))
  obj <- readRDS(files[j])
  if (is.null(prep)) prep <- obj$prep
  v <- variables(obj$draws)
  z_vars <- v[grepl("^z_person\\[", v)]
  z_list[[j]] <- unclass(as_draws_matrix(subset_draws(obj$draws, variable = z_vars)))
  rm(obj); invisible(gc())
}
z_mat <- do.call(rbind, z_list); rm(z_list)

N <- prep$stan_data$N; Q <- prep$Q
idx <- regmatches(colnames(z_mat), regexpr("\\[(\\d+),(\\d+)\\]", colnames(z_mat)))
i_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\1", idx))
q_idx <- as.integer(sub("\\[(\\d+),(\\d+)\\]", "\\2", idx))
stopifnot(length(i_idx) == N * Q)

z_hat <- matrix(0, N, Q); z_sd <- matrix(0, N, Q)
z_hat[cbind(i_idx, q_idx)] <- colMeans(z_mat)
z_sd[cbind(i_idx, q_idx)]  <- apply(z_mat, 2, sd)
z_sd <- pmax(z_sd, 1e-4)

# Keep the metadata Step 2 needs, not the (large) Step 1 data arrays.
prep_small <- prep
prep_small$stan_data <- list(N = N, V = prep$stan_data$V, K = prep$stan_data$K)

saveRDS(list(z_hat = z_hat, z_sd = z_sd, prep = prep_small,
             n_draws = nrow(z_mat), step1_files = basename(files)), out_file)
message("Saved bridge (", nrow(z_mat), " draws from ", length(files), " chains): ", out_file)
