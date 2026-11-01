# =============================================================================
# check_convergence.R
#
# Convergence diagnostics for the Step 1 fits (fitVAR.R output).
#
# Usage:
#   Rscript check_convergence.R <step1_folder> [output_folder]
#
#   step1_folder   folder holding dsem_step1_<n>.rds (one chain per file)
#   output_folder  default: <step1_folder>/../CONVERGENCE
#
# What it checks
#   1. Sampler diagnostics per chain: divergent transitions, iterations that
#      hit the maximum tree depth, E-BFMI (< 0.3 suggests poor exploration),
#      adapted step size.
#   2. R-hat and bulk/tail ESS for every model parameter (Vehtari et al.,
#      2021): R-hat < 1.01 and ESS > 400 are the usual targets.
#      z_person and A_pop_out are skipped: they are exact copies of
#      A/mu_person and A_pop.
#   3. Chain-by-chain comparison for the population-level parameters, to
#      find any single chain that disagrees with the rest (e.g. stuck at
#      sigma_mu ~ 0, as in the July run): each chain's mean lp__, how far its
#      means sit from the other chains', and R-hat recomputed without it.
#   4. Trace and rank plots for lp__, mu_pop, sigma_mu and the 12 worst
#      parameters by R-hat.
#
# Outputs (in output_folder)
#   convergence_report.txt              plain-text summary (also printed)
#   convergence_all_params.csv          rhat / ess for every parameter
#   convergence_population_level.csv    fuller summary for population-level
#   convergence_by_family.csv           counts per parameter family
#   convergence_by_chain.csv            sampler + chain-agreement diagnostics
#   convergence_plots.pdf               trace + rank plots
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript check_convergence.R <step1_folder> [output_folder]")
step1_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir   <- if (length(args) >= 2) args[2] else file.path(dirname(step1_dir), "CONVERGENCE")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(posterior))
have_bayesplot <- requireNamespace("bayesplot", quietly = TRUE) &&
                  requireNamespace("ggplot2", quietly = TRUE)

RHAT_OK <- 1.01
ESS_OK  <- 400

report_lines <- character(0)
say <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  cat(line, "\n", sep = "")
}

# -----------------------------------------------------------------------------
# 1. Read chains: keep draws (minus duplicated quantities), sampler diagnostics
#    and a little metadata; drop everything else straight away.
# -----------------------------------------------------------------------------
files <- list.files(step1_dir, pattern = "^dsem_step1_[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No dsem_step1_<n>.rds files in ", step1_dir)
chain_ids <- as.integer(sub("^dsem_step1_([0-9]+)\\.rds$", "\\1", basename(files)))
ord <- order(chain_ids); files <- files[ord]; chain_ids <- chain_ids[ord]
n_chains <- length(files)

draws_list <- vector("list", n_chains)
sampler_rows <- vector("list", n_chains)
keep_vars <- NULL

for (j in seq_len(n_chains)) {
  message("Reading chain file ", j, "/", n_chains, ": ", basename(files[j]))
  obj <- readRDS(files[j])
  d   <- obj$draws
  if (is.null(keep_vars)) {
    v <- variables(d)
    keep_vars <- v[!grepl("^(z_person|A_pop_out)\\[", v)]
  }
  draws_list[[j]] <- subset_draws(d, variable = keep_vars)

  sd_j <- tryCatch(obj$fit$sampler_diagnostics(format = "df"), error = function(e) NULL)
  md   <- tryCatch(obj$fit$metadata(), error = function(e) list())
  max_td <- if (!is.null(md$max_treedepth)) md$max_treedepth else 10
  if (!is.null(sd_j)) {
    e <- sd_j$energy__
    sampler_rows[[j]] <- data.frame(
      chain             = chain_ids[j],
      seed              = if (!is.null(md$seed)) md$seed else NA,
      iterations        = nrow(sd_j),
      divergences       = sum(sd_j$divergent__),
      max_treedepth_hit = sum(sd_j$treedepth__ >= max_td),
      ebfmi             = sum(diff(e)^2) / sum((e - mean(e))^2),
      step_size         = mean(sd_j$stepsize__),
      mean_accept       = mean(sd_j$accept_stat__),
      mean_leapfrog     = mean(sd_j$n_leapfrog__)
    )
  } else {
    sampler_rows[[j]] <- data.frame(chain = chain_ids[j], seed = NA, iterations = niterations(d),
      divergences = NA, max_treedepth_hit = NA, ebfmi = NA, step_size = NA,
      mean_accept = NA, mean_leapfrog = NA)
  }
  rm(obj, d, sd_j); invisible(gc())
}

iters <- vapply(draws_list, niterations, integer(1))
if (length(unique(iters)) > 1)
  stop("Chains have different numbers of draws: ", paste(unique(iters), collapse = ", "))

draws <- do.call(bind_draws, c(draws_list, along = "chain"))
rm(draws_list); invisible(gc())

say("Step 1 convergence check  —  ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Folder: ", step1_dir)
say("Chains: ", n_chains, " (files ", paste(chain_ids, collapse = ", "), ")")
say("Draws per chain: ", niterations(draws), "   total draws: ", ndraws(draws))
say("Parameters checked: ", nvariables(draws), " (z_person and A_pop_out skipped as duplicates)")
say("")

# -----------------------------------------------------------------------------
# 2. Sampler diagnostics
# -----------------------------------------------------------------------------
sampler_df <- do.call(rbind, sampler_rows)
say("== Sampler diagnostics ==")
if (all(is.na(sampler_df$divergences))) {
  say("Not available (sampler diagnostics were not saved in these files).")
} else {
  say("Divergent transitions: ", sum(sampler_df$divergences),
      " in total; chains with any: ", sum(sampler_df$divergences > 0))
  say("Iterations at max tree depth: ", sum(sampler_df$max_treedepth_hit),
      " in total; chains with any: ", sum(sampler_df$max_treedepth_hit > 0))
  say("E-BFMI below 0.3: ", sum(sampler_df$ebfmi < 0.3), " chain(s)",
      if (any(sampler_df$ebfmi < 0.3))
        paste0(" (", paste(sampler_df$chain[sampler_df$ebfmi < 0.3], collapse = ", "), ")") else "")
  say(sprintf("Step size: median %.4f, range %.4f - %.4f",
              median(sampler_df$step_size), min(sampler_df$step_size), max(sampler_df$step_size)))
}
say("")

# -----------------------------------------------------------------------------
# 3. R-hat / ESS for every parameter
# -----------------------------------------------------------------------------
message("Computing R-hat and ESS for ", nvariables(draws), " parameters (this takes a while)...")
summ <- summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
summ$family <- sub("\\[.*$", "", summ$variable)
write.csv(summ, file.path(out_dir, "convergence_all_params.csv"), row.names = FALSE)

fam <- do.call(rbind, lapply(split(summ, summ$family), function(g) data.frame(
  family          = g$family[1],
  n               = nrow(g),
  rhat_median     = median(g$rhat, na.rm = TRUE),
  rhat_max        = max(g$rhat, na.rm = TRUE),
  n_rhat_gt_1.01  = sum(g$rhat > 1.01, na.rm = TRUE),
  n_rhat_gt_1.05  = sum(g$rhat > 1.05, na.rm = TRUE),
  n_rhat_gt_1.10  = sum(g$rhat > 1.10, na.rm = TRUE),
  ess_bulk_min    = min(g$ess_bulk, na.rm = TRUE),
  ess_tail_min    = min(g$ess_tail, na.rm = TRUE),
  n_ess_lt_400    = sum(g$ess_bulk < ESS_OK | g$ess_tail < ESS_OK, na.rm = TRUE)
)))
fam_order <- c("lp__", "mu_pop", "sigma_mu", "A_pop", "sigma_A", "mu_z", "mu_person", "A_z", "A")
fam <- fam[order(match(fam$family, fam_order, nomatch = 99)), ]
write.csv(fam, file.path(out_dir, "convergence_by_family.csv"), row.names = FALSE)

say("== R-hat and ESS (targets: R-hat < ", RHAT_OK, ", bulk and tail ESS > ", ESS_OK, ") ==")
say(sprintf("%-10s %6s %8s %8s %8s %8s %8s %9s %9s %9s",
            "family", "n", "Rhat_med", "Rhat_max", ">1.01", ">1.05", ">1.10",
            "ESSb_min", "ESSt_min", "ESS<400"))
for (r in seq_len(nrow(fam))) with(fam[r, ], say(sprintf(
  "%-10s %6d %8.3f %8.3f %8d %8d %8d %9.0f %9.0f %9d",
  family, n, rhat_median, rhat_max, n_rhat_gt_1.01, n_rhat_gt_1.05, n_rhat_gt_1.10,
  ess_bulk_min, ess_tail_min, n_ess_lt_400)))
n_bad <- sum(summ$rhat > RHAT_OK | summ$ess_bulk < ESS_OK | summ$ess_tail < ESS_OK, na.rm = TRUE)
say("")
say("Parameters failing either target: ", n_bad, " of ", nrow(summ),
    sprintf(" (%.1f%%)", 100 * n_bad / nrow(summ)))
worst <- head(summ[order(-summ$rhat), c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")], 15)
say("Worst 15 by R-hat:")
for (r in seq_len(nrow(worst))) with(worst[r, ], say(sprintf(
  "  %-22s mean %9.3f  sd %8.3f  Rhat %.3f  ESS bulk %6.0f  tail %6.0f",
  variable, mean, sd, rhat, ess_bulk, ess_tail)))
say("")

# Fuller summary of the population-level parameters
pop_vars <- grep("^(lp__|mu_pop|sigma_mu|A_pop|sigma_A)(\\[|$)", variables(draws), value = TRUE)
summ_pop <- summarise_draws(subset_draws(draws, variable = pop_vars),
                            "mean", "median", "sd", ~quantile2(.x, c(0.05, 0.95)),
                            "rhat", "ess_bulk", "ess_tail")
write.csv(summ_pop, file.path(out_dir, "convergence_population_level.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Does any single chain disagree with the rest?
#    For lp__, mu_pop[] and sigma_mu[]:
#      - dev: |chain mean - median of all chains' means| / pooled
#        within-chain SD. The median is not pulled about by one odd chain,
#        so only that chain gets a large value. Above ~1 is notable.
#      - R-hat recomputed with that chain left out: if R-hat for the key
#        parameters drops to ~1.00 when one chain is removed, that chain is
#        the problem.
# -----------------------------------------------------------------------------
key_vars <- grep("^(lp__|mu_pop|sigma_mu)(\\[|$)", variables(draws), value = TRUE)
key <- as_draws_array(subset_draws(draws, variable = key_vars))
cm  <- apply(key, c(2, 3), mean)            # chains x variables
csd <- apply(key, c(2, 3), sd)
pooled_sd <- sqrt(colMeans(csd^2))
med_cm <- apply(cm, 2, median)
dev <- sapply(seq_len(n_chains), function(c)
  abs(cm[c, ] - med_cm) / pooled_sd)                              # variables x chains
rhat_key_all <- apply(key, 3, rhat)
rhat_drop <- if (n_chains > 2) sapply(seq_len(n_chains), function(c)
  max(apply(key[, -c, , drop = FALSE], 3, rhat))) else rep(NA_real_, n_chains)

chain_df <- cbind(sampler_df,
  mean_lp            = cm[, "lp__"],
  lp_dev             = dev["lp__", ],
  max_dev_mu_sigmamu = apply(dev[rownames(dev) != "lp__", , drop = FALSE], 2, max),
  n_dev_gt_1         = colSums(dev[rownames(dev) != "lp__", , drop = FALSE] > 1),
  max_rhat_key_without_chain = rhat_drop)
write.csv(chain_df, file.path(out_dir, "convergence_by_chain.csv"), row.names = FALSE)

say("== Chain-by-chain comparison (lp__, mu_pop, sigma_mu) ==")
say(sprintf("Max R-hat over these %d parameters with all chains: %.3f", length(key_vars), max(rhat_key_all)))
say("dev = distance of a chain's mean from the median chain mean, in within-chain SDs (> 1 is notable)")
say("maxRhat_w/o = max R-hat over these parameters with that chain left out")
say(sprintf("%6s %10s %7s %7s %10s %8s %6s %6s %7s %13s",
            "chain", "mean_lp", "lp_dev", "max_dev", "n_dev>1", "diverg", "maxTD", "EBFMI", "step", "maxRhat_w/o"))
for (r in seq_len(nrow(chain_df))) with(chain_df[r, ], say(sprintf(
  "%6d %10.1f %7.2f %7.2f %10d %8s %6s %6.2f %7.4f %13.3f",
  chain, mean_lp, lp_dev, max_dev_mu_sigmamu, n_dev_gt_1,
  ifelse(is.na(divergences), "NA", divergences), ifelse(is.na(max_treedepth_hit), "NA", max_treedepth_hit),
  ebfmi, step_size, max_rhat_key_without_chain)))
suspect <- chain_df$chain[chain_df$lp_dev > 1 | chain_df$max_dev_mu_sigmamu > 2]
if (n_chains > 2) {
  best <- which.min(rhat_drop)
  say(sprintf("Leaving out chain %d gives the lowest max R-hat (%.3f, vs %.3f with all chains).",
              chain_df$chain[best], rhat_drop[best], max(rhat_key_all)))
}
if (length(suspect) > 0) {
  say("Chains that look out of line with the others: ", paste(suspect, collapse = ", "),
      " — inspect their traces in convergence_plots.pdf")
} else {
  say("No single chain stands out from the others.")
}
say("")

# -----------------------------------------------------------------------------
# 5. Overall verdict (rule of thumb)
# -----------------------------------------------------------------------------
pop_fams <- fam[fam$family %in% c("lp__", "mu_pop", "sigma_mu", "A_pop", "sigma_A"), ]
pop_ok   <- all(pop_fams$rhat_max < RHAT_OK & pop_fams$n_ess_lt_400 == 0)
div_ok   <- isTRUE(sum(sampler_df$divergences, na.rm = TRUE) == 0)
say("== Summary ==")
say("Population-level parameters meet R-hat/ESS targets: ", ifelse(pop_ok, "YES", "NO"))
say("No divergent transitions: ", ifelse(div_ok, "YES", "NO"))
say("All parameters meet R-hat/ESS targets: ", ifelse(n_bad == 0, "YES", paste0("NO (", n_bad, " fail)")))

writeLines(report_lines, file.path(out_dir, "convergence_report.txt"))

# -----------------------------------------------------------------------------
# 6. Plots
# -----------------------------------------------------------------------------
if (have_bayesplot) {
  message("Writing plots...")
  worst12 <- head(summ$variable[order(-summ$rhat)], 12)
  groups <- list(
    "lp__ and mu_pop" = c("lp__", grep("^mu_pop\\[", variables(draws), value = TRUE)),
    "sigma_mu"        = grep("^sigma_mu\\[", variables(draws), value = TRUE),
    "Worst 12 by R-hat" = worst12)
  pdf(file.path(out_dir, "convergence_plots.pdf"), width = 11, height = 8.5)
  for (nm in names(groups)) {
    p <- groups[[nm]]
    print(bayesplot::mcmc_trace(draws, pars = p) + ggplot2::ggtitle(paste("Trace:", nm)))
    print(bayesplot::mcmc_rank_overlay(draws, pars = p) + ggplot2::ggtitle(paste("Rank plot:", nm)))
  }
  invisible(dev.off())
} else {
  message("bayesplot/ggplot2 not installed: skipping plots.")
}

message("Done. Outputs in ", out_dir)
