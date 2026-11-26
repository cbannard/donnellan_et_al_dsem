# =============================================================================
# offshot_check.R
#
# Checks that the Step 2 results are not driven by off-shot periods or by the
# padding that brings every recording to the same length.
#
# In each data file (data_full/p<id>_<11m|12m>.csv) rows 1-9 are the nine
# behaviours, row 10 marks samples when the dyad was off-shot and row 11 marks
# padding; columns are 200 ms samples. Off-shot and padded samples have all
# nine behaviours coded 0, so in the model they are indistinguishable from
# samples in which no behaviour occurred.
#
# The script reports
#   1. how much of the data is off-shot and padding, overall and per child,
#      and each child's coded on-shot time (both sessions together);
#   2. for each Step 2 predictor checked, the correlation across children of
#      its posterior mean with the child's share of off-shot samples and with
#      the child's on-shot time;
#   3. the predictor's association with vocabulary, log(1 + 18-month score),
#      as the change per SD of the predictor in a linear regression, without
#      and with off-shot share and on-shot time as covariates.
# A lag-block predictor x -> y is the sum of the child's five lag
# coefficients of x in the equation for y (posterior means); mean.y is the
# child's intercept for y.
#
# Usage (from the repository folder; reads ./CDI_data.csv):
#   Rscript offshot_check.R <run_dir> [data=<run_dir>/data_full]
#           [blocks=contingent.contingent,CV.contingent,non_CV.contingent,mean.contingent]
#   Needs <run_dir>/PERMUTATIONS/bridge.rds (make_bridge.R).
#
# Outputs in <run_dir>/OFFSHOT/: offshot_children.csv (one row per child),
#   offshot_predictors.csv (one row per predictor), offshot_report.txt
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript offshot_check.R <run_dir> [data=<folder>] [blocks=x.y,...]")
run_dir <- args[1]
keyed <- args[-1]
keyed_arg <- function(key, default) {
  hit <- keyed[startsWith(keyed, paste0(key, "="))]
  if (length(hit)) sub(paste0("^", key, "="), "", hit[length(hit)]) else default
}
data_dir   <- keyed_arg("data", file.path(run_dir, "data_full"))
blocks_arg <- keyed_arg("blocks", "contingent.contingent,CV.contingent,non_CV.contingent,mean.contingent")
blocks     <- strsplit(blocks_arg, ",", fixed = TRUE)[[1]]
out_dir    <- file.path(run_dir, "OFFSHOT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_s <- 0.2   # seconds per sample

# ---- 1. Off-shot and padding per child ---------------------------------------
files <- list.files(data_dir, pattern = "^p[0-9]+_1[12]m\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No data files in ", data_dir)
per_file <- do.call(rbind, lapply(files, function(f) {
  m <- as.matrix(read.csv(f, check.names = FALSE))
  if (nrow(m) != 11) stop(basename(f), ": expected 11 rows (9 behaviours, off-shot, padding), found ", nrow(m))
  any_beh <- colSums(m[1:9, , drop = FALSE]) > 0
  data.frame(child = sub("^p([0-9]+)_.*$", "\\1", basename(f)),
             samples = ncol(m), offshot = sum(m[10, ] == 1), padding = sum(m[11, ] == 1),
             overlap = sum(m[10, ] == 1 & m[11, ] == 1),
             behaviour_in_offshot = sum(any_beh & m[10, ] == 1),
             behaviour_in_padding = sum(any_beh & m[11, ] == 1))
}))
kids <- aggregate(. ~ child, data = per_file, FUN = sum)
kids$offshot_share <- kids$offshot / kids$samples
kids$padding_share <- kids$padding / kids$samples
kids$onshot_min    <- (kids$samples - kids$offshot - kids$padding) * sample_s / 60

# ---- 2-3. Predictors --------------------------------------------------------
bridge_file <- file.path(run_dir, "PERMUTATIONS", "bridge.rds")
if (!file.exists(bridge_file)) stop("Missing ", bridge_file, " (run make_bridge.R)")
bridge <- readRDS(bridge_file)
prep <- bridge$prep
V <- prep$V; K <- prep$K; beh <- prep$behaviours; per_pair <- as.integer(K / V)

cols_of <- function(tok) {
  p <- strsplit(tok, ".", fixed = TRUE)[[1]]
  yi <- match(p[2], beh)
  if (length(p) != 2 || is.na(yi)) stop("Predictors must look like x.y or mean.y: ", tok)
  if (p[1] == "mean") return(V * K + yi)
  xi <- match(p[1], beh)
  if (is.na(xi)) stop("Unknown behaviour in ", tok)
  ((xi - 1L) * per_pair + seq_len(per_pair) - 1L) * V + yi
}
label_of <- function(tok) {
  p <- strsplit(tok, ".", fixed = TRUE)[[1]]
  if (p[1] == "mean") paste0("mean ", p[2]) else paste0(p[1], " -> ", p[2])
}

cdi <- read.csv("CDI_data.csv", stringsAsFactors = FALSE)
d <- data.frame(child = as.character(prep$participants), stringsAsFactors = FALSE)
d <- cbind(d, kids[match(d$child, kids$child), c("offshot_share", "onshot_min")])
d$log_vocab <- log1p(cdi$us18_imp[match(d$child, as.character(cdi$participant))])
if (anyNA(d$offshot_share)) stop("No data file for children: ", paste(d$child[is.na(d$offshot_share)], collapse = ", "))
ok <- !is.na(d$log_vocab)

res <- do.call(rbind, lapply(blocks, function(tok) {
  x  <- rowSums(bridge$z_hat[, cols_of(tok), drop = FALSE])
  xs <- as.numeric(scale(x[ok]))
  f0 <- lm(d$log_vocab[ok] ~ xs)
  f1 <- lm(d$log_vocab[ok] ~ xs + scale(d$offshot_share[ok]) + scale(d$onshot_min[ok]))
  data.frame(predictor = label_of(tok),
             r_offshot_share = cor(x, d$offshot_share),
             r_onshot_time   = cor(x, d$onshot_min),
             b_per_sd        = unname(coef(f0)["xs"]),
             b_per_sd_controlled = unname(coef(f1)["xs"]),
             n_vocab = sum(ok))
}))

write.csv(kids, file.path(out_dir, "offshot_children.csv"), row.names = FALSE)
write.csv(res,  file.path(out_dir, "offshot_predictors.csv"), row.names = FALSE)

pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt <- function(df) paste(capture.output(print(format(df, digits = 2), row.names = FALSE)), collapse = "\n")
report <- paste0(
  "Off-shot periods and padding (", nrow(kids), " children, ", length(files), " recordings; ",
  sample_s * 1000, " ms samples)\n\n",
  "Off-shot samples: ", pct(sum(kids$offshot) / sum(kids$samples)), " of all samples; per child median ",
  pct(median(kids$offshot_share)), ", range ", pct(min(kids$offshot_share)), "-", pct(max(kids$offshot_share)), "\n",
  "Padding samples:  ", pct(sum(kids$padding) / sum(kids$samples)), " of all samples\n",
  "Coded on-shot time per child (both sessions): ", sprintf("%.1f-%.1f min", min(kids$onshot_min), max(kids$onshot_min)), "\n",
  "Samples marked both off-shot and padding: ", sum(kids$overlap), "\n",
  "Samples with a behaviour coded during off-shot / padding: ",
  sum(kids$behaviour_in_offshot), " / ", sum(kids$behaviour_in_padding), "\n\n",
  "Predictors (posterior means across children): correlation with the share of off-shot samples\n",
  "and with on-shot time; change in log(1 + vocabulary) per SD of the predictor, without and\n",
  "with both as covariates (", res$n_vocab[1], " children with a vocabulary score)\n\n",
  fmt(res[, c("predictor", "r_offshot_share", "r_onshot_time", "b_per_sd", "b_per_sd_controlled")]), "\n\n",
  "Largest |r| with off-shot share: ", sprintf("%.2f", max(abs(res$r_offshot_share))),
  "; with on-shot time: ", sprintf("%.2f", max(abs(res$r_onshot_time))), "\n")
writeLines(report, file.path(out_dir, "offshot_report.txt"))
cat(report)
