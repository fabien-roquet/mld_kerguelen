##### Random GLORYS subsampling experiment for PACE sensitivity tests.
###
### The script draws reproducible random space-time samples from the full
### GLORYS anomaly product, runs sparse PACE on each sample, and stores compact
### diagnostics used by appendix Figures A1-A2.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_arg <- args[startsWith(args, file_arg)]
  if (length(script_arg) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_arg[1]))))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }

  getwd()
}

get_arg <- function(args, flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

has_flag <- function(args, flag) {
  any(args == flag)
}

linear_trend <- function(y, x, min_points = 8) {
  valid <- is.finite(y) & is.finite(x)
  if (sum(valid) < min_points) {
    return(c(slope = NA_real_, pvalue = NA_real_))
  }
  fit <- lm(y[valid] ~ x[valid])
  c(
    slope = unname(coef(fit)[2]),
    pvalue = unname(coef(summary(fit))[2, 4])
  )
}

sample_space_time <- function(donmat, fraction, seed) {
  if (fraction >= 1) {
    return(donmat)
  }

  sampled <- matrix(NA_real_, nrow = nrow(donmat), ncol = ncol(donmat))
  available <- which(!is.na(donmat))
  set.seed(seed)
  n_keep <- max(1, round(length(available) * fraction))
  keep <- sample(available, size = n_keep, replace = FALSE)
  sampled[keep] <- donmat[keep]
  sampled
}

run_sparse_pace <- function(donmat, datemat, nmonth) {
  subject_has_data <- rowSums(!is.na(donmat)) > 0
  if (sum(subject_has_data) < 3) {
    stop("Too few spatial cells contain sampled observations for PACE.")
  }

  don <- t(donmat[subject_has_data, , drop = FALSE])
  dat <- t(datemat[subject_has_data, , drop = FALSE])
  max_k <- min(100, nmonth - 2, ncol(don) - 1)

  res <- FPCA(
    data.frame(don),
    data.frame(dat),
    optns = list(
      useBinnedData = "OFF",
      methodMuCovEst = "smooth",
      userBwCov = 0.25,
      userBwMu = 0.05,
      kernel = "epan",
      nRegGrid = nmonth,
      dataType = "Sparse",
      methodSelectK = "FVE",
      FVEthreshold = 1,
      methodXi = "CE",
      maxK = max_k,
      error = FALSE,
      plot = FALSE
    )
  )

  nbcp <- length(res$lambda)
  xest <- sweep(
    res$phi[, seq_len(nbcp), drop = FALSE] %*% t(res$xiEst[, seq_len(nbcp), drop = FALSE]),
    1,
    res$mu,
    "+"
  )

  full_xest <- matrix(NA_real_, nrow = nmonth, ncol = nrow(donmat))
  full_xest[, subject_has_data] <- xest
  full_xest
}

summarise_reconstruction <- function(xest, grid, time_table, percentage, replicate, seed, n_observations, source) {
  years <- sort(unique(time_table$year))
  annual_groups <- split(seq_len(nrow(time_table)), time_table$year)
  annual_by_cell <- sapply(annual_groups, function(idx) {
    colMeans(xest[idx, , drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(annual_by_cell))) {
    annual_by_cell <- matrix(annual_by_cell, ncol = 1)
  }
  annual_by_cell <- t(annual_by_cell)

  annual_global <- rowMeans(annual_by_cell, na.rm = TRUE)
  trend <- linear_trend(annual_global, years)
  map_slopes <- apply(annual_by_cell, 2, function(y) linear_trend(y, years)[["slope"]])

  series <- data.frame(
    percentage = percentage,
    replicate = replicate,
    seed = seed,
    year = years,
    mld = annual_global,
    source = source
  )
  global_trend <- data.frame(
    percentage = percentage,
    replicate = replicate,
    seed = seed,
    n_observations = n_observations,
    slope = trend[["slope"]],
    pvalue = trend[["pvalue"]],
    source = source
  )
  trend_map <- data.frame(
    percentage = percentage,
    replicate = replicate,
    seed = seed,
    long = grid$long,
    lat = grid$lat,
    slope = map_slopes,
    source = source
  )

  list(series = series, global_trend = global_trend, trend_map = trend_map)
}

read_dense_pace_or_full <- function(project_dir, nmonth, ncell, donmat) {
  dense_file <- file.path(project_dir, "processed", "2_fPCA", "GLORYS_masked_dense", "PCA_MLD.csv")
  if (file.exists(dense_file)) {
    dense <- as.matrix(read.csv(dense_file, row.names = 1, check.names = FALSE))
    if (nrow(dense) == nmonth && ncol(dense) == ncell) {
      return(list(xest = dense, source = "dense_pace"))
    }
    warning("Existing dense PACE reconstruction has unexpected dimensions; using full GLORYS anomalies for 100%.")
  } else {
    warning("Dense GLORYS PACE reconstruction is missing; using full GLORYS anomalies for 100%.")
  }
  list(xest = t(donmat), source = "full_glorys_anomaly")
}

write_replicate_outputs <- function(out_dir, tag, outputs) {
  write.csv(outputs$series, file.path(out_dir, paste0(tag, "_annual_series.csv")), row.names = FALSE)
  write.csv(outputs$global_trend, file.path(out_dir, paste0(tag, "_global_trend.csv")), row.names = FALSE)
  write.csv(outputs$trend_map, file.path(out_dir, paste0(tag, "_annual_trend_map.csv")), row.names = FALSE)
}

replicate_outputs_are_current <- function(expected_files, seed) {
  if (!all(file.exists(expected_files))) {
    return(FALSE)
  }
  trend_file <- expected_files[grepl("_global_trend[.]csv$", expected_files)]
  if (length(trend_file) != 1) {
    return(FALSE)
  }
  trend <- read.csv(trend_file)
  nrow(trend) == 1 && "seed" %in% names(trend) && trend$seed[1] == seed
}

aggregate_outputs <- function(replicate_dir, output_dir, percentages, replicates, base_seed) {
  read_group <- function(pattern) {
    files <- list.files(replicate_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      stop(paste("No replicate outputs found for", pattern))
    }
    do.call(rbind, lapply(files, read.csv))
  }

  annual_series <- read_group("_annual_series[.]csv$")
  global_trends <- read_group("_global_trend[.]csv$")
  trend_maps <- read_group("_annual_trend_map[.]csv$")

  requested <- expand.grid(percentage = percentages, replicate = seq_len(replicates))
  requested$seed <- base_seed + requested$percentage * 1000 + requested$replicate
  key_columns <- c("percentage", "replicate", "seed")

  annual_series <- merge(requested, annual_series, by = key_columns)
  global_trends <- merge(requested, global_trends, by = key_columns)
  trend_maps <- merge(requested, trend_maps, by = key_columns)

  if (nrow(global_trends) == 0) {
    stop("No current replicate outputs match the requested percentages, replicates, and seed.")
  }
  available_runs <- unique(global_trends[, key_columns])
  if (nrow(available_runs) < nrow(requested)) {
    warning("Some requested replicate outputs are missing from the aggregate diagnostics.")
  }

  trend_summary <- aggregate(
    cbind(slope, pvalue, n_observations) ~ percentage,
    data = global_trends,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), std = sd(x, na.rm = TRUE))
  )
  trend_summary <- do.call(data.frame, trend_summary)
  names(trend_summary) <- c(
    "percentage",
    "slope_mean",
    "slope_std",
    "pvalue_mean",
    "pvalue_std",
    "n_observations_mean",
    "n_observations_std"
  )

  trend_summary$n_replicates <- as.integer(aggregate(
    replicate ~ percentage,
    data = global_trends,
    FUN = function(x) length(unique(x))
  )$replicate)

  write.csv(annual_series, file.path(output_dir, "global_annual_series.csv"), row.names = FALSE)
  write.csv(global_trends, file.path(output_dir, "global_trends.csv"), row.names = FALSE)
  write.csv(trend_maps, file.path(output_dir, "annual_trend_maps.csv"), row.names = FALSE)
  write.csv(trend_summary, file.path(output_dir, "global_trend_summary.csv"), row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
script_dir <- get_script_dir()
project_dir <- normalizePath(get_arg(args, "--project-root", file.path(script_dir, "..")), mustWork = TRUE)
replicates <- as.integer(get_arg(args, "--replicates", "30"))
base_seed <- as.integer(get_arg(args, "--seed", "20260526"))
levels_text <- get_arg(args, "--levels", "5,10,20")
percentages <- as.integer(strsplit(levels_text, ",")[[1]])
force <- has_flag(args, "--force")

if (any(is.na(percentages)) || any(percentages <= 0) || any(percentages > 100)) {
  stop("--levels must contain comma-separated percentages in the interval 1..100.")
}
if (is.na(replicates) || replicates < 1) {
  stop("--replicates must be a positive integer.")
}

suppressPackageStartupMessages(library(fdapace))

input_file <- file.path(project_dir, "processed", "1_gridded_data", "r_input", "GLORYS_masked.txt")
output_dir <- file.path(project_dir, "processed", "2_fPCA", "GLORYS_random_sampling")
replicate_dir <- file.path(output_dir, "replicates")
dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.table(input_file, header = TRUE, na.strings = c("NA", ""))
rownames(data) <- seq_len(nrow(data))

time_table <- unique(data[, c("year", "mth")])
nmonth <- nrow(time_table)
first_time <- data$year == time_table$year[1] & data$mth == time_table$mth[1]
grid <- data[first_time, c("long", "lat")]
grid_size <- nrow(grid)

if (nrow(data) != grid_size * nmonth) {
  stop("GLORYS_masked.txt does not contain a complete regular grid for each month.")
}

st_cell <- matrix(data$mld, nrow = grid_size, ncol = nmonth)
date_values <- time_table$year + (time_table$mth - 1) / 12
datemat <- matrix(date_values, nrow = grid_size, ncol = nmonth, byrow = TRUE)

valid_cells <- rowSums(is.na(st_cell)) < nmonth
donmat <- st_cell[valid_cells, , drop = FALSE]
datemat <- datemat[valid_cells, , drop = FALSE]
grid <- grid[valid_cells, , drop = FALSE]
available_count <- sum(!is.na(donmat))
full_reconstruction <- NULL

for (percentage in percentages) {
  fraction <- percentage / 100
  for (replicate in seq_len(replicates)) {
    seed <- base_seed + percentage * 1000 + replicate
    tag <- sprintf("pct_%03d_rep_%02d", percentage, replicate)
    expected_files <- file.path(
      replicate_dir,
      paste0(tag, c("_annual_series.csv", "_global_trend.csv", "_annual_trend_map.csv"))
    )

    if (!force && replicate_outputs_are_current(expected_files, seed)) {
      message("Skipping existing ", tag)
      next
    }

    message("Running GLORYS random sampling PACE: ", percentage, "% replicate ", replicate)
    if (percentage == 100) {
      if (is.null(full_reconstruction)) {
        full_reconstruction <- read_dense_pace_or_full(project_dir, nmonth, nrow(donmat), donmat)
      }
      xest <- full_reconstruction$xest
      source <- full_reconstruction$source
      n_observations <- available_count
    } else {
      sampled <- sample_space_time(donmat, fraction, seed)
      xest <- run_sparse_pace(sampled, datemat, nmonth)
      source <- "sparse_pace"
      n_observations <- sum(!is.na(sampled))
    }

    outputs <- summarise_reconstruction(
      xest = xest,
      grid = grid,
      time_table = time_table,
      percentage = percentage,
      replicate = replicate,
      seed = seed,
      n_observations = n_observations,
      source = source
    )
    write_replicate_outputs(replicate_dir, tag, outputs)
  }
}

aggregate_outputs(replicate_dir, output_dir, percentages, replicates, base_seed)
message("Wrote GLORYS random sampling diagnostics to ", output_dir)
