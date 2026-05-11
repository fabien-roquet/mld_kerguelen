#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: setup_r_packages.R <project_dir> <repos> <package> [<package> ...]", call. = FALSE)
}

project_dir <- normalizePath(args[[1]], mustWork = TRUE)
repos <- args[[2]]
packages <- unique(args[-c(1, 2)])
r_lib <- Sys.getenv("MLD_R_LIB", file.path(project_dir, ".r-lib"))

dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(r_lib, .libPaths())))
options(repos = c(CRAN = repos))

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing R packages into ", r_lib, ": ", paste(missing, collapse = ", "))
  install.packages(missing, lib = r_lib, repos = repos, dependencies = TRUE)
} else {
  message("R packages already available: ", paste(packages, collapse = ", "))
}

still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Missing R packages after installation: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

message("R library path: ", r_lib)
