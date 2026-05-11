##### TRAITEMENT DATA COUCHE DE MELANGE KERGUELEN  // VINCENT DORIOT // avr. 2023
### Dense GLORYS fPCA/PCA export for the reorganized pipeline.

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

fill_missing_dense <- function(mat) {
  col_means <- colMeans(mat, na.rm = TRUE)
  global_mean <- mean(mat, na.rm = TRUE)
  col_means[is.nan(col_means)] <- global_mean

  missing_idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(missing_idx) > 0) {
    mat[missing_idx] <- col_means[missing_idx[, 2]]
  }

  mat
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))
setwd(project_dir) ### REPERTOIRE CONTENANT LES DATA + FILES ...
input_dir <- file.path(project_dir, "processed", "1_gridded_data", "r_input")
output_dir <- file.path(project_dir, "processed", "2_fPCA", "GLORYS_masked_dense")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.table(file.path(input_dir, "GLORYS_masked.txt"), header = TRUE, na.strings = c("NA", ""))
rownames(data) <- seq_len(nrow(data))

#### ORGANISATION DES DATA
#### DATES PAR an.mois
time_table <- unique(data[, c("year", "mth")])
nmonth <- nrow(time_table)
tps <- paste(time_table$year, sprintf("%02d", time_table$mth), sep = ".")

first_time <- data$year == time_table$year[1] & data$mth == time_table$mth[1]
grid_size <- sum(first_time)
if (nrow(data) != grid_size * nmonth) {
  stop("GLORYS_masked.txt does not contain a complete regular grid for each month.")
}

#### ON ORGANISE LES DONNEES DANS UNE MATRICE 1521 LIGNES X (12X17) COLONNES
### SERIE TEMPORELLE PAR CELLULE
st_cell <- t(matrix(data$mld, nrow = grid_size, ncol = nmonth))
st_cell <- t(st_cell)

############# PARTIE I : TRAVAIL SUR LA SERIE TEMPORELLE TOTALE / CELLULE
### DETECTION DES SERIE TEMPS VIDES = CELLULES VIDES
ligNA <- rowSums(is.na(st_cell))
indNA <- which(ligNA == nmonth)

#### ON CREE DEUX MATRICES AVEC DES NAs : LES DATES ET LES DATAS
dates <- seq(min(time_table$year), max(time_table$year), length = nmonth + 1)
dates <- dates[-(nmonth + 1)]
datemat <- matrix(dates, grid_size, nmonth, byrow = TRUE)

#### CASE indNA COMPLETEMENT VIDE
donmat <- st_cell
if (length(indNA) > 0) {
  datemat <- datemat[-indNA, , drop = FALSE]
  donmat <- donmat[-indNA, , drop = FALSE]
}

### ACP DENSE EXACTE
###
### The previous version used fdapace with smoothed covariance estimation even
### for the full GLORYS field. That is useful for sparse profiles, but it
### intentionally smooths the dense monthly fields and produced a large RMSE in
### Figure 6. Here we keep the original output schema but use an exact dense
### PCA/SVD reconstruction instead.
datemat <- t(datemat)
donmat <- t(donmat)
donmat_filled <- fill_missing_dense(donmat)

mu <- rowMeans(donmat_filled)
centered <- sweep(donmat_filled, 1, mu, "-")
pca <- svd(centered)
tol <- max(dim(centered)) * max(pca$d) * .Machine$double.eps
nbcp <- sum(pca$d > tol)
if (nbcp < 1) {
  stop("No non-zero GLORYS PCA modes found.")
}

phi <- pca$u[, seq_len(nbcp), drop = FALSE]
xiEst <- sweep(pca$v[, seq_len(nbcp), drop = FALSE], 2, pca$d[seq_len(nbcp)], "*")
lambda <- (pca$d[seq_len(nbcp)]^2) / (ncol(donmat_filled) - 1)
Xestot <- sweep(phi %*% t(xiEst), 1, mu, "+")
workGrid <- dates

write.csv(workGrid, file.path(output_dir, "PCA_GRID.csv"), row.names = TRUE)
write.csv(lambda, file.path(output_dir, "PCA_LAMBDA.csv"), row.names = TRUE)
write.csv(Xestot, file.path(output_dir, "PCA_MLD.csv"), row.names = TRUE)
write.csv(phi, file.path(output_dir, "PCA_PHI.csv"), row.names = TRUE)
write.csv(xiEst, file.path(output_dir, "PCA_XIEST.csv"), row.names = TRUE)
write.csv(mu, file.path(output_dir, "PCA_MU.csv"), row.names = TRUE)
