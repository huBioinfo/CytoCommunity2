suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(CCA)
  library(CCP)
})

# Config
LONG_CSV <- "../../data/Communication/config/EnrichScoreMatrix_long.csv"
OUT_DIR  <- "../../data/Communication/config"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NBOOT <- 999
MIN_COMMON_SAMPLES <- 3
SD_EPS <- 1e-8

# 1) Read long table
all_df <- read.csv(LONG_CSV, stringsAsFactors = FALSE)
all_df <- all_df %>%
  transmute(
    Sample   = as.character(Sample),
    CN       = as.character(CN),
    CellType = as.character(CellType),
    Score    = as.numeric(Score)
  )

# 2) Build Sample × CellType matrix for each CN
CN_list <- sort(unique(all_df$CN))
ScoreList <- lapply(CN_list, function(cn) {
  mat <- all_df %>%
    dplyr::filter(CN == cn) %>%
    dplyr::select(Sample, CellType, Score) %>%
    tidyr::pivot_wider(names_from = CellType, values_from = Score) %>%
    as.data.frame()
  
  rownames(mat) <- mat$Sample
  mat <- mat[, colnames(mat) != "Sample", drop = FALSE]
  mat[] <- lapply(mat, as.numeric)
  mat
})
names(ScoreList) <- CN_list

# 3) Data cleaning
clean_data <- function(X, Y) {
  common <- intersect(rownames(X), rownames(Y))
  if (length(common) < MIN_COMMON_SAMPLES) return(NULL)
  
  X <- X[common, , drop = FALSE]
  Y <- Y[common, , drop = FALSE]
  
  X <- X[, sapply(X, function(z) sd(z, na.rm = TRUE) > SD_EPS), drop = FALSE]
  Y <- Y[, sapply(Y, function(z) sd(z, na.rm = TRUE) > SD_EPS), drop = FALSE]
  
  Xs <- scale(as.matrix(X))
  Ys <- scale(as.matrix(Y))
  
  Xs <- Xs[, colSums(is.na(Xs)) == 0, drop = FALSE]
  Ys <- Ys[, colSums(is.na(Ys)) == 0, drop = FALSE]
  
  if (ncol(Xs) == 0 || ncol(Ys) == 0) return(NULL)
  if (nrow(Xs) < MIN_COMMON_SAMPLES) return(NULL)
  
  list(Xs = Xs, Ys = Ys)
}

# 4) Run CCA for all CN pairs
message("Run CCA for all CN pairs...")

CN_pairs <- combn(CN_list, 2, simplify = FALSE)
rho_results <- list()
cc_cache <- list()

for (pair in CN_pairs) {
  cnA <- pair[1]
  cnB <- pair[2]
  cleaned <- clean_data(ScoreList[[cnA]], ScoreList[[cnB]])
  if (is.null(cleaned)) next
  
  res <- tryCatch({
    cc_res <- CCA::cc(cleaned$Xs, cleaned$Ys)
    perm   <- CCP::p.perm(cleaned$Xs, cleaned$Ys, nboot = NBOOT)
    list(
      summary = data.frame(
        CN_A = cnA,
        CN_B = cnB,
        rho1 = cc_res$cor[1],
        pval = perm$p.value[1],
        n_samples = nrow(cleaned$Xs)
      ),
      cc_res = cc_res
    )
  }, error = function(e) NULL)
  
  if (!is.null(res)) {
    key <- paste(cnA, cnB, sep = "_vs_")
    rho_results[[key]] <- res$summary
    cc_cache[[key]] <- res$cc_res
  }
}

df_rho <- bind_rows(rho_results) %>% arrange(desc(abs(rho1)))

write.csv(df_rho, file.path(OUT_DIR, "CCA_config.csv"), row.names = FALSE)
print(df_rho[, c("CN_A", "CN_B", "rho1", "pval", "n_samples")])

# 5) Write coefficients for first canonical variate
for (i in seq_len(nrow(df_rho))) {
  cnA <- as.character(df_rho$CN_A[i])
  cnB <- as.character(df_rho$CN_B[i])
  key <- paste(cnA, cnB, sep = "_vs_")
  cc_res <- cc_cache[[key]]
  
  coefX <- cc_res$xcoef[, 1]
  coefY <- cc_res$ycoef[, 1]
  all_ct <- union(names(coefX), names(coefY))
  
  df_coef <- data.frame(
    CellType  = all_ct,
    CN_A_coef = unname(coefX[all_ct]),
    CN_B_coef = unname(coefY[all_ct]),
    rho1      = df_rho$rho1[i],
    pval      = df_rho$pval[i]
  )
  
  write.csv(
    df_coef,
    file.path(OUT_DIR, paste0("CCA_coefficients_CN", cnA, "_vs_CN", cnB, ".csv")),
    row.names = FALSE
  )
}

