options(stringsAsFactors = FALSE)

base_dir <- "D:/SX/TMFI_CPI_MR/01_stage1_instruments"

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  x[which(cumsum(w) >= 0.5)[1]]
}

run_mr <- function(path, outcome_name) {
  dat <- read.csv(path, check.names = FALSE)
  dat <- dat[dat$HARMONISE_STATUS %in% c("aligned", "swapped"), ]
  bx <- dat$BETA
  sx <- dat$SE
  by <- dat$BETA_OUT_HARMONISED
  sy <- dat$SE_OUT
  w <- 1 / sy^2
  k <- length(bx)

  ivw_beta <- sum(w * bx * by) / sum(w * bx^2)
  ivw_se_fixed <- sqrt(1 / sum(w * bx^2))
  q <- sum(w * (by - ivw_beta * bx)^2)
  q_df <- k - 1
  q_p <- pchisq(q, df = q_df, lower.tail = FALSE)
  phi <- max(1, q / q_df)
  ivw_se_mre <- ivw_se_fixed * sqrt(phi)
  ivw_p <- 2 * pnorm(abs(ivw_beta / ivw_se_mre), lower.tail = FALSE)

  egger <- lm(by ~ bx, weights = w)
  eg <- summary(egger)$coefficients
  egger_intercept <- eg[1, 1]
  egger_intercept_se <- eg[1, 2]
  egger_intercept_p <- eg[1, 4]
  egger_beta <- eg[2, 1]
  egger_se <- eg[2, 2]
  egger_p <- eg[2, 4]

  ratio <- by / bx
  ratio_w <- bx^2 / sy^2
  wm_beta <- weighted_median(ratio, ratio_w)
  set.seed(20260918)
  boot <- replicate(5000, {
    bx_b <- rnorm(k, bx, sx)
    by_b <- rnorm(k, by, sy)
    weighted_median(by_b / bx_b, bx_b^2 / sy^2)
  })
  wm_se <- sd(boot, na.rm = TRUE)
  wm_p <- 2 * pnorm(abs(wm_beta / wm_se), lower.tail = FALSE)

  loo <- sapply(seq_len(k), function(i) {
    keep <- seq_len(k) != i
    sum(w[keep] * bx[keep] * by[keep]) / sum(w[keep] * bx[keep]^2)
  })
  steiger_rx <- 2 * dat$AF1 * (1 - dat$AF1) * bx^2
  eaf_y <- dat$EAF_OUT_HARMONISED
  steiger_ry <- 2 * eaf_y * (1 - eaf_y) * by^2

  methods <- data.frame(
    outcome = outcome_name,
    method = c("IVW_multiplicative_random_effects", "IVW_fixed_effects", "Weighted_median_bootstrap", "MR_Egger"),
    nsnp = k,
    beta = c(ivw_beta, ivw_beta, wm_beta, egger_beta),
    se = c(ivw_se_mre, ivw_se_fixed, wm_se, egger_se),
    p = c(ivw_p, 2 * pnorm(abs(ivw_beta / ivw_se_fixed), lower.tail = FALSE), wm_p, egger_p),
    ci_low = c(ivw_beta - 1.96 * ivw_se_mre, ivw_beta - 1.96 * ivw_se_fixed, wm_beta - 1.96 * wm_se, egger_beta - 1.96 * egger_se),
    ci_high = c(ivw_beta + 1.96 * ivw_se_mre, ivw_beta + 1.96 * ivw_se_fixed, wm_beta + 1.96 * wm_se, egger_beta + 1.96 * egger_se)
  )
  diagnostics <- data.frame(
    outcome = outcome_name,
    nsnp = k,
    cochran_q = q,
    cochran_q_df = q_df,
    cochran_q_p = q_p,
    phi = phi,
    egger_intercept = egger_intercept,
    egger_intercept_se = egger_intercept_se,
    egger_intercept_p = egger_intercept_p,
    loo_min = min(loo),
    loo_max = max(loo),
    steiger_correct_count = sum(steiger_rx > steiger_ry, na.rm = TRUE),
    steiger_total = sum(is.finite(steiger_rx) & is.finite(steiger_ry))
  )
  list(methods = methods, diagnostics = diagnostics)
}

full <- run_mr(file.path(base_dir, "cpi_full_ldclumped_harmonised.csv"), "CPI_full")
noimg <- run_mr(file.path(base_dir, "cpi_no_imaging_ldclumped_harmonised.csv"), "CPI_no_imaging")

methods <- rbind(full$methods, noimg$methods)
diagnostics <- rbind(full$diagnostics, noimg$diagnostics)
write.csv(methods, file.path(base_dir, "formal_mr_methods.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(base_dir, "formal_mr_diagnostics.csv"), row.names = FALSE)
print(methods)
print(diagnostics)
