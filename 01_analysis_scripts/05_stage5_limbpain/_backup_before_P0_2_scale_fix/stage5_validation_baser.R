# =============================================================================
# Stage 5 - base-R replication of the site-specific pain MR estimates
# -----------------------------------------------------------------------------
# Purpose: arithmetic cross-check of the Python pipeline (stage5_mr.py) using
#          ONLY base R - no TwoSampleMR / MVMR / MRPRESSO required, so it runs
#          even if those packages cannot be installed.
#          Formulas mirror the Python implementation exactly:
#            * IVW-FE      : w = 1/seY^2 on (bX, bY)
#            * IVW-MRE     : same point estimate, SE = SE_FE * sqrt(max(1, Q/df))
#                            (multiplicative random effects)
#            * Weighted median : ratio-based, bootstrap SE (1000 reps)
#            * MR-Egger    : weighted lm(bY ~ bX, weights = 1/seY^2)
#            * MVMR (IVW)  : beta = (B'VB)^-1 B'Vy ; SE = sqrt(phi * diag((B'VB)^-1))
# Input : r_input/R_stage1_<OUTCOME>.tsv, r_input/R_stage2_<OUTCOME>.tsv
# Output: results/stage5_baser_all.csv
# =============================================================================

rin <- "D:/SX/TMFI_CPI_MR/05_stage5_limbpain/r_input"
out <- "D:/SX/TMFI_CPI_MR/05_stage5_limbpain/results"
dir.create(out, showWarnings = FALSE)

OUTCOMES <- c("BACK_PAIN_UKB", "KNEE_PAIN_UKB", "KNEE_PAIN_FIGSHARE")
SD_TMFI  <- 2.030884204160846
SD_BMI   <- 1.2307496342161395

# ---------------------------------------------------------------- IVW -------
ivw_fe <- function(bX, bY, seY) {
  w  <- 1 / seY^2
  b  <- sum(w * bX * bY) / sum(w * bX^2)
  se <- sqrt(1 / sum(w * bX^2))
  Q  <- sum(w * (bY - b * bX)^2)
  df <- length(bX) - 1
  phi <- max(1, Q / df)
  list(b = b, se = se, Q = Q, df = df, phi = phi,
       se_mre = se * sqrt(phi))
}

# ------------------------------------------------------- weighted median ----
wmedian <- function(bX, bY, seY, B = 1000) {
  ratio <- bY / bX
  seR   <- abs(seY / bX)
  w     <- 1 / seR^2
  # Bowden's weighted median WITH linear interpolation between the two
  # ratios that straddle the 0.5 cumulative-weight point
  wm <- function(r, ww) {
    o <- order(r); r <- r[o]; ww <- ww[o]
    cw <- cumsum(ww) / sum(ww)
    i  <- which(cw >= 0.5)[1]
    if (i == 1) return(r[1])
    BD <- r[i - 1]; BK <- r[i]
    WD <- cw[i - 1]; WK <- cw[i]
    ((0.5 - WD) * BK + (WK - 0.5) * BD) / (WK - WD)
  }
  est <- wm(ratio, w)
  set.seed(20260919)
  bs <- replicate(B, {
    s <- sample(length(ratio), replace = TRUE)
    wm(ratio[s], w[s])
  })
  list(b = est, se = sd(bs))
}

# -------------------------------------------------------------- MR-Egger ----
egger <- function(bX, bY, seY) {
  w   <- 1 / seY^2
  fit <- lm(bY ~ bX, weights = w)
  cf  <- summary(fit)$coefficients
  list(b = cf[2, 1], se = cf[2, 2],
       intercept = cf[1, 1], intercept_se = cf[1, 2],
       intercept_p = cf[1, 4])
}

# ------------------------------------------------------------ MVMR (IVW) ----
mvmr_ivw <- function(BX, BY, seBY) {
  W <- diag(1 / seBY^2)
  A <- t(BX) %*% W %*% BX
  b <- solve(A) %*% (t(BX) %*% W %*% BY)
  res <- BY - BX %*% b
  Q   <- as.numeric(sum(res^2 / seBY^2))
  df  <- length(BY) - ncol(BX)
  phi <- max(1, Q / df)
  se  <- sqrt(phi * diag(solve(A)))
  list(b = as.numeric(b), se = as.numeric(se), Q = Q, df = df, phi = phi)
}

or_ci <- function(b, se, sd) c(OR = exp(b * sd),
                               lo = exp((b - 1.96 * se) * sd),
                               hi = exp((b + 1.96 * se) * sd),
                               p  = 2 * pnorm(-abs(b / se)))

rows <- list()

for (oc in OUTCOMES) {

  # ------------------------------------------------------------ Stage 1 ----
  d1  <- read.delim(file.path(rin, paste0("R_stage1_", oc, ".tsv")),
                    stringsAsFactors = FALSE)
  bX  <- d1$tmfi_beta; seX <- d1$tmfi_se
  bY  <- d1$out_beta;  seY <- d1$out_se
  cat("\n========", oc, "========\n")
  cat("Stage 1 instruments:", nrow(d1), "\n")

  f  <- ivw_fe(bX, bY, seY)
  o1 <- or_ci(f$b, f$se, SD_TMFI)
  o2 <- or_ci(f$b, f$se_mre, SD_TMFI)
  cat(sprintf("IVW-FE   : OR=%.6f [%.6f, %.6f]  P=%.4g\n", o1[1], o1[2], o1[3], o1[4]))
  cat(sprintf("IVW-MRE  : OR=%.6f [%.6f, %.6f]  P=%.4g  (phi=%.4f, Q=%.2f, df=%d)\n",
              o2[1], o2[2], o2[3], o2[4], f$phi, f$Q, f$df))
  rows[[paste(oc, "IVW-FE")]]  <- data.frame(outcome = oc, method = "IVW-FE",
      beta = f$b, se = f$se, OR = o1[1], lo = o1[2], hi = o1[3], p = o1[4])
  rows[[paste(oc, "IVW-MRE")]] <- data.frame(outcome = oc, method = "IVW-MRE",
      beta = f$b, se = f$se_mre, OR = o2[1], lo = o2[2], hi = o2[3], p = o2[4])

  wm <- wmedian(bX, bY, seY)
  ow <- or_ci(wm$b, wm$se, SD_TMFI)
  cat(sprintf("WMedian  : OR=%.6f [%.6f, %.6f]  P=%.4g\n", ow[1], ow[2], ow[3], ow[4]))
  rows[[paste(oc, "WM")]] <- data.frame(outcome = oc, method = "weighted median",
      beta = wm$b, se = wm$se, OR = ow[1], lo = ow[2], hi = ow[3], p = ow[4])

  eg <- egger(bX, bY, seY)
  oe <- or_ci(eg$b, eg$se, SD_TMFI)
  cat(sprintf("MR-Egger : OR=%.6f [%.6f, %.6f]  P=%.4g  intercept=%.3g (P=%.3g)\n",
              oe[1], oe[2], oe[3], oe[4], eg$intercept, eg$intercept_p))
  rows[[paste(oc, "Egger")]] <- data.frame(outcome = oc, method = "MR-Egger",
      beta = eg$b, se = eg$se, OR = oe[1], lo = oe[2], hi = oe[3], p = oe[4])

  # ------------------------------------------------------------ Stage 2 ----
  d2 <- read.delim(file.path(rin, paste0("R_stage2_", oc, ".tsv")),
                   stringsAsFactors = FALSE)
  BX <- cbind(d2$tmfi_beta, d2$bmi_beta)
  BY <- d2$out_beta; seBY <- d2$out_se
  mv <- mvmr_ivw(BX, BY, seBY)
  cat("Stage 2 instruments:", nrow(d2), "\n")

  ot <- or_ci(mv$b[1], mv$se[1], SD_TMFI)
  ob <- or_ci(mv$b[2], mv$se[2], SD_BMI)
  cat(sprintf("MVMR TMFI: beta=%.6g se=%.6g OR=%.6f [%.6f, %.6f] P=%.4g\n",
              mv$b[1], mv$se[1], ot[1], ot[2], ot[3], ot[4]))
  cat(sprintf("MVMR BMI : beta=%.6g se=%.6g OR=%.6f [%.6f, %.6f] P=%.4g\n",
              mv$b[2], mv$se[2], ob[1], ob[2], ob[3], ob[4]))
  cat(sprintf("          (Q=%.2f, df=%d, phi=%.4f)\n", mv$Q, mv$df, mv$phi))

  rows[[paste(oc, "MVMR_TMFI")]] <- data.frame(outcome = oc, method = "MVMR:TMFI",
      beta = mv$b[1], se = mv$se[1], OR = ot[1], lo = ot[2], hi = ot[3], p = ot[4])
  rows[[paste(oc, "MVMR_BMI")]]  <- data.frame(outcome = oc, method = "MVMR:BMI",
      beta = mv$b[2], se = mv$se[2], OR = ob[1], lo = ob[2], hi = ob[3], p = ob[4])
}

all_rows <- do.call(rbind, rows)
write.csv(all_rows, file.path(out, "stage5_baser_all.csv"), row.names = FALSE)
cat("\n[done] wrote results/stage5_baser_all.csv\n")
print(all_rows)
