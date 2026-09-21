# =============================================================================
# Stage 5 - independent R validation of the site-specific pain analyses
# -----------------------------------------------------------------------------
# Input : r_input/R_stage1_<OUTCOME>.tsv   (42 TMFI instruments)
#         r_input/R_stage2_<OUTCOME>.tsv   (93 TMFI + BMI instruments)
#         Both are already harmonised against the 1000G EUR panel, so this
#         script does NOT re-harmonise; it feeds the estimates straight into
#         TwoSampleMR / MVMR.
#
# Run   : Rscript stage5_validation.R
# Output: results/stage5_R_stage1_<OUTCOME>.csv
#         results/stage5_R_stage2_<OUTCOME>.csv
#         results/stage5_R_all.csv
#
# Packages (uncomment the first run only):
#   install.packages(c("remotes","data.table"))
#   remotes::install_github("MRCIEU/TwoSampleMR")
#   remotes::install_github("WSpiller/RMVMR")        # MVMR
#   remotes::install_github("rondolab/MR-PRESSO")    # optional, MR-PRESSO
# =============================================================================

suppressMessages({
  library(data.table)
  library(TwoSampleMR)
})
# MVMR is the only package that may be unavailable (it comes from GitHub's
# RMVMR). Load it defensively so the Stage 1 (TwoSampleMR) results still run if
# MVMR failed to install. Stage 2 simply skips if MVMR is missing.
MVMR_OK <- suppressWarnings(require("MVMR", quietly = TRUE))
if (!MVMR_OK)
  cat("NOTE: MVMR not loaded -> Stage 2 (MVMR) will be skipped; Stage 1 still runs.\n")

rin  <- "D:/SX/TMFI_CPI_MR/05_stage5_limbpain/r_input"
out  <- "D:/SX/TMFI_CPI_MR/05_stage5_limbpain/results"
dir.create(out, showWarnings = FALSE)

OUTCOMES <- c("BACK_PAIN_UKB", "KNEE_PAIN_UKB", "KNEE_PAIN_FIGSHARE")
SD_TMFI  <- 2.030884204160846
SD_BMI   <- 1.2307496342161395

or_ci <- function(b, se, sd) {
  est <- b * sd
  c(OR = exp(est),
    lo = exp(est - 1.96 * se * sd),
    hi = exp(est + 1.96 * se * sd),
    p  = 2 * pnorm(-abs(b / se)))
}

all_rows <- list()

for (oc in OUTCOMES) {

  # ---------------------------------------------------------------- Stage 1 --
  d1 <- fread(file.path(rin, paste0("R_stage1_", oc, ".tsv")))
  dat <- data.frame(
    SNP                     = d1$SNP,
    beta.exposure           = d1$tmfi_beta,
    se.exposure             = d1$tmfi_se,
    beta.outcome            = d1$out_beta,
    se.outcome              = d1$out_se,
    eaf.exposure            = d1$tmfi_af,
    eaf.outcome             = d1$out_af,
    samplesize.exposure     = d1$tmfi_n,
    samplesize.outcome      = d1$out_n,
    exposure = "TMFI", outcome = oc,
    id.exposure = "TMFI", id.outcome = oc,
    units.exposure = "SD", units.outcome = "log odds",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )

  cat("\n================ ", oc, " ================\n")
  cat("Stage 1 instruments:", nrow(dat), "\n")

  res1 <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre",
                                  "mr_weighted_median",
                                  "mr_egger_regression",
                                  "mr_simple_mode", "mr_weighted_mode"))
  res1$OR_per_SD     <- exp(res1$b * SD_TMFI)
  res1$OR_CI_low     <- exp((res1$b - 1.96 * res1$se) * SD_TMFI)
  res1$OR_CI_high    <- exp((res1$b + 1.96 * res1$se) * SD_TMFI)
  res1$outcome <- oc
  print(res1[, c("method", "nsnp", "b", "se", "pval",
                 "OR_per_SD", "OR_CI_low", "OR_CI_high")])

  het <- try(mr_heterogeneity(dat), silent = TRUE)
  if (!inherits(het, "try-error")) {
    print(het[, c("method", "Q", "Q_df", "Q_pval")])
    fwrite(het, file.path(out, paste0("stage5_R_heterogeneity_", oc, ".csv")))
  }

  ple <- try(mr_pleiotropy_test(dat), silent = TRUE)
  if (!inherits(ple, "try-error")) {
    cat("MR-Egger intercept:", ple$egger_intercept, " P =", ple$pval, "\n")
  }

  ste <- try(directionality_test(dat), silent = TRUE)
  if (!inherits(ste, "try-error")) {
    cat("Steiger: correct causal direction =", ste$correct_causal_direction,
        " P =", ste$pval, "\n")
    ste$outcome <- oc
    fwrite(ste, file.path(out, paste0("stage5_R_steiger_", oc, ".csv")))
  }

  # MR-PRESSO (optional package)
  if (requireNamespace("MRPRESSO", quietly = TRUE)) {
    pr <- try(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                  BetaExposure = "beta.exposure",
                                  SdOutcome = "se.outcome",
                                  SdExposure = "se.exposure",
                                  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                  data = dat, NbDistribution = 1000,
                                  SignifThreshold = 0.05), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      cat("MR-PRESSO global test P =",
          pr$`Main Test`$`P-value`[1], "\n")
    }
  } else {
    cat("MRPRESSO not installed - skipped\n")
  }

  fwrite(res1, file.path(out, paste0("stage5_R_stage1_", oc, ".csv")))
  all_rows[[paste0("S1_", oc)]] <- res1

  # ---------------------------------------------------------------- Stage 2 --
  if (MVMR_OK) {
  d2 <- fread(file.path(rin, paste0("R_stage2_", oc, ".tsv")))
  mv <- format_mvmr(
    BXGs  = cbind(d2$tmfi_beta, d2$bmi_beta),
    BYG   = d2$out_beta,
    seBXGs = cbind(d2$tmfi_se, d2$bmi_se),
    seBYG = d2$out_se,
    RSID  = d2$SNP
  )
  cat("\nStage 2 instruments:", nrow(mv), "\n")

  ivw_mv <- mr_mvmr(mv)
  print(ivw_mv)

  str_mv <- try(strength_mvmr(mv), silent = TRUE)
  if (!inherits(str_mv, "try-error")) print(str_mv)

  pl_mv <- try(pleiotropy_mvmr(mv), silent = TRUE)
  if (!inherits(pl_mv, "try-error")) print(pl_mv)

  he_mv <- try(heterogeneity_mvmr(mv), silent = TRUE)
  if (!inherits(he_mv, "try-error")) print(he_mv)

  eg_mv <- try(mvmr_egger(mv), silent = TRUE)
  if (!inherits(eg_mv, "try-error")) print(eg_mv)

  # assemble a tidy table
  b   <- ivw_mv$Estimate
  se  <- ivw_mv$`Std. Error`
  p   <- ivw_mv$`P-value`
  tab <- data.frame(
    outcome   = oc,
    exposure  = c("TMFI", "BMI"),
    nsnp      = nrow(mv),
    beta      = b,
    se        = se,
    pval      = p,
    OR_per_SD = c(exp(b[1] * SD_TMFI), exp(b[2] * SD_BMI)),
    OR_lo     = c(exp((b[1] - 1.96 * se[1]) * SD_TMFI),
                  exp((b[2] - 1.96 * se[2]) * SD_BMI)),
    OR_hi     = c(exp((b[1] + 1.96 * se[1]) * SD_TMFI),
                  exp((b[2] + 1.96 * se[2]) * SD_BMI)),
    stringsAsFactors = FALSE
  )
  print(tab)
  fwrite(tab, file.path(out, paste0("stage5_R_stage2_", oc, ".csv")))
  all_rows[[paste0("S2_", oc)]] <- tab
  } else {
    cat("  MVMR not installed -> Stage 2 skipped for", oc, "\n")
  }
}

fwrite(rbindlist(all_rows, fill = TRUE),
       file.path(out, "stage5_R_all.csv"))
cat("\n[done] outputs in", out, "\n")
