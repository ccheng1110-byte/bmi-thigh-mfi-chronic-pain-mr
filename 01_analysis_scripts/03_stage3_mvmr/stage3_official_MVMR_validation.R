# Stage 3 independent validation with the official MVMR package
# Run in RStudio or with:
# Rscript D:/SX/TMFI_CPI_MR/03_stage3_mvmr/stage3_official_MVMR_validation.R

project_dir <- "D:/SX/TMFI_CPI_MR"
stage3_dir <- file.path(project_dir, "03_stage3_mvmr")
local_lib <- file.path(project_dir, "Rlib")

.libPaths(c(local_lib, .libPaths()))
if (!requireNamespace("MVMR", quietly = TRUE)) {
  stop("MVMR is not installed. Install it first, then rerun this script.")
}

run_one <- function(input_file, outcome_label) {
  dat <- read.csv(file.path(stage3_dir, input_file), check.names = FALSE)

  r_input <- MVMR::format_mvmr(
    BXGs = dat[, c("TMFI_BETA", "BMI_BETA_HARMONISED")],
    BYG = dat$OUTCOME_BETA_H,
    seBXGs = dat[, c("TMFI_SE", "BMI_SE")],
    seBYG = dat$OUTCOME_SE,
    RSID = dat$SNP
  )

  ivw <- MVMR::ivw_mvmr(r_input, gencov = 0)
  strength <- suppressWarnings(MVMR::strength_mvmr(r_input, gencov = 0))
  pleiotropy <- suppressWarnings(MVMR::pleiotropy_mvmr(r_input, gencov = 0))

  ivw_out <- data.frame(
    outcome = outcome_label,
    exposure = c("TMFI", "BMI"),
    beta = ivw[, "Estimate"],
    se = ivw[, "Std. Error"],
    t = ivw[, "t value"],
    p = ivw[, "Pr(>|t|)"],
    row.names = NULL
  )
  strength_out <- data.frame(
    outcome = outcome_label,
    exposure = c("TMFI", "BMI"),
    conditional_F = as.numeric(strength[1, ]),
    row.names = NULL
  )
  pleiotropy_out <- data.frame(
    outcome = outcome_label,
    Q = as.numeric(pleiotropy$Qstat),
    df = nrow(dat) - 2 - 1,
    p = as.numeric(pleiotropy$Qpval),
    row.names = NULL
  )

  list(ivw = ivw_out, strength = strength_out, pleiotropy = pleiotropy_out)
}

full <- run_one("cpi_full_mvmr_harmonised_usable.csv", "CPI_full")
no_imaging <- run_one(
  "cpi_no_imaging_mvmr_harmonised_usable.csv",
  "CPI_no_imaging"
)

write.csv(
  rbind(full$ivw, no_imaging$ivw),
  file.path(stage3_dir, "R_validation_ivw.csv"),
  row.names = FALSE
)
write.csv(
  rbind(full$strength, no_imaging$strength),
  file.path(stage3_dir, "R_validation_strength.csv"),
  row.names = FALSE
)
write.csv(
  rbind(full$pleiotropy, no_imaging$pleiotropy),
  file.path(stage3_dir, "R_validation_pleiotropy.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(stage3_dir, "R_validation_sessionInfo.txt")
)

message("Done. R validation files were written to: ", stage3_dir)
