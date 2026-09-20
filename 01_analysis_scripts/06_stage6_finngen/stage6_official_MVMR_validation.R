# Stage 6 (FinnGen external replication) validation with the official MVMR package.
# Run with:
# Rscript D:/SX/TMFI_CPI_MR/06_stage6_finngen/stage6_official_MVMR_validation.R
#
# Reads the SAME frozen tab-separated inputs that the Python (stage6_mr.py)
# and base-R (stage6_validation_baser.R) implementations consume, so the three
# implementations are compared on identical numbers.

project_dir <- "D:/SX/TMFI_CPI_MR"
stage6_dir <- file.path(project_dir, "06_stage6_finngen")
local_lib <- file.path(project_dir, "Rlib")

.libPaths(c(local_lib, .libPaths()))
if (!requireNamespace("MVMR", quietly = TRUE)) {
  stop("MVMR is not installed. Install it first, then rerun this script.")
}

endpoints <- c(
  "FG_PAIN" = "FG_PAIN",
  "FG_LIMBPAIN" = "FG_LIMBPAIN",
  "FG_FIBROMYALGIA" = "FG_FIBROMYALGIA"
)

run_one <- function(label) {
  f <- file.path(stage6_dir, "r_input", paste0("R_stage2_", label, ".tsv"))
  dat <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)

  r_input <- MVMR::format_mvmr(
    BXGs = dat[, c("tmfi_beta", "bmi_beta")],
    BYG = dat$out_beta,
    seBXGs = dat[, c("tmfi_se", "bmi_se")],
    seBYG = dat$out_se,
    RSID = dat$SNP
  )

  ivw <- MVMR::ivw_mvmr(r_input, gencov = 0)
  strength <- suppressWarnings(MVMR::strength_mvmr(r_input, gencov = 0))
  pleiotropy <- suppressWarnings(MVMR::pleiotropy_mvmr(r_input, gencov = 0))

  list(
    ivw = data.frame(
      outcome = label,
      exposure = c("TMFI", "BMI"),
      n_snp = nrow(dat),
      beta = as.numeric(ivw[, "Estimate"]),
      se = as.numeric(ivw[, "Std. Error"]),
      t = as.numeric(ivw[, "t value"]),
      p_t = as.numeric(ivw[, "Pr(>|t|)"]),
      row.names = NULL
    ),
    strength = data.frame(
      outcome = label,
      exposure = c("TMFI", "BMI"),
      conditional_F = as.numeric(strength[1, ]),
      row.names = NULL
    ),
    pleiotropy = data.frame(
      outcome = label,
      Q = as.numeric(pleiotropy$Qstat),
      df = nrow(dat) - 3,
      p = as.numeric(pleiotropy$Qpval),
      row.names = NULL
    )
  )
}

res <- lapply(names(endpoints), function(k) run_one(k))

write.csv(
  do.call(rbind, lapply(res, function(x) x$ivw)),
  file.path(stage6_dir, "results", "R_validation_ivw.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, lapply(res, function(x) x$strength)),
  file.path(stage6_dir, "results", "R_validation_strength.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, lapply(res, function(x) x$pleiotropy)),
  file.path(stage6_dir, "results", "R_validation_pleiotropy.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(stage6_dir, "results", "R_validation_sessionInfo.txt")
)

# Echo to stdout so the numbers can be read directly from the console.
for (x in res) {
  cat("\n=== ", x$ivw$outcome[1], " (n = ", x$ivw$n_snp[1], " SNPs) ===\n", sep = "")
  print(x$ivw, row.names = FALSE)
  print(x$strength, row.names = FALSE)
  cat("Q = ", x$pleiotropy$Q, " df = ", x$pleiotropy$df,
      " p = ", x$pleiotropy$p, "\n", sep = "")
}

message("Done. Files written to: ", file.path(stage6_dir, "results"))
