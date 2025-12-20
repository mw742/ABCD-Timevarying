#!/usr/bin/env Rscript

## - Keep ALL visits (baseline + follow-ups)
## - Merge them by FID/IID/eventname
## - Attach interview_age from abcd_y_lt.csv and compute time_months since baseline
## - Merge core covariates (sex, age, site, PCs) from covariates_final.tsv
## - Do basic cleaning:
##     * drop all-NA columns
##     * drop near-constant numeric columns
##     * drop highly correlated numeric columns
##   (NO high-missingness filter; keep baseline-only variables even if very sparse)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

message_ts <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = " "))
}

## -----------------------------
## Paths
## -----------------------------

## Core age table
age_file <- "/gpfs/group/ehlers/mwei/Tables/core/abcd-general/abcd_y_lt.csv"

## Core covariates (sex, age, site, PCs)
covariates_file <- "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/covariates_final.tsv"

## survival covariate directories
surv_dirs <- c(
  "/gpfs/group/ehlers/mwei/Tables/core/abcd-general/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/culture-environment/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/gender-identity-sexual-health/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/linked-external-data/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/mental-health/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/neurocognition/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/novel-technologies/surv",
  "/gpfs/group/ehlers/mwei/Tables/core/physical-health/surv"
)

## Output directory = current directory
out_dir <- getwd()
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

raw_out  <- file.path(out_dir, "cov_timevarying_raw.tsv")
filt_out <- file.path(out_dir, "cov_timevarying_filtered.tsv")
meta_out <- file.path(out_dir, "cov_timevarying_metadata.tsv")

## -----------------------------
## 1. Load age table, compute time_months
## -----------------------------

message_ts("Reading age table from: ", age_file)
age_dt <- fread(age_file)

if (!all(c("src_subject_id", "eventname") %in% names(age_dt))) {
  stop("age_file must contain 'src_subject_id' and 'eventname'.")
}

age_col <- grep("^interview_age", names(age_dt), value = TRUE)[1]
if (is.na(age_col)) {
  stop("No column starting with 'interview_age' found in abcd_y_lt.csv.")
}

setnames(age_dt, age_col, "interview_age")
message_ts("Using age column: ", age_col, " (renamed to 'interview_age').")

med_age <- median(age_dt$interview_age, na.rm = TRUE)
if (is.na(med_age)) stop("interview_age is all NA; cannot compute time.")

if (med_age < 25) {
  message_ts("interview_age seems to be in YEARS; converting to MONTHS.")
  age_dt[, interview_age := interview_age * 12]
} else {
  message_ts("interview_age seems to be in MONTHS; using as-is.")
}

age_dt <- age_dt[!is.na(interview_age)]

age_dt[, baseline_age := min(interview_age, na.rm = TRUE), by = src_subject_id]
age_dt[!is.finite(baseline_age), baseline_age := NA_real_]

age_dt[, time_months := interview_age - baseline_age]

age_dt <- age_dt[, .(src_subject_id, eventname, interview_age, baseline_age, time_months)]

message_ts("Age table dimensions (after processing): ",
           nrow(age_dt), " x ", ncol(age_dt))

## -----------------------------
## 2. Collect all files
## -----------------------------

all_long_files <- character(0)
for (d in surv_dirs) {
  if (!dir.exists(d)) {
    message_ts("Directory does not exist, skipping: ", d)
    next
  }
  files <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0L) {
    message_ts("No .tsv files in: ", d)
  } else {
    message_ts("Found ", length(files), " .tsv files in: ", d)
    all_long_files <- c(all_long_files, files)
  }
}
all_long_files <- sort(unique(all_long_files))

message_ts("Total .tsv covariate files found: ", length(all_long_files))

if (length(all_long_files) == 0L) {
  stop("No *_long.tsv files found in the given surv directories.")
}

## -----------------------------
## 3. Read and merge long covariates
## -----------------------------

read_long_file <- function(fpath) {
  dt <- tryCatch(
    fread(fpath),
    error = function(e) {
      message_ts("Failed to read: ", fpath, " | ", e$message)
      return(NULL)
    }
  )
  if (is.null(dt)) return(NULL)

  if (!"FID" %in% names(dt) & "IID" %in% names(dt)) {
    dt[, FID := IID]
  }
  if (!"IID" %in% names(dt) & "FID" %in% names(dt)) {
    dt[, IID := FID]
  }
  if (!"src_subject_id" %in% names(dt) & "IID" %in% names(dt)) {
    dt[, src_subject_id := IID]
  }

  required <- c("FID", "IID", "src_subject_id")
  if (!all(required %in% names(dt))) {
    message_ts("Skipping (missing FID/IID/src_subject_id): ", fpath)
    return(NULL)
  }

  if (!"eventname" %in% names(dt)) {
    message_ts("No eventname column in ", fpath,
               "; assuming all rows are baseline_year_1_arm_1.")
    dt[, eventname := "baseline_year_1_arm_1"]
  }

  dt <- unique(dt)

  stem <- basename(fpath)
  stem <- sub("\\.tsv$", "", stem)

  key_cols <- c("FID", "IID", "src_subject_id", "eventname")
  other_cols <- setdiff(names(dt), key_cols)

  if (length(other_cols) > 0L) {
    new_names <- paste0(stem, "__", other_cols)
    setnames(dt, other_cols, new_names)
  }

  return(dt[, c(key_cols, setdiff(names(dt), key_cols)), with = FALSE])
}

cov_list <- list()
for (f in all_long_files) {
  message_ts("Reading covariates from: ", f)
  dt <- read_long_file(f)
  if (!is.null(dt)) cov_list[[length(cov_list) + 1L]] <- dt
}

if (length(cov_list) == 0L) {
  stop("No usable covariate tables were loaded after reading *_long.tsv files.")
}

message_ts("Merging ", length(cov_list),
           " covariate tables by FID/IID/eventname...")

cov_long <- Reduce(function(x, y) {
  merge(x, y,
        by = c("FID", "IID", "eventname"),
        all = TRUE)
}, cov_list)

message_ts("Merged long covariates dimension: ",
           nrow(cov_long), " x ", ncol(cov_long))

## -----------------------------
## 4. Attach time_months from age table
## -----------------------------

message_ts("Merging age table to covariates by src_subject_id + eventname...")
cov_long <- merge(cov_long, age_dt,
                  by = c("src_subject_id", "eventname"),
                  all.x = TRUE)

message_ts("After merging age: ",
           nrow(cov_long), " x ", ncol(cov_long))

## -----------------------------
## 5. Merge core covariates (sex, age, site, PCs)
## -----------------------------

message_ts("Reading core covariates from: ", covariates_file)
cov_core <- fread(covariates_file)

if (!"FID" %in% names(cov_core) & "IID" %in% names(cov_core)) {
  cov_core[, FID := IID]
}
if (!"IID" %in% names(cov_core) & "FID" %in% names(cov_core)) {
  cov_core[, IID := FID]
}

pc_cols <- grep("^PC[0-9]+$", names(cov_core), value = TRUE)

keep_core <- c("FID", "IID", "sex", "age", "site", pc_cols)
keep_core <- intersect(keep_core, names(cov_core))

cov_core <- cov_core[, ..keep_core]

message_ts("Core covariates columns kept: ",
           paste(keep_core, collapse = ", "))

message_ts("Merging core covariates...")
cov_long <- merge(cov_long, cov_core, by = c("FID", "IID"), all.x = TRUE)

message_ts("Dimension after adding core covariates: ",
           nrow(cov_long), " x ", ncol(cov_long))

## -----------------------------
## 6. Cleaning: drop all-NA / near-constant / high-cor
## -----------------------------

key_cols <- c("FID", "IID", "eventname",
              "interview_age", "baseline_age", "time_months",
              "sex", "age", "site", pc_cols)
key_cols <- intersect(key_cols, names(cov_long))

non_key <- setdiff(names(cov_long), key_cols)

## Track drops for metadata
near_zero    <- character(0)
drop_missing <- character(0)   ## We KEEP all vars regardless of missingness
drop_high_cor <- character(0)

if (length(non_key) == 0L) {
  message_ts("No non-key columns; skipping cleaning filters.")
  cov_filt <- cov_long[, ..key_cols]
} else {

  ## Drop all-NA
  all_na <- sapply(cov_long[, ..non_key], function(x) all(is.na(x)))
  message_ts("Dropping ", sum(all_na), " non-key columns that are all NA.")
  non_key <- non_key[!all_na]
  cov_long <- cov_long[, c(key_cols, non_key), with = FALSE]

  ## Near-constant numeric
  if (length(non_key) > 0L) {
    is_num <- sapply(cov_long[, ..non_key], is.numeric)
    num_cols <- non_key[is_num]

    if (length(num_cols) > 0L) {
      vars <- sapply(cov_long[, ..num_cols], function(x) {
        x <- x[!is.na(x)]
        if (length(x) < 2L) return(0)
        var(x)
      })
      near_zero <- names(vars)[vars < 1e-6]
      message_ts("Dropping ", length(near_zero),
                 " near-constant numeric non-key columns.")
      num_cols <- setdiff(num_cols, near_zero)
      non_key <- setdiff(non_key, near_zero)
      cov_long <- cov_long[, c(key_cols, non_key), with = FALSE]
    }
  }


  ## High collinearity among numeric vars
  if (length(non_key) > 1L) {
    is_num2 <- sapply(cov_long[, ..non_key], is.numeric)
    num_cols2 <- non_key[is_num2]

    if (length(num_cols2) > 1L) {
      message_ts("Computing correlation matrix for numeric non-key columns...")
      cor_mat <- tryCatch(
        cor(cov_long[, ..num_cols2], use = "pairwise.complete.obs"),
        error = function(e) {
          message_ts("Correlation computation failed: ", e$message)
          NULL
        }
      )
      if (!is.null(cor_mat)) {
        cor_cutoff <- 0.9
        abs_cor <- abs(cor_mat)
        np <- ncol(abs_cor)
        cn <- colnames(abs_cor)
        to_drop <- rep(FALSE, np)

        for (i in seq_len(np - 1L)) {
          if (to_drop[i]) next
          high <- which(abs_cor[i, (i + 1L):np] > cor_cutoff)
          if (length(high) > 0L) {
            idx <- (i + 1L):np
            to_drop[idx[high]] <- TRUE
          }
        }
        drop_high_cor <- cn[to_drop]
        message_ts("Dropping ", length(drop_high_cor),
                   " highly correlated numeric non-key columns (|r|>", cor_cutoff, ").")

        non_key <- setdiff(non_key, drop_high_cor)
        cov_long <- cov_long[, c(key_cols, non_key), with = FALSE]
      }
    }
  }

  cov_filt <- cov_long[, c(key_cols, non_key), with = FALSE]
}

message_ts("Final filtered covariate dimensions: ",
           nrow(cov_filt), " x ", ncol(cov_filt))

## -----------------------------
## 7. Save outputs + metadata
## -----------------------------

message_ts("Writing raw covariates (after structural cleaning but before high-cor filter) to: ", raw_out)
fwrite(cov_long, raw_out, sep = "\t", quote = FALSE, na = "NA")

message_ts("Writing filtered covariates to: ", filt_out)
fwrite(cov_filt, filt_out, sep = "\t", quote = FALSE, na = "NA")

meta_dt <- data.table(
  variable = names(cov_long),
  is_key = names(cov_long) %in% key_cols,
  is_numeric = sapply(cov_long, is.numeric),
  missing_fraction = sapply(cov_long, function(x) mean(is.na(x))),
  dropped_near_zero = names(cov_long) %in% near_zero,
  dropped_high_missing = names(cov_long) %in% drop_missing,
  dropped_high_cor = names(cov_long) %in% drop_high_cor
)

message_ts("Writing metadata to: ", meta_out)
fwrite(meta_dt, meta_out, sep = "\t", quote = FALSE, na = "NA")

message_ts("0_build_data.R completed.")

