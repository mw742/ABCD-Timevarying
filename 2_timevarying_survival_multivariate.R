#!/usr/bin/env Rscript

## 1) Keep predictors that appear later (no back-fill/leakage):
##    - Use LOCF forward-fill only 
##    - For modeling, create X__obs indicator + median-impute X (so late-measured / one-time vars can be kept)
## 2) Drop only predictors with <50% subject coverage (ever observed in <50% individuals)
## 3) Univariate prefilter: keep ALL p <= 0.05 (no TOP_N limit)
## 4) LASSO: lambda.min, keep ALL non-zero coefficients (no TopK)
## 5) Forced covariates: sex/age/site/PCs always adjusted, but NOT shown in predictors-only output
##    - site handled as strata(site) so it never appears as coefficients

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(glmnet)
})

message_ts <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          " | ", paste0(..., collapse = " "))
}

## -----------------------------
## THRESHOLD SETTINGS
## -----------------------------
P_THRESH <- 0.05
MIN_SUBJECT_COVERAGE <- 0.50     # drop predictors observed in <50% subjects (ever observed at any time)
N_FOLDS <- 10
SEED <- 12345

## glmnet stability knobs 
GLMNET_MAXIT <- 200000L
GLMNET_STANDARDIZE <- TRUE

## -----------------------------
## Paths
## -----------------------------
base_dir <- getwd()
cov_file <- file.path(base_dir, "cov_timevarying_filtered.tsv")

phen_dir <- "/gpfs/group/ehlers/mwei/ABCD/init_phenos"
phen_files <- c(
  alcohol  = file.path(phen_dir, "alcohol_survival.tsv"),
  cannabis = file.path(phen_dir, "cannabis_survival.tsv"),
  nicotine = file.path(phen_dir, "nicotine_survival.tsv"),
  any      = file.path(phen_dir, "any_survival.tsv")
)

prs_paths <- list(
  PAU      = "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/GWAS_sud/PRS_scores_PRScsx/PAU_multi_PRS_two_col.tsv",
  CanUD    = "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/GWAS_sud/Levey2023_CanUD/prs_scores_prscsx/CanUD_EURAFR_META_PRS_two_col.tsv",
  nicotine = "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/GWAS_sud1/prs_scores/nicotine_PRS_two_col.tsv",
  anySUD   = "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/GWAS_sud1/prs_scores/anySUD_PRS_two_col.tsv"
)

## Outputs
interval_prefix <- "env_surv_intervals"
lasso_prefix    <- "lasso_lambdaMin"
sel_prefix      <- "selected_lambdaMin"
cox_prefix      <- "cox_tv_lambdaMin"

## Forced adjustment covariates
forced_covars <- c("sex", "age", "site", paste0("PC", 1:20))

## Non-predictor columns
non_predictor_cols <- c(
  "FID", "IID", "src_subject_id", "eventname",
  "interview_age", "baseline_age", "time_months"
)

## -----------------------------
## Helpers
## -----------------------------

## Forward fill within subject (LOCF)
ffill_vec <- function(x) {
  if (all(is.na(x))) return(x)
  last <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      last <- x[i]
    } else if (!is.na(last)) {
      x[i] <- last
    }
  }
  x
}

## Load PRS two-column file
load_prs_two_col <- function(path, prs_name) {
  if (!file.exists(path)) {
    message_ts("PRS file not found for ", prs_name, ": ", path, " (skipping)")
    return(NULL)
  }
  message_ts("Reading PRS file for ", prs_name, " from: ", path)
  dt <- fread(path)

  if (!"FID" %in% names(dt) && "IID" %in% names(dt)) dt[, FID := IID]
  if (!"IID" %in% names(dt) && "FID" %in% names(dt)) dt[, IID := FID]
  if (!all(c("FID", "IID") %in% names(dt))) stop("PRS must contain FID/IID.")

  score_cols <- setdiff(names(dt), c("FID", "IID"))
  if (length(score_cols) != 1L) stop("PRS file must have exactly 1 score column.")

  setnames(dt, score_cols[1], prs_name)
  dt <- dt[, .(FID, IID, score = get(prs_name))]
  setnames(dt, "score", prs_name)

  message_ts("PRS ", prs_name, " loaded with ", nrow(dt), " rows.")
  dt
}

## Build start-stop intervals for one subject
build_intervals_for_subject <- function(dt_sub) {
  TE <- unique(dt_sub$time)
  EV <- unique(dt_sub$event)
  if (length(TE) != 1L || length(EV) != 1L) {
    stop("Multiple time/event values for same subject; check phenotype merge.")
  }
  TE <- TE[1]; EV <- EV[1]

  dt_sub <- dt_sub[!is.na(time_months) & time_months <= TE]
  if (nrow(dt_sub) == 0L) return(NULL)

  setorder(dt_sub, time_months)

  ## Ensure start at 0
  if (dt_sub$time_months[1] > 0) {
    bl <- copy(dt_sub[1])
    bl$time_months <- 0
    dt_sub <- rbindlist(list(bl, dt_sub), use.names = TRUE, fill = TRUE)
    setorder(dt_sub, time_months)
  }

  n <- nrow(dt_sub)
  start <- numeric(n)
  stop  <- numeric(n)
  ev_iv <- integer(n)

  for (i in seq_len(n)) {
    start[i] <- ifelse(i == 1, 0, dt_sub$time_months[i])
    stop[i]  <- if (i < n) min(dt_sub$time_months[i + 1], TE) else TE
    ev_iv[i] <- ifelse(stop[i] == TE && EV == 1 && i == n, 1L, 0L)
  }

  keep <- which(stop > start)
  if (length(keep) == 0L) return(NULL)

  dt_iv <- dt_sub[keep]
  dt_iv[, start := start[keep]]
  dt_iv[, stop := stop[keep]]
  dt_iv[, event_interval := ev_iv[keep]]
  dt_iv
}

## Coerce forced covars to sensible types
coerce_forced_covars <- function(dt) {
  if ("sex" %in% names(dt)) dt[, sex := as.factor(sex)]
  if ("site" %in% names(dt)) dt[, site := as.factor(site)]
  dt
}

## Subject coverage: fraction of subjects who ever have non-missing values
subject_coverage <- function(dt, cols, id_col = "IID") {
  n_subj <- uniqueN(dt[[id_col]])
  covg <- numeric(length(cols))
  names(covg) <- cols
  for (v in cols) {
    if (!v %in% names(dt)) { covg[v] <- 0; next }
    covg[v] <- uniqueN(dt[[id_col]][!is.na(dt[[v]])]) / n_subj
  }
  covg
}

## Add __obs indicator + median impute for modeling (no leakage; __obs carries “available now” info)
add_obs_and_impute <- function(df, cols, suffix="__obs") {
  for (v in cols) {
    if (!v %in% names(df)) next
    if (!is.numeric(df[[v]]) && !is.integer(df[[v]])) next

    obs <- !is.na(df[[v]])
    df[[paste0(v, suffix)]] <- as.integer(obs)

    med <- suppressWarnings(median(df[[v]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    x <- df[[v]]
    x[is.na(x)] <- med
    df[[v]] <- as.numeric(x)
  }
  df
}

## Prepare LASSO matrix
prepare_lasso_matrix <- function(tv_dt, outcome_name) {
  tv_df <- as.data.frame(tv_dt)

  drop_for_model <- c(
    non_predictor_cols,
    "time", "event",
    "start", "stop", "event_interval"
  )

  cand <- setdiff(names(tv_df), drop_for_model)
  cand <- setdiff(cand, forced_covars)  # never select forced covars

  ## Only numeric predictors for glmnet
  num_pred <- cand[sapply(tv_df[, cand, drop = FALSE], is.numeric)]

  if (length(num_pred) == 0L) return(NULL)

  ## Univariate p<=0.05 filter (no TOP_N cap)
  univar_file <- file.path(base_dir, paste0("cox_tv_univar_", outcome_name, "_results.tsv"))
  if (file.exists(univar_file)) {
    message_ts("Reading univariate results for pre-filter: ", univar_file)
    uv <- fread(univar_file)
    if (all(c("term", "p") %in% names(uv))) {
      uv <- uv[!is.na(p)]
      keep_terms <- unique(uv[p <= P_THRESH]$term)
      num_pred <- intersect(num_pred, keep_terms)
      message_ts("Pre-filter keeps ", length(num_pred), " numeric predictors (p<=0.05).")
    }
  } else {
    message_ts("Univariate file not found for ", outcome_name, "; using all numeric predictors.")
  }

  if (length(num_pred) == 0L) return(NULL)

  ## Subject coverage filter (keep predictors observed in >=50% subjects at any time)
  covg <- subject_coverage(tv_dt, num_pred, id_col = "IID")
  keep <- names(covg)[covg >= MIN_SUBJECT_COVERAGE]
  dropped <- setdiff(num_pred, keep)
  num_pred <- keep
  message_ts("Coverage>= ", MIN_SUBJECT_COVERAGE, " keeps ", length(num_pred),
             " predictors; dropped ", length(dropped), ".")

  if (length(num_pred) == 0L) return(NULL)

  ## Add __obs + impute for modeling
  tv_df <- add_obs_and_impute(tv_df, num_pred, suffix="__obs")

  ## Model matrix X includes both X and X__obs
  x_cols <- c(num_pred, paste0(num_pred, "__obs"))
  X <- as.matrix(tv_df[, x_cols, drop = FALSE])

  ## glmnet Cox does not support true start-stop counting-process risk sets.
  ## Use interval-end time (stop) with event at interval end as a selection approximation.
  y <- Surv(tv_df$stop, tv_df$event_interval)

  list(X = X, y = y, predictors = num_pred, x_cols = x_cols)
}

## -----------------------------
## Load covariates + PRS + LOCF
## -----------------------------
if (!file.exists(cov_file)) stop("Cannot find: ", cov_file)
message_ts("Reading covariates: ", cov_file)
cov_dt <- fread(cov_file)

if (!all(c("FID", "IID", "time_months") %in% names(cov_dt))) {
  stop("cov_timevarying_filtered.tsv must contain FID, IID, time_months.")
}

## Merge PRS
prs_list <- list()
for (nm in names(prs_paths)) {
  prs_dt <- load_prs_two_col(prs_paths[[nm]], paste0("PRS_", nm))
  if (!is.null(prs_dt)) prs_list[[length(prs_list) + 1L]] <- prs_dt
}
if (length(prs_list) > 0L) {
  message_ts("Merging ", length(prs_list), " PRS tables...")
  prs_all <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE), prs_list)
  message_ts("PRS merged dims: ", nrow(prs_all), " x ", ncol(prs_all))
  cov_dt <- merge(cov_dt, prs_all, by = c("FID", "IID"), all.x = TRUE)
  message_ts("Cov dims after PRS: ", nrow(cov_dt), " x ", ncol(cov_dt))
}

## LOCF within subject (safe; no back-fill)
message_ts("LOCF within subject...")
setkey(cov_dt, FID, IID, time_months)
cov_cols_for_ffill <- setdiff(
  names(cov_dt),
  c("FID", "IID", "src_subject_id", "eventname", "interview_age", "baseline_age", "time_months")
)
if (length(cov_cols_for_ffill) > 0L) {
  cov_dt[, (cov_cols_for_ffill) := lapply(.SD, ffill_vec),
         by = .(FID, IID),
         .SDcols = cov_cols_for_ffill]
  message_ts("LOCF done for ", length(cov_cols_for_ffill), " columns.")
}

## -----------------------------
## Main loop
## -----------------------------
for (outcome in names(phen_files)) {

  fpath <- phen_files[[outcome]]
  if (!file.exists(fpath)) {
    message_ts("[skip] phenotype missing: ", fpath)
    next
  }

  message_ts("==============================================")
  message_ts("Outcome: ", outcome)
  message_ts("Reading phenotype: ", fpath)
  ph <- fread(fpath)

  if (!all(c("FID", "IID", "time", "event") %in% names(ph))) {
    message_ts("[skip] phenotype missing required columns for ", outcome)
    next
  }

  dat <- merge(cov_dt, ph, by = c("FID", "IID"), all.y = TRUE, allow.cartesian = TRUE)
  dat <- dat[!is.na(time) & !is.na(event)]
  if (nrow(dat) == 0L) {
    message_ts("[skip] no rows after merge for ", outcome)
    next
  }

  setkey(dat, FID, IID)

  ## Build intervals
  tv_dt <- dat[, build_intervals_for_subject(.SD), by = .(FID, IID)]
  if (is.null(tv_dt) || nrow(tv_dt) == 0L) {
    message_ts("[skip] no intervals for ", outcome)
    next
  }

  tv_dt <- coerce_forced_covars(tv_dt)

  message_ts("Intervals: ", nrow(tv_dt),
             " | subjects: ", length(unique(tv_dt$IID)),
             " | events: ", sum(tv_dt$event_interval == 1, na.rm = TRUE))

  ## Save intervals
  interval_file <- file.path(base_dir, paste0(interval_prefix, "_", outcome, ".tsv"))
  message_ts("Writing intervals: ", interval_file)
  fwrite(tv_dt, interval_file, sep = "\t", quote = FALSE, na = "NA")

  ## LASSO selection (lambda.min, keep ALL nonzero)
  lasso_in <- prepare_lasso_matrix(tv_dt, outcome)
  if (is.null(lasso_in)) {
    message_ts("[warn] No usable numeric predictors for LASSO in ", outcome)
    next
  }

  X <- lasso_in$X
  y <- lasso_in$y
  base_predictors <- lasso_in$predictors   # original names (without __obs)
  x_cols <- lasso_in$x_cols               # includes __obs cols

  message_ts("Running cv.glmnet with ", ncol(X), " columns (X + X__obs)...")
  set.seed(SEED)
  cvfit <- tryCatch(
    cv.glmnet(x = X, y = y, family = "cox", alpha = 1,
              nfolds = N_FOLDS, standardize = GLMNET_STANDARDIZE,
              maxit = GLMNET_MAXIT),
    error = function(e) {
      message_ts("[fatal] cv.glmnet failed: ", e$message)
      NULL
    }
  )
  if (is.null(cvfit)) next

  lambda_use <- cvfit$lambda.min
  message_ts("lambda.min = ", lambda_use)

  coef_vec <- as.matrix(coef(cvfit$glmnet.fit, s = lambda_use))[, 1]
  coef_dt <- data.table(term = names(coef_vec), coef = as.numeric(coef_vec))
  coef_dt <- coef_dt[!is.na(coef)]
  coef_dt[, abscoef := abs(coef)]
  setorder(coef_dt, -abscoef)

  lasso_file <- file.path(base_dir, paste0(lasso_prefix, "_", outcome, "_coef.tsv"))
  message_ts("Writing LASSO coefficients: ", lasso_file)
  fwrite(coef_dt[, .(term, coef)], lasso_file, sep = "\t", quote = FALSE, na = "NA")

  selected_terms <- coef_dt[coef != 0]$term
  message_ts("Selected non-zero terms at lambda.min: ", length(selected_terms))

  sel_file <- file.path(base_dir, paste0(sel_prefix, "_", outcome, "_predictors.tsv"))
  fwrite(data.table(outcome = outcome, term = selected_terms),
         sel_file, sep = "\t", quote = FALSE, na = "NA")

  if (length(selected_terms) == 0L) {
    message_ts("[skip] no non-zero terms selected for ", outcome)
    next
  }

  ## Final CoxPH: forced adjustment covars + selected terms
  forced_in_data <- intersect(forced_covars, names(tv_dt))

  ## Build modeling data with only necessary columns
  keep_cols <- unique(c("FID", "IID", "start", "stop", "event_interval",
                       forced_in_data,
                       base_predictors))
  keep_cols <- keep_cols[keep_cols %in% names(tv_dt)]
  cox_dat <- tv_dt[, ..keep_cols]

  ## Add __obs + impute ONLY for base predictors that appear in selected terms (not for everything)
  ## Determine which base predictors are needed (selected term may be X or X__obs)
  needed_base <- unique(gsub("__obs$", "", selected_terms))
  needed_base <- intersect(needed_base, names(cox_dat))

  ## IMPORTANT: __obs computed BEFORE imputation
  cox_df <- as.data.frame(cox_dat)
  cox_df <- add_obs_and_impute(cox_df, needed_base, suffix="__obs")
  cox_dat <- as.data.table(cox_df)

  ## Coerce forced covars
  cox_dat <- coerce_forced_covars(cox_dat)

  ## Build RHS:
  ## - forced covars (except site; handled via strata)
  ## - selected_terms (X and/or X__obs)
  rhs_forced <- setdiff(forced_in_data, "site")
  rhs_sel <- selected_terms[selected_terms %in% names(cox_dat)]
  rhs <- unique(c(rhs_forced, rhs_sel))

  ## QC: drop any constant terms
  rhs2 <- c()
  for (v in rhs) {
    x <- cox_dat[[v]]
    if (is.numeric(x)) {
      if (is.finite(sd(x)) && sd(x) > 0) rhs2 <- c(rhs2, v)
    } else if (is.factor(x)) {
      if (nlevels(x) >= 2) rhs2 <- c(rhs2, v)
    } else {
      if (length(unique(x)) >= 2) rhs2 <- c(rhs2, v)
    }
  }
  rhs <- unique(rhs2)

  if (length(rhs) == 0L) {
    message_ts("[skip] no usable terms after QC for ", outcome)
    next
  }

  ## Final formula 
  has_site <- "site" %in% names(cox_dat)
  rhs_str <- paste(rhs, collapse = " + ")
  if (has_site) {
    fml <- as.formula(paste0("Surv(start, stop, event_interval) ~ strata(site) + ", rhs_str))
  } else {
    fml <- as.formula(paste0("Surv(start, stop, event_interval) ~ ", rhs_str))
  }

  ## Fit final CoxPH
  message_ts("Fitting final CoxPH with ", length(rhs), " terms (+ strata(site) if present)...")
  fit_final <- tryCatch(
    coxph(fml, data = cox_dat, ties = "efron", robust = TRUE, cluster = IID),
    error = function(e) {
      message_ts("[fatal] final Cox failed for ", outcome, ": ", e$message)
      NULL
    }
  )
  if (is.null(fit_final)) next

  sfit <- summary(fit_final)

  ## Stats
  n_subj <- length(unique(cox_dat$IID))
  n_events <- sum(cox_dat$event_interval == 1, na.rm = TRUE)
  n_intervals <- nrow(cox_dat)

  res_full <- data.table(
    outcome     = outcome,
    term        = rownames(sfit$coefficients),
    beta        = as.numeric(sfit$coefficients[, "coef"]),
    HR          = exp(as.numeric(sfit$coefficients[, "coef"])),
    robust_se   = as.numeric(sfit$coefficients[, "robust se"]),
    z           = as.numeric(sfit$coefficients[, "z"]),
    p           = as.numeric(sfit$coefficients[, "Pr(>|z|)"]),
    CI_lower    = exp(as.numeric(sfit$coefficients[, "coef"]) - 1.96 * as.numeric(sfit$coefficients[, "robust se"])),
    CI_upper    = exp(as.numeric(sfit$coefficients[, "coef"]) + 1.96 * as.numeric(sfit$coefficients[, "robust se"])),
    n_subjects  = n_subj,
    n_events    = n_events,
    n_intervals = n_intervals,
    lambda_min  = lambda_use
  )

  ## Write FULL (includes forced covars except site, which is stratified)
  out_full <- file.path(base_dir, paste0(cox_prefix, "_", outcome, "_FULL.tsv"))
  message_ts("Writing FULL Cox results: ", out_full)
  fwrite(res_full, out_full, sep = "\t", quote = FALSE, na = "NA")

  ## Write PREDICTORS-only (remove forced covars from display)
  forced_patterns <- c("^sex", "^age$", "^PC[0-9]+$")
  is_forced_term <- rep(FALSE, nrow(res_full))
  for (pat in forced_patterns) is_forced_term <- is_forced_term | grepl(pat, res_full$term)

  res_pred <- res_full[!is_forced_term]
  out_pred <- file.path(base_dir, paste0(cox_prefix, "_", outcome, "_PREDICTORS.tsv"))
  message_ts("Writing PREDICTORS-only Cox results (forced removed): ", out_pred)
  fwrite(res_pred, out_pred, sep = "\t", quote = FALSE, na = "NA")
}

message_ts("2_timevarying_survival_multivariate.R complete.")

