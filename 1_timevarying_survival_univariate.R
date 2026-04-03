#!/usr/bin/env Rscript

## For each outcome (alcohol, cannabis, nicotine, any):
## - Build time-varying data with intervals
## - For EACH predictor, fit a separate Cox model:
##     Surv(start, stop, event_interval) ~ predictor (+ covariates)
## - Use robust SEs clustered on IID
## - Output per-outcome result files in current directory:
##     cox_tv_univar_<outcome>_results.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(stringr)
})

message_ts <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = " "))
}

## -----------------------------
## Paths
## -----------------------------

base_dir <- getwd()
message_ts("Current working directory (base_dir): ", base_dir)

cov_file <- file.path(base_dir, "cov_timevarying_filtered.tsv")

phen_dir <- "/gpfs/group/ehlers/mwei/ABCD/init_phenos"
phen_files <- c(
  alcohol  = file.path(phen_dir, "alcohol_survival.tsv"),
  cannabis = file.path(phen_dir, "cannabis_survival.tsv"),
  nicotine = file.path(phen_dir, "nicotine_survival.tsv"),
  any      = file.path(phen_dir, "any_survival.tsv")
)

prs_paths <- list(
  PAU      = "Path",
  CanUD    = "Path",
  nicotine = "Path",
  anySUD   = "Path"
)

out_prefix <- "cox_tv_univar"

## -----------------------------
## 1. Load covariates + PRS
## -----------------------------

message_ts("Reading time-varying covariates from: ", cov_file)
cov_dt <- fread(cov_file)

needed_keys <- c("FID", "IID", "time_months")
if (!all(needed_keys %in% names(cov_dt))) {
  stop("cov_timevarying_filtered.tsv must contain FID, IID and time_months.")
}

## Detect PC columns (PC1–PC20 etc.)
pc_cols <- grep("^PC[0-9]+$", names(cov_dt), value = TRUE)
message_ts("Detected PC columns: ", paste(pc_cols, collapse = ", "))

## Load PRS and merge
load_prs_two_col <- function(path, prs_name) {
  if (!file.exists(path)) {
    message_ts("PRS file not found for ", prs_name, ": ", path, " (skipping this PRS)")
    return(NULL)
  }
  message_ts("Reading PRS file for ", prs_name, " from: ", path)
  dt <- fread(path)

  if (!"FID" %in% names(dt) && "IID" %in% names(dt)) {
    dt[, FID := IID]
  }
  if (!"IID" %in% names(dt) && "FID" %in% names(dt)) {
    dt[, IID := FID]
  }

  if (!all(c("FID", "IID") %in% names(dt))) {
    stop("PRS file for ", prs_name, " must contain FID and/or IID columns.")
  }

  score_cols <- setdiff(names(dt), c("FID", "IID"))
  if (length(score_cols) != 1L) {
    stop("PRS file for ", prs_name,
         " should have exactly one score column besides FID/IID; found: ",
         paste(score_cols, collapse = ", "))
  }

  setnames(dt, score_cols[1], prs_name)
  dt <- dt[, .(FID, IID, score = get(prs_name))]
  setnames(dt, "score", prs_name)

  message_ts("PRS ", prs_name, " loaded with ",
             nrow(dt), " rows; score column renamed to ", prs_name)

  return(dt)
}

prs_list <- list()
for (nm in names(prs_paths)) {
  prs_dt <- load_prs_two_col(prs_paths[[nm]], paste0("PRS_", nm))
  if (!is.null(prs_dt)) prs_list[[length(prs_list) + 1L]] <- prs_dt
}

if (length(prs_list) > 0L) {
  message_ts("Merging ", length(prs_list), " PRS tables with each other...")
  prs_all <- Reduce(function(x, y) {
    merge(x, y, by = c("FID", "IID"), all = TRUE)
  }, prs_list)

  message_ts("PRS merged table dimensions: ",
             nrow(prs_all), " x ", ncol(prs_all))

  message_ts("Merging PRS into covariate table by FID/IID...")
  cov_dt <- merge(cov_dt, prs_all, by = c("FID", "IID"), all.x = TRUE)

  message_ts("Covariate table dimensions after adding PRS: ",
             nrow(cov_dt), " x ", ncol(cov_dt))
} else {
  message_ts("No PRS files were successfully loaded; proceeding without PRS.")
}

## Non-predictor columns = IDs, time variables, and adjustment covariates
non_predictor_cols <- c(
  "FID", "IID", "src_subject_id", "eventname",
  "interview_age", "baseline_age", "time_months",
  "sex", "age", "site",
  pc_cols                  
)

## Covariates adjust for in each model
adjust_covars <- c("sex", "age", "site", pc_cols)
adjust_covars <- intersect(adjust_covars, names(cov_dt)) 
message_ts("Adjustment covariates in univariate models: ",
           paste(adjust_covars, collapse = ", "))

## -----------------------------
## 2. Helper: build start–stop intervals
## -----------------------------

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

  if (dt_sub$time_months[1] > 0) {
    bl <- dt_sub[1]
    bl$time_months <- 0
    dt_sub <- rbindlist(list(bl, dt_sub))
    setorder(dt_sub, time_months)
  }

  n <- nrow(dt_sub)

  start <- numeric(n)
  stop  <- numeric(n)
  ev_iv <- integer(n)

  for (i in seq_len(n)) {
    start[i] <- ifelse(i == 1, 0, dt_sub$time_months[i])

    if (i < n) {
      stop[i] <- min(dt_sub$time_months[i + 1], TE)
    } else {
      stop[i] <- TE
    }

    ev_iv[i] <- ifelse(stop[i] == TE && EV == 1 && i == n, 1L, 0L)
  }

  keep <- which(stop > start)
  if (length(keep) == 0L) return(NULL)

  dt_iv <- dt_sub[keep]
  dt_iv[, start := start[keep]]
  dt_iv[, stop  := stop[keep]]
  dt_iv[, event_interval := ev_iv[keep]]

  return(dt_iv)
}

## -----------------------------
## 3. Helper: univariate Cox per predictor
## -----------------------------

run_univar_cox <- function(df, outcome_name) {

  message_ts("[", outcome_name, "] Initial nrow in run_univar_cox: ", nrow(df))

  df <- df %>%
    filter(!is.na(start), !is.na(stop), !is.na(event_interval),
           stop > start)

  message_ts("[", outcome_name, "] After filtering valid intervals (start/stop/event_interval, stop>start): nrow=", nrow(df))

  if (nrow(df) == 0L) {
    message_ts("[", outcome_name, "] No usable intervals after filtering; returning NULL.")
    return(NULL)
  }

  n_subjects <- length(unique(df$IID))
  n_events   <- sum(df$event_interval == 1, na.rm = TRUE)
  message_ts("[", outcome_name, "] subjects=", n_subjects, ", events=", n_events)

  drop_for_model <- c(non_predictor_cols,
                      "time", "event", "start", "stop", "event_interval")

  cand_pred <- setdiff(names(df), drop_for_model)
  message_ts("[", outcome_name, "] Candidate predictors BEFORE dropping constant: ", length(cand_pred))

  if (length(cand_pred) == 0L) {
    message_ts("[", outcome_name, "] No predictors left after excluding non_predictor/time/event columns; returning NULL.")
    return(NULL)
  }

  const_pred <- sapply(df[, cand_pred, drop = FALSE], function(x) {
    x2 <- x[!is.na(x)]
    length(unique(x2)) <= 1L
  })
  if (any(const_pred)) {
    message_ts("[", outcome_name, "] Dropping ", sum(const_pred),
               " predictors that are all-NA or constant in this outcome dataset.")
    cand_pred <- cand_pred[!const_pred]
  }

  message_ts("[", outcome_name, "] Candidate predictors AFTER dropping constant: ", length(cand_pred))

  if (length(cand_pred) == 0L) {
    message_ts("[", outcome_name, "] No varying predictors left; returning NULL.")
    return(NULL)
  }

  adj_here <- intersect(adjust_covars, names(df))
  adj_part <- if (length(adj_here) > 0L) {
    paste("+", paste(adj_here, collapse = " + "))
  } else {
    ""
  }

  message_ts("[", outcome_name, "] Adjusting for covariates: ",
             ifelse(adj_part == "", "none", paste(adj_here, collapse = ", ")))

  res_list <- vector("list", length(cand_pred))
  names(res_list) <- cand_pred

  for (term in cand_pred) {
    fml_str <- paste0("Surv(start, stop, event_interval) ~ ",
                      term, " ", adj_part)
    fml <- as.formula(fml_str)
    message_ts("[", outcome_name, "] Fitting univariate Cox for term: ", term)

    fit <- tryCatch(
      coxph(fml, data = df, ties = "efron",
            robust = TRUE, cluster = df$IID),
      error = function(e) {
        message_ts("[", outcome_name, "] Cox fit FAILED for term ", term, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(fit)) next

    cm <- coef(fit)
    vc <- tryCatch(vcov(fit), error = function(e) NULL)
    if (is.null(vc)) {
      message_ts("[", outcome_name, "] vcov FAILED for term ", term, "; skipping.")
      next
    }

    if (!term %in% names(cm)) {
      message_ts("[", outcome_name, "] Term ", term, " not found in coef names; skipping.")
      next
    }

    beta <- cm[term]
    se   <- sqrt(diag(vc))[term]
    if (is.na(beta) || is.na(se)) {
      message_ts("[", outcome_name, "] NA beta or se for term ", term, "; skipping.")
      next
    }

    z     <- beta / se
    p     <- 2 * pnorm(-abs(z))
    HR    <- exp(beta)
    lower <- exp(beta - 1.96 * se)
    upper <- exp(beta + 1.96 * se)

    res_list[[term]] <- data.frame(
      outcome    = outcome_name,
      term       = term,
      beta       = beta,
      HR         = HR,
      CI_lower   = lower,
      CI_upper   = upper,
      z          = z,
      p          = p,
      n_subjects = n_subjects,
      n_events   = n_events,
      stringsAsFactors = FALSE
    )
  }

  valid <- !sapply(res_list, is.null)
  message_ts("[", outcome_name, "] Number of terms with non-null results: ", sum(valid))

  if (!any(valid)) {
    message_ts("[", outcome_name, "] All terms failed or produced NULL; returning NULL.")
    return(NULL)
  }

  res <- do.call(rbind, res_list[valid])
  if (is.null(res) || nrow(res) == 0L) {
    message_ts("[", outcome_name, "] Combined results table empty; returning NULL.")
    return(NULL)
  }
  rownames(res) <- NULL
  return(res)
}

## -----------------------------
## 4. Loop over outcomes
## -----------------------------

for (nm in names(phen_files)) {
  fpath <- phen_files[[nm]]
  if (!file.exists(fpath)) {
    message_ts("Phenotype file not found for ", nm, ": ", fpath, " | Skipping.")
    next
  }

  message_ts("Processing outcome (univariate): ", nm)
  message_ts("Reading phenotype file: ", fpath)
  ph <- fread(fpath)

  required_ph <- c("FID", "IID", "time", "event")
  if (!all(required_ph %in% names(ph))) {
    warning("Phenotype file ", fpath, " missing required columns: ",
            paste(setdiff(required_ph, names(ph)), collapse = ", "),
            "; skipping.")
    next
  }

  message_ts("[", nm, "] ph dim: ", nrow(ph), " x ", ncol(ph))

  dat <- merge(cov_dt, ph,
               by = c("FID", "IID"),
               all.y = TRUE,
               allow.cartesian = TRUE)

  message_ts("[", nm, "] Merged dim (cov + phenotype): ",
             nrow(dat), " x ", ncol(dat))

  dat <- dat[!is.na(time) & !is.na(event)]
  message_ts("[", nm, "] After filtering non-missing time/event: nrow = ", nrow(dat))
  if (nrow(dat) == 0L) {
    warning("No rows with non-missing time/event for ", nm)
    next
  }

  setkey(dat, FID, IID)
  tv_list <- dat[, {
    dt_sub <- .SD
    build_intervals_for_subject(dt_sub)
  }, by = .(FID, IID)]

  message_ts("[", nm, "] After building intervals: nrow(tv_list) = ",
             nrow(tv_list), ", unique subjects = ", length(unique(tv_list$IID)))

  if (nrow(tv_list) == 0L) {
    warning("No intervals created for ", nm, "; skipping.")
    next
  }

  res <- run_univar_cox(as.data.frame(tv_list), nm)
  if (is.null(res)) {
    message_ts("[", nm, "] run_univar_cox returned NULL; skipping write for this outcome.")
    next
  }

  out_file <- file.path(base_dir,
                        paste0(out_prefix, "_", nm, "_results.tsv"))
  message_ts("Writing UNIVARIATE Cox results for ", nm, " to: ", out_file)
  fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")
}

message_ts("1_timevarying_survival_univariate.R completed.")

