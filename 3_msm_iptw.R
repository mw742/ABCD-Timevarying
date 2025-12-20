#!/usr/bin/env Rscript

## Robust MSM-IPTW pooled logistic regression over interval data.
## - Stabilized weights via propensity clipping (eps)
## - Weight truncation (trunc_q)
## - Strict weight QC gating (swp99_max, cap_max) before outcome model
## - Robust event handling
## - Reproducible manifest logging
##
##   Rscript 3_msm_iptw.R \
##     --tag=primary_eps0.01_q0.95_c10 \
##     --eps=0.01 --trunc_q=0.95 --max_confounders=10 \
##     --swp99_max=20 --cap_max=200 --seed=12345

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(sandwich)
  library(lmtest)
})

## -----------------------------
## CLI args (no extra deps)
## -----------------------------
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key, default=NULL) {
  hit <- grep(paste0("^--", key, "="), args, value=TRUE)
  if (length(hit)==0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

TAG             <- get_arg("tag", "main")
EPS             <- as.numeric(get_arg("eps", "0.01"))        # clip propensity to [eps, 1-eps]
TRUNC_Q         <- as.numeric(get_arg("trunc_q", "0.95"))    # weight truncation quantile
MAX_CONFOUNDERS <- as.integer(get_arg("max_confounders", "10"))
POS_LOW         <- as.numeric(get_arg("pos_low", "0.02"))
POS_HIGH        <- as.numeric(get_arg("pos_high", "0.98"))
SWP99_MAX       <- as.numeric(get_arg("swp99_max", "20"))
CAP_MAX         <- as.numeric(get_arg("cap_max", "200"))
SEED            <- as.integer(get_arg("seed", "12345"))

if (!is.finite(EPS) || EPS <= 0 || EPS >= 0.5) stop("Bad --eps: must be in (0,0.5). Got: ", EPS)
if (!is.finite(TRUNC_Q) || TRUNC_Q <= 0.5 || TRUNC_Q >= 1) stop("Bad --trunc_q: must be in (0.5,1). Got: ", TRUNC_Q)
if (!is.finite(POS_LOW) || !is.finite(POS_HIGH) || POS_LOW <= 0 || POS_HIGH >= 1 || POS_LOW >= POS_HIGH) {
  stop("Bad positivity thresholds: pos_low < pos_high within (0,1). Got: ", POS_LOW, " / ", POS_HIGH)
}
if (!is.finite(SWP99_MAX) || SWP99_MAX <= 0) stop("Bad --swp99_max: ", SWP99_MAX)
if (!is.finite(CAP_MAX) || CAP_MAX <= 0) stop("Bad --cap_max: ", CAP_MAX)

set.seed(SEED)

message_ts <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse=""))
`%||%` <- function(x, y) if (is.null(x)) y else x

## -----------------------------
## PATHS 
## -----------------------------
OUTCOMES <- c("alcohol","cannabis","nicotine","any")
INTERVAL_PREFIX <- "env_surv_intervals"
CAND_FILE <- "candidate_exposures_shared.tsv"

MAX_EXPOSURES_PER_OUTCOME <- 100

BASELINE_COVARS <- c("sex","age","site", paste0("PC",1:20))
INCLUDE_SITE <- TRUE

## output names include TAG
out_res_name <- function(outcome) paste0("msm_", TAG, "_iptw_pooledlogit_", outcome, "_results.tsv")
out_wqc_name <- function(outcome) paste0("msm_", TAG, "_iptw_weights_", outcome, "_qc.tsv")
out_manifest_name <- paste0("msm_", TAG, "_run_manifest.tsv")

## -----------------------------
## Manifest
## -----------------------------
manifest <- data.table::data.table(
  param = c("TAG","EPS","TRUNC_Q","MAX_CONFOUNDERS","POS_LOW","POS_HIGH","SWP99_MAX","CAP_MAX","SEED","R_version","data_table_version"),
  value = c(TAG, EPS, TRUNC_Q, MAX_CONFOUNDERS, POS_LOW, POS_HIGH, SWP99_MAX, CAP_MAX, SEED,
            paste(R.version$major, R.version$minor, sep="."),
            as.character(utils::packageVersion("data.table")))
)
data.table::fwrite(manifest, out_manifest_name, sep="\t")

## -----------------------------
## Helpers
## -----------------------------
append_fwrite <- function(dt, file) {
  if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
  data.table::fwrite(dt, file, sep="\t", quote=FALSE, na="NA",
                     append = file.exists(file),
                     col.names = !file.exists(file))
}

pick_exposures_for_outcome <- function(outcome, cand, nmax=100) {
  if ("outcome" %in% names(cand)) {
    dd <- cand[outcome == !!outcome & !is.na(term)]
    if (nrow(dd) == 0) return(character(0))
    if ("score" %in% names(dd)) setorder(dd, -score) else if ("p" %in% names(dd)) setorder(dd, p)
    terms <- unique(dd$term)
  } else {
    dd <- cand[!is.na(term)]
    if (all(c("n_outcomes","min_p") %in% names(dd))) setorder(dd, -n_outcomes, min_p)
    terms <- unique(dd$term)
  }
  if (length(terms) > nmax) terms <- terms[1:nmax]
  terms
}

add_obs_and_impute_numeric <- function(dt, v) {
  obs_col <- paste0(v, "__obs")
  x <- dt[[v]]
  dt[, (obs_col) := as.integer(!is.na(x))]
  med <- suppressWarnings(median(x, na.rm=TRUE))
  if (!is.finite(med)) med <- 0
  x[is.na(x)] <- med
  dt[, (v) := suppressWarnings(as.numeric(x))]
  dt
}

make_binary_treatment_adaptive <- function(dt, v) {
  obs_col <- paste0(v, "__obs")
  if (!obs_col %in% names(dt)) dt <- add_obs_and_impute_numeric(dt, v)

  x <- dt[[v]]
  obs <- (dt[[obs_col]] == 1)

  ux <- sort(unique(x[obs]))
  ux <- ux[is.finite(ux)]
  A <- paste0("A__", v)

  if (length(ux) == 0) {
    dt[, (A) := 0L]
    return(list(dt=dt, A=A, Aobs=obs_col, med=NA_real_, def="all_missing"))
  }

  if (all(ux %in% c(0,1))) {
    dt[, (A) := as.integer(obs & (x == 1))]
    return(list(dt=dt, A=A, Aobs=obs_col, med=1, def="binary_01"))
  }

  med_obs <- suppressWarnings(median(x[obs], na.rm=TRUE))
  if (!is.finite(med_obs)) med_obs <- suppressWarnings(median(x, na.rm=TRUE))
  if (!is.finite(med_obs)) med_obs <- 0

  dt[, (A) := as.integer(obs & (x > med_obs))]
  list(dt=dt, A=A, Aobs=obs_col, med=med_obs, def="median_split_obs")
}

truncate_weights <- function(w, q=0.95) {
  cap <- as.numeric(stats::quantile(w, probs=q, na.rm=TRUE))
  list(w=pmin(w, cap), cap=cap)
}

safe_ns <- function(x, df) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  # if too few finite values, return simple zeros
  if (sum(is.finite(x)) < (df + 2)) {
    m <- matrix(0, nrow=length(x), ncol=df)
  } else {
    m <- splines::ns(x, df=df)
  }
  as.matrix(m)
}

compute_sw <- function(dt, A, baseline_covars, Lt_covars, include_site=TRUE, eps=0.01) {
  dt <- copy(dt)
  setorder(dt, IID, start, stop)

  lagA <- paste0(A, "_lag1")
  dt[, (lagA) := shift(.SD[[1]], 1L, type="lag"), by=IID, .SDcols=A]
  dt[is.na(get(lagA)), (lagA) := 0L]

  dt[, t_stop := as.numeric(stop)]

  tt <- safe_ns(dt$t_stop, df=3)
  dt[, t_stop_spline1 := tt[,1]]
  dt[, t_stop_spline2 := tt[,2]]
  dt[, t_stop_spline3 := tt[,3]]

  base_pool <- baseline_covars
  if (include_site && "site" %in% names(dt)) base_pool <- c(base_pool, "site")
  base_in <- intersect(base_pool, names(dt))
  Lt_in   <- intersect(Lt_covars, names(dt))

  num_rhs <- paste(c(lagA, base_in, "t_stop_spline1","t_stop_spline2","t_stop_spline3"), collapse=" + ")
  den_rhs <- paste(c(lagA, base_in, "t_stop_spline1","t_stop_spline2","t_stop_spline3", Lt_in), collapse=" + ")

  f_num <- stats::as.formula(paste0(A, " ~ ", num_rhs))
  f_den <- stats::as.formula(paste0(A, " ~ ", den_rhs))

  m_num <- suppressWarnings(stats::glm(f_num, data=dt, family=stats::binomial(), control=stats::glm.control(maxit=200)))
  m_den <- suppressWarnings(stats::glm(f_den, data=dt, family=stats::binomial(), control=stats::glm.control(maxit=200)))

  p_num <- stats::predict(m_num, type="response")
  p_den <- stats::predict(m_den, type="response")

  # stabilize: clip propensities
  p_num <- pmin(pmax(p_num, eps), 1-eps)
  p_den <- pmin(pmax(p_den, eps), 1-eps)

  Avals <- dt[[A]]
  contrib <- ifelse(Avals==1, p_num/p_den, (1-p_num)/(1-p_den))

  dt[, contrib := as.numeric(contrib)]
  dt[, sw := cumprod(contrib), by=IID]
  dt[, contrib := NULL]
  dt[, t_stop := NULL]

  dt
}

fit_weighted_hazard <- function(dt, A, baseline_covars, extra_adj=NULL, include_site=TRUE) {
  dt <- copy(dt)

  # robust event
  dt[, event_interval := suppressWarnings(as.integer(as.numeric(event_interval)))]
  dt[is.na(event_interval), event_interval := 0L]
  dt[, y := as.integer(event_interval == 1L)]

  dt[, stop := as.numeric(stop)]
  dt <- dt[is.finite(stop) & is.finite(sw) & sw > 0]
  if (nrow(dt) == 0) return(NULL)

  dt[, t_stop := stop]
  tt <- safe_ns(dt$t_stop, df=4)
  dt[, t1 := tt[,1]]; dt[, t2 := tt[,2]]; dt[, t3 := tt[,3]]; dt[, t4 := tt[,4]]

  base_pool <- baseline_covars
  if (include_site && "site" %in% names(dt)) base_pool <- c(base_pool, "site")
  base_in  <- intersect(base_pool, names(dt))
  extra_in <- intersect(extra_adj %||% character(0), names(dt))

  rhs <- paste(c(A, base_in, extra_in, "t1","t2","t3","t4"), collapse=" + ")
  f <- stats::as.formula(paste0("y ~ ", rhs))

  m <- stats::glm(f, data=dt, family=stats::binomial(), weights=sw, control=stats::glm.control(maxit=200))
  V <- sandwich::vcovCL(m, cluster = dt$IID)
  ct <- lmtest::coeftest(m, vcov.=V)

  if (!A %in% rownames(ct)) return(NULL)
  est <- ct[A, "Estimate"]; se <- ct[A, "Std. Error"]
  z <- ct[A, "z value"]; p <- ct[A, "Pr(>|z|)"]

  if (!is.finite(est) || !is.finite(se) || se <= 0) return(NULL)

  data.table::data.table(
    term=A, beta=est, OR=exp(est), robust_se=se, z=z, p=p,
    CI_lower=exp(est - 1.96*se), CI_upper=exp(est + 1.96*se),
    n_subjects=data.table::uniqueN(dt$IID),
    n_events=sum(dt$y==1, na.rm=TRUE),
    n_intervals=nrow(dt)
  )
}

## -----------------------------
## Load candidates
## -----------------------------
if (!file.exists(CAND_FILE)) stop("Cannot find candidate file: ", CAND_FILE)
cand <- data.table::fread(CAND_FILE)
if (!"term" %in% names(cand)) stop("Candidate file must contain column 'term'.")

## -----------------------------
## Main loop
## -----------------------------
for (outcome in OUTCOMES) {
  message_ts("==============================================")
  message_ts("Outcome: ", outcome)

  int_file <- paste0(INTERVAL_PREFIX, "_", outcome, ".tsv")
  if (!file.exists(int_file)) { message_ts("[skip] missing intervals: ", int_file); next }

  dt <- data.table::fread(int_file)

  # required columns
  req <- c("IID","start","stop","event_interval")
  if (!all(req %in% names(dt))) stop("Intervals file lacks required columns: ", int_file)

  # types
  dt[, start := as.numeric(start)]
  dt[, stop  := as.numeric(stop)]
  if ("sex" %in% names(dt)) dt[, sex := as.factor(sex)]
  if ("site" %in% names(dt)) dt[, site := as.factor(site)]

  exposures <- pick_exposures_for_outcome(outcome, cand, nmax=MAX_EXPOSURES_PER_OUTCOME)
  exposures <- exposures[exposures %in% names(dt)]
  message_ts("Exposures to test (present in data): ", length(exposures))
  if (length(exposures)==0) next

  out_res <- out_res_name(outcome)
  out_wqc <- out_wqc_name(outcome)

  # overwrite on each run
  if (file.exists(out_res)) file.remove(out_res)
  if (file.exists(out_wqc)) file.remove(out_wqc)

  n_try <- 0; n_ok <- 0; n_posskip <- 0; n_swfail <- 0; n_outfail <- 0; n_wunstable <- 0

  for (v in exposures) {
    if (!(is.numeric(dt[[v]]) || is.integer(dt[[v]]))) next
    n_try <- n_try + 1
    message_ts("---- Exposure: ", v)

    tmp <- tryCatch(make_binary_treatment_adaptive(copy(dt), v),
                    error=function(e){ message_ts("[skip] make_treatment failed: ", e$message); NULL })
    if (is.null(tmp)) next
    dtx <- tmp$dt; A <- tmp$A; Aobs <- tmp$Aobs

    prev <- mean(dtx[[A]]==1, na.rm=TRUE)
    if (!is.finite(prev) || prev < POS_LOW || prev > POS_HIGH) {
      n_posskip <- n_posskip + 1
      message_ts("[skip] poor positivity, prev=", round(prev,4), " def=", tmp$def)
      next
    }

    # time-varying confounders = other exposures, limited
    Lt <- setdiff(exposures, v)
    Lt <- Lt[Lt %in% names(dtx)]
    if (length(Lt) > MAX_CONFOUNDERS) Lt <- Lt[1:MAX_CONFOUNDERS]

    for (cc in Lt) {
      if ((is.numeric(dtx[[cc]]) || is.integer(dtx[[cc]]))) dtx <- add_obs_and_impute_numeric(dtx, cc)
    }
    Lt2 <- unique(c(Lt, Aobs))

    dsw <- tryCatch(
      compute_sw(dtx, A=A, baseline_covars=BASELINE_COVARS, Lt_covars=Lt2, include_site=INCLUDE_SITE, eps=EPS),
      error=function(e){ n_swfail <<- n_swfail + 1; message_ts("[skip] SW failed: ", e$message); NULL }
    )
    if (is.null(dsw)) next

    tr <- truncate_weights(dsw$sw, q=TRUNC_Q)
    dsw[, sw := tr$w]
    cap <- tr$cap

    wqc <- data.table::data.table(
      outcome=outcome, exposure=v, A_name=A, A_def=tmp$def, A_median=tmp$med,
      prev=prev, EPS=EPS, TRUNC_Q=TRUNC_Q, cap=cap,
      sw_mean=mean(dsw$sw, na.rm=TRUE), sw_sd=sd(dsw$sw, na.rm=TRUE),
      sw_p95=as.numeric(stats::quantile(dsw$sw, 0.95, na.rm=TRUE)),
      sw_p99=as.numeric(stats::quantile(dsw$sw, 0.99, na.rm=TRUE)),
      sw_max=max(dsw$sw, na.rm=TRUE)
    )
    append_fwrite(wqc, out_wqc)

    # strict gating estimates
    if (is.finite(wqc$sw_p99) && is.finite(wqc$cap) && (wqc$sw_p99 > SWP99_MAX || wqc$cap > CAP_MAX)) {
      n_wunstable <- n_wunstable + 1
      message_ts("[skip] unstable weights: sw_p99=", signif(wqc$sw_p99,4), " cap=", signif(wqc$cap,4))
      next
    }

    r <- tryCatch(
      fit_weighted_hazard(dsw, A=A, baseline_covars=BASELINE_COVARS, extra_adj=c(Aobs), include_site=INCLUDE_SITE),
      error=function(e){ n_outfail <<- n_outfail + 1; message_ts("[skip] outcome model failed: ", e$message); NULL }
    )
    if (is.null(r)) next

    r[, `:=`(outcome=outcome, exposure=v, A_name=A, A_def=tmp$def, A_median=tmp$med,
             prev=prev, EPS=EPS, TRUNC_Q=TRUNC_Q, cap=cap)]
    append_fwrite(r, out_res)
    n_ok <- n_ok + 1
  }

  message_ts("OK results file: ", out_res)
  message_ts("OK weights file: ", out_wqc)
  message_ts("Summary ", outcome, ": tried=", n_try, " ok=", n_ok,
             " pos_skip=", n_posskip, " sw_fail=", n_swfail, " out_fail=", n_outfail,
             " w_unstable_skip=", n_wunstable)
}

message_ts("DONE.")

