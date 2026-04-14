# Time-Varying Environmental and Polygenic Predictors of Substance Use Initiation in Youth  
### A Survival and Causal Modeling Study in the ABCD Cohort

Mengman Wei\* and Qian Peng\* 


## Overview
This repository contains a lightweight, reproducible pipeline for:
1) building **time-varying (start–stop) interval data** from repeated ABCD visits,  
2) running **time-varying Cox survival models** (univariate screening and multivariable models), and  
3) estimating **causal effects** using **marginal structural models (MSM) with IPTW** and weighted pooled logistic regression.

The main outcomes are **time to first initiation** for:
- alcohol
- cannabis
- nicotine
- any substance initiation

> **Important note about data access:** ABCD data are controlled-access. This repo does **not** include ABCD data. To run the pipeline, you must have appropriate access and place the required input files in your environment.

---

## Repository contents
- `0_building_data.R`  
  Builds a merged longitudinal covariate dataset across visits (baseline + follow-ups) and computes **time since baseline (months)**.

- `1_timevarying_survival_univariate.R`  
  Creates **start–stop intervals** per participant per outcome and fits **one Cox model per predictor** (adjusted for standard covariates).

- `2_timevarying_survival_multivariate.R`  
  Fits a **multivariable time-varying survival model** per outcome using a two-stage approach:
  - screening/selection (e.g., LASSO)  
  - final **counting-process Cox model** on start–stop intervals  
  Includes safe handling of missing / late-measured variables (observation indicators + imputation) and robust SE clustered by participant.

- `3_msm_iptw.R`  
  Estimates causal effects using **MSM-IPTW**:
  - defines a binary exposure per interval (binary or median-split for continuous)
  - estimates **stabilized weights** from treatment models
  - applies positivity and weight-stability checks (clipping/truncation)
  - fits **weighted pooled logistic regression** with robust SE clustered by participant
  - outputs stable results and multiple-testing adjusted q-values

---

## Requirements
### Software
- R (4.2+)


---

## Citation

If you use or adapt this project's ideas and scripts, please cite the associated manuscript:

@article{wei2026timevarying,
  title={Time-Varying Environmental and Polygenic Predictors of Substance Use Initiation in Youth: A Survival and Causal Modeling Study in the ABCD Cohort},
  author={Wei, Mengman and Peng, Qian},
  journal={arXiv preprint arXiv:2604.07368},
  year={2026},
  doi={10.48550/arXiv.2604.07368}
}
