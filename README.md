# Time-Varying Environmental and Polygenic Predictors of Substance Use Initiation in Youth  
### A Survival and Causal Modeling Study in the ABCD Cohort

Mengman Wei and Qian Peng


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

If you use this work, please cite:

Wei, M., & Peng, Q. (2026). Time-Varying Environmental and Polygenic Predictors of Substance Use Initiation in Youth. arXiv:2604.07368. 
https://doi.org/10.48550/arXiv.2604.07368.


---

## Data access and responsible use

This repository provides analysis code and does not include participant-level ABCD Study data. Users must obtain authorized access to the required data through the appropriate official provider before running the analyses. Data use must comply with the applicable data-use agreements, approved research purposes, and institutional requirements.

The software license applies to the code and does not grant access to ABCD Study data or permission to redistribute external datasets or restricted derived outputs.

Generated files may contain participant identifiers, longitudinal measurements, genetic scores, substance-use outcomes, or other sensitive information. Do not upload these files to public repositories or share them with unauthorized individuals. Any sharing must comply with the applicable data-use terms and required approvals.

For reproducibility, document the cohort, data release/version, accession or DOI where available, preprocessing steps, software versions, and exact code commit used. Include acknowledgments and citations required by the data provider.

## Version status and updates

The scripts currently posted may reflect an earlier development version. Later revisions may exist outside this repository, and the correspondence between the posted scripts and the code used for the reported analyses has not yet been verified.

The review conducted to date involved source-code and repository-history inspection. It did not include end-to-end execution or independent reproduction of the associated study results. Identified issues apply to the reviewed repository version; their relevance to later versions and their effects on the reported results have not been established. Undetected errors or limitations may remain.

The posted scripts will be checked against available later versions and analysis records, with corrections made where needed. Updates will be documented in the commit history or release notes, including the validation performed and any known effects on outputs.

Users should record the exact commit used and independently assess the code’s suitability for their intended analyses.



