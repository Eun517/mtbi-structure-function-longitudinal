# mtbi-structure-function-longitudinal

Analysis code for a longitudinal study investigating structural and functional brain network changes in patients with mild Traumatic Brain Injury (mTBI).

## Overview
This repository contains MATLAB scripts used to analyze longitudinal rs-fMRI data. The pipeline includes computing functional network connectivity (FNC), performing edge-wise statistical modeling, and assessing correlations with clinical/cognitive measures.

## Key Analysis Features
- **FNC Computation**: Calculates subject-level Fisher-Z transformed connectivity matrices from brain network time-series.
- **Longitudinal Modeling**: Employs Linear Mixed-Effects (LME) models to identify Group × Time interactions:
  `FNC ~ Age + Sex + Group * Time + (1|Subject)`
- **Post-hoc Analysis**: Performs T-tests for group-wise and time-wise comparisons on significant edges.
- **Cognitive Association**: Assesses the relationship between FNC changes and cognitive scores (e.g., WCST, Digit Span) using LME and Pearson correlations.
- **Multiple Comparison Correction**: Uses Benjamini-Hochberg FDR correction.



## Project Structure
To run the script `brain_comm_fnc_clean.m`, please organize your files as follows:

```text
repo_root/
├──brain_comm_fnc_clean.m       # Main analysis script
├── inputs/                  # Place your data files here (not included in repo)
│   ├── age_patients.txt
│   ├── age_controls.txt
│   ├── sex_patients.txt
│   ├── sex_controls.txt
│   ├── NetInfo.txt          # Node/Network names
│   ├── data_T0.mat          # Time-series data at Baseline
│   ├── data_T1.mat          # Time-series data at Follow-up
│   ├── cognitive_patients.mat
│   └── cognitive_controls.mat
└── outputs/                 # Generated results will be saved here


External Function: This code requires fdr_bh.m for multiple comparison correction.
Download from: MATLAB Central File Exchange - fdr_bh
