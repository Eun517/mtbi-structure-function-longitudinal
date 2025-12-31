mtbi-structure-function-longitudinal
Analysis code for a longitudinal study investigating structural and functional brain network changes in patients with mild Traumatic Brain Injury (mTBI). This repository provides a complete pipeline to analyze how brain connectivity evolves over time and its relationship with cognitive recovery.

Overview
This repository contains MATLAB scripts for processing two types of brain connectivity data:

Functional Network Connectivity (FNC): Derived from resting-state fMRI (rs-fMRI) time-series.

Structural Connectivity (SC): Derived from Diffusion Tensor Imaging (DTI) streamline counts or tractography matrices.

The pipeline performs edge-wise statistical modeling to identify longitudinal changes and brain-behavior associations.

Key Analysis Features
1. Connectivity Computation
FNC (Functional): Calculates subject-level Fisher-Z transformed connectivity matrices from Schaefer-network time-series.

SC (Structural): Loads subject-wise DTI network matrices and extracts upper-triangular edge values for statistical analysis.

2. Statistical Modeling
Longitudinal LME: Employs Linear Mixed-Effects (LME) models to identify Group × Time interactions while controlling for Age and Sex:

Connectivity ~ Age + Sex + Group * Time + (1|Subject)

Post-hoc Analysis: Performs T-tests for group-wise (Patient vs. Control) and time-wise (T0 vs. T1) comparisons on significant edges.

Cognitive Association: Assesses the relationship between connectivity changes and cognitive scores (e.g., WCST, Digit Span) using LME and Pearson correlations.

Joint Modeling: Explores the triple interaction (Score ~ Connectivity * Group * Time) to identify group-specific connectivity-behavior relationships.

3. Multiple Comparison Correction
Uses Benjamini-Hochberg FDR correction to control for Type I errors across all network edges.

Project Structure
To run the analysis, organize your files as follows:

Plaintext

repo_root/
├── brain_comm_fnc_clean.m              # Main script for rs-fMRI FNC analysis
├── brain_comm_DTI_connectivity_clean.m # Main script for DTI Structural analysis
├── inputs/                             # Data directory (NOT included in repo)
│   ├── age_patients.txt                # Age list (nPat x 1)
│   ├── age_controls.txt                # Age list (nCon x 1)
│   ├── sex_patients.txt                # Sex coding (0/1)
│   ├── sex_controls.txt
│   ├── NetInfo.txt                     # List of Node/Network names
│   ├── data_T0.mat                     # rs-fMRI time-series at Baseline
│   ├── data_T1.mat                     # rs-fMRI time-series at Follow-up
│   ├── cognitive_patients.mat          # Cognitive scores [nPat x scores x time]
│   ├── cognitive_controls.mat          # Cognitive scores [nCon x scores x time]
│   └── dti_network/                    # Folder for DTI matrices (see below)
│       ├── Patients/
│       │   └── Sub01/
│       │       ├── T0/Sub01_T0_fdt_network_matrix.txt
│       │       └── T1/Sub01_T1_fdt_network_matrix.txt
│       └── Controls/
│           └── Con01/
│               ├── T0/Con01_T0_fdt_network_matrix.txt
│               └── T1/Con01_T1_fdt_network_matrix.txt
└── outputs/                            # Results (.mat files) saved here



External Function: This code requires fdr_bh.m for multiple comparison correction.
Download from: MATLAB Central File Exchange - fdr_bh


Usage
Prepare Data: Place your demographics and cognitive data in the inputs/ folder.

FNC Analysis: Run brain_comm_fnc_clean.m to analyze functional connectivity changes.

DTI Analysis: Run brain_comm_DTI_connectivity_clean.m to analyze structural connectivity changes.

Ensure your DTI text files follow the naming convention: <SubID>_<TimePoint>_fdt_network_matrix.txt.

Results: Check the outputs/ folder for generated summary tables and statistical results.

Citation
If you use this code in your research, please cite:

[Author Names], "Structural and Functional Brain Network Changes in mild Traumatic Brain Injury: A Longitudinal Study", Brain Communications, [Year].
