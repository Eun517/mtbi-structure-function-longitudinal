**mtbi-structure-function-longitudinal**
Analysis code for a longitudinal study investigating structural and functional brain network changes in patients with mild Traumatic Brain Injury (mTBI). 

This repository provides a complete pipeline to analyze connectivity evolution and its relationship with cognitive recovery.

📌 Overview
This repository contains MATLAB scripts for a multi-modal longitudinal analysis:

Functional Network Connectivity (FNC): Derived from resting-state fMRI (rs-fMRI).

Structural Connectivity (SC): Network matrices derived from DTI tractography.

DTI Diffusion Metrics: White matter integrity indices (FA, AD, RD, MD) extracted from network-related regions.

**Key Analysis Features**
**1. Functional & Structural Network Analysis**
Connectivity Computation: Fisher-Z transformation for FNC and upper-triangular edge extraction for DTI matrices.

Edge-wise LME Modeling: Identifies Group × Time interactions while controlling for Age and Sex:

Connectivity ~ Age + Sex + Group * Time + (1|Subject)

Post-hoc Tests: T-tests for group-wise (Patient vs. Control) and time-wise (T0 vs. T1) comparisons.

**2. DTI Metrics & Cognition Analysis**
Diffusion Indices: Analyzes FA (Fractional Anisotropy), AD (Axial), RD (Radial), and MD (Mean Diffusivity) for edges showing significant group effects in FNC.

Patient-Specific LME: Assesses how white matter integrity changes relate to cognitive recovery over time:

Cognition ~ DTI_Metric * Time + (1|Subject)

Brain-Behavior Correlation: Pearson correlations between diffusion metrics and cognitive scores (e.g., WCST, Digit Span) for both Patients and Controls at each time point.


**Project Structure**
To ensure the scripts run correctly, please organize your inputs/ directory as follows:

Plaintext

repo_root/
├── brain_comm_fnc_clean.m              # rs-fMRI FNC analysis script
├── brain_comm_DTI_connectivity_clean.m # DTI Network (SC) analysis script
├── brain_comm_DTI_metrics_analysis.m   # DTI metrics (FA/AD/RD/MD) & Cognition
├── inputs/                             # Data directory (NOT in repo)
│   ├── NetInfo.txt                     # Network node labels
│   ├── data_T0.mat / data_T1.mat       # fMRI time-series
│   ├── dti_network/                    # Subject-wise DTI matrices (.txt)
│   ├── WMint.mat                       # DTI metrics [nSub x Metric x Edge x Time]
│   ├── EdgeTable.mat                   # Master edge information table
│   ├── sig_edges_unc_G.mat             # Edges with significant group effects
│   ├── cognitive_patients.mat          # [nPat x Score x Time]
│   └── cognitive_controls.mat          # [nCon x Score x Time]
└── outputs/                            # Generated results


External Function: Requires fdr_bh.m (Download from MATLAB Central).

**Usage**
FNC/SC Analysis: Run brain_comm_fnc_clean.m and brain_comm_DTI_connectivity_clean.m to identify significant edges.
DTI Metrics Analysis: Ensure WMint.mat and sig_edges_unc_G.mat are prepared based on the FNC results, then run brain_comm_DTI_metrics_analysis.m.
