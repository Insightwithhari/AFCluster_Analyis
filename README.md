# Conformational State Analysis Scripts

This repository contains R scripts used to identify and analyze
conformational states from molecular dynamics simulations, as reported
in the accompanying JCIM manuscript.

## Workflow Overview

The analysis is performed in three steps:

1. Conformational state definition using PCA and Gaussian mixture modeling
2. State-resolved structural and contact analysis
3. Generation of publication-quality figures

## Repository Contents

### 01_state_definition_PCA_GMM.R
- Computes helix–pocket heavy-atom distances
- Performs global PCA on Cα atoms
- Orients PC1 for physical interpretability
- Identifies conformational states using GMM
- Writes centroid structures and session information

### 02_state_analysis_contacts.R
- Computes per-residue minimum heavy-atom distances
- Calculates state-resolved residue contact fractions

### 03_plotting_figures.R
- Generates all figures reported in the manuscript and Supporting Information

## Required Input Files

- Reference.pdb
- Trajectory.dcd

## Software Requirements

- R (≥ 4.2)
- bio3d, mclust, ggplot2, dplyr, tidyr, reshape2, viridis

Exact package versions are recorded in `sessionInfo.txt`.

## Execution Order

Run scripts in the following order:

1. 01_state_definition_PCA_GMM.R
2. 02_state_analysis_contacts.R
3. 03_plotting_figures.R

## Notes

Molecular dynamics trajectories are not deposited in accordance with
JCIM guidelines. Scripts are provided for transparency and reproducibility.
