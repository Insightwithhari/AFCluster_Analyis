# Conformational State Analysis from Molecular Dynamics Simulations

This repository contains R scripts used to identify, characterize, and visualize
conformational states from molecular dynamics (MD) simulations, as reported in the
accompanying manuscript.

The workflow combines global principal component analysis (PCA), Gaussian mixture
modeling (GMM), and state-resolved structural contact analysis.

---

## Workflow Overview

The analysis is performed in four sequential steps:

1. Definition of conformational states using PCA and Gaussian mixture modeling
2. State-resolved per-residue distance and contact analysis
3. Generation of publication-quality figures
4. State-resolved heavy-atom contact heatmap visualization

Each step is implemented as an independent script to ensure transparency,
reproducibility, and ease of review.

---

## Repository Contents

### `01_state_definition_PCA_GMM.R`
- Loads reference structure and MD trajectory
- Computes average heavy-atom distances between a regulatory helix and binding pocket
- Performs global PCA on Cα atoms
- Orients PC1 for physical interpretability using distance correlation
- Identifies conformational states using Gaussian mixture modeling
- Writes state assignments, centroid structures, and session information

**Outputs:**
- `global_PC1_GMM_<N>_states.csv`
- `Cluster_<N>_centroid.pdb`
- `sessionInfo.txt`

---

### `02_state_analysis_contacts.R`
- Computes per-residue minimum heavy-atom distances between helix and pocket
- Calculates state-resolved residue contact fractions using a fixed distance cutoff
- Uses GMM state assignments defined in Script 01

**Outputs:**
- `per_residue_contact_fractions.csv`

---

### `03_plotting_figures.R`
Generates publication-quality figures used in the main text and Supporting Information,
including:

- PC1 density distributions by state
- PCA conformational landscape (PC1 vs PC2)
- Helix–pocket distance distributions by state
- Correlation between PC1 and helix–pocket distance
- State-resolved per-residue contact fraction bar plots

**Outputs:**
- `Fig_PC1_density.png`
- `Fig_PCA_scatter.png`
- `Fig_distance_box.png`
- `Fig_PC1_vs_distance.png`
- `Fig_contact_fraction.png`

---

### `04_heavy_atom_contact_heatmap.R`
- Computes heavy-atom contacts between helix and pocket residues
- Uses fixed GMM state assignments from Script 01
- Calculates state-resolved contact frequencies
- Visualizes contacts as residue–residue heatmaps

**Outputs:**
- `Fig_heavy_atom_contact_heatmap.png`

---

## Required Input Files

The following input files are required but not included in the repository:

- `Reference.pdb` — reference structure used for trajectory alignment
- `Trajectory.dcd` — molecular dynamics trajectory

---

## Software Requirements

- R (≥ 4.2)
- Required R packages:
  - `bio3d`
  - `mclust`
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `reshape2`
  - `viridis`
  - `pbapply`

Exact package versions used are recorded in `sessionInfo.txt`.

---

## Execution Order

Scripts must be executed in the following order:

1. `01_state_definition_PCA_GMM.R`
2. `02_state_analysis_contacts.R`
3. `03_plotting_figures.R`
4. `04_heavy_atom_contact_heatmap.R`

---

## Notes

- All state definitions are derived exclusively from PCA and GMM in Script 01.
- Downstream analyses do not redefine or modify conformational states.

---
