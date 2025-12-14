This repository contains analysis scripts used in the study.

Scripts
-------
01_state_definition_PCA_GMM.R
  - Performs PCA, PC1 orientation, Gaussian mixture modeling
  - Identifies conformational states and centroid structures

02_state_analysis_contacts.R
  - Computes per-residue heavy-atom distances
  - Calculates state-resolved contact fractions

03_plotting_figures.R
  - Generates publication-quality figures

Inputs
------
Reference.pdb
Trajectory.dcd

Software
--------
R >= 4.2
bio3d, mclust, ggplot2, dplyr, tidyr
