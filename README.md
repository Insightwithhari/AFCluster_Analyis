# Integrative Modeling Identifies a Dual-Function Helix in Full-Length MDM2

This repository contains the R scripts required to reproduce the molecular dynamics (MD) analysis and figures reported in the manuscript.

## Repository Contents

* **`01_global_analysis.R`**: Performs Principal Component Analysis (PCA) and Gaussian Mixture Model (GMM) clustering on the MD trajectory. Calculates global geometric metrics (average distances, RMSD, etc.) and generates Figures A-E (Density, PCA Landscape, Distance Boxplots, Violin Plot).
* **`02_residue_analysis.R`**: Performs detailed per-residue contact analysis based on the states defined in script 01. Generates gate residue frequency plots, contact fraction bar charts, and residue-residue contact heatmaps.

## Instructions for Use

1.  **Dependencies**: Ensure the following R packages are installed:
    * `bio3d`
    * `ggplot2`
    * `mclust`
    * `dplyr`
    * `tidyr`
    * `reshape2`
    * `RColorBrewer`
    * `viridis`
    * `pbapply`

2.  **Input Data**:
    The scripts require the following input files in the working directory:
    * `Refrence.pdb`: The reference PDB structure.
    * `Trajectory.dcd`: The MD trajectory file.
    *(Note: Users wishing to reproduce the analysis on their own data should ensure their file names match these or update the variables in the "User Inputs" section of the scripts.)*

3.  **Execution Order**:
    * Run `01_global_analysis.R` first. This will generate the `MDM2_Analysis_gmm_states.csv` file.
    * Run `02_residue_analysis.R` second. This requires the CSV output from the first script.

## Data Availability
The raw MD trajectories and reference structure files used in the original study are available from the corresponding author upon request due to file size constraints.
