# ==============================================================================
# Script 01: Global Dynamics, GMM Clustering, and Conformational States
# ==============================================================================
# Description: 
#   1. Loads trajectory and reference structure.
#   2. Calculates global geometric metrics (Avg Heavy Atom Dist, Min C-alpha Dist).
#   3. Performs PCA and GMM clustering to identify metastable states.
#   4. Generates publication-quality figures (A-E) and summary statistics.
# ==============================================================================

# --- 1. Libraries ---
library(bio3d)
library(ggplot2)
library(mclust)       # For Gaussian Mixture Modeling
library(RColorBrewer)
library(dplyr)    

# --- 2. User Inputs & Parameters ---
# Ensure these filenames match exactly what is in your folder
topology_pdb  <- "Refrence.pdb"       
dcd_file      <- "Trajectory.dcd"     
helix_resnos  <- 190:205              # Residue range for Helix
pocket_resnos <- 21:109               # Residue range for Pocket
project_prefix <- "MDM2_Analysis"     # Prefix for all output files

# --- 3. Define Plotting Theme ---
theme_publication <- function(base_size=14) {
    theme_bw(base_size=base_size) +
        theme(
            plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(0,0,10,0)),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.text = element_text(size = rel(0.9), colour = "black"),
            legend.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_line(colour = "black")
        )
}

# --- 4. Load Data ---
cat("Loading structural data...\n")
if(!file.exists(topology_pdb)) stop(paste("File not found:", topology_pdb))
if(!file.exists(dcd_file)) stop(paste("File not found:", dcd_file))

topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

# --- 5. Calculate Distances ---
cat("Calculating geometric metrics...\n")

# A. Average Heavy Atom Distance (for Correlation/GMM)
heavy_helix  <- atom.select(topology, string = "noh", resno = helix_resnos)
heavy_pocket <- atom.select(topology, string = "noh", resno = pocket_resnos)

avg_heavy_distance <- function(i) {
    xyz <- traj_xyz[i, ]
    H <- matrix(xyz[heavy_helix$xyz],  ncol = 3, byrow = TRUE)
    P <- matrix(xyz[heavy_pocket$xyz], ncol = 3, byrow = TRUE)
    if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
    mean(bio3d::dist.xyz(H, P), na.rm = TRUE)
}

d_avg_heavy_all <- sapply(seq_len(nrow(traj_xyz)), avg_heavy_distance)
ok <- is.finite(d_avg_heavy_all)
d_avg_heavy <- d_avg_heavy_all[ok]

# B. Minimum C-alpha Distance (for Violin Plot)
ca_helix  <- atom.select(topology, elety = "CA", resno = helix_resnos)
ca_pocket <- atom.select(topology, elety = "CA", resno = pocket_resnos)

min_ca_distance <- function(i) {
    xyz <- traj_xyz[i, ]
    H <- matrix(xyz[ca_helix$xyz],  ncol = 3, byrow = TRUE)
    P <- matrix(xyz[ca_pocket$xyz], ncol = 3, byrow = TRUE)
    if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
    min(bio3d::dist.xyz(H, P), na.rm = TRUE)
}
d_min_ca_all <- sapply(seq_len(nrow(traj_xyz)), min_ca_distance)
d_min_ca <- d_min_ca_all[ok]

# --- 6. Global PCA & Correlation Check ---
cat("Performing PCA...\n")
sel_ca_global <- atom.select(topology, elety = "CA")
pc            <- pca.xyz(traj_xyz[, sel_ca_global$xyz], rm.gaps = TRUE)

# Invert PC1 if negatively correlated with distance (for consistent interpretation)
initial_corr <- suppressWarnings(cor.test(pc$z[ok, 1], d_avg_heavy, method = "spearman"))
if (!is.na(initial_corr$estimate) && initial_corr$estimate < 0) {
    cat("NOTE: Inverting PC1 for consistent interpretation.\n")
    pc$z <- -pc$z
}
final_corr_test <- suppressWarnings(cor.test(pc$z[ok, 1], d_avg_heavy, method = "spearman"))

# --- 7. GMM Clustering ---
cat("Running GMM Clustering on PC1...\n")
pc1_filtered <- pc$z[ok, 1]
gm_optimal <- Mclust(pc1_filtered)
num_states <- gm_optimal$G
state_final <- factor(paste0("Cluster ", gm_optimal$classification))

# Compile Results Dataframe
gmm_results <- data.frame(
    Frame    = which(ok),
    PC1      = pc1_filtered,
    PC2      = pc$z[ok, 2],
    Distance = d_avg_heavy, 
    State    = state_final
)

# IMPORTANT: Save this CSV as it is required for the next script
output_csv <- paste0(project_prefix, "_gmm_states.csv")
write.csv(gmm_results, output_csv, row.names = FALSE)
cat("GMM States saved to:", output_csv, "\n")

# --- 8. Centroid Extraction ---
centroid_frames <- gmm_results %>%
    group_by(State) %>%
    mutate(
        centroid_pc1 = mean(PC1),
        centroid_pc2 = mean(PC2),
        dist_to_centroid = sqrt((PC1 - centroid_pc1)^2 + (PC2 - centroid_pc2)^2)
    ) %>%
    slice_min(order_by = dist_to_centroid, n = 1) %>%
    ungroup()

# Write PDBs
for(i in 1:nrow(centroid_frames)) {
    frame_info <- centroid_frames[i,]
    cluster_name <- gsub(" ", "_", frame_info$State)
    file_name <- paste0(project_prefix, "_", cluster_name, "_centroid.pdb")
    write.pdb(pdb = topology, file = file_name, xyz = traj_xyz[frame_info$Frame, ])
}

# --- 9. Save Summary Statistics ---
sink(paste0(project_prefix, "_summary_stats.txt"))
cat("--- Summary of Analysis ---\n")
print(final_corr_test)
cat("\n--- State Occupancy ---\n")
print(table(gmm_results$State))
print(round(prop.table(table(gmm_results$State)) * 100, 1))
sink()

# --- 10. Generate Plots ---
cat("Generating Plots...\n")
cluster_colors <- brewer.pal(n = max(3, num_states), name = "Dark2")

# Plot A: Density
variance_explained <- (pc$L / sum(pc$L)) * 100
x_label_pc1 <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
plot_A <- ggplot(gmm_results, aes(x = PC1, fill = State)) +
    geom_density(alpha = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "State Distribution along PC1", x = x_label_pc1, y = "Probability Density") +
    theme_publication() +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5))

# Plot B: PCA Landscape
plot_B <- ggplot(gmm_results, aes(x = PC1, y = PC2, color = State)) +
    geom_point(alpha = 0.6, size = 1.5) +
    stat_ellipse(type = "norm", level = 0.80, linewidth = 0.8) +
    scale_color_manual(values = cluster_colors) +
    labs(title = "Conformational Landscape", x = x_label_pc1, y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
    theme_publication() +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5))

# Plot C: Distance Boxplot
plot_C <- ggplot(gmm_results, aes(x = State, y = Distance, fill = State)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "Helix-Pocket Distance by State", x = "Conformational State", y = "Avg. Heavy Atom Distance (Å)") +
    theme_publication() + theme(legend.position = "none")

# Plot D: Correlation
rho_val <- round(final_corr_test$estimate, 2)
p_val <- final_corr_test$p.value
p_txt <- if(p_val < 0.001) "p < 0.001" else paste("p =", round(p_val, 3))
plot_D <- ggplot(gmm_results, aes(x = PC1, y = Distance)) +
  geom_point(alpha = 0.4, size = 2.0, color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.7) +
  labs(title = "PC1 vs. Heavy Atom Distance", x = x_label_pc1, y = "Avg. Heavy Atom Distance (Å)") +
  theme_publication() +
  annotate("text", x = -Inf, y = Inf, label = paste("Rho =", rho_val, "\n", p_txt),
           hjust = -0.1, vjust = 1.2, size = 5, fontface = "italic")

# Plot E: Violin Plot (Min C-alpha)
df_dist_violin <- data.frame(D = d_min_ca)
p_violin <- ggplot(df_dist_violin, aes(x = "", y = D)) +
    geom_violin(fill = "skyblue", alpha = 0.6, color = "black") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.3, color = "red") +
    labs(title = "Distribution of Helix-Pocket Distances", x = "", y = "Minimum Cα–Cα distance (Å)") +
    theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save Plots
ggsave(paste0(project_prefix, "_FigA_Density.png"), plot_A, width = 7, height = 5, dpi = 600)
ggsave(paste0(project_prefix, "_FigB_PCA.png"), plot_B, width = 7, height = 5, dpi = 600)
ggsave(paste0(project_prefix, "_FigC_Box.png"), plot_C, width = 7, height = 5, dpi = 600)
ggsave(paste0(project_prefix, "_FigD_Corr.png"), plot_D, width = 7, height = 5, dpi = 600)
ggsave(paste0(project_prefix, "_FigE_Violin.png"), p_violin, width = 7, height = 5, dpi = 600)

cat("Script 01 Complete. Global plots saved.\n")
