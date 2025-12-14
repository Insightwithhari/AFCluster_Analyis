# --- 1. Load Libraries ---
library(bio3d)
library(ggplot2)
library(mclust) # For Gaussian Mixture Modeling
library(RColorBrewer)
library(dplyr)    # For data manipulation, used to find centroids

# --- 2. Define Publication Theme ---
# A custom ggplot theme for creating publication-quality plots.
theme_publication <- function(base_size=14) {
    (theme_bw(base_size=base_size)
     + theme(
         plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(0,0,10,0)),
         axis.title = element_text(face = "bold", size = rel(1)),
         axis.text = element_text(size = rel(0.9), colour = "black"),
         legend.title = element_text(face = "bold"),
         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(colour = "black")
     )
    )
}

# --- 3. Inputs and Data Loading ---
topology_pdb <- "mdm25121025_57a20_512_1024_ref.pdb"
dcd_file     <- "mdm25121025_57a20_512_1024_traj.dcd"
helix_resnos <- 170:185
pocket_resnos <- 1:89

# Load the topology and trajectory files
topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

# --- 4. MODIFIED: Calculate Average Helix-Pocket Heavy Atom Distance ---
# Select all non-hydrogen ("heavy") atoms for the two regions
heavy_helix  <- atom.select(topology, string = "noh", resno = helix_resnos)
heavy_pocket <- atom.select(topology, string = "noh", resno = pocket_resnos)

# Function to calculate the AVERAGE distance between all heavy atoms for a single frame
avg_heavy_distance <- function(i) {
    xyz <- traj_xyz[i, ]
    H <- matrix(xyz[heavy_helix$xyz],  ncol = 3, byrow = TRUE)
    P <- matrix(xyz[heavy_pocket$xyz], ncol = 3, byrow = TRUE)
    if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
    mean(bio3d::dist.xyz(H, P), na.rm = TRUE)
}

# Apply the function across all frames of the trajectory
d_avg_heavy_all <- sapply(seq_len(nrow(traj_xyz)), avg_heavy_distance)
ok <- is.finite(d_avg_heavy_all)
d_avg_heavy <- d_avg_heavy_all[ok]

# --- 5. Perform Global PCA on C-alpha Atoms and Orient PC1 ---
# MODIFIED: Performing PCA on all C-alpha atoms
sel_ca_global <- atom.select(topology, elety = "CA")
pc            <- pca.xyz(traj_xyz[, sel_ca_global$xyz], rm.gaps = TRUE)

# First, check the initial correlation to see if PC1 needs to be inverted
initial_corr <- suppressWarnings(cor.test(pc$z[ok, 1], d_avg_heavy, method = "spearman"))
rho <- initial_corr$estimate

# If correlation is negative, invert the entire pc$z object for consistency
if (!is.na(rho) && rho < 0) {
    cat("NOTE: Initial correlation was negative. Inverting PC1 for consistent interpretation.\n")
    pc$z <- -pc$z
}

# Now that PC1 is correctly oriented, run the final correlation test for reporting
final_corr_test <- suppressWarnings(cor.test(pc$z[ok, 1], d_avg_heavy, method = "spearman"))

print("--- Final Correlation between PC1 and Avg. Heavy Atom Distance ---")
print(final_corr_test)
print("------------------------------------------------------------------")


# --- 6. GMM Analysis on PC1 with Optimal States ---
# Use the correctly oriented PC1 for clustering
pc1_filtered <- pc$z[ok, 1]

gm_optimal <- Mclust(pc1_filtered)
num_states <- gm_optimal$G
cat("Optimal number of states found by BIC:", num_states, "\n\n")

state_final <- factor(paste0("Cluster ", gm_optimal$classification))

gmm_results_optimal <- data.frame(
    Frame    = which(ok),
    PC1      = pc1_filtered,
    PC2      = pc$z[ok, 2],
    Distance = d_avg_heavy, # Using the new heavy atom distance
    State    = state_final
)

write.csv(gmm_results_optimal, paste0("5121024_global_pc1_gmm_", num_states, "_states.csv"), row.names = FALSE)

state_counts <- table(gmm_results_optimal$State)
state_percentages <- round(prop.table(state_counts) * 100, 1)

print(paste("--- State Occupancy (Global PCA, G=", num_states, ") ---", sep=""))
print(state_counts)
print(state_percentages)
print("-----------------------------------------")


# --- 7. Find and Write Centroid Structures for Each Cluster ---
centroid_frames <- gmm_results_optimal %>%
    group_by(State) %>%
    mutate(
        centroid_pc1 = mean(PC1),
        centroid_pc2 = mean(PC2),
        dist_to_centroid = sqrt((PC1 - centroid_pc1)^2 + (PC2 - centroid_pc2)^2)
    ) %>%
    slice_min(order_by = dist_to_centroid, n = 1) %>%
    ungroup()

print("--- Centroid Frames for Each Cluster ---")
print(select(centroid_frames, State, Frame, PC1, PC2))
print("----------------------------------------")

# Write each centroid frame to a PDB file
for(i in 1:nrow(centroid_frames)) {
    frame_info <- centroid_frames[i,]
    cluster_name <- gsub(" ", "_", frame_info$State)
    file_name <- paste0("5121024_", cluster_name, "_centroid.pdb")
    
    write.pdb(
        pdb = topology,
        file = file_name,
        xyz = traj_xyz[frame_info$Frame, ]
    )
    cat("Saved centroid structure for", as.character(frame_info$State), "to", file_name, "\n")
}

# --- NEW SECTION: Save Summary Statistics to a Text File ---
summary_file_name <- paste0("5121024_summary_statistics_", num_states, "_states.txt")

# Open a connection to a file
sink(summary_file_name)

# Print all the desired information, which will be redirected to the file
cat("--- Summary of Global PCA and GMM Analysis ---\n\n")

cat("--- Final Correlation between PC1 and Avg. Heavy Atom Distance ---\n")
print(final_corr_test)
cat("\n")

cat(paste("--- State Occupancy (Global PCA, G=", num_states, ") ---\n", sep=""))
cat("Cluster Counts:\n")
print(state_counts)
cat("\nCluster Percentages (%):\n")
print(state_percentages)
cat("\n")

cat("--- Centroid Frames for Each Cluster ---\n")
# Need to print the data.frame part of the tibble for clean output
print.data.frame(select(centroid_frames, State, Frame, PC1, PC2))
cat("\n")

# Close the connection, redirecting output back to the console
sink()

cat("\nSummary statistics saved to", summary_file_name, "\n")


# --- 8. Generate and Save Separate Publication-Quality Plots ---
variance_explained <- (pc$L / sum(pc$L)) * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)
x_label_pc1 <- paste0("PC1 (", pc1_var, "%)")
y_label_pc2 <- paste0("PC2 (", pc2_var, "%)")

cluster_colors <- brewer.pal(n = max(3, num_states), name = "Dark2")

plot_A <- ggplot(gmm_results_optimal, aes(x = PC1, fill = State)) +
    geom_density(alpha = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "State Distribution along PC1", x = x_label_pc1, y = "Probability Density") +
    theme_publication() +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5))

plot_B <- ggplot(gmm_results_optimal, aes(x = PC1, y = PC2, color = State)) +
    geom_point(alpha = 0.6, size = 1.5) +
    stat_ellipse(type = "norm", level = 0.80, linewidth = 0.8) +
    scale_color_manual(values = cluster_colors) +
    labs(title = "Conformational Landscape", x = x_label_pc1, y = y_label_pc2) +
    theme_publication() +
    theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5))

plot_C <- ggplot(gmm_results_optimal, aes(x = State, y = Distance, fill = State)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "Helix-Pocket Distance by State", x = "Conformational State", y = "Avg. Heavy Atom Distance (Å)") + # Label updated
    theme_publication() +
    theme(legend.position = "none")

# --- Create correlation label for Plot D ---
rho_val <- round(final_corr_test$estimate, 2)
p_val <- final_corr_test$p.value
p_val_text <- if (p_val < 2.2e-16) {
    "p < 2.2e-16"
} else if (p_val < 0.001) {
    paste("p =", formatC(p_val, format = "e", digits = 2))
} else {
    paste("p =", round(p_val, 3))
}
corr_label <- paste("Spearman's rho =", rho_val, "\n", p_val_text)


#plotD
plot_D <- ggplot(gmm_results_optimal, aes(x = PC1, y = Distance)) +
  geom_point(alpha = 0.4, size = 2.0, color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.7) +
  labs(title = "PC1 vs. Average Heavy Atom Distance",
       x = x_label_pc1, y = "Avg. Heavy Atom Distance (Å)") +
  theme_publication() +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  annotate("text",
           x = -Inf, y = Inf,
           label = corr_label,  # plain text label already created above
           hjust = -0.05, vjust = 1.2,
           size = 4.5, fontface = "italic")


ggsave(paste0("5121024_figure_A_density_global_", num_states, "state.png"), plot_A, width = 7, height = 5, dpi = 800)
ggsave(paste0("5121024_figure_B_PCA_scatter_global_", num_states, "state.png"), plot_B, width = 7, height = 5, dpi = 800)
ggsave(paste0("5121024_figure_C_distance_box_global_", num_states, "state.png"), plot_C, width = 7, height = 5, dpi = 800)
ggsave(paste0("5121024_figure_D_correlation_", num_states, "state.png"), plot_D, width = 7, height = 5, dpi = 800)


ggsave(paste0("5121024_figure_A_density_global_", num_states, "state.pdf"), plot_A, width = 7, height = 5, device = cairo_pdf)
ggsave(paste0("5121024_figure_B_PCA_scatter_global_", num_states, "state.pdf"), plot_B, width = 7, height = 5, device = cairo_pdf)
ggsave(paste0("5121024_figure_C_distance_box_global_", num_states, "state.pdf"), plot_C, width = 7, height = 5, device = cairo_pdf)
ggsave(paste0("5121024_figure_D_correlation_", num_states, "state.pdf"), plot_D, width = 7, height = 5, device = cairo_pdf)

print(paste("All plots for the Global PCA", num_states, "-state model saved as separate, high-resolution files."))












# --- 2. Define Publication Theme ---
# A custom ggplot theme for creating publication-quality plots.
theme_publication <- function(base_size=14) {
    (theme_bw(base_size=base_size)
     + theme(
         plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(0,0,10,0)),
         axis.title = element_text(face = "bold", size = rel(1)),
         axis.text = element_text(size = rel(0.9), colour = "black"),
         legend.title = element_text(face = "bold"),
         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(colour = "black")
     )
    )
}



# Load the topology and trajectory files
topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

# --- 4. Calculate Minimum Helix-Pocket C-alpha Distance ---
# Select C-alpha atoms for the two regions
ca_helix  <- atom.select(topology, elety = "CA", resno = helix_resnos)
ca_pocket <- atom.select(topology, elety = "CA", resno = pocket_resnos)

# Function to calculate the MINIMUM distance for a single frame
min_ca_distance <- function(i) {
    xyz <- traj_xyz[i, ]
    H <- matrix(xyz[ca_helix$xyz],  ncol = 3, byrow = TRUE)
    P <- matrix(xyz[ca_pocket$xyz], ncol = 3, byrow = TRUE)
    if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
    min(bio3d::dist.xyz(H, P), na.rm = TRUE)
}

# Apply the function across all frames of the trajectory
d_min_ca_all <- sapply(seq_len(nrow(traj_xyz)), min_ca_distance)
ok <- is.finite(d_min_ca_all)
d_min_ca <- d_min_ca_all[ok]


# --- 5. Generate and Save Publication-Quality Violin Plot ---
df_dist_violin <- data.frame(D = d_min_ca)

p_violin_pub <- ggplot(df_dist_violin, aes(x = "", y = D)) +
    geom_violin(fill = "skyblue", alpha = 0.6, color = "black") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.3, color = "red") +
    labs(title = "Distribution of Helix-Pocket Distances", x = "", y = "Minimum Cα–Cα distance (Å)") +
    theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save the plot in both PNG and PDF formats
ggsave("minimum_ca_distance_violin_plot.png", p_violin_pub, width = 7, height = 5, dpi = 800)
ggsave("minimum_ca_distance_violin_plot.pdf", p_violin_pub, width = 7, height = 5, device = cairo_pdf)

print("Violin plot for the minimum C-alpha distance has been saved.")











#perResiduesAnalysis

# --- 1. Load Libraries ---
library(bio3d)
library(ggplot2)
library(RColorBrewer)
library(tidyr) # For data reshaping (pivot_longer)

# --- 2. Define Publication Theme ---
# A custom ggplot theme for creating publication-quality plots.
theme_publication <- function(base_size=14) {
    (theme_bw(base_size=base_size)
     + theme(
         plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, margin = margin(0,0,10,0)),
         axis.title = element_text(face = "bold", size = rel(1)),
         axis.text = element_text(size = rel(0.9), colour = "black"),
         legend.title = element_text(face = "bold"),
         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(colour = "black")
     )
    )
}

# --- 3. Inputs and Data Loading ---
# IMPORTANT: Update gmm_file to the name of your specific GMM results CSV.
topology_pdb <- "mdm25121025_57a20_512_1024_ref.pdb"
dcd_file     <- "mdm25121025_57a20_512_1024_traj.dcd"
gmm_file     <- "5121024_global_pc1_gmm_3_states.csv" #<-- EXAMPLE FILENAME, PLEASE UPDATE
helix_resnos <- 170:185
pocket_resnos <- 1:89

# Load structural data
topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

# --- 4. Load GMM States ---
# Load the state definitions from the CSV file generated by the previous script
gmm_data <- read.csv(gmm_file)

# Ensure 'State' is a factor for proper ordering in plots
gmm_data$State <- as.factor(gmm_data$State)


# --- 5. Per-Residue Heavy-Atom Distance Calculation ---
cat("Calculating per-residue minimum distances... (this may take a moment)\n")

# Get a sorted list of unique residue IDs in the helix
helix_resids <- sort(unique(topology$atom$resno[topology$atom$resno %in% helix_resnos]))

# Prepare a matrix to store the results
perres_min <- matrix(NA_real_, nrow = length(helix_resids), ncol = nrow(traj_xyz), dimnames = list(helix_resids, NULL))

# Get atom indices for the pocket (we only need to do this once)
pocket_heavy_inds <- atom.select(topology, resno = pocket_resnos, string = "noh")$xyz

# Loop through each residue in the helix
for (ri in seq_along(helix_resids)) {
    rno <- helix_resids[ri]
    helix_heavy_inds  <- atom.select(topology, resno = rno, string = "noh")$xyz
    
    # Calculate the minimum distance from this one helix residue to the whole pocket for every frame
    perres_min[ri, ] <- vapply(seq_len(nrow(traj_xyz)), function(i){
        H <- matrix(traj_xyz[i, helix_heavy_inds], ncol=3, byrow=TRUE)
        P <- matrix(traj_xyz[i, pocket_heavy_inds], ncol=3, byrow=TRUE)
        if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
        min(bio3d::dist.xyz(H, P), na.rm=TRUE)
    }, numeric(1))
}
cat("Per-residue distance calculation complete.\n")


# --- 6. Detailed Analysis and Plotting ---

# 6a. Gate Residue Frequency Plot
# Find which residue is closest in each frame
min_res_each_frame <- helix_resids[max.col(-t(perres_min), ties.method="first")]
df_freq <- data.frame(table(min_res_each_frame))
colnames(df_freq) <- c("Residue", "Count")

plot_gate_freq <- ggplot(df_freq, aes(x = reorder(Residue, -Count), y = Count)) +
    geom_col(fill = "steelblue", color = "black") +
 scale_x_discrete(labels = function(x) as.character(as.numeric(as.character(x)) + 20)) +
    labs(title = "Gate Residue Frequencies", x = "Helix Residue", y = "Frames as Closest Residue") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_gate_residue_freq.png", plot_gate_freq, width = 7, height = 5, dpi = 600)

# 6b. State-Resolved Contact Fraction Plot
cutoff <- 5 # Contact cutoff in Angstroms
resnos <- as.integer(rownames(perres_min))

# Get the frame numbers that correspond to each state
state_indices <- split(gmm_data$Frame, gmm_data$State)

# Calculate the fraction of time each residue is in contact for each state
contact_list <- lapply(resnos, function(r) {
    v <- perres_min[as.character(r), ]
    # Create a list to hold results for this residue
    res_fractions <- list(Residue = r)
    # DYNAMICALLY calculate fraction for each state found in the data
    for (state_name in names(state_indices)) {
        res_fractions[[state_name]] <- mean(v[state_indices[[state_name]]] <= cutoff, na.rm = TRUE)
    }
    as.data.frame(res_fractions)
})

contact_df <- do.call(rbind, contact_list)
# Reshape the data for plotting with ggplot
contact_df_long <- pivot_longer(contact_df, cols = -Residue, names_to = "State", values_to = "Fraction")

# FIX: Replace the dot ('.') introduced by R back to a space to match the original factor levels
contact_df_long$State <- gsub("\\.", " ", contact_df_long$State)

contact_df_long$State <- factor(contact_df_long$State, levels = levels(gmm_data$State))

plot_contact_frac <- ggplot(contact_df_long, aes(x = factor(Residue), y = Fraction, fill = State)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = function(x) as.character(as.numeric(as.character(x)) + 20)) +
    labs(title = "Per-Residue Contact Fraction by State", x = "Helix Residue", y = paste("Contact Fraction (<", cutoff, "Å)")) +
    theme_publication() +
theme(legend.position = c(0.20, 0.99), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_contact_fractions.png", plot_contact_frac, width = 8, height = 5, dpi = 600)

# 6c. Generate PDB Colored by Contact Fraction in a specific state
# Safely select the first state that is actually PRESENT in the data (e.g., "Cluster 1")
target_state_original <- names(state_indices)[1] 
# FIX: Convert the state name to a format R uses for column names (e.g., "Cluster 1" -> "Cluster.1")
target_state_safe <- gsub(" ", ".", target_state_original)

pdb_ca <- trim.pdb(topology, atom.select(topology, "calpha"))
# Use the "safe" name to access the correct column in the data frame
contact_frac_target <- contact_df[[target_state_safe]][match(pdb_ca$atom$resno, contact_df$Residue)]
contact_frac_target[is.na(contact_frac_target)] <- 0 # Residues not in helix get 0
pdb_ca$atom$b <- contact_frac_target * 100 # Scale fraction to 0-100 for B-factor

# Use the original name for a clean output filename
file_name_pdb <- paste0("helix_contact_map_", gsub(" ", "_", target_state_original), ".pdb")
write.pdb(pdb_ca, file = file_name_pdb)

print("Detailed analysis complete. Plots and PDB file have been saved.")



# --- 1. Load Libraries ---
library(bio3d)
library(ggplot2)
library(reshape2)
library(pbapply) # For the progress bar
library(viridis) # For better color palettes

# --- 2. Define Publication Theme ---
theme_publication <- function(base_size=14) {
    theme_bw(base_size=base_size) +
        theme(
            plot.title = element_text(face="bold", size=rel(1.2), hjust=0.5, margin=margin(0,0,10,0)),
            axis.title = element_text(face="bold", size=rel(1)),
            axis.text = element_text(size=rel(0.9), colour="black"),
            legend.title = element_text(face="bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_line(colour = "black")
        )
}

# --- 3. Inputs and Parameters ---
topology_pdb  <- "mdm25121025_57a20_512_1024_ref.pdb"
dcd_file      <- "mdm25121025_57a20_512_1024_traj.dcd"
gmm_file      <- "5121024_global_pc1_gmm_3_states.csv" # <-- IMPORTANT: File with pre-computed states
helix_resnos  <- 170:185
pocket_resnos <- 1:89
contact_cutoff <- 5

# --- 4. Load Data and Pre-calculate Atom Indices ---
topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

cat("Pre-calculating atom indices for speed...\n")

H_res_atom_indices <- lapply(helix_resnos, function(r) {
  atom.select(topology, resno = r, heavy = TRUE)$atom
})
P_res_atom_indices <- lapply(pocket_resnos, function(r) {
  atom.select(topology, resno = r, heavy = TRUE)$atom
})

to_xyz_triplets <- function(atom_indices) {
  as.vector(rbind(3 * atom_indices - 2, 3 * atom_indices - 1, 3 * atom_indices))
}

H_res_xyz_indices <- lapply(H_res_atom_indices, to_xyz_triplets)
P_res_xyz_indices <- lapply(P_res_atom_indices, to_xyz_triplets)

# --- 5. Contact Calculation Function ---
calc_contacts_fast <- function(frame_xyz) {
    nH <- length(H_res_xyz_indices)
    nP <- length(P_res_xyz_indices)
    out <- logical(nH * nP)
    k <- 1
    
    for (h in seq_len(nH)) {
        H_xyz <- matrix(frame_xyz[H_res_xyz_indices[[h]]], ncol = 3, byrow = TRUE)
        if (any(!is.finite(H_xyz))) return(rep(NA, length(out)))
        
        for (p in seq_len(nP)) {
            P_xyz <- matrix(frame_xyz[P_res_xyz_indices[[p]]], ncol = 3, byrow = TRUE)
            if (any(!is.finite(P_xyz))) {
                out[k] <- NA
            } else {
                min_dist <- min(bio3d::dist.xyz(H_xyz, P_xyz))
                out[k] <- (min_dist <= contact_cutoff)
            }
            k <- k + 1
        }
    }
    return(out)
}

# --- 6. Calculate Contact Matrix Over All Frames ---
cat("Calculating heavy-atom contacts for all frames... (This may take a moment)\n")
n_frames <- nrow(traj_xyz)

contact_matrix <- t(pbsapply(1:n_frames, function(i) calc_contacts_fast(traj_xyz[i, ])))

colnames(contact_matrix) <- paste0("H", rep(helix_resnos, each=length(pocket_resnos)),
                                     "_P", rep(pocket_resnos, times=length(helix_resnos)))

# --- 7. Load Pre-Computed GMM States ---
cat("\nLoading pre-computed cluster data from:", gmm_file, "\n")
gmm_data <- read.csv(gmm_file)

# Ensure 'State' is a factor for consistent plotting order
gmm_data$State <- as.factor(gmm_data$State)

# The 'Frame' column in the GMM file tells us which frames are valid and were clustered
# The 'State' column provides the cluster label for each of those frames
cluster_labels <- gmm_data$State
frames_to_analyze <- gmm_data$Frame

cat(sprintf("Analysis will proceed using the %d frames defined in the GMM file.\n", length(frames_to_analyze)))

# Filter the full contact matrix to only include the frames from the GMM analysis
contact_matrix_filtered <- contact_matrix[frames_to_analyze, , drop=FALSE]


# --- 8. Calculate and Plot Contact Frequencies per Cluster ---
cluster_ids <- levels(cluster_labels)
contact_freq_list <- lapply(cluster_ids, function(cl) {
    # Find which rows in our filtered data correspond to the current cluster
    rows <- which(cluster_labels == cl)
    colMeans(contact_matrix_filtered[rows, , drop=FALSE], na.rm = TRUE)
})
names(contact_freq_list) <- cluster_ids

contact_freq_df <- do.call(cbind, contact_freq_list)
contact_freq_df <- as.data.frame(contact_freq_df)
contact_freq_df$Contact <- colnames(contact_matrix_filtered)
contact_freq_df$HelixRes  <- as.numeric(sub("H(\\d+)_P\\d+", "\\1", contact_freq_df$Contact))
contact_freq_df$PocketRes <- as.numeric(sub("H\\d+_P(\\d+)", "\\1", contact_freq_df$Contact))

contact_melt <- melt(contact_freq_df, id.vars=c("Contact", "HelixRes", "PocketRes"),
                     variable.name="Cluster", value.name="Frequency")

# Plotting
p_heatmap <- ggplot(contact_melt, aes(x=PocketRes, y=HelixRes, fill=Frequency)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(option = "magma", limits=c(0,1), name="Contact\nFrequency") +
  # If PocketRes/HelixRes are numeric:
  scale_x_continuous(labels = function(x) x + 20) +
  scale_y_continuous(labels = function(y) y + 20) +
    facet_wrap(~Cluster) +
    labs(title=paste("Heavy-Atom Contact Frequencies by Pre-Computed State (Cutoff =", contact_cutoff, "Å)"),
         x="Pocket Residue Number", y="Helix Residue Number") +
    theme_publication() +
    theme(
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.background = element_rect(fill = "grey85", color = "black"),
        strip.text = element_text(face = "bold", size = rel(1.1))
    )

ggsave("heavy_atom_contact_heatmap_from_GMM.png", p_heatmap, width=12, height=8, dpi=300)

cat("\nHeatmap based on pre-computed GMM states saved to 'heavy_atom_contact_heatmap_from_GMM.png'\n")


# For your four main plots:
ggsave(paste0("5121024_figure_A_density_global_", num_states, "state.svg"), plot_A, width = 7, height = 5)
ggsave(paste0("5121024_figure_B_PCA_scatter_global_", num_states, "state.svg"), plot_B, width = 7, height = 5)
ggsave(paste0("5121024_figure_C_distance_box_global_", num_states, "state.svg"), plot_C, width = 7, height = 5)
ggsave(paste0("5121024_figure_D_correlation_", num_states, "state.svg"), plot_D, width = 7, height = 5)

# For violin plot:
ggsave("minimum_ca_distance_violin_plot.svg", p_violin_pub, width = 7, height = 5)

# For gate residue frequency and contact fraction plots:
ggsave("figure_gate_residue_freq.svg", plot_gate_freq, width = 7, height = 5)
ggsave("figure_contact_fractions.svg", plot_contact_frac, width = 8, height = 5)
