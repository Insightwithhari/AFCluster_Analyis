# ==============================================================================
# Script 02: Per-Residue Contacts and Molecular Interactions
# ==============================================================================
# Description:
#   1. Loads GMM states generated in Script 01.
#   2. Calculates per-residue minimum distances and contact frequencies.
#   3. Generates Gate Residue, Contact Fraction, and Contact Heatmap plots.
#   4. Maps contact data onto PDB B-factors.
# ==============================================================================

# --- 1. Libraries ---
library(bio3d)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(reshape2)
library(pbapply) 
library(viridis)

# --- 2. User Inputs ---
# MUST MATCH SCRIPT 01
topology_pdb  <- "Refrence.pdb"
dcd_file      <- "Trajectory.dcd"
# Load the output from Script 1
gmm_file      <- "MDM2_Analysis_gmm_states.csv" 
helix_resnos  <- 190:205
pocket_resnos <- 21:109
contact_cutoff <- 5  # Angstroms

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
cat("Loading Data and GMM States...\n")
if(!file.exists(topology_pdb)) stop(paste("File not found:", topology_pdb))
if(!file.exists(dcd_file)) stop(paste("File not found:", dcd_file))
if(!file.exists(gmm_file)) stop(paste("GMM file not found. Please run Script 01 first."))

topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)
gmm_data <- read.csv(gmm_file)
gmm_data$State <- as.factor(gmm_data$State) # Ensure factor

# --- 5. Per-Residue Minimum Distance Analysis ---
cat("Calculating per-residue minimum distances...\n")
helix_resids <- sort(unique(topology$atom$resno[topology$atom$resno %in% helix_resnos]))
pocket_heavy_inds <- atom.select(topology, resno = pocket_resnos, string = "noh")$xyz

perres_min <- matrix(NA_real_, nrow = length(helix_resids), ncol = nrow(traj_xyz), 
                     dimnames = list(helix_resids, NULL))

# Calculation Loop
for (ri in seq_along(helix_resids)) {
    rno <- helix_resids[ri]
    helix_heavy_inds  <- atom.select(topology, resno = rno, string = "noh")$xyz
    
    perres_min[ri, ] <- vapply(seq_len(nrow(traj_xyz)), function(i){
        H <- matrix(traj_xyz[i, helix_heavy_inds], ncol=3, byrow=TRUE)
        P <- matrix(traj_xyz[i, pocket_heavy_inds], ncol=3, byrow=TRUE)
        if (any(!is.finite(H)) || any(!is.finite(P))) return(NA_real_)
        min(bio3d::dist.xyz(H, P), na.rm=TRUE)
    }, numeric(1))
}

# --- 6. Plot: Gate Residue Frequencies ---
min_res_each_frame <- helix_resids[max.col(-t(perres_min), ties.method="first")]
df_freq <- data.frame(table(min_res_each_frame))
colnames(df_freq) <- c("Residue", "Count")

plot_gate <- ggplot(df_freq, aes(x = reorder(Residue, -Count), y = Count)) +
    geom_col(fill = "steelblue", color = "black") +
    # REMOVED +20 OFFSET
    labs(title = "Gate Residue Frequencies", x = "Helix Residue", y = "Frames as Closest Residue") +
    theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figure_gate_residue_freq.png", plot_gate, width = 7, height = 5, dpi = 600)

# --- 7. Plot: Contact Fractions by State ---
state_indices <- split(gmm_data$Frame, gmm_data$State)
resnos <- as.integer(rownames(perres_min))

contact_list <- lapply(resnos, function(r) {
    v <- perres_min[as.character(r), ]
    res_fractions <- list(Residue = r)
    for (state_name in names(state_indices)) {
        res_fractions[[state_name]] <- mean(v[state_indices[[state_name]]] <= contact_cutoff, na.rm = TRUE)
    }
    as.data.frame(res_fractions)
})

contact_df <- do.call(rbind, contact_list)
contact_long <- pivot_longer(contact_df, cols = -Residue, names_to = "State", values_to = "Fraction")
contact_long$State <- gsub("\\.", " ", contact_long$State) # Fix naming
contact_long$State <- factor(contact_long$State, levels = levels(gmm_data$State))

plot_frac <- ggplot(contact_long, aes(x = factor(Residue), y = Fraction, fill = State)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Dark2") +
    # REMOVED +20 OFFSET
    labs(title = "Per-Residue Contact Fraction", x = "Helix Residue", y = "Contact Fraction") +
    theme_publication() +
    theme(legend.position = c(0.20, 0.99), legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth=0.5))

ggsave("figure_contact_fractions.png", plot_frac, width = 8, height = 5, dpi = 600)

# --- 8. Heatmap Analysis (Heavy Atoms) ---
cat("Calculating full heavy-atom contact matrix (High Compute)...\n")
# Helper: pre-calculate indices for speed
H_indices <- lapply(helix_resnos, function(r) atom.select(topology, resno=r, heavy=TRUE)$atom)
P_indices <- lapply(pocket_resnos, function(r) atom.select(topology, resno=r, heavy=TRUE)$atom)
to_xyz <- function(idx) as.vector(rbind(3*idx-2, 3*idx-1, 3*idx))
H_xyz_idx <- lapply(H_indices, to_xyz)
P_xyz_idx <- lapply(P_indices, to_xyz)

calc_contacts_fast <- function(frame_xyz) {
    nH <- length(H_xyz_idx); nP <- length(P_xyz_idx)
    out <- logical(nH * nP); k <- 1
    for (h in seq_len(nH)) {
        H_xyz <- matrix(frame_xyz[H_xyz_idx[[h]]], ncol = 3, byrow = TRUE)
        if (any(!is.finite(H_xyz))) return(rep(NA, length(out)))
        for (p in seq_len(nP)) {
            P_xyz <- matrix(frame_xyz[P_xyz_idx[[p]]], ncol = 3, byrow = TRUE)
            min_dist <- min(bio3d::dist.xyz(H_xyz, P_xyz))
            out[k] <- (min_dist <= contact_cutoff)
            k <- k + 1
        }
    }
    return(out)
}

# Run calculation on frames identified in GMM
frames_to_ana <- gmm_data$Frame
contact_matrix <- t(pbsapply(frames_to_ana, function(i) calc_contacts_fast(traj_xyz[i, ])))
colnames(contact_matrix) <- paste0("H", rep(helix_resnos, each=length(pocket_resnos)),
                                   "_P", rep(pocket_resnos, times=length(helix_resnos)))

# Aggregate by Cluster
cluster_ids <- levels(gmm_data$State)
contact_freq_list <- lapply(cluster_ids, function(cl) {
    # Match GMM rows to Contact Matrix rows
    rows <- which(gmm_data$State == cl) 
    colMeans(contact_matrix[rows, , drop=FALSE], na.rm = TRUE)
})
names(contact_freq_list) <- cluster_ids

# Reshape for Heatmap
freq_df <- as.data.frame(do.call(cbind, contact_freq_list))
freq_df$Contact <- colnames(contact_matrix)
freq_df$HelixRes  <- as.numeric(sub("H(\\d+)_P\\d+", "\\1", freq_df$Contact))
freq_df$PocketRes <- as.numeric(sub("H\\d+_P(\\d+)", "\\1", freq_df$Contact))
melted <- melt(freq_df, id.vars=c("Contact", "HelixRes", "PocketRes"), variable.name="Cluster", value.name="Frequency")

# Plot Heatmap
p_heatmap <- ggplot(melted, aes(x=PocketRes, y=HelixRes, fill=Frequency)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(option = "magma", limits=c(0,1), name="Contact\nFrequency") +
    # REMOVED +20 OFFSET
    facet_wrap(~Cluster) +
    labs(title=paste("Contact Heatmap by State"), x="Pocket Residue", y="Helix Residue") +
    theme_publication() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

ggsave("figure_heatmap_GMM.png", p_heatmap, width=12, height=8, dpi=300)

cat("Script 02 Complete. All contact analysis and plots saved.\n")
