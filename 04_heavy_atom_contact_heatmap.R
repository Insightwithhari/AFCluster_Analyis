############################################################
# Script 04: State-resolved heavy-atom contact heatmap
############################################################

# --- Load libraries ---
library(bio3d)
library(ggplot2)
library(reshape2)
library(pbapply)
library(viridis)

# --- Inputs (MUST match other scripts) ---
pdb_file  <- "Reference.pdb"
traj_file <- "Trajectory.dcd"
gmm_file  <- "global_PC1_GMM_3_states.csv"

helix_resnos  <- 190:205
pocket_resnos <- 21:109
contact_cutoff <- 5  # Å

# --- Load data ---
pdb  <- read.pdb(pdb_file)
traj <- read.dcd(traj_file)
gmm  <- read.csv(gmm_file)
gmm$State <- factor(gmm$State)

frames_to_use <- gmm$Frame
states <- gmm$State

cat("Using", length(frames_to_use), "frames from GMM state definition\n")

# --- Precompute heavy-atom indices ---
helix_xyz <- lapply(helix_resnos, function(r)
  atom.select(pdb, resno = r, string = "noh")$xyz
)

pocket_xyz <- lapply(pocket_resnos, function(r)
  atom.select(pdb, resno = r, string = "noh")$xyz
)

# --- Heavy-atom contact calculation for one frame ---
calc_contacts <- function(xyz) {
  out <- matrix(FALSE,
                nrow = length(helix_xyz),
                ncol = length(pocket_xyz))
  
  for (i in seq_along(helix_xyz)) {
    H <- matrix(xyz[helix_xyz[[i]]], ncol = 3, byrow = TRUE)
    for (j in seq_along(pocket_xyz)) {
      P <- matrix(xyz[pocket_xyz[[j]]], ncol = 3, byrow = TRUE)
      out[i, j] <- min(bio3d::dist.xyz(H, P), na.rm = TRUE) <= contact_cutoff
    }
  }
  out
}

# --- Compute contacts for all clustered frames ---
cat("Computing heavy-atom contacts...\n")

contact_array <- pbsapply(
  frames_to_use,
  function(f) calc_contacts(traj[f, ]),
  simplify = "array"
)

# --- State-resolved contact frequencies ---
state_levels <- levels(states)

contact_freq <- lapply(state_levels, function(st) {
  idx <- which(states == st)
  apply(contact_array[, , idx, drop = FALSE], c(1, 2), mean)
})

names(contact_freq) <- state_levels

# --- Prepare dataframe for plotting ---
df_list <- lapply(names(contact_freq), function(st) {
  mat <- contact_freq[[st]]
  df <- melt(mat)
  colnames(df) <- c("HelixIndex", "PocketIndex", "Frequency")
  df$HelixRes  <- helix_resnos[df$HelixIndex]
  df$PocketRes <- pocket_resnos[df$PocketIndex]
  df$State <- st
  df
})

df_all <- do.call(rbind, df_list)

# --- Plot heatmap ---
p_heatmap <- ggplot(df_all,
                    aes(x = PocketRes,
                        y = HelixRes,
                        fill = Frequency)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1),
    name = "Contact\nFrequency"
  ) +
  facet_wrap(~State) +
  labs(
    title = paste(
      "State-Resolved Heavy-Atom Contact Frequencies",
      "(cutoff =", contact_cutoff, "Å)"
    ),
    x = "Pocket Residue",
    y = "Helix Residue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.background = element_rect(fill = "grey85", color = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = "Fig_heavy_atom_contact_heatmap.png",
  plot     = p_heatmap,
  width    = 12,
  height   = 8,
  dpi      = 300
)

cat("Heatmap generated: Fig_heavy_atom_contact_heatmap.png\n")
