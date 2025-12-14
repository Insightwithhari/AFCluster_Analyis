############################################################
# Script 01: Conformational state definition using PCA + GMM
############################################################

# --- Load libraries ---
library(bio3d)
library(mclust)
library(dplyr)

# --- Input files ---
topology_pdb  <- "Reference.pdb"
trajectory_dcd <- "Trajectory.dcd"

# --- Structural regions (defined from structural analysis) ---
helix_resnos  <- 190:205
pocket_resnos <- 21:109

# --- Load structure and trajectory ---
pdb  <- read.pdb(topology_pdb)
traj <- read.dcd(trajectory_dcd)

# --- Average heavy-atom helix–pocket distance ---
helix_sel  <- atom.select(pdb, resno = helix_resnos, string = "noh")
pocket_sel <- atom.select(pdb, resno = pocket_resnos, string = "noh")

avg_heavy_dist <- function(i) {
  xyz <- traj[i, ]
  H <- matrix(xyz[helix_sel$xyz],  ncol = 3, byrow = TRUE)
  P <- matrix(xyz[pocket_sel$xyz], ncol = 3, byrow = TRUE)
  mean(bio3d::dist.xyz(H, P), na.rm = TRUE)
}

dist_all <- sapply(seq_len(nrow(traj)), avg_heavy_dist)
ok <- is.finite(dist_all)

# --- Global PCA on C-alpha atoms ---
ca_sel <- atom.select(pdb, elety = "CA")
pc <- pca.xyz(traj[, ca_sel$xyz], rm.gaps = TRUE)

# --- Orient PC1 for physical interpretability ---
rho <- suppressWarnings(cor(pc$z[ok, 1], dist_all[ok], method = "spearman"))
if (!is.na(rho) && rho < 0) {
  pc$z <- -pc$z
}

# --- Gaussian Mixture Modeling on PC1 ---
pc1 <- pc$z[ok, 1]
gmm <- Mclust(pc1)

cat("Optimal number of states (BIC):", gmm$G, "\n")

states <- factor(paste0("Cluster ", gmm$classification))

state_table <- data.frame(
  Frame    = which(ok),
  PC1      = pc1,
  PC2      = pc$z[ok, 2],
  Distance = dist_all[ok],
  State    = states
)

write.csv(
  state_table,
  paste0("global_PC1_GMM_", gmm$G, "_states.csv"),
  row.names = FALSE
)

# --- Identify centroid frames ---
centroids <- state_table %>%
  group_by(State) %>%
  mutate(
    c1 = mean(PC1),
    c2 = mean(PC2),
    d  = sqrt((PC1 - c1)^2 + (PC2 - c2)^2)
  ) %>%
  slice_min(d, n = 1) %>%
  ungroup()

# --- Write centroid structures ---
for (i in seq_len(nrow(centroids))) {
  write.pdb(
    pdb = pdb,
    xyz = traj[centroids$Frame[i], ],
    file = paste0(gsub(" ", "_", centroids$State[i]), "_centroid.pdb")
  )
}

# --- Save session information ---
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
