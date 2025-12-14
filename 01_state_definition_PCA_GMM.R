############################################################
# Script 01: Reference PCA–GMM pipeline (DO NOT MODIFY)
############################################################

library(bio3d)
library(ggplot2)
library(mclust)
library(RColorBrewer)
library(dplyr)

theme_publication <- function(base_size=14) {
  theme_bw(base_size=base_size) +
    theme(
      plot.title = element_text(face="bold", hjust=0.5),
      axis.title = element_text(face="bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour="black", fill=NA)
    )
}

topology_pdb <- "Refrence.pdb"
dcd_file     <- "Trajectory.dcd"
helix_resnos <- 190:205
pocket_resnos <- 21:109

topology <- read.pdb(topology_pdb)
traj_xyz <- read.dcd(dcd_file)

# --- Average heavy-atom distance ---
heavy_helix  <- atom.select(topology, string="noh", resno=helix_resnos)
heavy_pocket <- atom.select(topology, string="noh", resno=pocket_resnos)

avg_heavy_distance <- function(i) {
  xyz <- traj_xyz[i,]
  H <- matrix(xyz[heavy_helix$xyz], ncol=3, byrow=TRUE)
  P <- matrix(xyz[heavy_pocket$xyz], ncol=3, byrow=TRUE)
  mean(bio3d::dist.xyz(H,P), na.rm=TRUE)
}

d_avg_heavy_all <- sapply(seq_len(nrow(traj_xyz)), avg_heavy_distance)
ok <- is.finite(d_avg_heavy_all)
d_avg_heavy <- d_avg_heavy_all[ok]

# --- PCA ---
sel_ca_global <- atom.select(topology, elety="CA")
pc <- pca.xyz(traj_xyz[, sel_ca_global$xyz], rm.gaps=TRUE)

rho <- suppressWarnings(
  cor.test(pc$z[ok,1], d_avg_heavy, method="spearman")
)$estimate

if (!is.na(rho) && rho < 0) pc$z <- -pc$z

final_corr_test <- suppressWarnings(
  cor.test(pc$z[ok,1], d_avg_heavy, method="spearman")
)

# --- GMM ---
pc1_filtered <- pc$z[ok,1]
gm <- Mclust(pc1_filtered)
num_states <- gm$G

gmm_results <- data.frame(
  Frame = which(ok),
  PC1 = pc1_filtered,
  PC2 = pc$z[ok,2],
  Distance = d_avg_heavy,
  State = factor(paste0("Cluster ", gm$classification))
)

write.csv(
  gmm_results,
  paste0("global_PC1_GMM_", num_states, "_states.csv"),
  row.names=FALSE
)

# --- Centroids ---
centroids <- gmm_results %>%
  group_by(State) %>%
  mutate(
    c1 = mean(PC1),
    c2 = mean(PC2),
    d = sqrt((PC1-c1)^2 + (PC2-c2)^2)
  ) %>%
  slice_min(d, n=1) %>%
  ungroup()

for (i in 1:nrow(centroids)) {
  write.pdb(
    topology,
    file=paste0(gsub(" ","_",centroids$State[i]),"_centroid.pdb"),
    xyz=traj_xyz[centroids$Frame[i],]
  )
}
