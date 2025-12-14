############################################################
# Script 03: Per-residue contact fraction
############################################################

library(bio3d)
library(tidyr)
library(ggplot2)

gmm <- read.csv("global_PC1_GMM_3_states.csv")
topology <- read.pdb("Refrence.pdb")
traj_xyz <- read.dcd("Trajectory.dcd")

helix_resnos <- 190:205
pocket_resnos <- 21:109
cutoff <- 5

state_indices <- split(gmm$Frame, gmm$State)

helix_resids <- helix_resnos
pocket_inds <- atom.select(topology, resno=pocket_resnos, string="noh")$xyz

perres_min <- sapply(helix_resids, function(r){
  h_inds <- atom.select(topology, resno=r, string="noh")$xyz
  sapply(seq_len(nrow(traj_xyz)), function(i){
    H <- matrix(traj_xyz[i,h_inds], ncol=3, byrow=TRUE)
    P <- matrix(traj_xyz[i,pocket_inds], ncol=3, byrow=TRUE)
    min(bio3d::dist.xyz(H,P))
  })
})

contact_frac <- lapply(state_indices, function(fr){
  colMeans(perres_min[fr,] <= cutoff)
})

df <- as.data.frame(contact_frac)
df$Residue <- helix_resids

df_long <- pivot_longer(df, -Residue, names_to="State", values_to="Fraction")

p <- ggplot(df_long, aes(factor(Residue), Fraction, fill=State)) +
  geom_col(position="dodge")

ggsave("per_residue_contact_fraction.png", p, width=8, height=5, dpi=600)
