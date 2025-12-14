############################################################
# Script 04: Heavy-atom contact heatmap
############################################################

library(bio3d)
library(reshape2)
library(ggplot2)
library(pbapply)
library(viridis)

gmm <- read.csv("global_PC1_GMM_3_states.csv")
topology <- read.pdb("Refrence.pdb")
traj_xyz <- read.dcd("Trajectory.dcd")

helix_resnos <- 190:205
pocket_resnos <- 21:109
cutoff <- 5

H_xyz <- lapply(helix_resnos, function(r)
  atom.select(topology, resno=r, string="noh")$xyz)

P_xyz <- lapply(pocket_resnos, function(r)
  atom.select(topology, resno=r, string="noh")$xyz)

contact <- pbsapply(gmm$Frame, function(i){
  sapply(seq_along(H_xyz), function(h){
    sapply(seq_along(P_xyz), function(p){
      min(bio3d::dist.xyz(
        matrix(traj_xyz[i,H_xyz[[h]]], ncol=3, byrow=TRUE),
        matrix(traj_xyz[i,P_xyz[[p]]], ncol=3, byrow=TRUE)
      )) <= cutoff
    })
  })
})

freq <- apply(contact, c(1,2), mean)
df <- melt(freq)

p <- ggplot(df, aes(Var2, Var1, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c()

ggsave("heavy_atom_contact_heatmap.png", p, width=10, height=8, dpi=300)
