############################################################
# Script 02: Violin plot of min CA–CA distance
############################################################

library(bio3d)
library(ggplot2)

topology <- read.pdb("Refrence.pdb")
traj_xyz <- read.dcd("Trajectory.dcd")

helix_resnos <- 190:205
pocket_resnos <- 21:109

ca_helix  <- atom.select(topology, elety="CA", resno=helix_resnos)
ca_pocket <- atom.select(topology, elety="CA", resno=pocket_resnos)

min_ca <- sapply(seq_len(nrow(traj_xyz)), function(i){
  H <- matrix(traj_xyz[i, ca_helix$xyz], ncol=3, byrow=TRUE)
  P <- matrix(traj_xyz[i, ca_pocket$xyz], ncol=3, byrow=TRUE)
  min(bio3d::dist.xyz(H,P), na.rm=TRUE)
})

df <- data.frame(Distance=min_ca)

p <- ggplot(df, aes(x="", y=Distance)) +
  geom_violin(fill="skyblue") +
  geom_jitter(width=0.1, alpha=0.3)

ggsave("violin_min_CA_distance.png", p, width=7, height=5, dpi=600)
