############################################################
# Script 02: State-resolved contact and distance analysis
############################################################

library(bio3d)
library(dplyr)
library(tidyr)

# --- Inputs ---
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

# --- Per-residue minimum heavy-atom distance ---
helix_resids <- sort(unique(pdb$atom$resno[pdb$atom$resno %in% helix_resnos]))
pocket_inds  <- atom.select(pdb, resno = pocket_resnos, string = "noh")$xyz

perres_min <- matrix(
  NA_real_,
  nrow = length(helix_resids),
  ncol = nrow(traj),
  dimnames = list(helix_resids, NULL)
)

for (r in seq_along(helix_resids)) {
  h_inds <- atom.select(pdb, resno = helix_resids[r], string = "noh")$xyz
  perres_min[r, ] <- sapply(seq_len(nrow(traj)), function(i) {
    min(
      bio3d::dist.xyz(
        matrix(traj[i, h_inds], ncol = 3, byrow = TRUE),
        matrix(traj[i, pocket_inds], ncol = 3, byrow = TRUE)
      ),
      na.rm = TRUE
    )
  })
}

# --- State-resolved contact fractions ---
state_frames <- split(gmm$Frame, gmm$State)

contact_df <- lapply(helix_resids, function(r) {
  v <- perres_min[as.character(r), ]
  frac <- sapply(state_frames, function(fr)
    mean(v[fr] <= contact_cutoff, na.rm = TRUE)
  )
  c(Residue = r, frac)
})

contact_df <- as.data.frame(do.call(rbind, contact_df))
write.csv(contact_df, "per_residue_contact_fractions.csv", row.names = FALSE)
