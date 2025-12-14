############################################################
# Script 03: Publication-quality figure generation
# (Density, PCA scatter, boxplot, correlation, violin,
#  gate residue frequency, contact fraction)
############################################################

# --- Load libraries ---
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)

# --- Publication theme ---
theme_publication <- function(base_size = 14) {
    theme_bw(base_size = base_size) +
        theme(
            plot.title   = element_text(face = "bold", hjust = 0.5),
            axis.title   = element_text(face = "bold"),
            axis.text    = element_text(colour = "black"),
            panel.grid   = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA)
        )
}

# ---------------------------------------------------------
# Load input data
# ---------------------------------------------------------
gmm <- read.csv("global_PC1_GMM_3_states.csv")
contacts <- read.csv("per_residue_contact_fractions.csv")

gmm$State <- factor(gmm$State)

cluster_colors <- brewer.pal(nlevels(gmm$State), "Dark2")

# ---------------------------------------------------------
# Figure 1: PC1 density plot
# ---------------------------------------------------------
p_density <- ggplot(gmm, aes(x = PC1, fill = State)) +
    geom_density(alpha = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "State Distribution along PC1",
         x = "PC1", y = "Probability Density") +
    theme_publication()

ggsave(
    filename = "Fig_PC1_density.png",
    plot     = p_density,
    width    = 7,
    height   = 5,
    dpi      = 600
)

# ---------------------------------------------------------
# Figure 2: PCA scatter plot
# ---------------------------------------------------------
p_pca <- ggplot(gmm, aes(x = PC1, y = PC2, color = State)) +
    geom_point(alpha = 0.6, size = 1.5) +
    stat_ellipse(type = "norm", level = 0.80, linewidth = 0.8) +
    scale_color_manual(values = cluster_colors) +
    labs(title = "Conformational Landscape",
         x = "PC1", y = "PC2") +
    theme_publication()

ggsave(
    filename = "Fig_PCA_scatter.png",
    plot     = p_pca,
    width    = 7,
    height   = 5,
    dpi      = 600
)

# ---------------------------------------------------------
# Figure 3: Distance by state (boxplot)
# ---------------------------------------------------------
p_box <- ggplot(gmm, aes(x = State, y = Distance, fill = State)) +
    geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "Helix–Pocket Distance by State",
         x = "Conformational State",
         y = "Average Heavy-Atom Distance (Å)") +
    theme_publication() +
    theme(legend.position = "none")

ggsave(
    filename = "Fig_distance_box.png",
    plot     = p_box,
    width    = 7,
    height   = 5,
    dpi      = 600
)

# ---------------------------------------------------------
# Figure 4: PC1 vs distance correlation
# ---------------------------------------------------------
p_corr <- ggplot(gmm, aes(x = PC1, y = Distance)) +
    geom_point(alpha = 0.4, size = 1.8, color = "red") +
    geom_smooth(method = "lm", se = FALSE,
                linetype = "dashed", color = "black") +
    labs(title = "PC1 vs. Helix–Pocket Distance",
         x = "PC1",
         y = "Average Heavy-Atom Distance (Å)") +
    theme_publication()

ggsave(
    filename = "Fig_PC1_vs_distance.png",
    plot     = p_corr,
    width    = 7,
    height   = 5,
    dpi      = 600
)

# ---------------------------------------------------------
# Figure 5: Per-residue contact fraction
# ---------------------------------------------------------
contacts_long <- pivot_longer(
    contacts,
    cols = -Residue,
    names_to = "State",
    values_to = "Fraction"
)

p_contact <- ggplot(contacts_long,
                    aes(x = factor(Residue),
                        y = Fraction,
                        fill = State)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_manual(values = cluster_colors) +
    labs(title = "Per-Residue Contact Fraction by State",
         x = "Helix Residue",
         y = "Contact Fraction") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
    filename = "Fig_contact_fraction.png",
    plot     = p_contact,
    width    = 8,
    height   = 5,
    dpi      = 600
)

cat("All figures generated successfully.\n")
