############################################################
# Script 03: Publication-quality figure generation
############################################################

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)

# --- Theme ---
theme_publication <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black")
    )
}

# --- Load data ---
gmm <- read.csv("global_PC1_GMM_3_states.csv")
contacts <- read.csv("per_residue_contact_fractions.csv")
gmm$State <- factor(gmm$State)

cols <- brewer.pal(nlevels(gmm$State), "Dark2")

# --- PC1 Density ---
p1 <- ggplot(gmm, aes(PC1, fill = State)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = cols) +
  theme_publication()

ggsave("Fig_PC1_density.png", p1, 7, 5, dpi = 600)

# --- PCA Scatter ---
p2 <- ggplot(gmm, aes(PC1, PC2, color = State)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = cols) +
  theme_publication()

ggsave("Fig_PCA_scatter.png", p2, 7, 5, dpi = 600)

# --- Distance by State ---
p3 <- ggplot(gmm, aes(State, Distance, fill = State)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = cols) +
  theme_publication() +
  theme(legend.position = "none")

ggsave("Fig_distance_box.png", p3, 7, 5, dpi = 600)

# --- PC1 vs Distance ---
p4 <- ggplot(gmm, aes(PC1, Distance)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
  theme_publication()

ggsave("Fig_PC1_vs_distance.png", p4, 7, 5, dpi = 600)

# --- Per-residue contact fraction ---
contacts_long <- pivot_longer(
  contacts, cols = -Residue,
  names_to = "State", values_to = "Fraction"
)

p5 <- ggplot(contacts_long, aes(factor(Residue), Fraction, fill = State)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig_contact_fraction.png", p5, 8, 5, dpi = 600)
