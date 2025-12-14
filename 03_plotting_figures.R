############################################################
# Script 03: Figure generation
############################################################

library(ggplot2)
library(RColorBrewer)

theme_pub <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

data <- read.csv("global_PC1_GMM_3_states.csv")

# Density plot
p1 <- ggplot(data, aes(PC1, fill = State)) +
  geom_density(alpha = 0.7) +
  scale_fill_brewer(palette = "Dark2") +
  theme_pub()

ggsave("Figure_PC1_density.png", p1, width = 7, height = 5, dpi = 600)
