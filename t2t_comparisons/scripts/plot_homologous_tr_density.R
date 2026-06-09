required_packages <- c("ggplot2", "dplyr", "tidyr")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(dplyr)
library(tidyr)

# get human density
human <- read.table("human_density_lookup.txt", header=F, col.names = c("humanID", "human_density"))

# get NHP density with human window ID
nhp <- read.table("lifted_vs_chimp_density.bed",header=F, col.names = c(
    "nhp_chr",
    "nhp_start",
    "nhp_end",
    "humanID",
    "count_TRs",
    "covered_bp",
    "lifted_length",
    "nhp_density"))

# join both using the human window ID
df <- human %>%
  inner_join(nhp %>% 
    select(humanID, nhp_density, lifted_length), by = "humanID")

# normalize so densities sum up to 1
df <- df %>%
  mutate(
    human_norm = human_density / sum(human_density),
    nhp_norm = nhp_density / sum(nhp_density))

r2 <- df %>%
  summarize(
    r2 = cor(human_norm, nhp_norm, use = "complete.obs")^2) %>% pull(r2)

maxval <- max(
  c(df$human_norm, df$nhp_norm),
  na.rm = TRUE)

# Heatmap comparing normalized densities across species
p <- ggplot(df, aes(x=human_norm, y=nhp_norm)) +
  geom_hex(bins = 30) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "red") +
  scale_fill_viridis_c(option = "F", direction=-1, name = "Window count") +
  annotate("text", x = max(df$human_norm) * 0.75, y = max(df$nhp_norm) * 0.95, label = paste0("R² = ", round(r2, 3))) +
  scale_x_continuous(limits = c(0, maxval), name = "Human normalized TR density") +
  scale_y_continuous(limits = c(0, maxval), name = "Chimpanzee normalized TR density") +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"),
    legend.position = c(0.9, 0.2),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.key.size = unit(1, "cm"))

ggsave("human_chimp_tr_density.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)
