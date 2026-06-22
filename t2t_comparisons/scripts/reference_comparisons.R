required_packages <- c("ggplot2", "dplyr", "tidyr", "broom")

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
library(broom)

# concatenated catalog with all spp-specific TRs
df <- read.table("all_spp_catalog.bed", header=F, col.names = c("motif_len", "tr_len", "motif_seq", "spp"))

color_mapping <- c(
  "bonobo" = "#ff8c69ff",
  "gorilla" = "#165016ff",
  "homo" = "#5c5c8aff",
  "pabelii" = "#08519c",
  "ppyg" = "#c83737ff",
  "pantro" = "#8a865dff",
  "symsyn" = "#c49a00",
  "macaca" = "black")

ggplot() +
geom_density(data=df%>%filter(motif_len<=400), aes(x=motif_len, y=after_stat(count), colour=spp), alpha=0.5, position="identity", show.legend=F) +
scale_x_continuous(name="Motif length (bp)") +
scale_y_continuous(trans="log1p", name= "Density", breaks = c(0, 100, 1000, 10000, 100000, 1000000, 3000000)) +
scale_color_manual(values=color_mapping) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

##Plot motif size proportion per species

# Bin motif lengths
df <- df %>%
  mutate(motif_class = case_when(
    motif_len %in% 1:6 ~ as.character(motif_len),
    motif_len >= 7 & motif_len <= 20 ~ "7-20",
    motif_len > 20 ~ ">20"))

# Count per species and motif class
df_summary <- df %>%
  group_by(spp, motif_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(spp) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

df_summary$motif_class <- factor(df_summary$motif_class, levels = rev(c("1", "2", "3", "4", "5", "6", "7-20", ">20")))
df_summary$spp <- factor(df_summary$spp, levels = c("homo", "pantro", "bonobo", "gorilla", "pabelii", "ppyg", "symsyn", "macaca"))

ggplot(df_summary, aes(x = spp, y = prop, fill = motif_class)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), name = "Proportion of TRs") +
  scale_fill_viridis_d(direction = 1) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

# homologous TR catalogs - count of TRs shared between human and each NHP
count<-read.table("homologous_tr_count_spp.txt", header=F, col.names = c("spp", "homologous_count"))

ggplot(count, aes(x = homologous_count, y = spp, fill = spp)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_x_continuous(name = "TR count") +
  scale_y_discrete(name = "Species") +
  scale_fill_manual(values = color_mapping) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

# Comparison of human and NHP TR length - considering only TRs shared between human and all NHP species

df<-read.table("homologous_tr_catalog_all_spp.bed", header=F,
  col.names=c("chr_homo", "start_homo", "end_homo", "motif_len_homo", "tr_len_homo", "copy_n_homo", "motif_seq_homo", "tr_seq_homo", "spp_homo",
              "chr_nhp", "start_nhp", "end_nhp", "motif_len_nhp", "tr_len_nhp", "copy_n_nhp", "motif_seq_nhp", "tr_seq_nhp", "spp_nhp", "TRID", "genomic_feature"))

color_mapping <- c(
  "cds" = "darkgoldenrod1",
  "promoter" = "dodgerblue",
  "5utr" = "firebrick",
  "3utr" = "coral1",
  "intronic"= "mediumseagreen",
  "intergenic" = "darkviolet")

df <- df %>% arrange(factor(genomic_feature, levels = c("intergenic", "intronic", "promoter", "cds", "3utr", "5utr")))

ggplot(df, aes(x = tr_len_homo, y = tr_len_nhp, colour = genomic_feature)) + 
  geom_abline(linewidth=1, slope = 1, intercept = 0, color = "black") + 
  facet_wrap(~ spp_nhp, scales="free") + 
  geom_smooth(method="lm", linewidth=1, alpha = 1, show.legend=F) +
  ggrastr::geom_point_rast(size = 2, alpha = 0.8, shape=16, show.legend=F) + 
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(name = "Human TR length") +
  scale_y_continuous(name = "NHP TR length") +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

# Calculate r2 per feature - also for filtered top most divergent TRs 
df <- df %>%
  mutate(relative_diff = abs(tr_len_nhp - tr_len_homo) / ((tr_len_nhp + tr_len_homo) / 2))

# Filter more divergent TRs based on percentage per genomic feature
dff<- df %>%                                                                          
  group_by(spp_nhp, genomic_feature) %>%                                               
  filter(relative_diff <= quantile(relative_diff, probs = 0.999, na.rm = TRUE)) %>%
  ungroup()

r_squared_list <- df %>%
  group_by(spp_nhp, genomic_feature) %>%
  summarise(
    r_squared = round(
      summary(lm(tr_len_homo ~ tr_len_nhp, data = pick(tr_len_homo, tr_len_nhp)))$r.squared,3)) %>% arrange(genomic_feature)

r_squared_list_f <- dff %>%
  group_by(spp_nhp, genomic_feature) %>%
  summarise(
    r_squared = round(
      summary(lm(tr_len_homo ~ tr_len_nhp, data = pick(tr_len_homo, tr_len_nhp)))$r.squared,3)) %>% arrange(genomic_feature)

g<-cbind(r_squared_list, r_squared_list_f$r_squared)

write.table(g,"r_squared_per_feature_01.txt", col.names=T, row.names=T, sep="\t", quote=F) # added divergence time

# Plot correlation and divergence time after concatenating across all species
ggplot(g, aes(x = divergence_time, y = r2_01, colour = spp, shape = feature)) +
  geom_point(size = 2, alpha = 0.7, show.legend=F) +
  scale_colour_manual(values=color_mapping) +
  scale_x_continuous(name = "Divergence time (Mya)") +
  scale_y_continuous(name = "R2 (TR length correlation with human)") +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=7, family="Helvetica"),
axis.title.x=element_text(size=7, family="Helvetica"),
axis.title.y=element_text(size=7, family="Helvetica"),
axis.text.x=element_text(colour="black", size=7, family="Helvetica"),
axis.text.y=element_text(colour="black", size=7, family="Helvetica"))

# calculate slopes per genomic feature

feature_slopes <- df %>%
  group_by(feature) %>%
  do(
    tidy(lm(r2 ~ divergence_time, data = .))) %>%
  filter(term == "divergence_time") %>%
  select(feature, slope = estimate, std.error, statistic, p.value)

# Compare TR length per motif class between human and each NHP species

df <- df %>%
  mutate(
    tr_len_homo = as.numeric(tr_len_homo),
    tr_len_nhp = as.numeric(tr_len_nhp))

df <- df %>%
  mutate(tr_len_nhp = ifelse(is.na(tr_len_nhp), motif_len_nhp * copy_n_nhp, tr_len_nhp))

results <- df %>%
  group_by(spp_nhp, motif_len_homo) %>%
  summarize(
    n = sum(!is.na(tr_len_homo) & !is.na(tr_len_nhp)),
    mean_human = sprintf("%.2f", mean(tr_len_homo, na.rm = TRUE)),
    mean_nhp   = sprintf("%.2f", mean(tr_len_nhp, na.rm = TRUE)),
    p_value = if (n >= 2) {
      wilcox.test(tr_len_homo, tr_len_nhp, paired = TRUE)$p.value
    } else {
      NA_real_
    }, .groups = "drop")

write.table(results,"tr_length_per_motif_all_spp.txt", col.names=T, row.names=T, sep="\t", quote=F)

color_mapping <- c(
  "bonobo" = "#ff8c69ff",
  "gorilla" = "#165016ff",
  "pabelii" = "#08519c",
  "ppyg" = "#c83737ff",
  "pantro" = "#8a865dff",
  "symsyn" = "#c49a00",
  "macaca" = "black")

# for plotting, use only STR lengths
df<-df%>%filter(motif_len_homo<=6)

results <- results %>%
  mutate(
    mean_human = as.numeric(mean_human),
    mean_nhp = as.numeric(mean_nhp))

plot_df <- results %>%
  mutate(
    diff = mean_human-mean_nhp,
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""))

ggplot(plot_df, aes(x = motif_len_homo, y = diff, color = spp_nhp)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  geom_point(size = 2, show.legend=F) +
  geom_line(show.legend=F, linewidth=0.6) +
  geom_text(aes(label = sig), vjust = -0.7, size = 3, show.legend=F) +
  scale_x_continuous(name = "Motif length", breaks = 1:6, labels = c("Mono", "Di", "Tri", "Tetra", "Penta", "Hexa")) +
  scale_y_continuous(name = "Human - NHP TR length") +
  scale_color_manual(values = color_mapping) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=10, family="Helvetica"),
axis.title.x=element_text(size=10, family="Helvetica"),
axis.title.y=element_text(size=10, family="Helvetica"),
axis.text.x=element_text(colour="black", size=10, family="Helvetica"),
axis.text.y=element_text(colour="black", size=10, family="Helvetica"))
