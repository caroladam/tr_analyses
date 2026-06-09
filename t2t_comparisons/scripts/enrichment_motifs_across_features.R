########################################
### Plot motif sizes for all species ###
########################################

library(ggplot2)
library(ggdark)
library(dplyr)
library(viridis)

color_mapping <- c("bonobo" = "#ff8c69ff", "gorilla" = "#165016ff", "homo"= "#5c5c8aff", "pabelli" = "#08519c", "ppyg" = "#c83737ff", "pantro" = "#8a865dff", "symsyn" = "black")

ggplot() +
geom_density(data=df, aes(x=V1, y=after_stat(count), colour=V5), alpha=0.5, position="identity", show.legend=F) +
scale_x_continuous(name="Motif length") +
scale_y_continuous(trans="log1p", name= "Density (log)", breaks = c(0, 100, 1000, 10000, 100000, 1000000, 3000000)) +
#scale_y_continuous(name= "Density") +
scale_color_manual(values=color_mapping) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=20),
axis.title.x=element_text(size= 20),
axis.title.y=element_text(size = 20),
axis.text.x=element_text(colour="black", size = 20),
axis.text.y=element_text(colour="black", size = 20))

#or as points
# Aggregate the data: count of V4 per V8
df_summary <- df %>%
  group_by(V8, V4) %>%
  summarise(count = n(), .groups = "drop")

# Plot with lines and points
ggplot(df_summary, aes(x = V4, y = count, colour = V8, group = V8)) +
  geom_line(linewidth = 2, show.legend=F) +
  scale_x_continuous(name = "Motif length") +
  scale_y_continuous(trans = "log1p", name = "Count (log)") +
scale_color_manual(values=color_mapping) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 21),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20))

###############################################################
########## Plot motif size proportion per species##############
###############################################################

#Bin motif lengths
df <- df %>%
  mutate(motif_class = case_when(
    V1 %in% 1:6 ~ as.character(V1),
    V1 >= 7 & V1 <= 20 ~ "7-20",
    V1 > 20 ~ ">20"))

#Count per species and motif class
df_summary <- df %>%
  group_by(V4, motif_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(V4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# order motif classes
df_summary$motif_class <- factor(df_summary$motif_class,
                                 levels = rev(c("1", "2", "3", "4", "5", "6", "7-20", ">20")))

df_summary$V4 <- factor(df_summary$V4, 
                         levels = c("homo", "pantro", "bonobo", "gorilla", 
                                    "pabelii", "ppyg", "symsyn"))

# barplot
ggplot(df_summary, aes(x = V4, y = prop, fill = motif_class)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), name = "Proportion of TRs") +
  scale_fill_viridis_d(direction = 1) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=21),
axis.title.x=element_text(size= 20),
axis.title.y=element_text(size = 20),
axis.text.x=element_text(colour="black", size = 20),
axis.text.y=element_text(colour="black", size = 20))

#### Test if motif classes are enriched or depleted for each species
df_summary <- df %>% group_by(V4) %>% count(motifs)

chisq_enrichment <- function(data) {
  
  results <- data %>%
    group_by(V4) %>%
    summarise(
      test = list(chisq.test(n)),
      counts = list(n),
      classes = list(motif_class),
      .groups = "drop") %>%
    rowwise() %>%
    mutate(
      observed = list(test$observed),
      expected = list(test$expected),
      residuals = list(test$residuals)) %>%
    select(V4, classes, observed, expected, residuals)
  
  # Expand into long format for motif classes
  per_class <- results %>%
    unnest(c(classes, observed, expected, residuals)) %>%
    mutate(
      per_class_p = 2 * (1 - pnorm(abs(residuals))),
      per_class_fdr = p.adjust(per_class_p, method = "fdr"),
      fold_enrichment = observed / expected,
      status = case_when(
        residuals > 2 & per_class_fdr < 0.05 ~ "enriched",
        residuals < -2 & per_class_fdr < 0.05 ~ "depleted",
        TRUE ~ "ns")) 
  
  return(per_class)}

# df_summary should have columns: V4 (species), motif_class, n (counts)
results <- chisq_enrichment(df_summary)

######## barplot of shared TRs ########

df<-trs_count_per_species

color_mapping <- c(
  "bonobo" = "#ff8c69ff",
  "gorilla" = "#165016ff",
  "homo_total" = "#5c5c8aff",
  "pabelii" = "#08519c",
  "ppyg" = "#c83737ff",
  "pantro" = "#8a865dff",
  "symsyn" = "#c49a00",
  "macaca" = "black")

ggplot(df, aes(x = V2, y = V1, fill = V1)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_x_continuous(name = "TR count") +
  scale_y_discrete(name = "Species") +
  scale_fill_manual(values = color_mapping) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(size = 21),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20))


