required_packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "stringr", "ggrastr")

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
library(viridis)
library(ggrepel)
library(tidyr)
library(stringr)
library(ggrastr)

df <- read.table("homo_catalog.no_overlaps_censat_annotated.bed", header=F, col.names = c(
    "chr",
    "start",
    "end",
    "motif_len",
    "tr_len",
    "copy_n",
    "constancy",
    "motif_seq",
    "spp",
    "feature"))

# Get proportion of chromosomes makeup by TRs
kar<-read.table("karyotipe_chm13.txt", header=F, col.names = c(
    "chr",
    "start",
    "end",
    "centromere_start",
    "centromere_end"))

# Get total length of TR per chromosome
length_per_chr<-aggregate(tr_len~chr, data=df, sum)

# merge sum data with the chromosome sizes
merged_df<-merge(length_per_chr, kar[,c("chr", "end")], by="chr")

# Calculate the proportion of the chromosome that is make up by TRs
merged_df$proportion<-(merged_df$tr_len/merged_df$end)*100

merged_df$chr <- factor(merged_df$chr, levels = arrange(merged_df$chr))

ggplot(merged_df, aes(y=chr, x=proportion)) +
  geom_bar(stat="identity", fill="lightgray") +
scale_y_discrete(name = "Chromosome") +
scale_x_continuous(name = "Proportion (%)") +
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

color_mapping <- c("cds" = "darkgoldenrod1", "promoter" = "dodgerblue", "5utr" = "firebrick", "3utr" = "coral1", "intron"= "mediumseagreen", "intergenic" = "darkviolet")

df <- df %>% arrange(factor(feature, levels = c("intergenic","intron","promoter", "cds", "3utr", "5utr")))

df <- df %>%
  mutate(motif_class = case_when(
    motif_len %in% 1:12 ~ as.character(motif_len),
    motif_len >= 13 & motif_len <= 50  ~ "13-50", 
    motif_len > 50 ~ ">50")) %>%
  mutate(motif_class = factor(motif_class, levels = c(
    as.character(1:12), "13-50", ">50")))

#Count per feature per motif class
mt_c <- df %>% count(feature, motif_class) %>%
  group_by(feature) %>%
  mutate(prop = n / sum(n))

mt_c <- mt_c %>% mutate(feature = factor(feature, levels = c("intergenic","intronic","promoter", "cds", "5utr", "3utr")))
mt_c <- mt_c %>% mutate(motif_class = factor(motif_class, levels = c(">50","13-50","12","11","10", "9", "8", "7", "6", "5", "4", "3", "2", "1")))

# stacked barplot
ggplot(mt_c, aes(x = feature, y = prop, fill = motif_class)) +
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

# Faceted barplots?
mt_c <- mt_c %>% mutate(motif_class = factor(motif_class, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13-50", ">50")))

ggplot(mt_c, aes(x = motif_class, y = prop, fill = feature)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ feature, ncol = 3, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), name = "Proportion of TRs") +
  scale_x_discrete(name="Motif length") + 
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

# motif enrichement per feature

chisq_enrichment <- function(df) {
  
  # Compute genome-wide proportions of motifs
  genome_props <- df %>%
    group_by(motif_class) %>%
    summarise(total = n(), .groups = "drop") %>%
    mutate(prop = total / sum(total))
  
  genome_props_vec <- setNames(genome_props$prop, genome_props$motif_class)
  
  # Count per feature per motif
  df_counts <- df %>%
    group_by(feature, motif_class) %>%
    summarise(count = n(), .groups = "drop")
  
# Run chi-square per feature

chisq_enrichment <- function(df) {
  
  # Count per feature per motif
  df_counts <- df %>%
    group_by(feature, motif_class) %>%
    summarise(count = n(), .groups = "drop")
  
  # Run chi-square per feature using leave-one-out background
  results <- df_counts %>%
    group_by(feature) %>%
    group_map(~ {
      focal <- .x
      focal_feature <- .y$feature  # get feature being used in that iteration
      
      # all other features minus the focal
      others <- df_counts %>% filter(feature != focal_feature)
      
      # proportions of motifs in "others"
      others_props <- others %>%
        group_by(motif_class) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        mutate(prop = total / sum(total))
      
      test <- chisq.test(
        x = focal$count,
        p = setNames(others_props$prop, others_props$motif_count))
      
      tibble(
        feature = focal_feature,
        classes = focal$motif_class,
        observed = test$observed,
        expected = test$expected,
        residuals = test$residuals
      )}) %>%
    bind_rows()
  
  # Per-class expansion
  per_class <- results %>%
    tidyr::unnest(c(classes, observed, expected, residuals)) %>%
    mutate(
      per_class_p = 2 * (1 - pnorm(abs(residuals))),
      per_class_fdr = p.adjust(per_class_p, method = "fdr"),
      fold_enrichment = observed / expected,
      status = case_when(
        residuals > 2 & per_class_fdr < 0.05 ~ "enriched",
        residuals < -2 & per_class_fdr < 0.05 ~ "depleted",
        TRUE ~ "ns"))
  
  return(per_class)}

results <- chisq_enrichment(df)

write.table(results,"chistest_motif_per_features.txt", col.names=T, row.names=F, quote=F, sep="\t")

# proportion of TR bp per genomic feature
len_chr <- df %>%
  group_by(feature, chr) %>%
  summarize(total_length = sum(tr_len, na.rm = TRUE), .groups = "drop")

df_chr<-read.table("genome_tr_percentages_per_chr.txt", header=T)

df_chr_enrich <- df_chr %>%
  group_by(chr) %>%
  mutate(
    # totals per chromosome
    total_tr_chr = sum(total_tr_length),
    total_len_chr = sum(total_genome_len),
    
    # expected TR length per feature within chromosome
    expected = total_tr_chr * (total_genome_len / total_len_chr),
    
    # standardized residual (Z-score, Poisson approx)
    residual = (total_tr_length - expected) / sqrt(expected),
    
    # two-tailed p-value
    per_feature_p = 2 * (1 - pnorm(abs(residual))),
    
    # fold enrichment
    fold_enrichment = total_tr_length / expected) %>%
  ungroup() %>%
  mutate(
    # FDR correction across all rows
    per_feature_p_fdr = p.adjust(per_feature_p, method = "BH"),
    
    status = case_when(
      residual > 3.3 & per_feature_p_fdr < 0.001 ~ "enriched",
      residual < -3.3 & per_feature_p_fdr < 0.001 ~ "depleted",
      TRUE ~ "ns"))

df_chr_enrich %>%
  ggplot(aes(x = fold_enrichment, y = feature, color = feature)) + 
  geom_boxplot(alpha = 0.8, outlier.shape = NA, show.legend = FALSE, position = position_dodge(width = 0.75)) +
  geom_point(data=df_chr_enrich, aes(x = fold_enrichment, y=feature, color=feature), 
    shape = 16, size = 2, alpha = 0.5, show.legend = FALSE, position = position_jitter(width = 0, height = 0.2)) +
  stat_summary(fun = mean, geom = "point", shape = 23, aes(fill=feature),
    size = 3, show.legend = FALSE, position = position_dodge(width = 0.75)) +
  geom_vline(linetype="dashed", xintercept=1, color="black") +
  scale_x_continuous(name = "Fold enrichment") +
  scale_y_discrete(name = "Genomic feature") +
  scale_color_manual(values=color_mapping) +
  scale_fill_manual(values=color_mapping) +
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

df <- df %>%
  mutate(gc_content = str_count(toupper(motif_seq), "[GC]") / str_length(motif_seq) * 100)

# per feature
df %>%
  ggplot(aes(x = factor(feature), y = gc_content)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, show.legend = FALSE,
               position = position_dodge(width = 2)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1, show.legend = FALSE,
               position = position_dodge(width = 2)) +
  scale_color_manual(values=color_mapping) +
  scale_x_discrete(name = "Feature") +
  scale_y_continuous(name = "GC content") +
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
