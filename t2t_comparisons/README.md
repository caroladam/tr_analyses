## Analyzing TR data in T2T primate genomes

### scripts:
- `get_chm13_appris_features.sh`: Determines CHM13 genomic features from an annotation file: [`chm13.v2 annotation`](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3).
Outputs are BED files containing genomic coordinates for CDS, 5' and 3' UTRs, promoters (1 kb upstream of the TSS), and introns. We filter annotations to include only APPRIS principal annotations.
- `chm13_plots.R`: R code used to generate all CHM13 plots: proportion of chromosomes made up by TRs, distribution of TRs per motif length for each genomic feature, and fold-enrichment of TRs per genomic feature compared to genome-wise composition.
- `reference_comparisons.R`: R code used to generate all T2T reference comparisons between human and each NHP species: motif length distribution per species, count of homologous TRs between human and NHPs, scatterplots of human x NHP TR lengths for shared TRs, scatterplot of TR length correlation between human and NHP and the divergence time, and comparison of TR length per motif length between human and NHP.
- `homologous_tr_density.sh`: Estimate TR density in human x NHP homologous regions in 1Mb windows
- `plot_homologous_tr_density.R`: R code to generate heatmaps of normalized TR densities in homologous 1Mb windows of human and NHP species. 
