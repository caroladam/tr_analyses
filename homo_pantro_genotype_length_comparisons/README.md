## Comparisons of mean allele length for homologous TRs between humans and chimpanzees

### files:
- `homo_pantro_homolog_trs.txt`: contains genomic location and summary data for homologous TRs between humans and chimpanzees, including mean, variance and the 5th-95th quantile of allele lengths for each species.

### scripts:
- `scripts/get_percentiles.py`: calculates summary stats from TR genotype data (mean, variance, 5th and 95th percentiles of allele length per locus)
- `scripts/plot_shared_tr_lens.R`: generate plots for mean allele length comparison of homologous TRs between humans and chimps (input `homo_pantro_homolog_trs.txt`)
