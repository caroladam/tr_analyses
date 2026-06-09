## Expression divergence x TR divergence

### files:
- `NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt`: normalized expression values from a primate-specific set of orthologs obtained from [Brawand et al. 2011](https://doi.org/10.1038/s41586-022-04510-w)
- `pollen_expr_results.txt`: differential gene expression results obtained from [Pollen et al. 2019](https://doi.org/10.1016/j.cell.2019.01.017)
- `pollen_derived_genes.txt`: list of genes with human-specific regulatory changes during cerebral development obtained from [Pollen et al. 2019](https://doi.org/10.1016/j.cell.2019.01.017)

### scripts:
- `scripts/get_logfc_expr.sh`: calculates log2 fold-change expression divergence from normalized expression values between humans and chimps.
- `scripts/plot_expr_div.R`: generate plots for expression divergence x TR divergence
