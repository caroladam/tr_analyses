## Mutation rate estimates from human-chimp TR divergence data

### script:
- `estimate_mutation_rates.R`: Calculate TR mutation rate using genetic distance-based methods and compare to de novo pedigree estimates.

**Inputs:**
- Tab-delimited files containing de novo TR mutation rate estimates from pedigree data from Porubsky et al. (2025) - doi:10.1038/s41586-025-08922-2.
-  A tab-delimited file with TR copy number from humans and chimpanzees.
  
**Outputs:**
- Genetic distance-based TR mutation rate estimates per motif category.
- Line plot comparing distance and pedigree-based TR mutation rate estimates across motif lengths.
