#!/bin/bash

# ---
# Title: homologous_tr_density.sh
# Date: 2026
# Author: Adam, Carolina L.
# Purpose: Estimate TR density in homo x NHP homologous regions
# ---

# generate 1Mb non-overlapping windows
bedops --chop 1000000 karyotipe_chm13.bed > chm13_window_1mb.bed
bedtools sort -i chm13_window_1mb.bed > chm13_window_1mb_sorted.bed

# get TR density
bedtools coverage -a chm13_window_1mb_sorted.bed -b homo_catalog.no_overlaps.bed > homo_1mb_coveraged.bed
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' homo_1mb_coveraged.bed > homo_density.bed

# add window ID before lifting over
awk 'BEGIN{OFS="\t"}
{
 id=$1":"$2"-"$3
 print $1,$2,$3,id
}' chm13_window_1mb_sorted.bed > homo_windows_ID.bed

# Lift homo windows to NHP coordinates
liftOver -minMatch=0.1 homo_windows_ID.bed homo_to_pantro.chain homo_windows_ID_lifted.bed unmapped.bed

# Filter out weird mapping - 0.8Mb for all spp except macaque and gibbon (0.2Mb)
awk 'BEGIN{OFS="\t"}
{
 len=$3-$2
 if(len>=800000 && len<=2000000)
     print
}' homo_windows_ID_lifted.bed > homo_windows_ID_lifted.filtered.bed


# Get coverage of NHP TRs in the lifted windows
bedtools coverage -a homo_windows_ID_lifted.filtered.bed -b pantro/pantro_catalog.no_overlaps.bed > lifted_vs_pantro_density.bed

# make look up table of homo density with the window IDs
awk 'BEGIN{OFS="\t"}
{
 id=$1":"$2"-"$3
 print id,$4
}' homo_density.bed > homo_density_lookup.txt
