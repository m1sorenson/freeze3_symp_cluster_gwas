#!/bin/bash

## ------------------------------------------------------------------------
# this loop iterates through each study with symptom clusters data and plots
# the symptom cluster GWAS p value quantiles vs the broad p value quantiles
## ------------------------------------------------------------------------

## MUST LOAD THESE MODULES AND INSTALL PACKAGES FOR THIS TO WORK:
# module load 2021
# module load matplotlib/3.4.2-foss-2021a
# pip3 install -r requirements.txt


# Script
mkdir tmp
broaddir=/home/maihofer/freeze3_gwas/results_cat

IFS=$'\n'
for line in $(tail -n+2 dosage_locations_symp_clusters.csv | head -n 1)
do
  study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $4}')
  study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')
  ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $6}')
  outkey=$(echo $line | awk 'BEGIN{FS=","}  {print $9}')
  exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $8}')
  if [[ $exclude -eq 0 ]]; then
    file1=results_cat/${outkey}_pcs.assoc.gz
    file2=$(ls ${broaddir}/${study_1}_${study_2}_${ancgroup}*Continuous.assoc.gz)
    echo $file1 $file2
    SNP_COL_1=$(zcat $file1 | head -n 1 | \
      awk '{ for(i=1; i<=NF; i++){ \
        if($i ~ "([Ss][Nn][Pp])|([Rr][sS])|([Ii][Dd])|([Mm][Aa][Rr][Kk][Ee][Rr])"){ print i } \
      } }')
    SNP_COL_2=$(zcat $file1 | head -n 1 | \
      awk '{ for(i=1; i<=NF; i++){ \
        if($i ~ "([Ss][Nn][Pp])|([Rr][sS])|([Ii][Dd])|([Mm][Aa][Rr][Kk][Ee][Rr])"){ print i } \
      } }')
    P_COL_1=$(zcat $file1 | head -n 1 | \
      awk '{ for(i=1; i<=NF; i++){ \
        if($i ~ "(^[Pp])(($)|([_-]?[Vv][Aa][Ll]))"){ print i } \
      } }')
    P_COL_2=$(zcat $file1 | head -n 1 | \
      awk '{ for(i=1; i<=NF; i++){ \
        if($i ~ "(^[Pp])(($)|([_-]?[Vv][Aa][Ll]))"){ print i } \
      } }')
    zcat $file1 | awk -v snp_col=$SNP_COL_1 -v p_col=$P_COL_1 'BEGIN{OFS=","}{if(NR==1){print "#SNP,P"}else{print $snp_col,$p_col}}' | LC_ALL=C sort -k 1 | gzip -c > tmp/${outkey}_broad.csv.gz
    zcat $file2 | awk -v snp_col=$SNP_COL_2 -v p_col=$P_COL_2 'BEGIN{OFS=","}{if(NR==1){print "#SNP,P"}else{print $snp_col,$p_col}}' | LC_ALL=C sort -k 1 | gzip -c > tmp/${outkey}_sympclust.csv.gz
    #make a similar loop to this locally

    #plot_qq.py
    python3 plot_qq.py ${outkey} tmp/${outkey}_broad.csv.gz tmp/${outkey}_sympclust.csv.gz
  fi
done
