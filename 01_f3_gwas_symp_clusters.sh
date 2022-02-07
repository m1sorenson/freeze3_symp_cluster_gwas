#!/bin/bash

### 1)Study Level Analysis steps:

##Group 1 and 2 GWAS only: GWAS
cov=pcs

#Make sure time codes are correct
IFS=$'\n'
for line in $(tail -n+2 dosage_locations_symp_clusters.csv)
do
  study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $4}')
  study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')
  phenofile=$(echo $line | awk 'BEGIN{FS=","}  {print $3}')

  ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $6}')
  timecode=$(echo $line | awk 'BEGIN{FS=","} {print $7}')
  exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $8}')

  outkey=$(echo $line | awk 'BEGIN{FS=","}  {print $9}')

  if [[ $exclude -ne 1 ]]
  then
    echo gwas for $outkey
    sbatch --time=$timecode --error errandout/$outkey.e --output errandout/$outkey.o \
    --export=ALL,phenofile="$phenofile",study="$study_1",study_2="$study_2",ancgroup="$ancgroup",outkey="$outkey",cov="$cov" run_trauma_gwas_symp_clusters.slurm -D $WORKING_DIR
  fi

done


### 2)QQ-plots:
#First, load python and install python packages
module load 2021
module load matplotlib/3.4.2-foss-2021a
pip3 install -r requirements.txt

#Now run plots
sbatch --time=01:00:00 --error errandout/plot_qq.e --output errandout/plot_qq.o run_qq.sh


### 3)Meta-Analysis step:
for chr in {1..22} X; do
  # Add no extension to genotyped results and add _broad extension to summary stat files
  sed s/{CHR_NUM}/${chr}/g f3_symp_PHENOB_jan26_2022.mi > metal_scripts/f3_symp_PHENOB_jan26_2022_chr${chr}.mi
  sed s/{CHR_NUM}/${chr}/g f3_symp_PHENOD_jan26_2022.mi > metal_scripts/f3_symp_PHENOD_jan26_2022_chr${chr}.mi
done

#PHENO B and PHENO D
ls metal_scripts/*PHENOB* > metafilelist.txt
ls metal_scripts/*PHENOD* >> metafilelist.txt
dataset=f3_symp_clusters_jan26_2022
sbatch -t 01:00:00  --error errandout/"$dataset".e --output errandout/"$dataset".o   --export=ALL,metafile=metafilelist.txt -D /home/pgca1pts/freeze3_symp_cluster_gwas run_meta_v2_loo_v2.slurm

### 3) Combine METAL results and generate final output

PHENO_B_N=8371
PHENO_D_N=8311
percentN=0.8

cat metal_results/eur_ptsd_symp_cluster_PHENO_B_jan26_2022_*.tbl | awk -v totalN=$PHENO_B_N -v percentN=$percentN 'BEGIN{OFS="\t"}{ if(NR==1){ print "#"$1,$2,$3,$4,$5,$6,$10,$11,$12} else if ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN) print $1,$2,$3,$4,$5,$6,$10,$11,$12}' | grep -v : | LC_ALL=C sort -g -k 9 > metal_results/eur_ptsd_symp_cluster_PHENO_B_jan26_2022.tbl
cat metal_results/eur_ptsd_symp_cluster_PHENO_D_jan26_2022_*.tbl | awk -v totalN=$PHENO_D_N -v percentN=$percentN 'BEGIN{OFS="\t"}{ if(NR==1){ print "#"$1,$2,$3,$4,$5,$6,$10,$11,$12} else if ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN) print $1,$2,$3,$4,$5,$6,$10,$11,$12}' | grep -v : | LC_ALL=C sort -g -k 9 > metal_results/eur_ptsd_symp_cluster_PHENO_D_jan26_2022.tbl
