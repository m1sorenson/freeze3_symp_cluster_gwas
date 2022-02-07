# Freeze3 symptom cluster GWAS meta-analysis
This repository contains code to run the symptom cluster (re-experiencing, hyperarousal) gwas and meta-analysis on PGC PTSD freeze 3. This is different from the overall meta-analysis in a few ways:
1. This is only done on studies we have genotype access to, we don't have the summary statistics data from other studies
2. The number of samples is varied depending on what symptom data we have (there is some missing data)
3. We don't have a meta-analysis yet for symptom C and E because in some studies these are combined into one symptom, and in others they are two separate symptoms

## Files
1. *01_f3_gwas_symp_clusters.sh* - this is the main bash script with the blocks of code to run the GWAS's and the meta-analysis
2. *dosage_locations_symp_clusters.csv* - this has the data for which studies to run the symp cluster gwas for, and the pheno files to use
3. *f3_symp_PHENOB_jan26_2022.mi* - this is the METAL script to run the meta-analysis for symptom B (re-experiencing)
4. *f3_symp_PHENOD_jan26_2022.mi* - this is the METAL script to run the meta-analysis for symptom D (hyperarousal)
5. *run_meta_v2_loo_v2.slurm* - this is the SLURM bash script used to run METAL (meta-analysis)
6. *run_trauma_gwas_symp_clusters.slurm* - this is the SLURM bash script used to run the study-level GWAS's
7. *symp_pheno_fixing_and_counting.ipynb* - this is a jupyter notebook (open the file in GitHub or with `jupyter notebook` to see it in nice formatting) used to transform the study-level phenotype txt files into pheno-specific .pheno files to be passed into the study-level GWAS step, and count the number of samples on a per-study per-pheno basis. This also fixes the issues in the EHVP and ONGB FID and IID's.

## Usage
### Prerequisites
1. Must have a folder named `pheno`, with phenotype .txt files with the columns `FID IID pheno1 pheno2 pheno3`, and the file name should be in the form [STUDY_NAME]\_[CAPS4/CAPS5/PCL4/PCL5]\_[BCD/BCDE].txt, e.g. *1_MRSC_C_MAX_PCL4_BCD.txt*
2. Must have metal and plink2 available in your $PATH

### Run it
#### Set up & GWAS's
1. To start, make a few directories: `mkdir errandout metal_results metal_scripts results_cat`
2. Next run the block of code in *01_f3_gwas_symp_clusters.sh* that begins with `### 1)Study Level Analysis steps:`
3. Once each of the studies has finished running, check the `[study].e` files in the *errandout* folder for any errors

#### P-value distribution plotting versus broad
4. Now, you can check that the GWAS's ran correctly and make sense by plotting the P-value distributions for each study
5. First, load python and install the required packages:
```
module load 2021
module load matplotlib/3.4.2-foss-2021a
pip3 install -r requirements.txt
```

6. Now, run run_qq.sh to produce plots in the `plots` folder for each study:
```
sbatch --time=01:00:00 --error errandout/plot_qq.e --output errandout/plot_qq.o run_qq.sh
```

#### Meta-Analysis and results
7. Next, run the block of code in *01_f3_gwas_symp_clusters.sh* that begins with `### 2)Meta-Analysis step:`
8. Check the `f3_symp_PHENOB...mi_errorlogs` and `f3_symp_PHENOD...mi_errorlogs` for any errors
9. If all the above steps ran without errors, the final metal results will be in the metal_results folder, run the block starting with `### 3) Combine METAL results and generate final output` to get the final concatenated METAL files

### Troubleshooting
Some things to check if there are any errors:
1. Check that the .pheno files are the correct format, i.e. `FID IID pheno`
2. Check that the .pheno IDs match those that are in whatever dosage files you are using
3. Check that the *dosage_locations_symp_clusters.csv* has the correct pheno files in the PHENO column, check for typos
4. If there is an error that says something like `CANCELLED DUE TO TIME LIMIT`, increase the value of the **timecode** column in *dosage_locations_symp_clusters.csv* for whatever study failed (the advantage of having lower time codes is that the SLURM scheduler prioritizes lower time-limit jobs first)
5. If you have any `Segmentation Fault` errors, they are likely due to incorrect formatting or values in the .pheno files, dropping rows with NA values is an easy first step to solve this.
