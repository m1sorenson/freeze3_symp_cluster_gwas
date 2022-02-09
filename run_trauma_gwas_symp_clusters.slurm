#!/bin/bash
# module load plink2

#Only has minor modifications to be run for f3, namely changing datadir, reslength -le, and looking for chr10_000_020 to account for chunk size
#have also added sex chromosome. Sex chromosomes are imputed separately by sex so I do some stuff to account for this, namely dosage list command
#Don't really worry about PAR, ricopili doesn't really impute it (some overlap, have to worry about this in filtering).
# see this for https://www.cog-genomics.org/plink/1.9/assoc x chromosome handling. may have to reduce to plink 1.9 and produce output separately!!
# https://www.cog-genomics.org/plink/2.0/assoc#glm ill only hae to do this if I can't figure out how to combine male and female x imputation files!!
#may just want to combine them  prior to any of this! if plink2 doesnt take dosage-lsit that sucks!
#I think this may be the best option...

echo $study
echo "--export=phenofile=${phenofile},study=${study_1},study_2=${study_2},ancgroup=${ancgroup},outkey=${outkey},cov=${cov}"

datadir=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/${study}/qc1

mkdir ${TMPDIR}/${study}_${study_2}

#for file in $( ls $datadir | grep .gz$ | sed 's/.gz//g' | grep  chr10_000_020) # TEST CASE
for file in $( ls $datadir | grep .gz$ | sed 's/.gz//g' | grep -v "chrX.mal" | grep -v "chrX.fem")
do
  outname=$(echo $file | awk 'BEGIN{FS="ch.fl.chr"}{print $2}'  | sed 's/.out.dosage//g') #get output filename...
  echo $outname

  plink2 --import-dosage "$datadir"/"$file".gz --fam "$datadir"/"$file".fam --map "$datadir"/"$file".map --pheno pheno/$phenofile --pheno-col-nums 3 \
  --covar pheno/p2_"$study"_"$ancgroup"_"$study_2"_"$cov".cov --maf 0.01 \
  --glm cols=chrom,pos,ref,alt,ax,a1freq,machr2,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p   hide-covar firth-fallback  \
  --covar-variance-standardize --threads 11 --ci 0.95 $plink_flag \
  --out "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_chr"$outname"_"$cov"
done

#Auto detect pheno levels evaluated and loop over these...
for phenoset in $(ls "$TMPDIR"/"$study"_"$study_2" | awk 'BEGIN{FS="."}{print $2}' | grep -v log | grep -v psam | grep -v pvar | grep -v pgen | sort | uniq)
do
  reslength=$(ls "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov"."$phenoset".glm.* | wc -l) #Get linecount to make sure thats good!!
  echo "Working on "$study"_"$study_2"_"$ancgroup"_*_"$cov"."$phenoset""
  if [ $reslength -le 130 ]
  then
    echo "Fewer than 130 results files --  $reslength -- something maybe be wrong "
  fi

  if [ $( echo $phenoset | grep -c "Case") -eq 1 ]
  # does this result in header being added to file multiple times
  then
    cat "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov"."$phenoset".glm.* | awk '{if (NR == 1 || $3 != "ID") print}'  | sed 's/#//g'  | awk '{ if (NR==1 || ($8 >= 0.01 && $8 <= 0.99 && $11 > 0.6)) print }' | LC_ALL=C sort -k 1,2n | gzip > $TMPDIR/"$study"_"$study_2"_"$ancgroup"_"$cov"."$phenoset".assoc.gz
  else
    cat "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov"."$phenoset".glm.* | awk '{if (NR == 1 || $3 != "ID") print}'  | sed 's/#//g'  | awk '{ if (NR==1 || ($8 >= 0.01 && $8 <= 0.99 && $9 > 0.6)) print }' | LC_ALL=C sort -k 1,2n | gzip > $TMPDIR/"$study"_"$study_2"_"$ancgroup"_"$cov"."$phenoset".assoc.gz
  fi

  #Now scale the BETA based on whether the study is CAPS4, CAPS5, PCL4, PCL5, or PSS4
  factor=1
  if [[ $phenoset =~ [cC][aA][pP][sS]-?_?4 ]]; then
    if [[ $phenoset =~ _b_ ]]; then
      factor=40
    elif [[ $phenoset =~ _d_ ]]; then
      factor=40
    fi
  elif [[ $phenoset =~ [cC][aA][pP][sS]-?_?5 ]]; then
    if [[ $phenoset =~ _b_ ]]; then
      factor=20
    elif [[ $phenoset =~ _d_ ]]; then
      factor=28
    fi
  elif [[ $phenoset =~ [pP][cC][lL]-?_?4 ]]; then
    if [[ $phenoset =~ _b_ ]]; then
      factor=25
    elif [[ $phenoset =~ _d_ ]]; then
      factor=25
    fi
  elif [[ $phenoset =~ [pP][cC][lL]-?_?5 ]]; then
    if [[ $phenoset =~ _b_ ]]; then
      factor=20
    elif [[ $phenoset =~ _d_ ]]; then
      factor=28
    fi
  elif [[ $phenoset =~ [pP][sS][sS]-?_?4 ]]; then
    if [[ $phenoset =~ _b_ ]]; then
      factor=15
    elif [[ $phenoset =~ _d_ ]]; then
      factor=15
    fi
  fi
  zcat "$TMPDIR"/"$study"_"$study_2"_"$ancgroup"_"$cov"."$phenoset".assoc.gz | \
    awk -v factor=$factor 'BEGIN{OFS="\t"}{if (NR > 1){ $12 = $12 / factor } print $0}' | \
    gzip > results_cat/"$outkey"_"$cov"_scaled.assoc.gz

  #And make a split one to make this easy for us...
  for chr in {1..22} X
  do
    zcat results_cat/"$outkey"_"$cov"_scaled.assoc.gz | awk -v chr=$chr '{if (NR == 1 || $1==chr) print}'  | gzip > results_cat/"$outkey"_"$cov"_scaled.assoc.chr${chr}.gz
  done

  echo "completed"
  #rm -f results/"$study"_*_"$phenolevel"_"$ancgroup"_"$study_2"_"$cov"."$pheno".glm.linear
done
