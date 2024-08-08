#!/bin/bash

## script to run permutation tests to determine a threshold for statistical significance
## June 2024

work_dir='/home/ivm/Documents/WORK_RECESSIVE'
TAG="recessive_gnh_flagship"
pheno=$2
consq=$3

# for each trait, permute the .samples file and test for association for all masks.
# these will result in 439+55 permutations, based on the phenotypes for the G&H flagship analysis.

# permute the samples file - once per phenotype
samples_permuted="$work_dir/bgen_genotypes/permuted.$pheno.genotypes.sample"
if [ ! -f $samples_permuted ]; then
    echo "Will generate a permuted sample index."
    head -2 $work_dir/bgen_genotypes/genotypes.sample > $samples_permuted
    tail -n+3 $work_dir/bgen_genotypes/genotypes.sample | shuf >> $samples_permuted
fi
mkdir -p regenie_results_permuted

if [ $1 = 'step2qt' ]; then
    echo "Run step-2 of regenie testing $pheno with $consq..."

    regenie  \
        --step 2 \
        --bgen $work_dir/bgen_genotypes/chrALL.$TAG.$consq.bgen \
        --sample $samples_permuted \
        --covarFile $work_dir/brava_inputs/covariates.tab \
        --phenoFile $work_dir/brava_inputs/phenotypes.55cur_qt.tab \
        --phenoCol $pheno \
        --pred $work_dir/brava_inputs/regenie_models/55cur_qt_pred.list \
        --out $work_dir/regenie_results_permuted/55cur_qt.permuted.$TAG.$consq \
        --qt \
        --apply-rint \
        --minMAC 3 \
        --bsize 500 \
        --gz \
        --verbose \
        --threads 4 \

elif [ $1 = 'step2bt' ]; then
    batch=$4
    echo "Run step-2 of regenie testing $batch-$pheno with $consq..."
    outfile="$work_dir/regenie_results_permuted/icd10.permuted.$TAG.${consq}_${pheno}.regenie.gz"
    if [ ! -f $outfile ]; then
        regenie  \
            --step 2 \
            --bgen $work_dir/bgen_genotypes/chrALL.$TAG.$consq.bgen \
            --sample $samples_permuted \
            --covarFile $work_dir/brava_inputs/covariates.tab \
            --phenoFile $work_dir/brava_inputs/phenotypes.ICD10_$batch.tab \
            --phenoCol $pheno \
            --pred $work_dir/brava_inputs/regenie_models/icd10_${batch}_pred.list \
            --out $work_dir/regenie_results_permuted/icd10.permuted.$TAG.$consq \
            --bt \
            --minCaseCount 50 \
            --firth \
            --approx \
            --pThresh 0.01 \
            --bsize 1000 \
            --gz \
            --threads 4 \
            # --af-cc \   --verbose \
    else
        echo "File $outfile already exists, moving on."
    fi
else
    echo "Wrong arguments!"
fi

## submit as 
# awk 'NR<25' brava_inputs/phenotypes.55cur_qt.list | while read pheno; do for consq in pLOF pLOF_deleterious synonymous; do bash run_regenie_perm.sh step2qt $pheno $consq; done; done
# awk 'NR==1' brava_inputs/phenotypes.ICD10_top439.list | while read batch pheno; do for consq in pLOF pLOF_deleterious synonymous; do bash run_regenie_perm.sh step2bt $pheno $consq $batch; done; done
