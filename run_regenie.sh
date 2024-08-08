#!/bin/bash

### run as:
## for consq in pLoF pLoF_damaging damaging_missense synonymous; do bash run_regenie.sh step2bt $consq; done
## for consq in pLOF pLOF_deleterious synonymous; do bash run_regenie.sh step2qt $consq; done

work_dir='/home/ivm/Documents/WORK_RECESSIVE'
TAG="recessive_gnh_flagship"
# TAG="additive_gnh_flagship"

# phenoFile="$work_dir/input/phenotypes.ICD10_3of4.tab"; pheno_tag="icd10_3of4"
phenoFile=$work_dir/input/phenotypes.55cur_qt.tab; pheno_tag="55cur_qt"

if [ $1 = 'step1_single' ]; then
    echo "Run step-1 of regenie for $2 alone."

    regenie \
        --step 1 \
        --bed $work_dir/input/chr1-22_38640-samples_with_GSA_indep-pairwise_500_50_0.2 \
        --covarFile $work_dir/input/covariates.tab \
        --phenoFile $phenoFile \
        --phenoCol $2 \
        --bt \
        --minCaseCount 50 \
        --bsize 500 \
        --lowmem \
        --lowmem-prefix $work_dir/input/regenie_models/$pheno_tag.$2 \
        --gz \
        --verbose \
        --threads 6 \
        --out $work_dir/input/regenie_models/$pheno_tag.$2 \
        # --apply-rint \

elif [ $1 = 'step1qt' ]; then
    echo "Run step-1 of regenie for quantitative traits"
    regenie \
        --step 1 \
        --bed $work_dir/input/chr1-22_38640-samples_with_GSA_indep-pairwise_500_50_0.2 \
        --covarFile $work_dir/input/covariates.tab \
        --phenoFile $phenoFile \
        --apply-rint \
        --bsize 500 \
        --lowmem \
        --lowmem-prefix $work_dir/input/regenie_models/$pheno_tag \
        --gz \
        --verbose \
        --threads 6 \
        --out $work_dir/input/regenie_models/$pheno_tag   \

elif [ $1 = 'step2qt' ]; then
    echo "Run step-2 of regenie testing $2 at $TAG..."
    regenie  \
        --step 2 \
        --bgen $work_dir/bgen_genotypes/chrALL.$TAG.$2.bgen \
        --sample $work_dir/bgen_genotypes/chrALL.$TAG.$2.sample \
        --covarFile $work_dir/input/covariates.tab \
        --phenoFile $phenoFile \
        --pred $work_dir/input/regenie_models/${pheno_tag}_pred.list \
        --out $work_dir/regenie_results_quant/${pheno_tag}.$TAG.$2 \
        --apply-rint \
        --minMAC 3 \
        --bsize 500 \
        --gz \
        --verbose \
        --threads 6 \

elif [ $1 = 'step1bt' ]; then
    echo "Run step-1 of regenie for binary traits"
    regenie \
        --step 1 \
        --bed $work_dir/input/chr1-22_38640-samples_with_GSA_indep-pairwise_500_50_0.2 \
        --covarFile $work_dir/input/covariates.tab \
        --phenoFile $phenoFile\
        --bsize 500 \
        --bt \
        --minCaseCount 50 \
        --gz \
        --verbose \
        --threads 8 \
        --lowmem \
        --lowmem-prefix $work_dir/input/regenie_models/test.tmp.$pheno_tag \
        --out $work_dir/input/regenie_models/$pheno_tag
  
elif [ $1 = 'step2bt' ]; then
    echo "Run step-2 of regenie testing $2 at $TAG..."
    regenie  \
        --step 2 \
        --bgen $work_dir/bgen_genotypes/chrALL.$TAG.$2.bgen \
        --sample $work_dir/bgen_genotypes/chrALL.$TAG.$2.sample \
        --phenoFile $phenoFile \
        --covarFile input/covariates.tab \
        --bsize 1000 \
        --bt \
        --minCaseCount 50 \
        --firth \
        --approx \
        --pThresh 0.01 \
        --af-cc \
        --pred $work_dir/input/regenie_models/${pheno_tag}_pred.list \
        --gz \
        --verbose \
        --threads 8 \
        --out $work_dir/regenie_results_binary/${pheno_tag}.$TAG.$2

elif [ $1 = 'step2female' ]; then
    echo "Run step-2 of regenie testing $2 at $TAG for FEMALEs..."
    awk '{print "1",$1}' /genesandhealth/red/SamHodgson/input/GRMs/FEMALE_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt > $work_dir/samples.females.txt
    regenie  \
        --step 2 \
        --bgen $work_dir/bgen_genotypes/chrALL.$TAG.$2.bgen \
        --sample $work_dir/bgen_genotypes/genotypes.sample \
        --keep $work_dir/samples.females.txt \
        --phenoFile input/phenotypes_20220615.tab \
        --phenoCol $pheno_all \
        --covarFile input/covariates.tab \
        --bsize 1000 \
        --bt \
        --minCaseCount 50 \
        --firth \
        --approx \
        --pThresh 0.01 \
        --af-cc \
        --pred $work_dir/input/regenie_models/test_brava.allAs_pred.list \
        --gz \
        --verbose \
        --threads 4 \
        --out $work_dir/regenie_results_binary/test_females.$TAG.$2

else
    echo "Wrong arguments!"
fi

# end of script