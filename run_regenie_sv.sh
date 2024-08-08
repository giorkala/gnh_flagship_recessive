#!/bin/bash

### NOTES ###
# 1. This analysis should be similar to Teng's (recessive effects using imputed genotypes), but with higher accuracy due to WES and more stringent QC.
# 2. This is different than the Burden analysis as here we dont need to subset to GSA samples, and so we can use up to 44k individuals
# 3. Remember that sample labels MUST contain underscores here, so thats another difference to the burden analysis

dir_base='/home/ivm/Documents/WORK_RECESSIVE'
dir_out='/home/ivm/Documents/WORK_RECESSIVE/SINGLE_VARIANT'
# PGEN='/genesandhealth/consortiumpriorityperiod-library-red/Callsets/2023_02_44kCallset/release_2023-JUL-07/chr1..22_hard_filters_stringent_CR.95_44028samples.tidied.vep105'
PGEN='/home/ivm/Documents/WORK_RECESSIVE/SINGLE_VARIANT/WES_44kCallset.mac5'
# plink2 --pfile chr1..22_hard_filters_stringent_CR.95_44028samples.tidied.vep105 --mac 5 --make-pgen --out ~/Documents/WORK_RECESSIVE/SINGLE_VARIANT/WES_44kCallset.mac5
# bcftools view $VCF -S samples.phenotyped.txt -R include.txt -Ob -o WES_44kCallset.chr18.pruned.bcf
covariates="/home/ivm/Documents/GNH_44k_WES/covariates_PC1-20_WES_44k_MAF1pc_44028_dvh.txt"

TAG="WES_44k"
pheno_tag="55cur_qt"
Ncpu=5

mkdir -p $dir_out/regenie_models $dir_out/regenie_results_quant $dir_out/regenie_results_binary

if [ $1 = 'step1qt' ]; then
    echo "Run step-1 of regenie for quantitative traits"
    regenie \
        --step 1 \
        --bed $dir_out/chr1-22_44028-samples_indep-pairwise_500_50_0.2 \
        --covarFile $dir_base/brava_inputs/covariates.tab \
        --phenoFile $dir_base/brava_inputs/phenotypes.$pheno_tag.tab \
        --apply-rint \
        --bsize 500 \
        --lowmem \
        --lowmem-prefix $dir_out/regenie_models/$pheno_tag.tmp \
        --gz \
        --verbose \
        --threads $Ncpu \
        --out $dir_out/regenie_models/$pheno_tag

elif [ $1 = 'step1qt_single' ]; then
    echo "Run step-1 of regenie for trait $2"
    regenie \
        --step 1 \
        --bed $dir_out/chr1-22_44028-samples_indep-pairwise_500_50_0.2 \
        --covarFile $covariates \
        --catCovarList S1QST_Gender \
        --covarExcludeList OrageneID,PseudoNHS_2023_07_14 \
        --phenoFile $dir_base/brava_inputs/phenotypes.$pheno_tag.tab \
        --phenoCol $2 \
        --apply-rint \
        --bsize 500 \
        --lowmem \
        --lowmem-prefix $dir_out/regenie_models/$pheno_tag.tmp.$2 \
        --gz \
        --verbose \
        --threads $Ncpu \
        --out $dir_out/regenie_models/$pheno_tag.$2
    
    # then run the following
    # brava_inputs/phenotypes.66quant.list | while read x; do cat SINGLE_VARIANT/regenie_models/66quant.$x\_pred.list; done > 66quant_pred.list

elif [ $1 = 'step2qt' ]; then
    echo "Run step-2 of regenie for trait $2 and $TAG..."
    regenie  \
        --step 2 \
        --pgen $PGEN \
        --extract $dir_out/WES_44kCallset.stringent_AC_Hom_gt3.chrALL.vars \
        --phenoFile $dir_base/brava_inputs/phenotypes.$pheno_tag.tab \
        --phenoCol $2 \
        --pred $dir_out/regenie.WES_44k.55cur_qt_pred.list \
        --covarFile $covariates \
        --catCovarList S1QST_Gender \
        --covarExcludeList OrageneID,PseudoNHS_2023_07_14 \
        --test recessive \
        --out $dir_out/regenie_results_quant/$TAG.chrALL_pruned.REC \
        --apply-rint \
        --bsize 500 \
        --gz \
        --verbose \
        --threads $Ncpu
            # --chr $2 \
        # --pred $dir_base/brava_inputs/regenie_models/${pheno_tag}_pred.list \

elif [ $1 = 'step1bt' ]; then
    batch=$2
    echo "Run step-1 of regenie for binary traits and batch $batch"
    regenie \
        --step 1 \
        --bed $dir_out/chr1-22_44028-samples_indep-pairwise_500_50_0.2 \
        --covarFile /home/ivm/Documents/GNH_44k_WES/covariates_PC1-20_WES_44k_MAF1pc_44028_dvh.txt \
        --catCovarList S1QST_Gender \
        --covarExcludeList OrageneID,PseudoNHS_2023_07_14 \
        --phenoFile $dir_base/brava_inputs/phenotypes.ICD10_$batch.tab  \
        --out $dir_out/regenie_models/test_ICD10_$batch \
        --lowmem-prefix $dir_out/regenie_models/test_ICD10_$batch.tmp \
        --lowmem \
        --bsize 500 \
        --bt \
        --minCaseCount 50 \
        --gz \
        --verbose \
        --threads $Ncpu
        # --phenoCol AL.amyloidosis,Abdominal.aortic.aneurysm..AAA.,Acute.appendicitis..AcApp.,Acute.lymphoid.leukemia,Adenomyosis,Age.related.macular.degeneration,Alcohol.use.disorder,Alopecia.Areata,Anorexia.nervosa,Aortic.stenosis,Asthma..Asthma.,Atopic.diseases,Atrial.Fibrillation,Attention_deficit_hyperactivity_disorder.ADHD.,Autism.spectrum.disorder..ASD. \
        # --phenoCol Type_2_Diabetes \
        # --covarExcludeList PseudoNHS_2022_07_20 \
        # --covarFile /genesandhealth/red/KlaudiaWalter/data/pca/PC1-20_WES_44k_MAF1pc_38640-samples.txt \
        # --phenoFile /genesandhealth/library-red/genesandhealth/phenotypes_curated/version005_2022_06/2022-06-15_big_regenie_phenoFile.BroadID.tab \
        # --phenoCol Enthesopathies__synovial_disorders,Unspecified_or_Rare_Diabetes,GNH0242_Type_2_Diabetes_narrow,GNH0018_EssentialHypertension,GNH0007_Type2Diabetes,Type_2_Diabetes,GNH0249_NormalUncomplicatedTermDelivery,Hypertension,Allergic_and_chronic_rhinitis,Dermatitis_atopc_contact_other_unspecified \

elif [ $1 = 'step2bt' ]; then
    batch=$2
    echo "Run step-2 of regenie testing for ICD10 batch $batch"
    regenie  \
        --step 2 \
        --pgen $PGEN \
        --extract $dir_out/WES_44kCallset.stringent_AC_Hom_gt3.chrALL.vars \
        --covarFile $covariates \
        --catCovarList S1QST_Gender \
        --covarExcludeList OrageneID,PseudoNHS_2023_07_14 \
        --phenoFile $dir_out/phenotypes.ICD10_$batch.tab\
        --pred $dir_out/regenie.WES_44k.icd10_${batch}_pred.list \
        --out $dir_out/regenie_results_binary/$TAG.ADD.icd10_${batch} \
        --minMAC 10 \
        --bsize 1000 \
        --bt \
        --minCaseCount 50 \
        --firth \
        --approx \
        --pThresh 0.01 \
        --af-cc \
        --gz \
        --verbose \
        --threads $Ncpu \
        # --test recessive \

elif [ $1 = 'prune' ]; then
    echo "Identify variants with at least 5 alt-homozygous carriers on chrom-$2."
    VCF="/genesandhealth/consortiumpriorityperiod-library-red/Callsets/2023_02_44kCallset/release_2023-JUL-07/chr$2_hard_filters.tidied.vep105.vcf.gz"
    # bcftools view $VCF -S $dir_out/samples.phenotyped.txt | bcftools query -i'GT="1/1"' -f'[%CHROM:%POS:%REF:%ALT\n]' | uniq -c | awk '($1>4){print $2}' | tr ':' '_' > $dir_out/chr$2.homz.pruned
    bcftools query $VCF -i'stringent_AC_Hom>4' -f'%CHROM:%POS:%REF:%ALT\t%AC_Hom\n' > $dir_out/chr$2.AC_Hom.txt
    # the output can then be used with "--extract" in regenie.

elif [ $1 = 'merge' ]; then
    phenotypes=$(head -1 brava_inputs/phenotypes.bottom20quant.tab | tr '\t' '\n' | awk '(NR>2)')
    # phenotypes=$(head -1 brava_inputs/phenotypes.tengsICD10.tab | tr '\t' '\n' | awk '(NR>2)')
    echo "Will merge any sumstats with prefix $2. Phenotypes found:"
    echo $phenotypes
    for pheno in $phenotypes; do
        newFile="$2.$pheno.chrALL.regenie"
        zcat $2.chr1_$pheno.regenie.gz > $newFile
        for C in {2..22}; do zcat $2.chr${C}_${pheno}.regenie.gz | awk '(NR>1)' >> $newFile; done
        gzip $newFile
    done
else
    echo "Wrong arguments!"
fi