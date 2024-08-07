#!bin/bash
module load common-apps/bcftools/1.18

if [[ -z $1 ]]; then
    echo "Error: no file-tag was given!"
else
    TAG=$1
    echo "Will processes files with $TAG."
fi

prefix='/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k/phasing/recessive_encoding/vcf_ready'

# concatenate several chromosome-based files to one
#for consq in pLoF pLoF_damaging damaging_missense_or_protein_altering synonymous; do
for consq in pLOF pLOF_deleterious deleterious_missense synonymous; do
    for mode in additive recessive; do
        echo $prefix/chr{1..22}.${mode}_${TAG}.$consq.vcf.gz | tr ' ' '\n' > file_to_concat.txt 
        bcftools concat -f file_to_concat.txt -Oz -o $prefix/chrALL.${mode}_${TAG}.$consq.vcf.gz
    done
done

tar -czvf GNH.recessive_and_additive.chrALL.$TAG.tar.gz vcf_ready/chrALL.additive_$TAG.* vcf_ready/chrALL.recessive_$TAG.*
rm file_to_concat.txt
