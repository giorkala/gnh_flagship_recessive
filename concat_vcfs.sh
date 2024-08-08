#!bin/bash
module load common-apps/bcftools/1.18

if [[ -z $1 ]]; then
    echo "Error: no file-tag was given!"
else
    TAG=$1
    echo "Will processes files with $TAG."
fi

work_dir='/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k/phasing/recessive_encoding/'

# concatenate several chromosome-based files to one
for consq in pLOF pLOF_deleterious synonymous; do
    for mode in additive recessive; do
        echo $work_dir/vcf_ready/chr{1..22}.${mode}_${TAG}.$consq.vcf.gz | tr ' ' '\n' > $work_dir/file_to_concat.txt 
        bcftools concat -f file_to_concat.txt -Oz -o $work_dir/vcf_ready/chrALL.${mode}_${TAG}.$consq.vcf.gz

        echo "Generating a BGEN file for $vcf_prefix.vcf.gz"
        vcf_prefix=$work_dir/vcf_ready/chrALL.${mode}_${TAG}.$consq
        bgen_prefix=$work_dir/bgen_genotypes/chrALL.${mode}_${TAG}.$consq
        
        plink2 --export bgen-1.2 bits=8 \
            --vcf $vcf_prefix.vcf.gz dosage=DS \
            --out $bgen_prefix
        bgenix -index -g $bgen_prefix.bgen

        # fix the samples file
        head -2 $bgen_prefix.sample > $bgen_prefix.tmp
        awk '(NR>2){print 1,$2,$3,$4}' $bgen_prefix.sample >> $bgen_prefix.tmp
        mv $bgen_prefix.tmp $bgen_prefix.sample

    done
done

tar -czvf GNH.recessive_and_additive.chrALL.$TAG.tar.gz vcf_ready/chrALL.additive_$TAG.* vcf_ready/chrALL.recessive_$TAG.*
rm file_to_concat.txt