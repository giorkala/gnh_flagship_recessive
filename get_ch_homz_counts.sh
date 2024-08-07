#!bin/bash

if [[ -z $1 ]]; then
    echo "Error: no file-tag was given!"
else
    TAG=$1
    echo "Will processes files with $TAG."
fi

work_dir="/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k/phasing/recessive_encoding"

# print summaries:
rm -f $work_dir/called_chets.$TAG.summary.csv
for C in {1..22}; do
    for bial in chet hom cis; do
        #for consq in pLoF pLoF_damaging damaging_missense_or_protein_altering synonymous; do
	for consq in pLOF pLOF_deleterious deleterious_missense synonymous; do
            tmp=$(grep $bial $work_dir/sandbox/chr$C.gen_$TAG.$consq.txt | wc -l | cut -d' ' -f1)
            echo "$C,$bial,$consq,$tmp" >> $work_dir/called_chets.$TAG.summary.csv
        done
    done
done

