#!/bin/bash

# requirements: call_chets binary, gene_map file, genotypes file
# output: biallelic carriers for each gene, for each consequence type

# no need to pull out genotypes from the BCFs, as this is done already.
# see https://github.com/giorkala/GNH_BRaVa_dev_misc/blob/main/recessive_encoding/run_call_chets.sh
# remember to change the gene_map file to the correct one

cpp_dir="/software/team281/gk18/call_chets/"
# this is a binary built from https://github.com/frhl/call_chets

if [[ -z $1 ]]; then
    CHR=$LSB_JOBINDEX
else
    CHR=$1
fi

AFthres="05"
TAG="gnh_flagship"

echo "Will identify biallelic carriers for chr-$CHR and $TAG."

work_dir="/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k"
tmp_dir="$work_dir/phasing/recessive_encoding/sandbox"
genotypes="$tmp_dir/chr$CHR.$TAG.txt.gz"

# lists of variants are provided already, eg:
# zcat GNH_LoF_HC_LC_variants.txt.gz | tail -n+2 | sed 's/_/:/g' | awk '($2>0)' > sandbox/gene_map.chrALL.LoF_HC_LC.txt
# zcat GNH_ClinVar_patho_variants.txt.gz | tail -n+2 | sed 's/_/:/g' | awk '($2>0)' > sandbox/gene_map.chrALL.ClinVar.txt
# cat sandbox/gene_map.chrALL.LoF_HC_LC.txt sandbox/gene_map.chrALL.ClinVar.txt | sort -k1 > sandbox/gene_map.chrALL.ClinVar_pLOF.txt

for consq in ClinVar LoF_HC_LC ClinVar_pLOF; do
    $cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chrALL.$consq.txt --show-variants > $tmp_dir/gnh.for_kousik.chr$CHR.$consq.txt
done
  
# print summaries:
for bial in chet hom cis; do
    for consq in ClinVar LoF_HC_LC ClinVar_pLOF; do
        tmp=$(grep $bial $tmp_dir/gnh.for_kousik.chr$CHR.$consq.txt | wc -l | cut -d' ' -f1)
        echo "$bial-$consq events found: $tmp"
    done
done

### how to LSF:
### bsub -J biall_encoding.[1-22] -o biallenc.%I.%J -q normal -R"select[mem>1000] rusage[mem=1000]" -M1000 -n 1 bash run_call_chets.for_kousik.sh