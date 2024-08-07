#!/bin/python3

# Make tables with ADD/DOM/REC genotypes from the VCFs - for G&H
# The script reads the VCFs and creates a table with the genotypes for each gene

# First, run the following accordingly to get the list of all available genes:
# bcftools view -HG chrALL.recessive_PP90af05.s50.$consq.vcf.gz | cut -f1-3 > gene_list.$consq.txt

import pandas as pd
import numpy as np

work_dir='some/dir'
out_prefix=work_dir+'/domdev_genotypes'
samples_file=work_dir+'/samples.vcf.txt'

# genes for G&H flagship analysis
# remember to exclude multiple occurrences here
# ie, the association ('CLYBL-pLOF', 'Vit_B12') is surpassed by the pLOF_deleterious homologue
genes_tmp = [
    ['PIEZO1', 'pLOF_deleterious'],
    ['UGT1A6', 'pLOF_deleterious'],
    ['UGT1A9', 'pLOF_deleterious'],
    ['DPH1', 'pLOF_deleterious'],
    ['NDUFB3', 'pLOF_deleterious'],
    ['PIEZO1', 'pLOF_deleterious'],
    ['IL7', 'pLOF_deleterious'],
    ['MIA2', 'pLOF_deleterious'],
    ['PSMB11', 'pLOF_deleterious'],
    ['APOE', 'pLOF_deleterious'],
    ['NCAPD2', 'pLOF_deleterious'],
    ['FUT8', 'pLOF_deleterious'],
    ['CLYBL', 'pLOF_deleterious'],
    ['MATN3', 'pLOF_deleterious'],
    ['UGT2B15', 'pLOF_deleterious'],
    ['ALG10B', 'pLOF_deleterious'],
    ['FAM208B', 'pLOF_deleterious'],
    ['SKIV2L', 'pLOF_deleterious'],
    ['NLRP10', 'pLOF_deleterious'],
    ['MAPKBP1', 'pLOF_deleterious'],
    ['HSD17B14', 'pLOF_deleterious'],
    ['MED1', 'pLOF_deleterious'],
    ['ZBP2', 'pLOF_deleterious'],
    ['MET', 'pLOF_deleterious'],
    ['POLD3', 'pLOF_deleterious'],
    ['TLDC2', 'pLOF_deleterious'],
    ['GPX2', 'pLOF_deleterious'],
    # ['UGT1A6', 'synonymous'],
    # ['RSPH1', 'synonymous'],
    # ['UGT1A3', 'synonymous'],
    # ['ALG10', 'pLOF_deleterious'],
    # ['UGT1A9', 'synonymous'],
    # ['PIEZO1', 'synonymous'],
    # ['NDUFB3', 'synonymous'],
    # ['PIEZO1', 'synonymous'],
    # ['DPH1', 'synonymous'],
    # ['IL7', 'synonymous'],
    # ['MIA2', 'synonymous'],
]

assert len(set([x[0] for x in genes_tmp])) == len(genes_tmp), "wrong input"

genes_of_interest = { 'pLOF_deleterious':[], 'synonymous':[] }
for x in genes_tmp:
    if x[1]=='pLOF_deleterious':
        genes_of_interest['pLOF_deleterious'].append(x[0])
    else:
        genes_of_interest['synonymous'].append(x[0])

# start reading genotypes from disc
df_add = []
df_rec = []
for consq in ['pLOF_deleterious','synonymous']:
    vcf0=f'{work_dir}/chrALL.recessive_gnh_flagship.{consq}.vcf.gz'
    vcf1=f'{work_dir}/chrALL.additive_gnh_flagship.{consq}.vcf.gz'
    
    gene_list = pd.read_csv(f'{work_dir}/gene_list.gnh_flagship.{consq}.txt', header=None, 
                            sep='\s+', names=['CHR','POS','GENE','other'], index_col='GENE')
    
    gene_list['ID'] = np.arange(len(gene_list))
    for gene in genes_of_interest[consq]:
        print(gene, gene_list.loc[gene].CHR)
        lines_to_skip = gene_list.loc[gene].ID + 34
        df_rec.append( pd.read_csv(vcf0, skiprows=lines_to_skip, sep='\t', nrows=1, header=None) )
        df_add.append( pd.read_csv(vcf1, skiprows=lines_to_skip, sep='\t', nrows=1, header=None) )

df_rec = pd.concat(df_rec).set_index(2)
df_add = pd.concat(df_add).set_index(2)

# create new files with REC/DOM/ADD/HetsOnly encodings
for gene in [x[0] for x in genes_tmp ]:
    df_geno = pd.read_csv(samples_file, header=None)

    g_rec = df_rec.loc[gene,9:].to_numpy().astype(int)
    g_add = df_add.loc[gene,9:].to_numpy().astype(int)
    g_hetsOnly = g_add - g_rec

    df_geno['REC'] = g_rec
    df_geno['ADD'] = g_add
    df_geno['HetsOnly'] = g_hetsOnly

    df_geno.set_index(0, inplace=True)
    hets = df_geno.query('ADD>=1').index
    
    df_geno['DOM'] = 0
    df_geno.loc[hets, 'DOM'] = 1
    print(gene)
    print(df_geno.sum(axis=0) )

    df_geno.to_csv(f'{out_prefix}.{gene}.txt.gz', sep='\t', compression='gzip')

# end of script