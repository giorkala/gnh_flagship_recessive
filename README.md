## Rare Recessive Burden Analysis
Scripts used for the recessive association analysis, within the **G&amp;H Flagship paper**.  Prior to this, we phased our WES dataset using a Snakemake pipeline developed [here](https://github.com/BRaVa-genetics/snakemake_pipeline_for_phasing).

This repository contains:
1. `get_biallelic_carriers.sh` to identify biallelic carriers, given phased genotypes. Also, see [here](https://github.com/giorkala/GNH_BRaVa_dev_misc/blob/main/recessive_encoding/run_call_chets.sh) for a more involved script that starts from the phased genotypes and creates VCF files ready for association.
2. `get_ch_homz_counts.sh` generates a tally with the numbers of compound-het and homozygous genotypes.
3. `concat_vcfs.sh` concatenates chromosome-based VCFs to one, for a given consequence, and creates a BGEN for more efficient association testing.
4. `run_regenie_{burden/sv/perm}.sh` are scripts to invoke REGENIE for gene-burden testing, variant-level, or permutation (ie FDR analysis), respectively. Each script has several modules, e.g. step1 (when needed) or step2, for quant or binary traits, which are selected based on the arguments, e.g. `bash run_regenie.sh step2qt pLOF` will run step2 for quant traits and test for the pLOF burdens.
5. `dominance_deviation_{linear/logistic}.R`: R scripts to test for dominance deviation, for a candidate gene-trait pair. This requires a table with ADD/DOM/REC encodes, per gene, which can be obtained with `dominance_deviation_prep.py`.

Georgios Kalantzis, gk18@sanger.ac.uk  
August 2024
