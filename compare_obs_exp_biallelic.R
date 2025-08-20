#!/usr/bin/env Rscript

# Developed by F. Lassen for UKB and BRaVa
# Adjusted by G. Kalantzis for G&H 

library(data.table)
library(argparse)

get_maf_mac_from_frqx_file <- function(AC_file){
  dt_AC <- fread(AC_file)
  dt_AC[, MAF := ALT_CTS/OBS_CT]
  dt_AC[, MAC := ALT_CTS]
  dt_AC <- dt_AC[MAC > 0]
  cols <- c("ID","MAF", "MAC")
  dt_AC <- dt_AC[,..cols]
  setkey(dt_AC,"ID")
  return(dt_AC)
}

sim_biallelic_genotypes <- function(gene_id, dt_AC, seed_id, n_samples=39148){
  MAF <- dt_AC$MAF[dt_AC$gene_id == gene_id]
  n_variants_in_gene <- length(MAF)
  # write(paste0(cur_gene_id, " with ", n_variants_in_gene," ", annotation,"s."), stderr())

  H1 <- matrix(rbinom(n_samples * length(MAF), 1, rep(MAF, each = n_samples)), nrow = n_samples)
  H2 <- matrix(rbinom(n_samples * length(MAF), 1, rep(MAF, each = n_samples)), nrow = n_samples)
  
  # homs <- which(apply((X_1 + X_2) == 2, 1, any))
  homs <- which(rowSums((H1 + H2) == 2) > 0)
  chets <- ((rowSums(H1) >= 1) & (rowSums(H2) >= 1))
  chets <- length(setdiff(which(chets), homs))
  homs <- length(homs)
  data.table(gene_id, chets, homs, seed_id)
}

# parameters - check argparser for others
annotation<-c("LOFTEE-HC-LoF","LOFTEE-LC-LoF", "deleterious_missense_CADDgt20_AND_PPgt0.446_AND_SIFT-CLASSeqDeleterious")
annotation_label<-'pLOF_pDM'
annotation<-c("LOFTEE-HC-LoF")
annotation_label<-'pLOF'

main <- function(args){

    print(args)
    
    AC_file <- args$AC_file
    anno_file <- args$anno_file
    geno_file <- args$geno_file
    samples_file <- args$samples_file    
    seeds <- 1:args$N_reps
    max_maf <- args$max_maf 
    # annotation <- args$annotation
    print(paste("Using annotation:", annotation_label))
    
    # read AC file
    dt_AC <- get_maf_mac_from_frqx_file(AC_file)    
    stopifnot(nrow(dt_AC)>0)

    # read annotation file 
    vep <- read.table(anno_file)
    colnames(vep) <- c("ID", "gene_id", "consq")
    vep <- vep[vep$consq %in% annotation,]
    vep <- setDT(vep[!duplicated(vep),])
    setkey(vep, "ID")    
    # vep$symbol <- sub("\\(.*", "", vep$Gene)
    # NOTE: some genes might be missing a symbol and will be skipped!
    stopifnot(nrow(vep)>0)

    # merge tables
    dt_AC <- merge(dt_AC, vep)
    dt_AC <- dt_AC[ dt_AC$MAF <= max_maf]
    stopifnot(nrow(dt_AC)>0)
    valid_variants <- dt_AC$ID
    valid_genes <- unique(dt_AC$gene_id)
    rm(vep)

    cat("Valid genes:", length(valid_genes), '\n')
    cat("Valid variants:", length(valid_variants), '\n')

    is_valid_variant <- function(x){
      all(unlist(strsplit(x, split="\\|")) %in% valid_variants)
    }
    # function to check if a variant is valid
    # split based on '|' or ';'
    # and check if all parts are in valid_variants
    is_valid_variant <- function(x){
      all(unlist(strsplit(x, "[|;]")) %in% valid_variants)
    }

    # read samples file
    samples <- fread(samples_file, header=FALSE)$V1
    n_samples <- length(samples)
    stopifnot(n_samples>0)

    # read genotypes and filter to selected samples
    dt_geno <- read.table(geno_file)
    colnames(dt_geno) <- c("eid", "chrom", "gene_id", "knockout", "ds", "variant")
    dt_geno <- dt_geno[dt_geno$eid %in% samples, ]
    # subset to valid genes and variants that we have in VEP/AC file
    dt_geno <- dt_geno[dt_geno$knockout %in% c("hom", "chet"),]
    cat("Number of biallelic genotypes loaded:", nrow(dt_geno), '\n')
    dt_geno <- dt_geno[dt_geno$gene_id %in% valid_genes,]
    dt_geno$is_valid_variant <- unlist(lapply(dt_geno$variant, is_valid_variant))
    dt_geno <- dt_geno[dt_geno$is_valid_variant == TRUE,]
    cat("Number of genotypes under consideration:", nrow(dt_geno), '\n')

    # count Homs
    hom_counts <- data.table(table(dt_geno$gene_id[dt_geno$knockout %in% "hom"]))
    colnames(hom_counts) <- c("gene_id", "obs_homs")
    setkeyv(hom_counts, "gene_id")

    # count CHs
    ch_counts <- data.table(table(dt_geno$gene_id[dt_geno$knockout %in% "chet"]))
    colnames(ch_counts) <- c("gene_id", "obs_ch")
    setkeyv(ch_counts, "gene_id")

    # merge CH and homs
    obs_counts <- merge(hom_counts, ch_counts, all=TRUE)
    obs_counts[is.na(obs_counts)] <- 0
    obs_counts <- obs_counts[rev(order(obs_counts$obs_homs, obs_counts$obs_ch)),]
    # obs_counts$obs_homs_per_ch <- obs_counts$obs_homs / obs_counts$obs_ch
    
    genes_to_sim <- unique(dt_geno$gene_id[dt_geno$knockout %in% c("hom", "chet")])
    write(paste(length(genes_to_sim), "genes to simulate.."), stderr())
    
    # perform simulation across seeds
    rbinom_tmp <- function(MAF) {
      rbinom(n_samples, 1, MAF)
    }
    cat("Starting simulations...\n")
    combined <- rbindlist(lapply(seeds, function(seed_id){
        # cat(paste0("Seed=",seed_id,'\n')) #, stderr())
        set.seed(seed_id)
        t0 <- Sys.time()
        sim <- rbindlist(lapply(genes_to_sim, function(x) sim_biallelic_genotypes(x, dt_AC, seed_id, n_samples)), use.names = TRUE, fill = TRUE)
        cat(paste0("Time taken for seed ", seed_id, ": ", round(difftime(Sys.time(), t0, units = "secs"), 2), " seconds.\n"))
        return(sim)
    }), use.names = TRUE, fill = TRUE)

    # almost done, now merge OBS with EXP ans save to disk
    results <- aggregate(
      cbind(chets, homs) ~ gene_id,
      data = combined,
      FUN = mean
    )
    
    results <- merge(obs_counts, results, by.x='gene_id', by.y='gene_id')
    setDT(results)
    colnames(results) <- c("gene_id", "o_homs", "o_chet", 'e_chet', 'e_homs')

    results <- melt(results, 
                    id.vars = c("gene_id"),
                    measure.vars = list(
                      c("e_homs","e_chet"),
                      c("o_homs","o_chet")
                    ),
                    variable.name = "gtype",
                    value.name = c("expected", "observed"))
    results[, gtype := factor(gtype, 
                              levels = c(1,2), 
                              labels = c("Homozygotes","Compound Heterozygotes"))]
    # colnames(results)[3:4] <- c("expected", "observed")
    results<-results[, consq:=annotation_label]
    write(paste("Writing results to", args$output), stderr())
    fwrite(results, args$output, sep="\t")
    
    # check if total observed is correct
    if (nrow(dt_geno) != sum(results$observed)) {
        stop("Warning: Total observed genotypes do not match the sum of observed values in results.")
    }
    
    # Group by type and consq and compute stats
    correlation_stats <- results[, .(
      correlation = cor(expected, observed, method = "pearson", use = "complete.obs"),
      r_squared = cor(expected, observed, use = "complete.obs")^2,
      n = length(observed)
    ), by = .(gtype)]
    
    # Add formatted label column
    correlation_stats[, label := sprintf("R² = %.2f", r_squared)]
    print(correlation_stats)
}

parser <- ArgumentParser()
parser$add_argument("--AC_file", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--anno_file", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--geno_file", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--samples_file", default=NULL, required = TRUE, help = "path to samples")
parser$add_argument("--N_reps", type = "integer", default=10, help = "Number of repetitions for simulation")
parser$add_argument("--max_maf", type = "numeric", default=0.05, help = "Maximum MAF to consider")
# parser$add_argument("--annotation", default=NULL, help = "e.g. pLoF")
parser$add_argument("--output", default=NULL, required = TRUE, help = "Path to output file")
args <- parser$parse_args()

main(args)
# end of file