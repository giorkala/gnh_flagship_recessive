setwd("~/Documents/WORK_RECESSIVE")
library(RNOmni)

### candidates for the ICD10 codes
geno_prefix="recessive_encoding/genotypes_gnh_flagship"
associations = c( 'MATN3'='N20',
                  'UGT2B15'='N20',
                  'FAM200B'='L02',
                  'SKIV2L'='J15',
                  'NLRP10'='J12',
                  'MAPKBP1'='L90',
                  'HSD17B14'='M85',
                  'NLRP10'='B97',
                  'MED1'='F81',
                  'ZPBP2'='F81',
                  'MET'='K75',
                  'TLDC2'='N48',
                  'GPX2'='N40'
)

# ca didates for BoG'24 analysis
# geno_prefix="recessive_encoding/genotypes_analysis_bog24"
# associations = c("HBB"="D58",
#                  "TMEM106C"="N20", 
#                  "CNTN5"="I48",
#                  "SLC7A6"="M05",
#                  "AMBRA1"="M20",
#                  "SPTBN4"="E10",
#                  "AMBN"="E10",
#                  "UBXN1"="J10",
#                  "VPS51"="J15",
#                  "RIPOR1"="J47",
#                  "SLC27A3"="M53",
#                  "KIF7"="I31",
#                  "C1QTNF8"="O46",
#                  "DPP8"="C79",
#                  "LRP6"="I80",
#                  "ITIH3"="I87",
#                  "GPX2"="N40",
#                  "LSR"="J30"
# ) #

genes_all = names(associations)
traits_all = c()
for (x in genes_all){ traits_all = c(traits_all, associations[[x]])}

# read the files with the ICD10 data and keep only the required codes
df_tmp = read.csv( 'input/phenotypes.ICD10_2of4.tab', sep='\t')
keep_these = c('IID',intersect(traits_all, colnames(df_tmp)))
df_pheno = df_tmp[, keep_these]
for ( x in c(3:4) ){
  df_tmp = read.csv( paste('input/phenotypes.ICD10_',x,'of4.tab',sep=''), sep='\t')
  keep_these = c('FID','IID',intersect(traits_all, colnames(df_tmp)))
  df_pheno = merge( df_pheno, df_tmp[, keep_these], by='IID' )
}

covariates  = read.csv('input/covariates.tab', sep='')
df_pheno = merge( df_pheno, covariates, by='IID' )
rownames(covariates) <- covariates$IID

### end 

samples = read.csv('/genesandhealth/library-red/genesandhealth/GSAv3EAMD/KING_relationships_GSA-51k/v2.3.2/KING_degree3_unrelated_final.fam', sep='\t', header=FALSE)
keep_these <- paste('GNH-',samples$V1, sep='')
# select a set of unrelated samples
df_pheno <- df_pheno[ df_pheno$IID %in% keep_these, ]
rownames(df_pheno) <- df_pheno$IID

results=matrix(0,nrow=length(genes_all), ncol=12)
i=1
# genes_all = c('HBB')
for (gene in genes_all ){
  trait = associations[[gene]]
  print('')
  print(paste('Reading genotypes and initializing tables for',gene,'and',trait))
  genotypes = read.table(paste(geno_prefix,'/genotypes.add_rec_dom.',gene,'.txt.gz', sep=''), header=TRUE, row.names ='X0' )
  # genotypes <- genotypes[ keep_these, ]
  
  df_pruned <- df_pheno[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10',
                        'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19',
                        'PC20', 'age', 'age2', 'sex', 'age_sex', 'age2_sex',trait), ]

  # df_pruned <- df_pheno[, c('PC1', 'PC2', 'PC3', 'PC4', 'age', 'age2', 'sex',
  #                       'age_sex', 'age2_sex',trait), ]
  
  df_pruned <- merge( df_pruned, genotypes[c('ADD','DOM','REC','HetsOnly')], by='row.names', all=TRUE )
  # df_pruned <- merge( df_pruned, genotypes[c('ADD')], by='row.names', all=TRUE )
  df_pruned <- df_pruned[ df_pruned$Row.names %in% keep_these, ]
  
  # drop NAs based on the focal trait
  df_pruned <- df_pruned[ complete.cases(df_pruned[[trait]]), ]
  # make a new column and drop the original trait (to simplify coding)
  df_pruned['y'] = df_pruned[[trait]]
  df_pruned <- df_pruned[, !colnames(df_pruned) %in% c('Row.names',trait)] #Row.names',
  
  print('Done, proceeding with testing...')
  print(sprintf('Effective sample size: %d',nrow(df_pruned)))
  nvars = ncol(df_pruned)-3
  
  # additive [only add]
  model1 = glm( formula = y ~ . -REC -DOM -HetsOnly, data=df_pruned, family = binomial(link = "logit"))
  print( sprintf("Additive: %.4f %.2e", summary(model1)$coefficients[nvars,1], summary(model1)$coefficients[nvars,4] ) )
  results[i,1:3] <- c(summary(model1)$coefficients[27,1], summary(model1)$coefficients[27,4], AIC(model1)) #$r.squared
  
  # dominant 
  model2 = glm( formula = y ~ . -ADD -REC -HetsOnly, data=df_pruned, family = binomial(link = "logit"))
  print( sprintf("Dominance: %.4f %.2e", summary(model2)$coefficients[nvars,1], summary(model2)$coefficients[nvars,4] ) )
  results[i,4:6] <- c(summary(model2)$coefficients[27,1], summary(model2)$coefficients[27,4], AIC(model2))
  
  # recessive [only rec]
  if (sum(df_pruned$REC)>0){
    model3 = glm( formula = y ~ . -ADD -DOM -HetsOnly, data=df_pruned, family = binomial(link = "logit"))
    # summary(model3)
    print( sprintf("Recessive: %.4f %.2e",summary(model3)$coefficients[nvars,1], summary(model3)$coefficients[nvars,4] ) )
    results[i,7:9] <- c(summary(model3)$coefficients[27,1], summary(model3)$coefficients[27,4], AIC(model3))
  } else {
    print("Warning: No biallelic carriers found!")
    results[i,7:9] <- c(0,1,9999)
  }

  # keep add and dom
  model_dd = glm( formula = y ~ . -REC -DOM, data=df_pruned, family = binomial(link = "logit"))
  
  if ( nrow(summary(model_dd)$coefficients) > nvars ){
    print( sprintf("Dom-deviation: %.4f %.2e",summary(model_dd)$coefficients[nvars+1,1], summary(model_dd)$coefficients[nvars+1,4] ) )
    results[i,10:12] <- c(summary(model_dd)$coefficients[28,1], summary(model_dd)$coefficients[28,4],AIC(model_dd))
  } else {
    print("Warning: NAs for HetsOnly")
    results[i,10:12] <- c(0,1,9999)
  }
   
  # # only hets 
  # model4 = glm( formula = y ~ . -ADD -REC -DOM, data=df_pruned, family = binomial(link = "logit"))
  # 
  # # hets vs homz 
  # model5 = glm( formula = y ~ . -REC -DOM -HetsOnly, data=df_pruned[df_pruned$ADD != 0,], family = binomial(link = "logit"))
  # 
  ###

  i <- i+1
}
results <- results[, c(1,4,7,10,3,6,9,12,2,5,8,11)]

