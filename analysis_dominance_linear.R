setwd("~/Documents/WORK_RECESSIVE")
library(RNOmni)

### associations for BoG'24 analysis
# geno_prefix="recessive_encoding/genotypes_analysis_bog24"
# associations <- c( 'HBB'='HDL_C', 
#                    'CATSPER2'='Immature_granulocytes',
#                    'CLYBL'='Vitamin_B12',
#                    'UGT1A6'='Bilirubin',
#                    'NBPF3'='ALP'
# )

geno_prefix="/genesandhealth/red/GeorgiosKalantzis/WORK_RECESSIVE_new/recessive_encoding/genotypes_gnh_flagship"
# associations <- c('NBPF3'='ALP',
#                   'UGT1A3'='Bilirubin',
#                   'UGT1A6'='Bilirubin',
#                   'ALG10'='Eosinophils',
#                   'RSPH1'='Eosinophils',
#                   'PIEZO1'='Hb',
#                   'PIEZO1'='HbA1c',
#                   'APOE'='LDL.C',
#                   'APOE'='Non_HDL_cholesterol',
#                   'DIO1'='T4',
#                   'CLYBL'='Vitamin_B12',
#                   'FUT6'='Vitamin_B12'
# )

# significant for G&H flagship
genes_all = c( "NBPF3",  "UGT1A3", "UGT1A6", "ALG10",  "RSPH1"  ,"PIEZO1", "PIEZO1", "APOE",   "APOE","DIO1",   "CLYBL",  "FUT6"  )
traits_all = c("ALP", "Bilirubin","Bilirubin","Eosinophils", "Eosinophils", "Hb", "HbA1c","LDL.C", "Non_HDL_cholesterol", "T4","Vitamin_B12", "Vitamin_B12")
# traits_all = c(); for (x in genes_all){ traits_all = c(traits_all, associations[[x]])}

# p_rec < 5e-6 and p_rec>5e-8 for G&H flagship
# genes_all = c('PIEZO1', 'UGT1A6', 'UGT1A9', 'DPH1', 'NUDT13', 'IL7', 'MIA2', 'PIEZO1','PSMB11', 'APOE', 'NCAPD2', 'MIA2')
# traits_all = c('Bilirubin', 'Bilirubin', 'Bilirubin', 'CRP', 'Hb', 'Lymphocytes',
               # 'Neutrophils', 'RDW', 'Systolic_BP', 'Total_cholesterol', 'Vitamin_B12', 'WBC')

covariates  = read.csv('input/covariates.tab', sep='')
rownames(covariates) <- covariates$IID

pheno_quant = read.csv('input/phenotypes.55cur_qt.tab', sep='\t')
rownames(pheno_quant) <- pheno_quant$IID
# pheno_quant <- pheno_quant[,c('HDL_C','MCH','MCV','Ferritin','Bilirubin','Hb','Haemoglobin_F','Haemotocrit','RDW','ALP','Vitamin_B12','Immature_granulocytes')]
keep_these = c('IID',intersect(traits_all, colnames(pheno_quant)))
pheno_quant <- pheno_quant[,keep_these]

pheno_quant <- merge( pheno_quant, covariates, by='row.names', all=TRUE )
rownames(pheno_quant) <- pheno_quant$Row.names

# select a set of unrelated samples
samples = read.csv('/genesandhealth/library-red/genesandhealth/GSAv3EAMD/KING_relationships_GSA-51k/v2.3.2/KING_degree3_unrelated_final.fam', sep='\t', header=FALSE)
keep_these <- paste('GNH-',samples$V1, sep='')
# work with the more stringent set of unrelated samples:
# samples = read.csv('/genesandhealth/green/forGeorgios/WES-44k_ban_pak_unrelated_from_GSA-51k_stringent.fam', sep='\t', header=FALSE)
# keep_these <- samples$V1

results=matrix(0,nrow=length(traits_all), ncol=13)
for (i in c(1:length(genes_all)) ){
# for (i in c(8:9) ){
  gene = genes_all[i]
  trait = traits_all[i]
  print(paste('Reading genotypes and initializing tables for',gene,'and',trait))
  
  genotypes = read.table(paste(geno_prefix,'/genotypes.add_rec_dom.',gene,'.txt.gz', sep=''), header=TRUE, row.names ='X0' )
  # genotypes <- genotypes[ keep_these, ]
  
  df_pruned <- pheno_quant[,c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10',
                              'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19',
                              'PC20', 'age', 'age2', 'sex', 'age_sex', 'age2_sex',trait)]
  
  df_pruned <- merge( df_pruned, genotypes[c('ADD','DOM','REC','HetsOnly')], by='row.names', all=TRUE )
  df_pruned <- df_pruned[ df_pruned$Row.names %in% keep_these, ]

  # drop NAs based on the focal trait and REC
  df_pruned <- df_pruned[ complete.cases(df_pruned[[trait]]), ]
  df_pruned <- df_pruned[ complete.cases(df_pruned$REC), ]
  
  # perform rank-inverse transformation
  df_pruned['y_rint'] = RankNorm(df_pruned[[trait]])    
  # drop the original trait (to simplify testing)
  df_pruned <- df_pruned[, !colnames(df_pruned) %in% c('Row.names',trait)]
  #y=df_pruned[[trait]]

  print('Done, proceeding with testing...')
  print(sprintf('Effective sample size: %d',nrow(df_pruned)))
  nvars = ncol(df_pruned)-3
  
  # additive [only add]
  model1 = lm( formula = y_rint ~ . -REC -DOM -HetsOnly, data=df_pruned)
  print( sprintf("Additive: %.4f %.2e", summary(model1)$coefficients[nvars,1], summary(model1)$coefficients[nvars,4] ) )
  results[i,1:3] <- c(summary(model1)$coefficients[nvars,1], summary(model1)$coefficients[nvars,4], summary(model1)$r.squared)
  
  # dominant 
  model2 = lm( formula = y_rint ~ . -ADD -REC -HetsOnly, data=df_pruned, )
  print( sprintf("Dominance: %.4f %.2e", summary(model2)$coefficients[nvars,1], summary(model2)$coefficients[nvars,4] ) )
  results[i,4:6] <- c(summary(model2)$coefficients[nvars,1], summary(model2)$coefficients[nvars,4], summary(model2)$r.squared)
  
  # recessive [only rec]
  if (sum(df_pruned$REC)>0){
    model3 = lm( formula = y_rint ~ . -ADD -DOM -HetsOnly, data=df_pruned )
    # summary(model3)
    results[i,7:9] <- c(summary(model3)$coefficients[nvars,1], summary(model3)$coefficients[nvars,4], summary(model3)$r.squared)
    print( sprintf("Recessive: %.4f %.2e",summary(model3)$coefficients[nvars,1], summary(model3)$coefficients[nvars,4] ) )
  } else {
    print("Warning: No biallelic carriers found!")
    results[i,7:9] <- c(0,1,9999)
  }
  
  # keep add and dom
  model_dd = lm( formula = y_rint ~ . -REC -DOM, data=df_pruned )
  if ( nrow(summary(model_dd)$coefficients) > nvars ){
    print( sprintf("Dom-deviation: %.4f %.2e",summary(model_dd)$coefficients[nvars+1,1], summary(model_dd)$coefficients[nvars+1,4] ) )
    results[i,10:12] <- c(summary(model_dd)$coefficients[nvars+1,1], summary(model_dd)$coefficients[nvars+1,4], summary(model_dd)$r.squared )
  } else {
    print("Warning: NAs for HetsOnly")
    results[i,10:12] <- c(0,1,9999)
  }
  # # only hets 
  # model4 = lm( formula = y_rint ~ . -ADD -REC -DOM, data=df_pruned)
  # 
  # # hets vs homz 
  # model5 = lm( formula = y_rint ~ . -REC -DOM -HetsOnly, data=df_pruned[df_pruned$ADD != 0,])
  # #summary(model2) RankNorm(y)
  
  ###
  #print(c('AIC:',extractAIC(model1)[2],extractAIC(model2)[2],extractAIC(model3)[2],extractAIC(model_dd)[2]))
  print(c('AIC add-dom vs add:',extractAIC(model_dd)[2] - extractAIC(model1)[2]))
  
  results[i,13] = sum(df_pruned$REC)
}
results <- results[, c(1,4,7,10,3,6,9,12,2,5,8,11,13)]
