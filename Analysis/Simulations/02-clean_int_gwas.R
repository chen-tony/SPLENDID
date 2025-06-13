# re-format GxPC interaction GWAS from PLINK for SPLENDID

library(dplyr)
library(data.table)
library(tidyr)

args = commandArgs(trailingOnly = T)
train = args[[1]]
rho = as.numeric(args[[2]])
gr = as.numeric(args[[3]])
GA = as.numeric(args[[4]])
model = args[[5]]

for (rep in 1:10) {
  print(rep)
  
  # interaction GWAS
  int = fread(paste0('GWAS/gwas_INT', train, '_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.PHENO', rep, '.glm.linear')) %>%
    filter(TEST %in% c('ADD', 'USER_5DF')) %>%
    select(ID, TEST, P) %>%
    mutate(P = as.numeric(P)) %>%
    spread(key=TEST, value=P) %>%
    select(ID, P_ADD=ADD, P_HET=USER_5DF)
  
  int$P_ADD[is.na(int$P_ADD)] = 1
  int$P_HET[is.na(int$P_HET)] = 1
  
  cat(train, rho, rep, gr, GA, nrow(int), ' ')
  
  fwrite(int, 
         paste0('SUMMARY/interaction_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'),
         sep='\t')
  
  rm(int); gc()
}