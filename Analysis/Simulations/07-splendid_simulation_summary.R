suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  # library(glmnet)
  
  library(bigstatsr)
  library(bigsnpr)
  library(bigassertr)
  library(bigparallelr)
  library(snpStats)
  library(params)
  
  library(foreach)
  library(dplyr)
  library(data.table)
  library(tidyr)
})

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]]
GA = params[[4]]
model = params[[5]]

print('Rcpp')
sourceCpp('Package/big_gglassoUtils.cpp')

##########
## DATA ##
##########
print('data')

train_id = fread(paste0('ID/train', train, '.id'), header=F,
                 col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

plink = snp_attach(paste0('PLINK/sim_train', train, '.rds'))
snps = read.table('Generate/aou_snp.txt', header=F)$V1

pheno_full = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt'))

registerDoParallel(cores=10)

out = foreach(rep = 1:10, .combine='c') %dopar% {
  pheno = pheno_full %>%
    select(FID, IID, Y=paste0('PHENO', rep)) 
  
  data = plink$fam %>%
    select(FID=family.ID, IID=sample.ID) %>%
    inner_join(train_id) %>%
    left_join(pheno) %>%
    left_join(pc)
  
  ind_row = which(plink$fam$fam %in% data$FID)
  ind_col = which(plink$map$marker.ID %in% snps)
  
  all.equal(plink$fam$fam[ind_row], data$FID)
  print(length(ind_col))
  
  Y_train = data$Y
  
  Z_train = data %>%
    select(contains('PC')) %>% 
    data.matrix()
  
  Z_mean = colMeans(Z_train)
  Z_sd = apply(Z_train, 2, sd)
  
  for (k in 1:ncol(Z_train)) {
    Z_train[,k] = (Z_train[,k] - Z_mean[k]) / Z_sd[k]
  }
  
  Z_train = cbind(intercept=1, Z_train)
  
  E_train = Z_train[,paste0('PC', 1:5)]
  
  print('continuous')
  Y_mean = mean(Y_train)
  Y_sd = sd(Y_train)
  Y_train_std = (Y_train - Y_mean) / Y_sd
  
  res_train = residuals(lm(Y_train_std ~ 0 + Z_train))
  
  rm(data, pheno, pc); gc()
  
  test_sum = bigsummaries(plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
                          rowInd=ind_row, colInd=ind_col,
                          E=E_train, Y=res_train)
  
  test_sum$snps = plink$map$marker.ID[ind_col]
  
  saveRDS(data.frame(test_sum), paste0('SUMMARY/summary_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep,  '_', model, '.RDS'))
  
  cat(rep, 'done!\n')
  rep
}

q('no')

####################
### PC summaries ###
####################
pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))
for (train_id in c('2a', '4a')) {
  train = fread(paste0(home, 'ID/train', train_id, '.id'), header=F,
                col.names=c('FID', 'IID'))
  
  pc_train = pc %>% filter(FID %in% train$FID) %>% 
    select(paste0('PC', 1:5)) %>% 
    data.matrix()
  Esum = list(mean=colMeans(pc_train), sd=apply(pc_train, 2, sd))
  saveRDS(Esum, paste0('SUMMARY/Esum_train', train_id, '.RDS'))
  
}




