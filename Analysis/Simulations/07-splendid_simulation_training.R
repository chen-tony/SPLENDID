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
rep = as.numeric(params[[5]])
model = params[[6]]

print('Rcpp')
sourceCpp('Package/big_gglassoUtils.cpp')
sourceCpp('Package/big_gglassoLin.cpp')

##########
## DATA ##
##########
print('data')
plink = snp_attach(paste0('PLINK/sim_train', train, '_tun.rds'))

full_summary = readRDS(paste0('SUMMARY/summary_train', train, 
                              '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS')) %>%
  data.frame() %>%
  filter(X_scale > 0)

train_id = fread(paste0('ID/train', train, '.id'), header=F,
                 col.names=c('FID', 'IID'))
tun_id = fread('ID/subtuning.id', header=F,
               col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) %>%
  select(FID, IID, Y=paste0('PHENO', rep)) 

data = plink$fam %>%
  select(FID=family.ID, IID=sample.ID) %>%
  left_join(pheno) %>%
  left_join(pc)

train_data = data %>%
  filter(FID %in% train_id$FID)
tun_data = data %>%
  filter(FID %in% tun_id$FID)

rm(pc, pheno); gc()

# match SNPs
print('match')
map_match = plink$map %>%
  inner_join(full_summary, by=c('marker.ID'='snps'))

ind_train = which((plink$fam$fam %in% train_data$FID) & (plink$fam$samp %in% train_data$IID))
ind_tun = which((plink$fam$fam %in% tun_data$FID) & (plink$fam$samp %in% tun_data$IID))
ind_col = which(plink$map$marker.ID %in% map_match$marker.ID)

all.equal(plink$fam$fam[ind_train], train_data$FID)
all.equal(plink$fam$fam[ind_tun], tun_data$FID)
all.equal(plink$map$marker.ID[ind_col], map_match$marker.ID)

length(ind_train)
length(ind_tun)
length(ind_col)

print('outcome')
Y_train = train_data$Y
Y_mean = mean(Y_train)
Y_sd = sd(Y_train)
Y_train_std = (Y_train - Y_mean) / Y_sd

Y_tun = tun_data$Y
Y_tun_std = (Y_tun - Y_mean) / Y_sd

print('covariates')
Z_train = train_data %>%
  select(contains('PC')) %>% 
  data.matrix()

Z_tun = tun_data %>%
  select(contains('PC')) %>% 
  data.matrix()

Z_mean = colMeans(Z_train)
Z_sd = apply(Z_train, 2, sd)

for (k in 1:ncol(Z_train)) {
  Z_train[,k] = (Z_train[,k] - Z_mean[k]) / Z_sd[k]
  Z_tun[,k] = (Z_tun[,k] - Z_mean[k]) / Z_sd[k]
}

Z_train = cbind(intercept=1, Z_train)
Z_tun = cbind(intercept=1, Z_tun)

E_train = Z_train[,paste0('PC', 1:5)]
E_tun = Z_tun[,paste0('PC', 1:5)]

Z2 = colSums(Z_train^2)

############
### L0L1 ###
############

# Parameters
print('meta-analysis')
grouped_list = list()

meta_analysis = fread(paste0('SUMMARY/interaction_gwas_train', train, 
                             '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'))

print(' meta-analysis match')
meta_match = map_match %>%
  left_join(meta_analysis, by=c('marker.ID'='ID'))
meta_match[is.na(meta_match$P_HET)] = 1

grouped_list[['5e-6']] = as.numeric(meta_match$P_HET) <= 5e-6
grouped_list[['5e-8']] = as.numeric(meta_match$P_HET) <= 5e-8
grouped_list[['0']] = rep(F, length(ind_col))
lapply(grouped_list, sum)

grid = expand.grid(pval = c('0', '5e-8', '5e-6'), lambda_ix = 1:5) %>%
  mutate(file_exists = file.exists(paste0('L0L1/OUTPUT/results_train', train, 
                                          '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', 
                                          model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))) 
grid$n_grouped = sapply(grid$pval, function(x) sum(grouped_list[[as.character(x)]])) 

grid = grid %>%
  filter(!file_exists) %>%
  filter(pval=='0' | n_grouped>0) %>%
  select(-file_exists, -n_grouped)

print(grid)

rm(meta_analysis); gc()

lambda1_vec = exp(seq(log(1e-2), log(1e-4), length.out=5))

if (rho == 1) {
  dfmax = 1000
} else if (rho == 1.5) {
  dfmax = 2500
} else if (rho == 2) {
  dfmax = 7500
} else if (rho == 3) {
  dfmax = 15000
} else if (rho == 4) {
  dfmax = 25000
}

registerDoParallel(cores=15)

out = foreach(i = 1:nrow(grid), .combine='c') %dopar% {
  pval = as.character(grid[i,1])
  lambda_ix = grid[i,2]
  lambda1 = lambda1_vec[lambda_ix]
  
  grouped = grouped_list[[pval]]
  print(sum(grouped))
  
  XE = meta_match$XE_scale * grouped + meta_match$X_scale * (1 - grouped)
  XEY = meta_match$XEY_norm * grouped + meta_match$XY_norm * (1 - grouped)
  lambda0_max = max(XE * pmax(abs(XEY/XE) - 0.5*lambda1/XE, 0))^2
  rm(XE, XEY); gc()
  
  lambda_grid = expand.grid(
    lambda0_G = exp(seq(log(lambda0_max*0.95), log(lambda0_max*5e-4), length.out=50)),
    lambda1_G = lambda1, 
    lambda2_G = 0
  ) %>%
    arrange(desc(lambda1_G), desc(lambda0_G)) %>%
    data.matrix()
  
  test_L0 = COPY_cdfit_gaussian_hsr_par(
    BM = plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
    y = Y_train_std,
    rowInd = ind_train,
    colInd = ind_col,
    Z = Z_train, 
    E = E_train,
    y_val = Y_tun_std,
    rowInd_val = ind_tun,
    Z_val = Z_tun,
    E_val = E_tun,
    lambda = lambda_grid,
    Z2 = Z2,
    X_scale = map_match$X_scale,
    XE_scale = map_match$XE_scale,
    XY_norm = map_match$XY_norm,
    XEY_norm = map_match$XEY_norm,
    grouped = grouped,
    lambda_ix = as.character(lambda_ix),
    pval = pval,
    tol = 1e-3,
    maxit = 100,
    dfmax = dfmax,
    n_abort = 5,
    nlam_min = 10,
    ncheck = 1
  )
  
  # Results
  colnames(test_L0$beta) = outer(c('ADD', paste0('INT', 1:ncol(E_train))), 1:nrow(lambda_grid), FUN=function(k,l) {
    paste0('lambda', l, '_', k)
  }) %>% c()
  
  test_L0$snps = plink$map$marker.ID[ind_col]
  
  results = data.frame(lambda_grid,
                       nb_active = test_L0$nb_active,
                       nb_interact = test_L0$nb_interact,
                       nb_candidate = test_L0$nb_candidate,
                       nb_checks = test_L0$nb_checks,
                       loss = test_L0$loss,
                       iter = test_L0$iter)
  
  message = paste0(pval, '_', lambda_ix, '_done!')
  saveRDS(test_L0, paste0('L0L1/OUTPUT/model_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))
  saveRDS(results, paste0('L0L1/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))
  
  print(message)
  
  i
}

print(paste0('all ', nrow(grid), 'done!'))
