suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  
  library(bigstatsr)
  library(bigsnpr)
  library(bigassertr)
  library(bigparallelr)
  library(snpStats)
  library(params)
  
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

full_summary = readRDS(paste0('SUMMARY/summary_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS')) %>%
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

rm(pc, pheno, data); gc()

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

E_train = matrix(nrow=0, ncol=0)
E_tun = matrix(nrow=0, ncol=0)

Z2 = colSums(Z_train^2)

############
### L0L1 ###
############
# Parameters
print('meta-analysis')
grouped = rep(F, length(ind_col))
cat(sum(grouped), ' of ', length(grouped), ' grouped \n')

print('lambda')
lambda_max = 2*max(full_summary$XY_norm)

lambda_grid = expand.grid(
  lambda0_G = 0,
  lambda1_G = exp(seq(log(lambda_max), log(lambda_max*1e-3), length.out=200)),
  lambda2_G = 0
) %>%
  arrange(desc(lambda1_G), desc(lambda0_G)) %>%
  data.matrix()

# Regression
print('regression')
if (rho == 1) {
  dfmax = 5000
} else if (rho == 1.5) {
  dfmax = 10000
} else if (rho == 2) {
  dfmax = 15000
} else if (rho == 3) {
  dfmax = 30000
} else if (rho == 4) {
  dfmax = 40000
}

test_L0 = COPY_cdfit_gaussian_hsr(
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
  tol = 1e-3,
  maxit = 100,
  dfmax = dfmax,
  n_abort = 5,
  nlam_min = 30,
  ncheck = 1
)

# Results
print('results')
test_L0$snps = plink$map$marker.ID[ind_col]

results = data.frame(lambda_grid,
                     nb_active = test_L0$nb_active,
                     nb_interact = test_L0$nb_interact,
                     nb_candidate = test_L0$nb_candidate,
                     nb_checks = test_L0$nb_checks,
                     loss = test_L0$loss,
                     iter = test_L0$iter)

print('save model')
saveRDS(test_L0, paste0('Lasso/OUTPUT/model_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
saveRDS(results, paste0('Lasso/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

