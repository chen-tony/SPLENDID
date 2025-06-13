#!/usr/bin/env Rscript

## PARAMETERS ##
# --env LAMBDA_IX
# --env PVAL
# --env DFMAX
# --input-recursive PACKAGE: C++ directory
# --input UTILS: C++
# --input LIN: C++
# --input PLINK_RDS: bigsnpr .rds file
# --input PLINK_BK: bigsnpr .bk file
# --input MAP: AOU-UKB mapping
# --input PHENO
# --input COV
# --input RELATED
# --input ANCESTRY
# --input PCA
# --input SUMMARY: summary statistics for analysis
# --input META: GWAS meta-analysis
# --output-recursive OUT

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Rcpp)
  library(RcppArmadillo)
  library(bigsnpr)
})

lambda_ix = Sys.getenv('LAMBDA_IX')
pval = Sys.getenv('PVAL')
print(lambda_ix)
print(pval)

if (Sys.getenv('DFMAX') == '') {
  DFMAX == 10000
} else {
  DFMAX = as.numeric(Sys.getenv('DFMAX'))
}
print(DFMAX)

print('Utils'); sourceCpp(Sys.getenv('UTILS'))
print('Lin'); sourceCpp(Sys.getenv('LIN'))

############
### DATA ###
############
print('genotype')
plink = snp_attach(Sys.getenv('PLINK_RDS'))

# print('map')
aou_map = fread(Sys.getenv('MAP'))

print('data')
print(' pheno')
phenotype = fread(Sys.getenv('PHENO'))
names(phenotype)[1:3] = c('FID', 'IID', 'Y')
print(' cov')
covariates = fread(Sys.getenv('COV'))
print(' related')
related = fread(Sys.getenv('RELATED'))
print(' ancestry')
ancestry = fread(Sys.getenv('ANCESTRY')) %>% select(FID, IID, ancestry)
print(' pca')
pc_scores = fread(Sys.getenv('PCA')) %>% select(FID=`#FID`, IID, PC=contains('SCORE'))

data = plink$fam %>%
  select(FID=family.ID, IID=sample.ID) %>%
  left_join(phenotype) %>% 
  left_join(covariates) %>%
  left_join(ancestry) %>%
  left_join(pc_scores) %>%
  filter(!(IID %in% related$sample_id)) %>%
  na.omit()
rm(phenotype, covariates, ancestry, pc_scores, related); gc()

full_summary0 = readRDS(Sys.getenv('SUMMARY'))
full_summary = data.frame(
  marker.ID=full_summary0$snps,
  XY_norm=full_summary0$XY_norm,
  XEY_norm=full_summary0$XEY_norm,
  X_scale=full_summary0$X_scale,
  XE_scale=full_summary0$XE_scale
)

map_match = plink$map %>%
  inner_join(full_summary) %>%
  filter(marker.ID %in% aou_map$aou_code)

aou_ukb_snps = map_match$marker.ID

ind_row = which((plink$fam$fam %in% data$FID) & (plink$fam$samp %in% data$IID))
ind_col = which(plink$map$marker.ID %in% aou_ukb_snps)
ind_train = which((data$FID %in% plink$fam$fam) & (data$IID %in% plink$fam$samp))

cat('ind_row: ', length(ind_row), '\n')
cat('ind_col: ', length(ind_col), '\n')
cat('ind_train: ', length(ind_train), '\n')

cat(nrow(plink$map), 'SNPs in bigsnpr \n')
cat(length(full_summary$snps), 'SNPs in summaries \n')
cat(length(ind_col), 'SNPs kept \n')
cat(sum(plink$map$marker.ID[ind_col] %in% plink$map$marker.ID), 'SNPs overlap with bigsnpr \n')
cat(sum(plink$map$marker.ID[ind_col] %in% full_summary$snps), 'SNPs overlap with summaries \n')

cat('plink and data IDs: ', all.equal(plink$fam$samp[ind_row], data$IID[ind_train]), '\n')
cat('plink and data SNPs: ', all.equal(plink$map$marker.ID[ind_col], map_match$marker.ID), '\n')
# if (!all.equal(plink$fam$samp[ind_row], data$IID[ind_train])) quit('no')

## standardize
print('standardize')
Y = data$Y[ind_train]
Z = data[ind_train,] %>% select(age, female, paste0('PC', 1:20)) %>% data.matrix()
E = Z[,paste0('PC', 1:5)]

print(' Y')
Y_mu = mean(Y)
Y_sd = sd(Y)
Y_std = c(data.matrix(scale(Y)))

print(' Z')
Z_mu = colMeans(Z)
Z_sd = apply(Z, 2, sd)
Z_std = data.matrix(scale(Z))

print(' E')
E_mu = colMeans(E)
E_sd = apply(E, 2, sd)
E_std = data.matrix(scale(E))

# Esum = list(
#   mean = E_mu,
#   sd = E_sd
# )
# saveRDS(Esum, paste0(Sys.getenv('OUT'), '/Esum.RDS'))

## summaries
print('summaries: ')
Z2 = colSums(Z^2)

############
### L0L1 ###
############

# Parameters
print('meta-analysis')
meta_analysis = fread(Sys.getenv('META'))

meta_match = map_match %>%
  left_join(meta_analysis, by=c('marker.ID'='rsid'))
meta_match$P_het[is.na(meta_match$P_het)] = 1
  
if (pval != '0') {
  grouped = meta_match$P_het < as.numeric(pval)
} else {
  grouped = rep(F, length(ind_col))
}
cat(sum(grouped), ' of ', length(grouped), ' grouped \n')

print('lambda')

lambda1 = exp(seq(log(1e-2), log(1e-4), length.out=5))[as.numeric(lambda_ix)]

XE = meta_match$XE_scale * grouped + meta_match$X_scale * (1-grouped)
XEY = meta_match$XEY_norm * grouped + meta_match$XY_norm * (1-grouped)
lambda0_max = max(XE * pmax(abs(XEY/XE) - 0.5*lambda1/XE, 0))^2
print(lambda0_max)
rm(XE, XEY); gc()

lambda_grid = expand.grid(lambda0_G = exp(seq(log(lambda0_max*0.95), log(lambda0_max*5e-4), length.out=50)),
                          lambda1_G = lambda1,
                          lambda2_G = 0) %>%
  data.matrix()

# lambda_grid = expand.grid(
#   lambda0_G = exp(seq(log(1e-3), log(5e-7), length.out=50)),
#   lambda1_G = lambda1,
#   lambda2_G = 0
# ) %>%
#   arrange(desc(lambda1_G), desc(lambda0_G)) %>%
#   data.matrix()
 
# Regression
print('regression')

test_L0 = COPY_cdfit_gaussian_hsr0(
  BM = plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
  y = Y_std,
  rowInd = ind_row,
  colInd = ind_col,
  Z = Z_std, 
  E = E_std,
  lambda = lambda_grid,
  Z2 = Z2,
  X_scale = meta_match$X_scale,
  XE_scale = meta_match$XE_scale,
  XY_norm = meta_match$XY_norm,
  XEY_norm = meta_match$XEY_norm,
  grouped = grouped,
  tol = 1e-3,
  maxit = 100,
  dfmax = DFMAX,
  ncheck = 1
)

# Results
print('results')
colnames(test_L0$beta) = outer(c('ADD', paste0('INT', 1:ncol(E_std))), 1:nrow(lambda_grid), FUN=function(k,l) {
  paste0('lambda', l, '_', k)
}) %>% c()

test_L0$snp = plink$map$marker.ID[ind_col]

results = data.frame(lambda_grid,
                     nb_active = test_L0$nb_active,
                     nb_interact = test_L0$nb_interact,
                     nb_candidate = test_L0$nb_candidate,
                     loss = test_L0$loss,
                     iter = test_L0$iter)

print('save model')
saveRDS(test_L0, paste0(Sys.getenv('OUT'), '/model_pval_', pval, '_lambda', lambda_ix, '.RDS'))

print('save results')
saveRDS(results, paste0(Sys.getenv('OUT'), '/results_pval_', pval, '_lambda', lambda_ix, '.RDS'))

