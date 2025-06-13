#!/usr/bin/env Rscript

## PARAMETERS ##
# --input-recursive PACKAGE: C++ directory
# --input UTILS: C++
# --input PLINK_RDS: bigsnpr .rds file
# --input PLINK_BK: bigsnpr .bk file
# --input MAP: AOU-UKB mapping
# --input PHENO
# --input COV
# --input RELATED
# --input ANCESTRY
# --input PCA
# --output OUT: summary data
# --output Esum: covariate summary data

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Rcpp)
  library(RcppArmadillo)
  library(bigsnpr)
})

print('Utils'); sourceCpp(Sys.getenv('UTILS'))

############
### DATA ###
############
print('genotype')
plink = snp_attach(Sys.getenv('PLINK_RDS'))
cat('fam: ', dim(plink$fam), '\n')
cat('geno: ', dim(plink$genotypes), '\n')
cat('map: ', dim(plink$map), '\n')

print('map')
map = fread(Sys.getenv('MAP'))
aou_ukb_snps = map %>% pull(aou_code)
cat(length(aou_ukb_snps), 'in MAP \n')

print('data')
print(' pheno')
phenotype = fread(Sys.getenv('PHENO'))
colnames(phenotype)[1:3] = c('FID', 'IID', 'Y')
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
  left_join(phenotype %>% select(FID, IID, Y)) %>% 
  left_join(covariates) %>%
  left_join(ancestry) %>%
  left_join(pc_scores) %>%
  filter(!(IID %in% related$sample_id)) %>%
  na.omit()
cat('data: ', dim(data), '\n')
head(data)

rm(phenotype, covariates, related, ancestry, pc_scores, map); gc()

# subsets for bigsnpr data
ind_row = which((plink$fam$fam %in% data$FID) & (plink$fam$samp %in% data$IID))
ind_col = which(plink$map$marker.ID %in% aou_ukb_snps)
cat('ind_row: ', length(ind_row), '\n')
cat('ind_col: ', length(ind_col), '\n')

# subsets for AOU data
ind_train = which((data$FID %in% plink$fam$fam) & (data$IID %in% plink$fam$samp))
cat('ind_train: ', length(ind_train), '\n')

cat('plink and data IDs: ', all.equal(plink$fam$samp[ind_row], data$IID[ind_train]), '\n')

## standardize
print('standardize')

print(' Y')
Y = data$Y[ind_train]
# Y_mu = mean(Y)
# Y_sd = sd(Y)
Y_std = c(data.matrix(scale(Y)))

print(' Z')
Z = data[ind_train,] %>% select(age, female, paste0('PC', 1:20)) %>% data.matrix()
# Z_mu = colMeans(Z)
# Z_sd = apply(Z, 2, sd)
Z_std = data.matrix(scale(Z))

print(' E')
E = Z[,paste0('PC', 1:5)]
E_mu = colMeans(E)
E_sd = apply(E, 2, sd)
E_std = data.matrix(scale(E))

saveRDS(list(mean=E_mu, sd=E_sd), Sys.getenv('Esum'))

## summaries
print('residuals')
res = residuals(lm(Y_std ~ 0 + Z_std))
rm(Y, Z, E, Y_std, data); gc()

print('start summaries')
test_summary = bigsummaries(BM=plink$genotypes$copy(code=c(0,1,2,rep(0,253))),
                            rowInd=ind_row, colInd=ind_col,
                            E=E_std,
                            Y=res)
gc()

test_summary$snps = plink$map$marker.ID[ind_col]

test_summary$E_mu = E_mu
test_summary$E_sd = E_sd

saveRDS(test_summary, Sys.getenv('OUT'))

lapply(test_summary, summary)