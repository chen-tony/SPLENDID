##############
### inputs ###
##############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c('-b', '--bigsnpr'), type='character',
              default=NULL,
              help='bigsnpr file (with .rds extension)',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='slingshot',
              help='output file name',
              metavar='character'),
  make_option('--path', type='character',
              default='~/Slingshot',
              help='path to Slingshot package'),
  make_option('--keep', type='character',
              default=NULL,
              help='text file with FID/IID of training data',
              metavar='character'),
  make_option('--extract', type='character',
              default=NULL,
              help='text file with SNP IDs',
              metavar='character'),
  make_option('--covar', type='character',
              default=NULL,
              help='text file with covariates'),
  make_option('--covar-names', type='character',
              default=NULL,
              help='covariate names, separated by comma'),
  make_option('--interact-names', type='character',
              default=NULL,
              help='covariate names for interactions, separated by comma'),
  make_option('--pheno', type='character',
              default=NULL,
              help='text file with phenotype'),
  make_option('--pheno-name', type='character',
              default=NULL,
              help='phenotype name'),
  make_option('--ncores', type='integer',
              default=NULL,
              help='number of cores for parallel computing'),
  make_option('--maxit', type='integer',
              default=1e6,
              help='maximum iterations for SVD'),
  make_option('--tol', type='numeric',
              default=1e-11,
              help='approx. tolerance for SVD'),
  make_option('--verbose', type='integer',
              default=1,
              help='1 for messages, 0 for no messages'),
  make_option('--imputed', type='integer',
              default=0,
              help='0 to replace NA with 0, 1 if bigsnpr already imputed')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

#######################
### clean up inputs ###
#######################
opt$covar_names = strsplit(opt$covar_names, split=',')[[1]] # turn comma-separated into a character vector
opt$interact_names = strsplit(opt$interact_names, split=',')[[1]] # turn comma-separated into a character vector

#####################
### read in files ###
#####################
suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
  library(dplyr)
  library(Rcpp)
  library(RcppArmadillo)
})

invisible(sourceCpp(paste0(opt$path, 'utils.cpp')))

# read in files
plink = snp_attach(opt$bigsnpr)

cov = fread(opt$covar) %>%
  select(FID, IID, all_of(opt$covar_names))
pheno = fread(opt$pheno) %>%
  select(FID, IID, Y=all_of(opt$pheno_name))

ids = fread(opt$keep, header=F, col.names=c('FID', 'IID'))
snps = fread(opt$extract, header=F)$V1

# get row/col IDs
data = plink$fam %>%
  select(FID=family.ID, IID=sample.ID) %>%
  inner_join(ids) %>%
  left_join(pheno) %>%
  left_join(cov) %>%
  na.omit()

ind_row = which(plink$fam$fam %in% data$FID)
ind_col = which(plink$map$marker.ID %in% snps)

####################
### prepare data ###
####################
# covariates
Z_train = data %>%
  select(all_of(opt$covar_names)) %>%
  data.matrix()

Z_mean = colMeans(Z_train)
Z_sd = apply(Z_train, 2, sd)

for (k in 1:ncol(Z_train)) {
  Z_train[,k] = (Z_train[,k] - Z_mean[k]) / Z_sd[k]
}

Z_train = cbind(intercept=1, Z_train)

# interaction covariates
E_train = data %>%
  select(all_of(opt$interact_names)) %>%
  data.matrix()

E_mean = colMeans(E_train)
E_sd = apply(E_train, 2, sd)

for (k in 1:ncol(E_train)) {
  E_train[,k] = (E_train[,k] - E_mean[k]) / E_sd[k]
}

Esum = list(mean = E_mean, sd = E_sd)

# outcome
Y_train = data$Y
Y_mean = mean(Y_train)
Y_sd = sd(Y_train)
Y_train_std = (Y_train - Y_mean) / Y_sd

res_train = residuals(lm(Y_train_std ~ 0 + Z_train))

rm(data, pheno, cov, ids, snps, k,
   interact_names,
   Y_train, Y_train_std,
   Z_mean, Z_sd, E_mean, E_sd, Y_mean, Y_sd)
invisible(gc())

#####################
### get summaries ###
#####################
if (opt$ncores > 1) {
  if (opt$imputed == 1) {
    big_summary = bigsummaries(plink$genotypes,
                               rowInd=ind_row, colInd=ind_col,
                               E=E_train, Y=res_train, maxit=opt$maxit, tol=opt$tol)
  } else {
    big_summary = bigsummaries(plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
                               rowInd=ind_row, colInd=ind_col,
                               E=E_train, Y=res_train, maxit=opt$maxit, tol=opt$tol)
  }
  
} else {
  if (opt$imputed == 1) {
    big_summary = bigsummaries_par(plink$genotypes,
                                   rowInd=ind_row, colInd=ind_col,
                                   E=E_train, Y=res_train, 
                                   ncores=opt$ncores, maxit=opt$maxit, tol=opt$tol)
  } else {
    big_summary = bigsummaries_par(plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
                                   rowInd=ind_row, colInd=ind_col,
                                   E=E_train, Y=res_train, 
                                   ncores=opt$ncores, maxit=opt$maxit, tol=opt$tol)
  }
}

big_summary$marker.ID = plink$map$marker.ID[ind_col]

###############
### outputs ###
###############
saveRDS(data.frame(big_summary), 
        paste0(opt$out, '_summaries.RDS'))

saveRDS(Esum, 
        paste0(opt$out, '_Esum.RDS'))

cat(paste0('Summaries written to ', opt$out, '_summaries.RDS \n'))
cat(paste0('Interaction summaries written to ', opt$out, '_Esum.RDS \n'))


