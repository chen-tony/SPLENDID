##############
### inputs ###
##############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option('--bigsnpr', type='character',
              default=NULL,
              help='bigsnpr file (with .rds extension)',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='splendid',
              help='output file name',
              metavar='character'),
  make_option('--path', type='character',
              default='~/SPLENDID',
              help='path to SPLENDID package'),
  make_option('--keep', type='character',
              default=NULL,
              help='text file with FID/IID of training data',
              metavar='character'),
  make_option('--keep-tun', type='character',
              default=NULL,
              help='text file with FID/IID of tuning data',
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
  make_option('--summaries', type='character',
              default=NULL,
              help='summary file from splendid_summaries.R'),
  make_option('--het-gwas', type='character',
              default=NULL,
              help='heterogeneity GWAS file (e.g. from splendid_gwas.R)'),
  make_option('--ncores', type='integer',
              default=1,
              help='number of cores for parallel computing (recommended = # lambda1 values x # p-value thresholds)'),
  make_option('--lambda', type='character',
              default=NULL,
              help='lambda1 values'),
  make_option('--pval', type='character',
              default='5e-8',
              help='p-value threshold(s) for heterogeneity'),
  make_option('--nlambda0', type='integer',
              default=50,
              help='number of lambda0 for each lambda1/p-value'),
  make_option('--nlambda0_min', type='integer',
              default=10,
              help='minimum number of lambda0 to test'),
  make_option('--weights', type='file',
              default=NULL,
              help='adaptive penalty weights [SNP name, weight]'),
  make_option('--ncheck', type='integer',
              default=1,
              help='number of active set checks'),
  make_option('--nabort', type='integer',
              default=10,
              help='number of lambda0 after minimum tuning error'),
  make_option('--maxsupp', type='integer',
              default=5000,
              help='maximum support for L0L1 regression'),
  make_option('--verbose', type='integer',
              default=1,
              help='1 for messages, 0 for no messages'),
  make_option('--maxit', type='integer',
              default=100,
              help='maximum iterations for coordinate descent'),
  make_option('--tol', type='numeric',
              default=1e-3,
              help='tolerance for coordinate descent'),
  make_option('--imputed', type='integer',
              default=0,
              help='0 to replace NA with 0, 1 if bigsnpr already imputed')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

#####################
### read in files ###
#####################
suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
  library(dplyr)
  library(Rcpp)
  library(RcppArmadillo)
  library(foreach)
  library(bigparallelr)
})

invisible(sourceCpp(paste0(opt$path, 'utils.cpp')))
invisible(sourceCpp(paste0(opt$path, 'linear.cpp')))

opt$covar_names = strsplit(opt$covar_names, split=',')[[1]] # turn comma-separated into a character vector
opt$interact_names = strsplit(opt$interact_names, split=',')[[1]] # turn comma-separated into a character vector

# read in files
plink = snp_attach(opt$bigsnpr)

cov = fread(opt$covar) %>%
  select(FID, IID, all_of(opt$covar_names))
pheno = fread(opt$pheno) %>%
  select(FID, IID, Y=all_of(opt$pheno_name))

full_summary = readRDS(opt$summaries) %>%
  data.frame() %>%
  filter(X_scale > 0)
meta_analysis = fread(opt$het_gwas)

train_ids = fread(opt$keep, header=F, col.names=c('FID', 'IID'))

if (!is.null(opt$keep_tun)) {
  tun_ids = fread(opt$keep_tun, header=F, col.names=c('FID', 'IID'))
  ids = rbind(train_ids, tun_ids) %>% arrange(FID, IID)
} else {
  ids = train_ids %>% arrange(FID, IID)
}
snps = fread(opt$extract, header=F)$V1

# get row/col IDs
data = plink$fam %>%
  select(FID=family.ID, IID=sample.ID) %>%
  inner_join(ids) %>%
  left_join(pheno) %>%
  left_join(cov) %>%
  na.omit()

train_data = data %>%
  filter(FID %in% train_ids$FID)
ind_row = which(plink$fam$fam %in% train_data$FID)

if (!is.null(opt$keep_tun)) {
  tun_data = data %>%
    filter(FID %in% tun_ids$FID)
  ind_row_tun = which(plink$fam$fam %in% tun_data$FID)
}

meta_match = plink$map %>%
  filter(marker.ID %in% snps) %>%
  inner_join(full_summary, by='marker.ID') %>%
  left_join(meta_analysis, by=c('marker.ID'))
meta_match[is.na(meta_match$P_HET)] = 1

# penalty weights
if (!is.null(opt$weights)) {
  penalty_weights = fread(opt$weights)
  names(penalty_weights) = c('marker.ID', 'weights')
  
  meta_match = meta_match %>%
    left_join(penalty_weights)
  rm(penalty_weights)
  meta_match$weights[is.na(meta_match$weights)] = 1
  
} else {
  meta_match$weights = 1
}

ind_col = which(plink$map$marker.ID %in% meta_match$marker.ID)

rm(full_summary, meta_analysis)
invisible(gc())

####################
### prepare data ###
####################
# training data
Y_train = train_data$Y
Y_mean = mean(Y_train)
Y_sd = sd(Y_train)
Y_train_std = (Y_train - Y_mean) / Y_sd

Z_train = train_data %>%
  select(all_of(opt$covar_names)) %>% 
  data.matrix()
Z_mean = colMeans(Z_train)
Z_sd = apply(Z_train, 2, sd)

for (k in 1:ncol(Z_train)) {
  Z_train[,k] = (Z_train[,k] - Z_mean[k]) / Z_sd[k]
}
Z_train = cbind(intercept=1, Z_train)

E_train = train_data %>%
  select(all_of(opt$interact_names)) %>%
  data.matrix()
E_mean = colMeans(E_train)
E_sd = apply(E_train, 2, sd)

for (k in 1:ncol(E_train)) {
  E_train[,k] = (E_train[,k] - E_mean[k]) / E_sd[k]
}

Z2 = colSums(Z_train^2)

# tuning data
if (!is.null(opt$keep_tun)) {
  Y_tun = tun_data$Y
  Y_tun_std = (Y_tun - Y_mean) / Y_sd
  
  Z_tun = tun_data %>%
    select(all_of(opt$covar_names)) %>% 
    data.matrix()
  
  E_tun = tun_data %>%
    select(all_of(opt$interact_names)) %>%
    data.matrix()
  
  for (k in 1:ncol(Z_tun)) {
    Z_tun[,k] = (Z_tun[,k] - Z_mean[k]) / Z_sd[k]
  }
  Z_tun = cbind(intercept=1, Z_tun)
  
  for (k in 1:ncol(E_tun)) {
    E_tun[,k] = (E_tun[,k] - E_mean[k]) / E_sd[k]
  }
}

rm(data, cov, pheno, k,
   snps, ids, train_ids, tun_ids,
   Y_train, Y_tun, Y_mean, Y_sd, 
   E_mean, E_sd, Z_mean, Z_sd,
   train_data, tun_data)
invisible(gc())

##################
### regression ###
##################
# parameters
if (is.null(opt$lambda)) {
  lambda1_vec = exp(seq(log(1e-2), log(1e-4), length.out=5))
} else {
  lambda1_vec = strsplit(opt$lambda, split=',')[[1]]
}
pval = strsplit(opt$pval, split=',')[[1]]
names(pval) = pval
grid = expand.grid(pval = pval, lambda_ix = seq_along(lambda1_vec))

# grouped penalty
grouped_list = lapply(pval, FUN=function(p) as.numeric(meta_match$P_HET) <= as.numeric(p))
grid$n_grouped = sapply(grid$pval, function(x) sum(grouped_list[[as.character(x)]])) 

# regression
registerDoParallel(cores=opt$ncores)

out = foreach(i = 1:nrow(grid), .combine='c') %dopar% {
  pval = as.character(grid[i,1])
  lambda_ix = grid[i,2]
  lambda1 = as.numeric(lambda1_vec[lambda_ix])
  
  grouped = grouped_list[[pval]]
  
  XE = meta_match$XE_scale * grouped + meta_match$X_scale * (1 - grouped)
  XEY = meta_match$XEY_norm * grouped + meta_match$XY_norm * (1 - grouped)
  lambda0_max = max(XE * pmax(abs(XEY/XE) - 0.5*lambda1/XE, 0))^2
  rm(XE, XEY); invisible(gc())
  
  if (lambda0_max > 0) {
    lambda_grid = expand.grid(
      lambda0_G = exp(seq(log(lambda0_max*0.95), log(lambda0_max*5e-4), length.out=opt$nlambda0)),
      lambda1_G = lambda1, 
      lambda2_G = 0
    ) %>%
      arrange(desc(lambda1_G), desc(lambda0_G)) %>%
      data.matrix()
    
    if (!is.null(opt$keep_tun)) {
      L0_reg = cdfit_gaussian_par(
        BM = plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
        y = Y_train_std,
        rowInd = ind_row,
        colInd = ind_col,
        Z = Z_train, 
        E = E_train,
        y_tun = Y_tun_std,
        rowInd_val = ind_row_tun,
        Z_tun = Z_tun,
        E_tun = E_tun,
        weights = meta_match$weights,
        lambda = lambda_grid,
        Z2 = Z2,
        X_scale = meta_match$X_scale,
        XE_scale = meta_match$XE_scale,
        XY_norm = meta_match$XY_norm,
        XEY_norm = meta_match$XEY_norm,
        grouped = grouped,
        lambda_ix = as.character(lambda_ix),
        pval = pval,
        tol = opt$tol,
        maxit = opt$maxit,
        dfmax = opt$maxsupp,
        n_abort = opt$nabort,
        nlam_min = opt$nlambda0_min,
        ncheck = opt$ncheck
      )
    } else {
      L0_reg = cdfit_gaussian_par0(
        BM = plink$genotypes$copy(code = c(0, 1, 2, rep(0, 253))),
        y = Y_train_std,
        rowInd = ind_row,
        colInd = ind_col,
        Z = Z_train, 
        E = E_train,
        weights = meta_match$weights,
        lambda = lambda_grid,
        Z2 = Z2,
        X_scale = meta_match$X_scale,
        XE_scale = meta_match$XE_scale,
        XY_norm = meta_match$XY_norm,
        XEY_norm = meta_match$XEY_norm,
        grouped = grouped,
        lambda_ix = as.character(lambda_ix),
        pval = pval,
        tol = opt$tol,
        maxit = opt$maxit,
        dfmax = opt$maxsupp,
        ncheck = opt$ncheck
      )
    }
    
    # Results
    colnames(L0_reg$beta) = outer(c('ADD', paste0('INT', 1:ncol(E_train))), 
                                  1:nrow(lambda_grid), FUN=function(k,l) {
                                    paste0('lambda', l, '_', k)
                                  }) %>% c()
    
    L0_reg$snps = plink$map$marker.ID[ind_col]
    
    results = data.frame(lambda_grid,
                         pval = pval,
                         nb_active = L0_reg$nb_active,
                         nb_interact = L0_reg$nb_interact,
                         nb_candidate = L0_reg$nb_candidate,
                         nb_checks = L0_reg$nb_checks,
                         loss = L0_reg$loss,
                         iter = L0_reg$iter)
    
    
    saveRDS(L0_reg, paste0(opt$out, '_pval_', pval, '_lambda_', lambda_ix, '_beta.RDS'))
    saveRDS(results, paste0(opt$out, '_pval_', pval, '_lambda_', lambda_ix, '_results.RDS'))
    
    if (opt$verbose > 0) {
      message = paste0(pval, '_', lambda_ix, '_done!')
      message
    }
  }
  
  i
}

cat('all', nrow(grid), 'done! \n')

###############
### outputs ###
###############
cat(paste0('Betas written to ', opt$out, '*_beta.RDS \n'))
cat(paste0('Grids written to ', opt$out, '*_results.RDS \n'))

