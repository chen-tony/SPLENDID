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
  make_option('--bfile', type='character',
              default=NULL,
              help='plink file (without .bed/.bim/.fam extension)',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='slingshot',
              help='output file name',
              metavar='character'),
  make_option('--path', type='character',
              default='~/Slingshot',
              help='path to Slingshot package'),
  make_option('--plink2', type='character',
              default='~/plink2',
              help='path to plink2 software'),
  make_option('--keep-tun', type='character',
              default=NULL,
              help='text file with FID/IID of tuning data',
              metavar='character'),
  make_option('--keep-val', type='character',
              default=NULL,
              help='text file with FID/IID of validation data',
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
              help='summary file from slingshot_summaries.R'),
  make_option('--Esum', type='character',
              default=NULL,
              help='interaction summary file from slingshot_summaries.R'),
  make_option('--het-gwas', type='character',
              default=NULL,
              help='heterogeneity GWAS file (e.g. from slingshot_gwas.R)'),
  make_option('--ncores', type='integer',
              default=NULL,
              help='number of cores for parallel computing (recommended = # lambda1 values x # p-value thresholds)'),
  make_option('--maxsupp', type='integer',
              default=NULL,
              help='maximum support for PRS'),
  make_option('--maxfactor', type='integer',
              default=NULL,
              help='maximum factor larger than best grid PRS'),
  make_option('--cleansup', type='integer',
              default=1,
              help='1 for cleanup of intermediate files, 0 for no actions'),
  make_option('--verbose', type='integer',
              default=1,
              help='1 for messages, 0 for no messages')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

opt = list(
  bigsnpr = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/chr21_ALL_sub.rds',
  bfile = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/chr21_ALL',
  out = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/chr21_ALL_sub',
  path = '/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/Package/',
  plink2 = '~/software/plink2',
  keep = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/keep.txt',
  keep_tun = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/keep_tun.txt',
  extract = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/extract.txt', 
  covar = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/covar.txt', 
  covar_names = 'PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10', 
  pheno = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/pheno.txt',
  pheno_name = 'PHENO1',
  interact_names = 'PC1,PC2,PC3,PC4,PC5',
  Esum = '/n/holyscratch01/xlin/tonychen/MultiAncestry/TEST/chr21_ALL_sub_Esum.RDS', 
  maxsupp = NULL,
  maxfactor=NULL,
  cleanup = 1,
  verbose=1
)


#####################
### read in files ###
#####################
suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
  library(dplyr)
  library(glmnet)
  library(Matrix)
})

opt$covar_names = strsplit(opt$covar_names, split=',')[[1]] # turn comma-separated into a character vector
opt$interact_names = strsplit(opt$interact_names, split=',')[[1]] # turn comma-separated into a character vector

K = length(opt$interact_names) + 1

plink = snp_attach(opt$bigsnpr)

tun_ids = fread(opt$keep_tun, header=F, col.names=c('FID', 'IID'))
if (!is.null(opt$keep_val)) val_ids = fread(opt$keep_val, header=F, col.names=c('FID', 'IID'))

cov = fread(opt$covar) %>%
  select(FID, IID, all_of(opt$covar_names))
pheno = fread(opt$pheno) %>%
  select(FID, IID, Y=all_of(opt$pheno_name))

data = pheno %>%
  left_join(cov) %>%
  filter(FID %in% tun_ids$FID & IID %in% tun_ids$IID)
data$res = residuals(lm(Y~., data %>% select(-FID, -IID)))

#############
### betas ###
#############
beta_files = Sys.glob(paste0(opt$out, '*beta.RDS'))
grid_files = Sys.glob(paste0(opt$out, '*results.RDS'))

beta = NULL
for (beta_file in beta_files) {
  model0 = readRDS(beta_file)
  beta0 = model0$beta
  rownames(beta0) = model0$snps
  
  beta = cbind(beta, beta0)
  rm(beta0); gc()
}

grid = NULL
for (grid_file in grid_files) {
  grid0 = readRDS(grid_file)
  
  grid = rbind(grid, grid0)
  rm(grid0); gc()
}

# write plink scoring files
ix_keep = which(grid$nb_active > 0)
ix_int = which(grid$nb_interact > 0)

beta_main_col = sapply(ix_keep, FUN=function(i) (i-1)*K + 1)
beta_main_row = which(rowSums(beta[,beta_main_col] != 0) > 0)
beta_main = beta[beta_main_row, beta_main_col]
colnames(beta_main) = paste0('ADD_', which(grid$nb_active > 0))

beta_main_map = plink$map %>%
  filter(marker.ID %in% rownames(beta_main)) %>%
  select(chromosome, marker.ID, allele1, allele2) %>%
  cbind(data.matrix(beta_main))
write.table(beta_main_map,
            paste0(opt$out, '_beta_main.txt'),
            row.names=F, quote=F)
write.table(beta_main_map$marker.ID,
            paste0(opt$out, '_beta_main_snp.txt'),
            row.names=F, quote=F, col.names=F)


if (length(ix_int) > 0) {
  beta_int_col = sort(c(sapply(ix_int, FUN=function(i) (i-1)*K + 2:K)))
  beta_int_row = which(rowSums(beta[,beta_int_col] != 0) > 0)
  
  beta_int = beta[beta_int_row, beta_int_col]
  
  if (length(beta_int_row) == 1) {
    beta_int_map = plink$map %>%
      filter(marker.ID %in% rownames(beta)[beta_int_row]) %>%
      select(chromosome, marker.ID, allele1, allele2) %>%
      cbind(matrix(beta_int, nrow=1))
  } else {
    beta_int_map = plink$map %>%
      filter(marker.ID %in% rownames(beta_int)) %>%
      select(chromosome, marker.ID, allele1, allele2) %>%
      cbind(data.matrix(beta_int))
  }
  
  write.table(beta_int_map,
              paste0(opt$out, '_beta_int.txt'),
              row.names=F, quote=F)
  write.table(beta_int_map$marker.ID,
              paste0(opt$out, '_beta_int_snp.txt'),
              row.names=F, quote=F)
}

##################
### tuning PRS ###
##################
silent_command = ifelse(opt$verbose > 1, '', ' --silent')

system(paste0(opt$plink2,
              silent_command,
              ' --bfile ', opt$bfile,
              ' --keep ', opt$keep_tun,
              ' --extract ', opt$out, '_beta_main_snp.txt ',
              ' --score ', opt$out, '_beta_main.txt 2 3 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 5-', ncol(beta_main_map), ' ',
              '--out ', opt$out, '_main'))

system(paste0(opt$plink2,
              silent_command,
              ' --bfile ', opt$bfile,
              ' --keep ', opt$keep_tun,
              ' --extract ', opt$out, '_beta_int_snp.txt ',
              ' --score ', opt$out, '_beta_int.txt 2 3 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 5-', ncol(beta_int_map), ' ',
              '--out ', opt$out, '_int'))

prs = fread(paste0(opt$out, '_main.sscore'),
            col.names=c('FID', 'IID', paste0('ADD_', which(grid$nb_active > 0)))) %>%
  data.frame()

# integrate interactions
if (length(ix_int) > 0) {
  prs_int = fread(paste0(opt$out, '_int.sscore'),
                  col.names = c('FID', 'IID', c(outer(1:(K-1),which(grid$nb_interact > 0),  
                                                      FUN=function(k,i) paste0('INT', k, '_', i))))) %>%
    data.frame()
  
  Esum = readRDS(opt$Esum)
  E_tun = prs %>%
    select(FID, IID) %>%
    left_join(data %>% select(FID, IID, all_of(opt$interact_names))) %>%
    select(-FID, -IID) %>%
    data.frame()
  for (k in 1:(K-1)) E_tun[,k] = (E_tun[,k] - Esum$mean[k]) / Esum$sd[k]
  
  for (i in which(grid$nb_interact > 0)) {
    for (k in 1:(K-1)) {
      prs[,paste0('ADD_', i)] = prs[,paste0('ADD_', i)] + prs_int[,paste0('INT1_', i)] * E_tun[,1]
    }
  }
}

###################
### grid search ###
###################
grid$tun_r2 = 0
grid$tun_r2[which(grid$nb_active > 0)] = apply(prs %>% select(-FID, -IID), 2, FUN=function(x) 
  ifelse(sd(x)>0, cor(data$res, x)^2, 0))

ix_grid = which.max(grid$tun_r2)

data = data %>%
  left_join(prs %>% select(FID, IID, PRS_grid=paste0('ADD_', ix_grid)))

##################
### ensembling ###
##################
if (is.null(opt$maxsupp)) opt$maxsupp = Inf
if (is.null(opt$maxfactor)) opt$maxfactor = Inf

if (is.infinite(opt$maxsupp) & is.infinite(opt$maxfactor)) {
  # keep everything
  ix_grid_keep = which(grid$nb_active > 0)
} else {
  # limit support size
  grid0 = grid %>%
    mutate(tun_r2 = ifelse(nb_active > min(opt$maxsupp, opt$maxfactor*grid$nb_active[which.max(grid$tun_r2)]),
                           tun_r2, NA))
  ix_grid_order = intersect(order(-grid$tun_r2), which(grid$nb_active > 0))
  beta_grid = beta[beta_main_row,(ix_grid_order[1]-1)*6 + 1]
  ix_grid_keep = ix_grid_order[1]
  for (i in 2:length(ix_grid_order)) {
    # for (i in 2:10) {
    if (is.na(grid0$tun_r2[ix_grid_order[i]] > 0)) {
      beta_grid0 = cbind(beta_grid, beta[beta_main_row,(ix_grid_order[i]-1)*6 + 1])
      if (sum(rowSums(beta_grid0 != 0) > 0) < min(opt$maxsupp, opt$maxfactor*grid$nb_active[which.max(grid$tun_r2)])) {
        beta_grid = beta_grid0
        ix_grid_keep = c(ix_grid_keep, ix_grid_order[i])
      }
    }
  }
  ix_grid_keep = sort(ix_grid_keep)
}

prs_ens = prs %>% select(paste0('ADD_', ix_grid_keep)) %>% data.matrix()
cv = cv.glmnet(prs_ens, data$res, alpha=0)

data$PRS_ens = as.numeric(predict(cv, prs_ens, s='lambda.min'))

if (opt$verbose > 1) {
  cat('grid R2: ', with(data, cor(PRS_grid, res)^2), '\n')
  cat('ensemble R2: ', with(data, cor(PRS_ens, res)^2), '\n')
}

#######################
### save final beta ###
#######################
# grid
beta_grid = beta[,(ix_grid-1)*6 + 1:6]
colnames(beta_grid) = c('ADD', paste0('INT', 1:(K-1)))
ix_beta_grid = which(rowSums(beta_grid != 0) > 0)
beta_grid_map = plink$map %>%
  filter(marker.ID %in% rownames(beta_grid)[ix_beta_grid]) %>%
  select(chromosome, marker.ID, allele1, allele2) %>%
  cbind(data.matrix(beta_grid[ix_beta_grid,]))
rownames(beta_grid_map) = NULL

write.table(beta_grid_map, 
            paste0(opt$out, '_grid_beta.txt'),
            row.names=F, quote=F)

# ensemble
ix_grid_all_add = (ix_grid_keep-1)*6 + 1
ix_grid_all_int = (ix_grid_keep-1)*6 + 2
ix_ens_row = which(rowSums(beta[,ix_grid_all_add] != 0) > 0)
ix_ens_col= sort(c(sapply(ix_grid_keep, FUN=function(i) (i-1)*K + 1:K)))

coef_cv = coef(cv, s='lambda.min')[-1]

beta_ens = matrix(nrow=length(ix_ens_row), ncol=K)
rownames(beta_ens) = rownames(beta)[ix_ens_row]
colnames(beta_ens) = c('ADD', paste0('INT', 1:(K-1)))
for (k in 1:K)  {
  beta_ens[,k] = data.matrix(beta[ix_ens_row,(ix_grid_keep-1)*K + k]) %*% coef_cv
}

beta_ens_map = plink$map %>%
  filter(marker.ID %in% rownames(beta_ens)) %>%
  select(chromosome, marker.ID, allele1, allele2) %>%
  cbind(beta_ens)

supp_ens = sum(beta_ens[,1] != 0)
het_ens = max(colSums(beta_ens[,2:K] != 0))

write.table(beta_ens_map, 
            paste0(opt$out, '_ensemble_beta.txt'),
            row.names=F, quote=F)

if (opt$cleanup == 1) {
  if (opt$verbose > 1) cat('cleaning up files \n')
  system(paste0('rm ', opt$out, '*_beta.RDS'))
  system(paste0('rm ', opt$out, '*_results.RDS'))
  system(paste0('rm ', opt$out, '*.sscore'))
  system(paste0('rm ', opt$out, '*_beta_int*'))
  system(paste0('rm ', opt$out, '*_beta_main*'))
  system(paste0('rm ', opt$out, '*.log'))
}

