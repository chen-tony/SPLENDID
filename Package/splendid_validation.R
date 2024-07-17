##############
### inputs ###
##############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option('--bfile', type='character',
              default=NULL,
              help='plink file (without .bed/.bim/.fam extension)',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='splendid',
              help='output file name',
              metavar='character'),
  make_option('--plink2', type='character',
              default='~/plink2',
              help='path to plink2 software'),
  make_option('--keep-val', type='character',
              default=NULL,
              help='text file with FID/IID of validation data',
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
  make_option('--Esum', type='character',
              default=NULL,
              help='interaction summary file from splendid_summaries.R'),
  make_option('--ncores', type='integer',
              default=NULL,
              help='number of cores for parallel computing (recommended = # lambda1 values x # p-value thresholds)'),
  make_option('--cleanup', type='integer',
              default=1,
              help='1 for cleanup of intermediate files, 0 for no actions'),
  make_option('--verbose', type='integer',
              default=1,
              help='1 for messages, 0 for no messages')
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
  library(glmnet)
  library(Matrix)
})

opt$covar_names = strsplit(opt$covar_names, split=',')[[1]] # turn comma-separated into a character vector
opt$interact_names = strsplit(opt$interact_names, split=',')[[1]] # turn comma-separated into a character vector

K = length(opt$interact_names) + 1

val_ids = fread(opt$keep_val, header=F, col.names=c('FID', 'IID'))

cov = fread(opt$covar) %>%
  select(FID, IID, all_of(opt$covar_names))
pheno = fread(opt$pheno) %>%
  select(FID, IID, Y=all_of(opt$pheno_name))

data = pheno %>%
  left_join(cov) %>%
  filter(FID %in% val_ids$FID & IID %in% val_ids$IID)
data$res = residuals(lm(Y~., data %>% select(-FID, -IID)))

######################
### validation PRS ###
######################
silent_command = ifelse(opt$verbose > 1, '', ' --silent')

system(paste0(opt$plink2,
              silent_command,
              ' --bfile ', opt$bfile,
              ' --keep ', opt$keep_val,
              ' --extract ', opt$out, '_grid_beta_snp.txt ',
              ' --score ', opt$out, '_grid_beta.txt 2 3 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 5-', K+4, ' ',
              '--out ', opt$out, '_grid_val'))

prs_grid = fread(paste0(opt$out, '_grid_val.sscore'),
                 col.names=c('FID', 'IID', 'ADD', paste0('INT_', 1:(K-1)))) %>%
  data.frame()

system(paste0(opt$plink2,
              silent_command,
              ' --bfile ', opt$bfile,
              ' --keep ', opt$keep_val,
              ' --extract ', opt$out, '_ensemble_beta_snp.txt ',
              ' --score ', opt$out, '_ensemble_beta.txt 2 3 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 5-', K+4, ' ',
              '--out ', opt$out, '_ensemble_val'))

prs_ens = fread(paste0(opt$out, '_ensemble_val.sscore'),
                col.names=c('FID', 'IID', 'ADD', paste0('INT_', 1:(K-1)))) %>%
  data.frame()

# integrate interactions
Esum = readRDS(opt$Esum)
E_val = prs_ens %>%
  select(FID, IID) %>%
  left_join(data %>% select(FID, IID, all_of(opt$interact_names))) %>%
  select(-FID, -IID) %>%
  data.frame()
for (k in 1:(K-1)) E_val[,k] = (E_val[,k] - Esum$mean[k]) / Esum$sd[k]

prs_grid$INT = prs_grid$ADD
prs_ens$INT = prs_grid$ADD
for (k in 1:(K-1)) {
  prs_grid$INT = prs_grid$INT + prs_grid[,paste0('INT_', k)] * E_val[,k]
  prs_ens$INT = prs_ens$INT + prs_ens[,paste0('INT_', k)] * E_val[,k]
}

data = data %>%
  left_join(prs_grid %>% select(FID, IID, PRS_grid=INT)) %>%
  left_join(prs_ens %>% select(FID, IID, PRS_ens=INT))

###############
### RESULTS ###
###############
if (opt$verbose > 1) {
  cat('grid R2: ', with(data, cor(PRS_grid, res)^2), '\n')
  cat('ensemble R2: ', with(data, cor(PRS_ens, res)^2), '\n')
}

if (opt$cleanup == 1) {
  if (opt$verbose > 1) cat('cleaning up files \n')
  system(paste0('rm ', opt$out, '*.log'))
}
