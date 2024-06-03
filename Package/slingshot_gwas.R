##############
### inputs ###
##############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option('--bfile', type='character',
              default=NULL,
              help='plink file without .bed/.bim/.fam extensions',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='slingshot',
              help='output file name',
              metavar='character'),
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
  make_option('--plink2', type='character',
              default='~/plink2',
              help='path for plink2'),
  make_option('--ncores', type='integer',
              default=NULL,
              help='number of cores for parallel computing'),
  make_option('--verbose', type='integer',
              default=1,
              help='0 for no messages, 1 for messages, 2 for more messages')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

if (opt$verbose > 1) {
  cat('inputs: \n')
  print(opt)
}

if (is.null(opt$interact_names)) {
  cat('ERROR: must specify interactions')
  quit('no')
}

######################
### clean commands ###
######################
# sample IDs
keep_command = ifelse(
  is.null(opt$keep),
  '', 
  paste0(' --keep ', opt$keep)
)

# SNP IDs
extract_command = ifelse(
  is.null(opt$extract),
  '', 
  paste0(' --extract ', opt$extract)
)

# interaction numbering
all_covar_names = strsplit(opt$covar_names, split=',')[[1]]
all_interact_names = paste0('SNPx', strsplit(opt$interact_names, split=',')[[1]])

all_reg_names = c('SNP', all_covar_names, paste0('SNPx', all_covar_names))

params_idx = which(all_reg_names %in% c('SNP', all_covar_names, all_interact_names))

tests_idx = which(all_reg_names[params_idx] %in% all_interact_names)

params_command = paste0(
  ' --parameters ', 
  paste(params_idx, collapse=',')
)

tests_command = paste0(
  ' --tests ',
  paste(tests_idx, collapse=',')
)

# multi-threading
threads_command = ifelse(
  is.null(opt$ncores),
  '',
  paste0(' --threads ', opt$ncores)
)

# messaging
silent_command = ifelse(
  opt$verbose == 0,
  ' --silent',
  ''
)

################
### run GWAS ###
################
if (opt$verbose > 0) cat('RUNNING GWAS \n')
system(paste0(
  opt$plink2, 
  ' --bfile ', opt$bfile, 
  silent_command, 
  keep_command, 
  extract_command, 
  ' --linear interaction hide-covar ',
  threads_command,
  ' --pheno ', opt$pheno, 
  ' --pheno-name ', opt$pheno_name,
  ' --covar ', opt$covar, 
  ' --covar-name ', opt$covar_name,
  ' --covar-variance-standardize ',
  params_command, 
  tests_command, 
  ' --out ', opt$out
))

##################
### clean GWAS ###
##################
if (opt$verbose > 0) cat('CLEANING GWAS \n')

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

gwas = fread(paste0(opt$out, '.', opt$pheno_name, '.glm.linear'))

clean_gwas = gwas %>%
  filter(TEST == 'ADD' | grepl('USER_', TEST)) %>%
  select(ID, TEST, P) %>%
  mutate(TEST = ifelse(TEST=='ADD', 'ADD', 'HET')) %>%
  mutate(P = as.numeric(P)) %>%
  spread(key=TEST, value=P) %>%
  select(marker.ID=ID, P_ADD=ADD, P_HET=HET)
clean_gwas[is.na(clean_gwas)] = 1

###############
### outputs ###
###############
fwrite(clean_gwas, paste0(opt$out, '_interaction_gwas.txt'), sep='\t')

rm = file.remove(paste0(paste0(opt$out, '.', opt$pheno_name, '.glm.linear')))

cat(paste0('GWAS written to ', opt$out, '_interaction_gwas.txt\n'))
