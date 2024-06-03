##############
### inputs ###
##############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c('-b', '--bfile'), type='character',
              default=NULL,
              help='plink file (without file .bed/.bim/.fam extensions)',
              metavar='character'),
  make_option(c('-o', '--out'), type='character',
              default='slingshot',
              help='output file name',
              metavar='character'),
  make_option('--keep', type='character',
              default=NULL,
              help='text file with FID/IID of training (and tuning) data',
              metavar='character'),
  make_option('--extract', type='character',
              default=NULL,
              help='text file with SNP IDs',
              metavar='character'),
  make_option('--ncores', type='integer',
              default=1,
              help='number of cores for parallel computing'),
  make_option(c('-v', '--verbose'), type='integer',
              default=1,
              help='0 for no messages, 1 for messages, 2 for more messages')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

if (opt$verbose > 1) {
  cat('inputs: \n')
  print(opt)
}

###################
### create file ###
###################
suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
})

# subset samples
if (!is.null(opt$keep)) {
  fam = fread(paste0(opt$bfile, '.fam'))
  all_ids = fread(paste0(opt$keep), header=F)
  ind.row = which(fam$V1 %in% all_ids$V1 & fam$V2 %in% all_ids$V2)
  
  if (opt$verbose > 0) {
    cat(length(ind.row), 'samples kept\n')
  }
} 

# subset SNPs
if (!is.null(opt$extract)) {
  bim = fread(paste0(opt$bfile, '.bim'))
  all_snps = fread(paste0(opt$extract), header=F)$V1
  ind.col = which(bim$V2 %in% all_snps)
  
  if (opt$verbose > 0) {
    cat(length(ind.col), 'SNPs kept\n')
  }
  
} 

if (file.exists(paste0(opt$out, '.bk'))) {
  if (opt$verbose > 1) cat('removing existing ', paste0(opt$out, '.bk'), '\n')
  rm = file.remove(paste0(opt$out, '.bk'))
}

# create bigsnpr object
if (is.null(opt$keep) & is.null(opt$extract)) {
  # keep all SNPs and IDs
  snp_readBed2(
    paste0(opt$bfile, '.bed'),
    backingfile = paste0(opt$out),
    ncores = opt$ncores
  )
} else if (is.null(opt$keep)) {
  # only subset on SNPs
  snp_readBed2(
    paste0(opt$bfile, '.bed'),
    backingfile = paste0(opt$out),
    ind.col = ind.col,
    ncores = opt$ncores
  )
} else if (is.null(opt$extract)) {
  # only subset on IDs
  snp_readBed2(
    paste0(opt$bfile, '.bed'),
    backingfile = paste0(opt$out),
    ind.col = ind.col,
    ncores = opt$ncores
  )
} else {
  # subset both SNPs and IDs
  snp_readBed2(
    paste0(opt$bfile, '.bed'),
    backingfile = paste0(opt$out),
    ind.row = ind.row,
    ind.col = ind.col,
    ncores = opt$ncores
  )
}

if (opt$verbose > 0) {
  cat('bigsnpr file written to ', paste0(opt$out, '[.bk/.rds] \n'))
}
