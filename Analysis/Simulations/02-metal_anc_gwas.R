# meta-analysis of ancestry-specific GWAS for iPGS+refit

# DOWNLOAD: https://csg.sph.umich.edu/abecasis/metal/download/
# DOCUMENTATION: https://genome.sph.umich.edu/wiki/METAL_Documentation

####################
## CREATE SCRIPTS ##
####################
# train=1; rho=3; gr=0.6; GA=3; rep=1

grGA = data.frame(
  GA = c(1, 1, 2, 2),
  gr = c(0.6, 0.8, 0.4, 0.6)
)
model = 'aou'

# create metal script files
for (train in c('2a', '4a')) {
  for (rho in 1:4) {
    for (ix in 1:nrow(grGA)) {
      gr = grGA$gr[ix]
      GA = grGA$GA[ix]
      for (rep in 1:10) {
        script = paste0('SCHEME STDERR \n',
                        # EUR
                        'MARKER rsid0 \n', 
                        'ALLELE a0 a1 \n',
                        'PVALUE p \n', 
                        'EFFECT beta \n',
                        'STDERR beta_se \n',
                        'PROCESS PROSPER/GWAS/EUR_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt \n',
                        # AFR
                        'MARKER rsid0 \n', 
                        'ALLELE a0 a1 \n',
                        'PVALUE p \n', 
                        'EFFECT beta \n',
                        'STDERR beta_se \n',
                        'PROCESS PROSPER/GWAS/AFR_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt \n',
                        # EAS
                        'MARKER rsid0 \n', 
                        'ALLELE a0 a1 \n',
                        'PVALUE p \n', 
                        'EFFECT beta \n',
                        'STDERR beta_se \n',
                        'PROCESS PROSPER/GWAS/EAS_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt \n',
                        # SAS
                        'MARKER rsid0 \n', 
                        'ALLELE a0 a1 \n',
                        'PVALUE p \n', 
                        'EFFECT beta \n',
                        'STDERR beta_se \n',
                        'PROCESS PROSPER/GWAS/SAS_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt \n',
                        # AMR
                        'MARKER rsid0 \n', 
                        'ALLELE a0 a1 \n',
                        'PVALUE p \n', 
                        'EFFECT beta \n',
                        'STDERR beta_se \n',
                        'PROCESS PROSPER/GWAS/AMR_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt \n',
                        
                        'OUTFILE METAL/Meta/metal_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '_ .tbl \n',
                        'ANALYZE HETEROGENEITY \n',
                        'QUIT')
        
        writeLines(script, paste0('METAL/Script/metal_script_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'))
      }
    }
  }
}

###############
## RUN METAL ##
###############
args = commandArgs(trailingOnly = T)
train = args[[1]]
rho = as.numeric(args[[2]])
gr = as.numeric(args[[3]])
GA = as.numeric(args[[4]])
rep = args[[5]]
model = args[[6]]

system(paste0('METAL-master/build/metal/metal ',
              'METAL/Script/metal_script_train', train, '_rho', rho, 
              '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'))




