library(dplyr)
library(data.table)
library(tidyr)
library(mvtnorm)

#######################################
### use GCTA to generate phenotypes ###
#######################################
# model='hm3'; rho=1; anc='EUR'; gr=1; GA=1
# model='hm3'; rho=1; anc='EUR'; gr=0.6; GA=2

#generate phenotypes given different genetic architecture
args = commandArgs(trailingOnly = T)

# mixed or shared causal SNPs
model = args[[1]]

#rho represent causal snps proportion
rho = as.numeric(args[[2]])

#gr represent genetic relatedness
gr = as.numeric(args[[3]])

#GA represent genetic architecture (amount of heterogeneity)
GA = as.numeric(args[[4]])

# model='hm3'
model='aou'
grGA = data.frame(
  GA = c(1, 1, 2, 2),
  gr = c(0.6, 0.8, 0.4, 0.6)
)
# rho=1; ix=1; anc='EUR'

for (rho in 1:4) {
  for (ix in 1:nrow(grGA)) {
    gr = grGA$gr[ix]
    GA = grGA$GA[ix]
    cat(rho, gr, GA, '\n')
    
    all_pheno = data.frame()
    for (anc in c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')) {
      select.cau = read.table(paste0('Generate/Beta/select.cau_', anc, '_rho',rho,'_gr',gr, '_GA',GA, '_', model), header=F)
      colnames(select.cau) = c('snpid','effect_size')
      
      #plink format used minor allele as coding allele
      #the fifth column is minor allele
      #need to match the minor allele with the coding allele
      snp.infor = read.table(paste0('PLINK/', model, '_', ceiling(rho), '_causal.bim'))
      colnames(snp.infor) = c('chr','snpid','posg','pos','minor','major')
      
      select.cau.infor = left_join(select.cau,snp.infor,by='snpid')
      
      select.cau.infor.split = select.cau.infor %>% separate(snpid,into=c('rsid','pos2','noncoding','coding'),sep=':')
      
      idx = which(select.cau.infor.split$coding!=select.cau.infor.split$minor)
      coding_effect_size = select.cau$effect_size
      minor_effect_size = coding_effect_size
      minor_effect_size[idx] = -coding_effect_size[idx]
      herit = nrow(select.cau)*var(minor_effect_size)
      write.table(herit, paste0('Generate/Heritability/h2_', anc, '_rho',rho,'_gr',gr, '_GA',GA, '_', model, '.txt'),
                  row.names=F, quote=F, col.names=F)
      
      causal_map = select.cau.infor %>% 
        mutate(SNP=snpid) %>%
        separate(snpid,into=c('rsid','pos2','noncoding','coding'),sep=':') %>%
        select(SNP, rsid, chr, pos2, a1=coding, a0=noncoding, beta=effect_size)
      write.table(causal_map, paste0('PHENO/causal_', 
                                     anc, '_rho',rho,'_gr',gr,'_GA',GA,'_', model, '.txt'),
                  row.names=F, quote=F)
      
      system(paste0('~/software/plink2 ',
                    '--silent ',
                    '--bfile PLINK/', model, '_', ceiling(rho), '_causal ',
                    '--keep ID/', anc, '.ID ',
                    '--score PHENO/causal_', anc, '_rho',rho,'_gr',gr,'_GA',GA,'_', model, '.txt ',
                    ' 1 5 7 header cols=+scoresums,-scoreavgs ',
                    '--out PHENO/causal_pheno_', anc, '_rho',rho,'_gr',gr,'_GA',GA,'_', model))
      
      pheno = fread(paste0('PHENO/causal_pheno_', anc, '_rho',rho,'_gr',gr,'_GA',GA,'_', model, '.sscore'))
      
      sigma_g = var(pheno$SCORE1_SUM)
      sigma_e = sigma_g * (1/herit - 1)
      
      set.seed(123)
      pheno_anc = pheno %>% select(FID=`#FID`, IID) %>% data.frame()
      for (rep in 1:10) {
        pheno_anc[,paste0('PHENO', rep)] = pheno$SCORE1_SUM + rnorm(nrow(pheno), 0, sqrt(sigma_e))
      }
      
      all_pheno = rbind(all_pheno, pheno_anc)
      
      fwrite(all_pheno, paste0('PHENO/phenotypes_rho',rho,
                               '_gr',gr,'_GA',GA,'_', model, '.txt'), sep='\t')
    }
  }
}

