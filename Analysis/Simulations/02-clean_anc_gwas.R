# re-format ancestry-specific GWAS from PLINK for PRS methods

library(dplyr)
library(data.table)
library(tidyr)
library(doParallel)

train = args[[1]]
rho = as.numeric(args[[2]])
gr = as.numeric(args[[3]])
GA = as.numeric(args[[4]])
model = args[[5]]

registerDoParallel(cores=10)# Registers 10 cores for parallel computing
sim.out = foreach(rep =1:10, .combine=rbind) %dopar% {
  print(rep)
  for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    gwas = fread(paste0('GWAS/gwas_', anc, train, '_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.PHENO', rep, '.glm.linear'))  %>%
      filter(TEST == 'ADD') %>%
      mutate(A0 = case_when(
        A1==ALT ~ REF,
        T ~ ALT
      )) %>%
      mutate(V1 = ID) %>%
      separate(ID,into=c('rsid','pos2','noncoding','coding'), sep=':') %>%
      select(rsid, chr=`#CHROM`, 
             a1=A1, a0=A0, beta=BETA, 
             beta_se=SE, n_eff=OBS_CT, 
             eaf=A1_FREQ, p=P, pos=POS, rsid0=V1) %>%
      na.omit()
    cat(nrow(gwas), ' ')
    
    # PROSPER: rsid, chr, a1, a0, beta, beta_se, n_eff
    fwrite(gwas, 
           paste0('PROSPER/GWAS/', anc, '_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'),
           sep='\t')
    
    # PRSCSX: SNP, A1, A2, BETA, SE
    fwrite(gwas %>% select(SNP=rsid, A1=a1, A2=a0, BETA=beta, SE=beta_se), 
           paste0('PRSCSX/GWAS/', anc, '_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'),
           sep='\t')
    
    # CTSLEB
    fwrite(gwas %>% select(CHR=chr, SNP=rsid, BP=pos, A1=a1, 
                           BETA=beta, SE=beta_se, P=p, rs_id=rsid0), 
           paste0('CTSLEB/GWAS/', anc, '_gwas_train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.txt'),
           sep='\t')
    
    rm(gwas); gc()
  }
  
  rep
}