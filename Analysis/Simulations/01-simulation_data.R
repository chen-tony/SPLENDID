# data from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP
# full data with all samples and SNPs is named all_chr.[bed/bim/fam]
## first row of all_chr.bim: [1 rs3107975:55326:T:C 0 55326 C T]
## first row of all_chr.fam: [AFR_1 AFR_1 0 0 0 -9]

library(dplyr)
library(data.table)
library(tidyr)

#################
## DATA SPLITS ##
#################
fam = fread('all_chr.fam') %>%
  mutate(ancestry = substr(V1, 1, 3))

tunval = fam %>%
  group_by(ancestry) %>%
  slice_tail(n=20000) %>% 
  ungroup()

tuning = tunval %>%
  group_by(ancestry) %>%
  slice_head(n=10000) %>%
  ungroup()

validation = tunval %>%
  group_by(ancestry) %>%
  slice_tail(n=10000) %>%
  ungroup()

training = fam %>%
  group_by(ancestry) %>%
  slice_head(n=100000) %>%
  ungroup()

fam_map = fam %>%
  select(FID=V1, IID=V2) %>%
  mutate(anc = substr(FID, 1, 3)) %>%
  mutate(split = case_when(
    FID %in% training$V1 ~ 'Training',
    FID %in% tuning$V1 ~ 'Tuning',
    FID %in% validation$V1 ~ 'Validation'
  ))

with(fam_map, table(anc, split))

write.table(fam_map, paste0(home, 'simulation_data_splits.txt'),
            row.names=F, quote=F)

fam = fread(paste0(home, 'simulation_data_splits.txt'), header=T)


# tuning samples (subsetted from full 50k)
tun_id = fam %>% filter(split == 'Tuning') %>% select(FID, IID)

subtun_id = rbind(
  tun_id %>% filter(substr(FID, 1, 3) == 'AFR') %>% slice_head(n=2000),
  tun_id %>% filter(substr(FID, 1, 3) == 'AMR') %>% slice_head(n=2000),
  tun_id %>% filter(substr(FID, 1, 3) == 'EAS') %>% slice_head(n=2000),
  tun_id %>% filter(substr(FID, 1, 3) == 'EUR') %>% slice_head(n=10000),
  tun_id %>% filter(substr(FID, 1, 3) == 'SAS') %>% slice_head(n=2000)
)
table(substr(subtun_id$FID, 1, 3))
write.table(subtun_id,
            paste0(home, 'ID/subtuning.id'),
            row.names=F, quote=F, col.names=F)

# validation samples
val_id = fam %>% filter(split == 'Validation') %>% select(FID, IID)
write.table(val_id,
            paste0(home, 'ID/validation.id'),
            row.names=F, quote=F, col.names=F)

# tuning and validation samples combined
subtunval_id = rbind(subtun_id, val_id)
table(substr(subtunval_id$FID, 1, 3))
write.table(subtunval_id,
            paste0(home, 'ID/subtunval.id'),
            row.names=F, quote=F, col.names=F)

# training samples
## 70k EUR + 30k nonEUR
train_id2a = c(
  fam$FID[substr(fam$FID, 1, 3) == 'EUR'][1:70000],
  fam$FID[substr(fam$FID, 1, 3) == 'AFR'][1:7500],
  fam$FID[substr(fam$FID, 1, 3) == 'EAS'][1:7500],
  fam$FID[substr(fam$FID, 1, 3) == 'SAS'][1:7500],
  fam$FID[substr(fam$FID, 1, 3) == 'AMR'][1:7500]
)
table(substr(train_id2a, 1, 3))
length(train_id2a)

write.table(cbind(train_id2a, train_id2a), 'ID/train2a.id',
            row.names=F, quote=F, col.names=F)

## 20k EUR + 80k nonEUR
train_id4a = c(
  fam$FID[substr(fam$FID, 1, 3) == 'EUR'][1:20000],
  fam$FID[substr(fam$FID, 1, 3) == 'AFR'][1:20000],
  fam$FID[substr(fam$FID, 1, 3) == 'EAS'][1:20000],
  fam$FID[substr(fam$FID, 1, 3) == 'SAS'][1:20000],
  fam$FID[substr(fam$FID, 1, 3) == 'AMR'][1:20000]
)
table(substr(train_id4a, 1, 3))
length(train_id4a)

write.table(cbind(train_id4a, train_id4a), 'ID/train4a.id',
            row.names=F, quote=F, col.names=F)

#####################
## SUBSET VARIANTS ##
#####################
bim = fread('all_chr.bim')

rsid_code = bim %>%
  select(V2) %>%
  separate(V2, into=c('rsid', 'pos', 'alt', 'ref'), sep=':') %>% 
  pull(rsid)

write.table('sim_snps.txt', col.names=F, row.names=F, quote=F)

bim_map = data.frame(bim, rsid=rsid_code) %>%
  select(CHR=V1, rsid, snp_code=V2, pos=V4, alt=V5, ref=V6)

write.table('sim_mapping.txt', row.names=F, quote=F)

##
bim_map = fread('sim_mapping.txt')

write.table(bim_map$rsid, 'Generate/sim_mapping_snp.txt',
            row.names=F, quote=F, col.names=F)

aou = fread('AOU/aou_ukb_map_final.txt') # SNPs from All of Us genotyping array

sim_map = bim_map %>%
  filter(rsid %in% aou$rsid)

write.table(sim_map$rsid, 'Generate/aou_rsid.txt',
            row.names=F, quote=F, col.names=F)
write.table(sim_map$snp_code, 'Generate/aou_snp.txt',
            row.names=F, quote=F, col.names=F)


#############
## RUN PCA ##
#############
# data from 1000G: https://www.cog-genomics.org/plink/2.0/resources
# processing data from 1000G: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html

system(paste0('~/software/plink --bfile all_phase3 ',
              '--extract Generate/aou_rsid.txt ',
              '--indep-pairwise 500 5 0.05 --out PCA/sim'))

pruned = fread('PCA/sim.prune.in', header=F)$V1
length(pruned)

pruned_rsid = bim_map %>%
  filter(rsid %in% pruned) %>%
  pull(snp_code)
write.table(pruned_rsid, 'PCA/sim.prune.in_rsid', row.names=F, quote=F, col.names=F)
length(pruned_rsid)

# PCA in 1000G
system(paste0('~/software/plink2 --bfile all_phase3 ',
              '--extract PCA/sim.prune.in ',
              '--pca 20 allele-wts ',
              '--out PCA/sim_1000G_pca'))

pca_loadings = fread('PCA/sim_1000G_pca.eigenvec.allele') %>%
  rename('CHR'=`#CHROM`) %>%
  left_join(bim_map %>% select('rsid', 'snp_code'), by=c('ID'='rsid')) %>%
  select(CHR, ID, SNP=snp_code, REF, ALT, A1, paste0('PC', 1:20)) %>%
  filter(ALT==A1)

write.table(pca_loadings, 'PCA/sim_1000G_pca.eigenvec.allele.bim',
            row.names=F, quote=F)

# PC scoring in 1000G
fam = fread('ALL/all_chr.fam')

set.seed(123)
pca_test_id = sample(1:nrow(fam), 1000)
fam[pca_test_id,] %>%
  select(FID=V1, IID=V2) %>%
  fwrite('PCA/sim_test.id', sep='\t')

system(paste0('plink2 --bfile all_phase3 ',
              '--extract PCA/sim.prune.in ',
              '--score PCA/sim_1000G_pca.eigenvec.allele.bim 2 6 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 7-26 ',
              '--out PCA/sim_1000G_pca'))

system(paste0('plink2 --bfile ', sim, 'ALL/all_chr ',
              '--extract PCA/sim.prune.in_rsid ',
              '--keep PCA/sim_test.id ',
              '--score PCA/sim_1000G_pca.eigenvec.allele.bim 3 6 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 7-26 ',
              '--out PCA/sim_test_pca'))

# PC scoring in simulation data
system(paste0('plink2 all_chr ', 
              '--extract PCA/sim.prune.in_rsid ',
              '--score PCA/sim_1000G_pca.eigenvec.allele.bim 3 6 header ', 
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 7-26 ',
              '--out /PCA/all_sim_pca'))


##############################
## PREPARE DATA FOR PRS-CSx ##
##############################
# tuning / validation data file
# ~/software/plink2 \
# --bfile all_chr \
# --keep ID/subtunval.id \
# --make-bed \
# --freq \
# --out PLINK/subtunval

##############################
## PREPARE DATA FOR PROSPER ##
##############################
tuning_id = fread('ID/subtuning.id', header=F, col.names=c('FID', 'IID'))
validation_id = fread('ID/validation.id', header=F, col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore')

pc_tun = pc %>% rename('FID'='#FID') %>% 
  filter(FID %in% tuning_id$FID) %>%
  select(FID, IID, PC=contains('SCORE'))

pc_val = pc %>% rename('FID'='#FID') %>% 
  filter(FID %in% validation_id$FID) %>%
  select(FID, IID, PC=contains('SCORE'))

# IDs
for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
  if (anc == 'ALL') {
    tun_anc = tuning_id
    val_anc = validation_id
  } else {
    tun_anc = tuning_id %>%
      filter(substr(FID, 1, 3) == anc)
    val_anc = validation_id %>%
      filter(substr(FID, 1, 3) == anc)
  }
  
  fwrite(tun_anc, paste0('ID/subtuning_', anc, '.id'), sep='\t')
  fwrite(val_anc, paste0('ID/validation_', anc, '.id'), sep='\t')
}

# change plink formatting
for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'ALL')) {
  print(anc)
  tun_bim = fread(paste0('PROSPER/PLINK/', anc, '_tuning.bim')) %>%
    separate(V2,into=c('rsid','pos2','noncoding','coding'), sep=':') %>%
    select(V1, rsid, V3, V4, V5, V6) %>%
    data.frame()
  write.table(tun_bim, paste0('PROSPER/PLINK/', anc, '_tuning.bim'),
              row.names=F, quote=F, col.names=F)
  
  val_bim = fread(paste0('PROSPER/PLINK/', anc, '_validation.bim')) %>%
    separate(V2,into=c('rsid','pos2','noncoding','coding'), sep=':') %>%
    select(V1, rsid, V3, V4, V5, V6) %>%
    data.frame()
  write.table(val_bim, paste0('PROSPER/PLINK/', anc, '_validation.bim'),
              row.names=F, quote=F, col.names=F)
}


grGA = data.frame(
  gr = c(0.6, 0.8, 0.4, 0.6),
  GA = c(1, 1, 2, 2)
) 

# phenotypes
model = 'aou'
for (rho in 1:4) {
  for (grGA_ix in 1:nrow(grGA)) {
    gr = grGA$gr[grGA_ix]
    GA = grGA$GA[grGA_ix]
    if (file.exists(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt'))) {
      cat(rho, gr, GA, model, '\n')
      pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt'))
      
      for (rep in 1:10) {
        pheno_tunval = pheno %>%
          filter(FID %in% tuning_id$FID | FID %in% validation_id$FID) %>%
          select(FID, IID, Y=paste0('PHENO', rep)) %>%
          left_join(pc)
        
        # regress PC main effects from outcome
        fit_tun = lm(paste0('Y~', paste(paste0('SCORE', 1:20, '_SUM'), collapse='+')), pheno_tunval %>% filter(FID %in% tuning_id$FID))
        fit_val = lm(paste0('Y~', paste(paste0('SCORE', 1:20, '_SUM'), collapse='+')), pheno_tunval %>% filter(FID %in% validation_id$FID))
        
        pheno_tunval$res = 0
        pheno_tunval$res[pheno_tunval$FID %in% tuning_id$FID] = residuals(fit_tun)
        pheno_tunval$res[pheno_tunval$FID %in% validation_id$FID] = residuals(fit_val)
        
        pheno_tun = pheno_tunval %>%
          filter(FID %in% tuning_id$FID) %>%
          select(FID, IID, Y=res) %>%
          data.frame()
        
        pheno_val = pheno_tunval %>%
          filter(FID %in% validation_id$FID) %>%
          select(FID, IID, Y=res) %>%
          data.frame()
        
        for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'ALL')) {
          if (anc != 'ALL') {
            pheno_tun_anc = pheno_tun %>%
              filter(substr(FID, 1, 3) == anc)
            pheno_val_anc = pheno_val %>%
              filter(substr(FID, 1, 3) == anc)
          } else {
            pheno_tun_anc = pheno_tun
            pheno_val_anc = pheno_val
            
          }
          
          write.table(pheno_tun_anc, 
                      paste0('PROSPER/Pheno/', anc, '_tun_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', model, '_rep', rep, '.txt'),
                      row.names=F, col.names=F, quote=F)
          write.table(pheno_val_anc, 
                      paste0('PROSPER/Pheno/', anc, '_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', model, '_rep', rep, '.txt'),
                      row.names=F, col.names=F, quote=F)
          
          rm(pheno_tun_anc, pheno_val_anc); gc()
        }
      }
      
      rm(pheno, pheno_tun, pheno_val); gc()
    }
  }
}

##########################
## 1000G LD FOR CT-SLEB ##
##########################
for (anc in c('EUR', 'AFR', 'EAS', 'AMR', 'SAS')) {
  system(paste0('~/software/plink2 ',
                '--bfile all_phase3 ',
                '--extract Generate/aou_rsid.txt ',
                '--keep 1000G/Continental_IDs/', anc, '.id ',
                '--make-bed ',
                '--out PLINK/', anc, '_1000G'))
}

bim = fread('PLINK/subtunval.bim') %>%
  tidytable::separate(V2,into=c('rsid','pos2','noncoding','coding'), sep=':') %>%
  select(V1, rsid, V3, V4, V5, V6) %>%
  data.frame()

write.table(bim, 'PLINK/subtunval.bim',
            row.names=F, quote=F, col.names=F)


