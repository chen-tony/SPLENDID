# CT-SLEB: https://github.com/andrewhaoyu/CTSLEB

# run 2-way CT and EB between EUR and non-EUR

##################
## PREPARE DATA ##
##################
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(CTSLEB)
})

sumstat=paste0(scratch, 'CTSLEB/GWAS/')

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]]
GA = params[[4]]
rep = as.numeric(params[[5]])
model = params[[6]]
anc = params[[7]]
# train='1a'; rho=1; gr=0.6; GA=1; rep=1; model='hm3'; anc='AFR'
# train='4a'; rho=1; gr=0.6; GA=1; rep=1; model='aou'; anc='AFR'

tun_id = fread('ID/subtuning.id', header=F, col.names=c('FID', 'IID'))
val_id = fread('ID/validation.id', header=F, col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore', col.names=c('FID', 'IID', paste0('PC', 1:20)))

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) %>%
  select(FID, IID, Y=paste0('PHENO', rep)) 

data = pheno %>%
  left_join(pc)

tunval_data = data %>%
  filter(FID %in% c(tun_id$FID, val_id$FID))

# compute residuals for tuning
fit_tun = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% tun_id$FID))
fit_val = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% val_id$FID))

tunval_data$res = 0
tunval_data$res[tunval_data$FID %in% tuning_id$FID] = residuals(fit_tun)
tunval_data$res[tunval_data$FID %in% validation_id$FID] = residuals(fit_val)

y_tun = tunval_data %>%
  select(FID, IID, Y=res) %>%
  filter(FID %in% tun_id$FID) %>%
  filter(substr(FID, 1, 3) == anc)

rm(pheno, pc, data, tunval_data); gc()

out_dir = paste0('CTSLEB/Analysis/', anc, '_', train, '_', rho, '_', gr, '_', GA, '_', rep, '/')
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

##############
## CLUMPING ##
##############
print('CT')
PRS_farm <- SetParamsFarm(
  threads = 2,
  pthres = c(5e-08, 5e-07, 5e-06, 5e-05, 5e-04, 0.005, 0.05, 0.1),
  plink19_exec = '~/software/plink',
  plink2_exec = '~/software/plink2'
)

sum_EUR = fread(paste0('CTSLEB/GWAS/EUR_gwas_train', train, '_rho', rho, 
                       '_gr', gr, '_rep', rep, 
                       '_GA', GA, '_', model, '.txt')) 
sum_AFR = fread(paste0('CTSLEB/GWAS/', anc, '_gwas_train', train, '_rho', rho, 
                       '_gr', gr, '_rep', rep, 
                       '_GA', GA, '_', model, '.txt')) 

prs_mat = dimCT(
  results_dir = out_dir, 
  sum_ref = sum_EUR,
  sum_target = sum_AFR,
  ref_plink = 'PLINK/EUR_1000G', 
  target_plink = paste0('PLINK/', anc, '_1000G'), 
  test_target_plink = paste0('PROSPER/PLINK/', anc, '_tuning'), 
  out_prefix = anc,
  params_farm = PRS_farm
)
saveRDS(prs_mat, paste0('CTSLEB/PRS/prsCT_', anc, '_train', train, '_rho', rho, 
                        '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS'))
dim(prs_mat)

score = fread(paste0(out_dir, 'temp/score_file'))
saveRDS(score$V1, paste0('CTSLEB/PRS/snpCT_', anc, '_train', train, '_rho', rho, 
                         '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS'))
rm(score)

n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test <- rep(0,n.total.prs)

prs_tun = prs_mat %>%
  filter(`#FID` %in% y_tun$FID)

for(p_ind in 1:n.total.prs){
  prs_r2_vec_test[p_ind] = cor(y_tun$Y, prs_tun[,(2+p_ind)])^2
}

max_ind <- which.max(prs_r2_vec_test)
print(prs_r2_vec_test[max_ind])

print(colnames(prs_tun)[max_ind+2])

##########
### EB ###
##########
print('EB')
# find the best snp set
best_snps <- colnames(prs_tun)[max_ind+2]

#calculate eb effect using EB coefficients
prs_mat_eb <- CalculateEBEffectSize(bfile = 'PLINK/subtunval',
                                    snp_ind = best_snps,
                                    plink_list = plink_list,
                                    out_prefix = anc,
                                    results_dir = out_dir,
                                    params_farm = PRS_farm)
saveRDS(prs_mat_eb, paste0('CTSLEB/PRS/prsEB_', anc, '_train', train, '_rho', rho, 
                           '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS'))

dim(prs_mat_eb)

score = fread(paste0(out_dir, 'temp/score_eb_file'))
saveRDS(score$V1, paste0('CTSLEB/PRS/snpEB_', anc, '_train', train, '_rho', rho, 
                         '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS'))

# clean up files
system(paste0('rm ', out_dir, 'temp/*sscore'))
system(paste0('rm ', out_dir, 'temp/*clumped'))
system(paste0('rm ', out_dir, 'temp/*log'))

print('DONE!')
q('no')
