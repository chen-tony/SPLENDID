suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
}) 

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]]
GA = params[[4]]
rep = as.numeric(params[[5]])
model = params[[6]]
# train='1a'; rho=1; gr=1; GA=1; rep=1; model='hm3'

cov = fread('PROSPER/Cov/ALL_cov_val.txt')

dir = paste0('PROSPER/Analysis/train_', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '/')

validation_id = fread('ID/validation.id', header=F)$V1

pheno_val = rbind(
  fread(paste0('PROSPER/Pheno/AFR_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', 
               model, '_rep', rep, '.txt')),
  fread(paste0('PROSPER/Pheno/AMR_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', 
               model, '_rep', rep, '.txt')),
  fread(paste0('PROSPER/Pheno/EAS_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', 
               model, '_rep', rep, '.txt')),
  fread(paste0('PROSPER/Pheno/EUR_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', 
               model, '_rep', rep, '.txt')),
  fread(paste0('PROSPER/Pheno/SAS_val_pheno_rho', rho, '_gr', gr, '_GA', GA, '_', 
               model, '_rep', rep, '.txt'))
) %>%
  select(FID = V1, IID=V2, res=V3)

# beta
afr_beta = fread(paste0(dir, 'PROSPER/after_ensemble_AFR/PROSPER_prs_file.txt'))
amr_beta = fread(paste0(dir, 'PROSPER/after_ensemble_AMR/PROSPER_prs_file.txt'))
eas_beta = fread(paste0(dir, 'PROSPER/after_ensemble_EAS/PROSPER_prs_file.txt'))
eur_beta = fread(paste0(dir, 'PROSPER/after_ensemble_EUR/PROSPER_prs_file.txt'))
sas_beta = fread(paste0(dir, 'PROSPER/after_ensemble_SAS/PROSPER_prs_file.txt'))
all_beta = fread(paste0(dir, 'PROSPER/after_ensemble_ALL/PROSPER_prs_file.txt'))

saveRDS(list(AFR=afr_beta, AMR=amr_beta, EAS=eas_beta, EUR=eur_beta, SAS=sas_beta, ALL=all_beta),
        paste0('PROSPER/Results/beta_train_', train,
               '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

snp = unique(c(afr_beta$rsid, amr_beta$rsid, eas_beta$rsid, eur_beta$rsid, sas_beta$rsid, all_beta$rsid))
saveRDS(snp, paste0('PROSPER/Results/snps_train_', train, '_rho', rho,
                    '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

rm(afr_beta, amr_beta, eas_beta, eur_beta, sas_beta, all_beta)

# PROSPER scores
prosper_val = rbind(
  fread(paste0(dir, '/PROSPER/tmp/sample_scores_AFR/after_ensemble_testing.txt')),
  fread(paste0(dir, '/PROSPER/tmp/sample_scores_AMR/after_ensemble_testing.txt')),
  fread(paste0(dir, '/PROSPER/tmp/sample_scores_EAS/after_ensemble_testing.txt')),
  fread(paste0(dir, '/PROSPER/tmp/sample_scores_EUR/after_ensemble_testing.txt')),
  fread(paste0(dir, '/PROSPER/tmp/sample_scores_SAS/after_ensemble_testing.txt'))
) %>%
  rename(PROSPER=ensemble_score) %>%
  left_join(
    fread(paste0(dir, '/PROSPER/tmp/sample_scores_ALL/after_ensemble_testing.txt'),
          col.names=c('V1', 'V2', 'PROSPER2'))
  ) %>%
  select(FID=V1, IID=V2, PROSPER, PROSPER2)

data = pheno_val %>%
  left_join(cov) %>%
  left_join(prosper_val) %>%
  mutate(Ancestry = substr(FID, 1, 3))

data$PROSPER_res = residuals(lm(PROSPER ~ ., data=data %>% select(PROSPER, paste0('PC', 1:20))))
data$PROSPER2_res = residuals(lm(PROSPER2 ~ ., data=data %>% select(PROSPER2, paste0('PC', 1:20))))

prosper_result = data %>%
  bind_rows(., data %>% mutate(Ancestry = 'ALL')) %>%
  bind_rows(., data %>% filter(Ancestry != 'EUR') %>% mutate(Ancestry = 'nonEUR')) %>%
  mutate(Ancestry = factor(Ancestry, 
                           levels=c('ALL', 'EUR', 'nonEUR', 'AFR', 'EAS', 'SAS', 'AMR'))) %>%
  group_by(Ancestry) %>%
  summarize(n=n(),
            R2_PROSPER = cor(PROSPER, res)^2,
            R2_PROSPER2 = cor(PROSPER2, res)^2,
            R2_PROSPER_res = cor(PROSPER_res, res)^2,
            R2_PROSPER2_res = cor(PROSPER2_res, res)^2)

prosper_result

saveRDS(data, paste0('PROSPER/Results/val_train_', train, 
                     '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
saveRDS(prosper_result, paste0('PROSPER/Results/results_train_', train, 
                               '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

# clear up everything
unlink(dir, recursive=T)

