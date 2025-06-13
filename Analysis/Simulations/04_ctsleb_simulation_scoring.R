# CT-SLEB: https://github.com/andrewhaoyu/CTSLEB

# run SuperLearning between all 5 populations

##################
## PREPARE DATA ##
##################
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(CTSLEB)
  library(SuperLearner)
  library(glmnet)
})

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]]
GA = params[[4]]
rep = as.numeric(params[[5]])
model = params[[6]]

# SUPER LEARNING
tun_id = fread('ID/subtuning.id', header=F,
               col.names=c('FID', 'IID'))
val_id = fread('ID/validation.id', header=F,
               col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

cov = fread('PROSPER/Cov/ALL_cov_val.txt')

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) %>%
  select(FID, IID, Y=paste0('PHENO', rep)) 

tunval_data = data %>%
  filter(FID %in% c(tun_id$FID, val_id$FID))

# compute residuals for tuning
fit_tun = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% tun_id$FID))
fit_val = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% val_id$FID))

tunval_data$res = 0
tunval_data$res[tunval_data$FID %in% tuning_id$FID] = residuals(fit_tun)
tunval_data$res[tunval_data$FID %in% validation_id$FID] = residuals(fit_val)

val_data = tunval_data %>%
  filter(FID %in% val_id$FID) %>%
  select(FID, IID, Y=res) %>%
  mutate(ancestry = substr(FID, 1 ,3))
val_data$CTSLEB = 0

prs_tunval = cbind(
  readRDS(paste0('CTSLEB/PRS/prsEB_AFR_train', train, '_rho', rho, 
                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')),
  readRDS(paste0('CTSLEB/PRS/prsEB_AMR_train', train, '_rho', rho, 
                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')) %>%
    select(-contains('ID')),
  readRDS(paste0('CTSLEB/PRS/prsEB_EAS_train', train, '_rho', rho, 
                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')) %>%
    select(-contains('ID')),
  readRDS(paste0('CTSLEB/PRS/prsEB_SAS_train', train, '_rho', rho, 
                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')) %>%
    select(-contains('ID'))
)
names(prs_tunval)[-c(1:2)] = paste0('PRS', 1:(ncol(prs_tunval) - 2))

####################
## SUPER LEARNING ##
####################
SL.library <- c(
  "SL.glmnet",
  "SL.ridge"
)

# ancestry-specific SL
for (anc in c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')) {
  cat(anc, ': ')
  set.seed(123)
  
  y_tun = tunval_data %>%
    select(FID, IID, Y=res) %>%
    filter(FID %in% tun_id$FID) %>%
    filter(substr(FID, 1, 3) == anc)
  
  y_val = tunval_data %>%
    select(FID, IID, Y=res) %>%
    filter(FID %in% val_id$FID) %>%
    filter(substr(FID, 1, 3) == anc)
  
  prs_tun = prs_tunval %>%
    filter(IID %in% y_tun$IID) %>%
    select(-contains('ID')) %>%
    data.matrix()
  
  prs_val = prs_tunval %>%
    filter(IID %in% y_val$IID) %>%
    select(-contains('ID')) %>%
    data.matrix()
  
  ix_keep = which(colSums(prs_tun != 0) > 0)
  
  prs_tun = prs_tun[,ix_keep]
  prs_val = prs_val[,ix_keep]
  cat(ncol(prs_tun), ', ')
  
  prs_cor = cor(prs_tun)
  cat('COR, ')
  
  ix_cor = caret::findCorrelation(prs_cor, 0.98)
  # print(length(ix_cor))
  
  prs_tun_sl = data.frame(prs_tun[,-ix_cor])
  cat(ncol(prs_tun_sl), ', ')
  
  prs_val_sl = data.frame(prs_val[,-ix_cor])
  rm(prs_tun, prs_val, prs_cor); gc()
  
  cat(' SL')
  sl <- SuperLearner(Y = as.numeric(y_tun$Y),
                     X = prs_tun_sl,
                     family = gaussian(),
                     SL.library = SL.library)
  pred_tun <- predict(sl, prs_tun_sl, onlySL = TRUE)[[1]]
  pred_val <- predict(sl, prs_val_sl, onlySL = TRUE)[[1]]
  
  val_data$CTSLEB[val_data$ancestry==anc] = pred_val
  # cor(pred_val, y_val$Y)^2
  
  cat('\n')
  
  rm(prs_tun_sl, prs_val_sl, 
     sl, pred_tun, pred_val); gc()
}

# altogether SL
{
  print('ALL')
  set.seed(123)
  
  y_tun = tunval_data %>%
    select(FID, IID, Y=res) %>%
    filter(FID %in% tun_id$FID)
  
  y_val = tunval_data %>%
    select(FID, IID, Y=res) %>%
    filter(FID %in% val_id$FID)
  
  prs_tun = prs_tunval %>%
    filter(IID %in% y_tun$IID) %>%
    select(-contains('ID')) %>%
    data.matrix()
  
  prs_val = prs_tunval %>%
    filter(IID %in% y_val$IID) %>%
    select(-contains('ID')) %>%
    data.matrix()
  
  ix_keep = which(colSums(prs_tun != 0) > 0)
  
  prs_tun = prs_tun[,ix_keep]
  prs_val = prs_val[,ix_keep]
  
  prs_cor = cor(prs_tun)
  ix_cor = caret::findCorrelation(prs_cor, 0.98)
  print(length(ix_cor))
  
  prs_tun_sl = data.frame(prs_tun[,-ix_cor])
  prs_val_sl = data.frame(prs_val[,-ix_cor])
  rm(prs_tun, prs_val, prs_cor); gc()
  
  print(' SL')
  sl <- SuperLearner(Y = as.numeric(y_tun$Y),
                     X = prs_tun_sl,
                     family = gaussian(),
                     SL.library = SL.library)
  pred_tun <- predict(sl, prs_tun_sl, onlySL = TRUE)[[1]]
  pred_val <- predict(sl, prs_val_sl, onlySL = TRUE)[[1]]
  
  val_data$CTSLEB2 = pred_val
  
  rm(prs_tun_sl, prs_val_sl, 
     sl, pred_tun, pred_val); gc()
}


val_data = val_data %>% left_join(cov)
val_data$CTSLEB_res = residuals(lm(CTSLEB ~ ., data=val_data %>% select(CTSLEB, paste0('PC', 1:20))))
val_data$CTSLEB2_res = residuals(lm(CTSLEB2 ~ ., data=val_data %>% select(CTSLEB2, paste0('PC', 1:20))))


saveRDS(val_data, paste0('CTSLEB/Results/val_train_', train, 
                         '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))


ctsleb_result = val_data %>%
  bind_rows(., val_data %>% mutate(ancestry = 'ALL')) %>%
  bind_rows(., val_data %>% filter(ancestry != 'EUR') %>% mutate(ancestry = 'nonEUR')) %>%
  mutate(Ancestry = factor(ancestry, 
                           levels=c('ALL', 'EUR', 'nonEUR', 'AFR', 'EAS', 'SAS', 'AMR'))) %>%
  group_by(Ancestry) %>%
  summarize(n=n(),
            R2_CTSLEB = cor(CTSLEB, Y)^2,
            R2_CTSLEB2 = cor(CTSLEB2, Y)^2,
            R2_CTSLEB_res = cor(CTSLEB_res, Y)^2,
            R2_CTSLEB2_res = cor(CTSLEB2_res, Y)^2)

saveRDS(ctsleb_result, paste0('CTSLEB/Results/results_train_', train, 
                              '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))



print('DONE')
