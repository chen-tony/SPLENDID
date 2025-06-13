suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(Matrix)
  library(glmnet)
})

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = as.numeric(params[[2]])
gr = as.numeric(params[[3]])
GA = as.numeric(params[[4]])
rep = as.numeric(params[[5]])
model = params[[6]]


map = readRDS(paste0('PLINK/sim_train', train, '_tun_map.rds'))

##########
## data ##
##########
# data
tuning_id = fread('ID/subtuning.id', header=F, col.names=c('FID', 'IID'))
validation_id = fread('ID/validation.id', header=F, col.names=c('FID', 'IID'))
tunval_id = fread('ID/subtunval.id', header=F, col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) %>%
  select(FID, IID, Y=paste0('PHENO', rep)) 

data = pheno %>%
  left_join(pc) %>%
  filter(FID %in% c(tuning_id$FID, validation_id$FID))

fit_tun = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% tuning_id$FID))
fit_val = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% validation_id$FID))

data$res = 0
data$res[data$FID %in% tuning_id$FID] = residuals(fit_tun)
data$res[data$FID %in% validation_id$FID] = residuals(fit_val)

rm(pheno); #gc()

tun_data = data %>%
  filter(FID %in% tuning_id$FID)

val_data = data %>%
  filter(FID %in% validation_id$FID)

rm(data); #gc()

#############################
## clean beta / grid files ##
#############################
print('beta')
beta = NULL
grid = NULL

if (file.exists(paste0('Lasso/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
  # read in
  model0 = readRDS(paste0('Lasso/OUTPUT/model_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  beta0 = model0$beta
  rownames(beta0) = model0$snps
  
  grid0 = readRDS(paste0('Lasso/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  # aggregate
  beta = cbind(beta, beta0)
  
  grid = rbind(grid, data.frame(grid0))
  
  rm(beta0, grid0, model0); #gc()
} else {
  quit('no')
}

ix_keep = which(grid$nb_active > 0)

beta_main_col = ix_keep
beta_main_row = which(rowSums(beta[,beta_main_col] != 0) > 0)
beta_main = beta[beta_main_row, beta_main_col]

beta_main_map = map %>%
  filter(marker.ID %in% rownames(beta_main)) %>%
  select(chromosome, marker.ID, allele1, allele2) %>%
  cbind(data.matrix(beta_main))
write.table(beta_main_map, 
            paste0('Lasso/PRS/beta_add_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
            row.names=F, quote=F)
write.table(beta_main_map$marker.ID, 
            paste0('Lasso/PRS/snps_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
            row.names=F, quote=F, col.names=F)

#################
## compute PRS ##
#################
# compute PRS
print('plink')
system(paste0('~/software/plink2 ',
              '--silent ',
              '--bfile ALL/all_chr ',
              '--keep ', home, 'ID/subtunval.id ',
              '--extract Lasso/PRS/snps_train', train, '_rho', rho,
              '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt ',
              '--score Lasso/PRS/beta_add_train', train, '_rho', rho,
              '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 2 3 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 5-', ncol(beta_main_map), ' ',
              '--out Lasso/PRS/prs_add_train', train, '_rho', rho,
              '_gr', gr, '_GA', GA, '_rep', rep, '_', model))

# read PRS
prs_add = fread(paste0('Lasso/PRS/prs_add_train', train, '_rho', rho, 
                       '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                col.names=c('FID', 'IID', paste0('ADD_', which(grid$nb_active > 0)))) %>%
  data.frame()

prs_tun = prs_add %>% filter(FID %in% tuning_id$FID) 
prs_val = prs_add %>% filter(FID %in% validation_id$FID)

#################
## grid search ##
#################
print('val')
# R2
grid$tun_r2 = 0
grid$tun_r2[which(grid$nb_active > 0)] = apply(prs_tun %>% select(-FID, -IID), 2, FUN=function(x) 
  ifelse(sd(x)>0, cor(tun_data$res, x)^2, 0))

grid$val_r2 = 0
grid$val_r2[which(grid$nb_active > 0)] = apply(prs_val %>% select(-FID, -IID), 2, FUN=function(x) 
  ifelse(sd(x)>0, cor(val_data$res, x)^2, 0))

grid %>% 
  mutate(ix=row_number()) %>% 
  slice_max(tun_r2)

ix_add = grid %>% 
  mutate(ix=row_number()) %>% 
  # filter(nb_active < 30000) %>%
  slice_max(tun_r2) %>% 
  slice_min(ix) %>%
  pull(ix)

val_data$PRS_add = prs_val[,paste0('ADD_', ix_add)]

saveRDS(grid, paste0('Lasso/RESULTS/grid_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

ix_ipgs = which.max(grid %>% filter(!is.na(tun_r2) & tun_r2 > 0) %>% pull(tun_r2))
beta_out = beta_main_map[,c(1:4, ix_ipgs + 4)] 
names(beta_out)[5] = 'BETA'
beta_out = beta_out %>% filter(BETA != 0)
write.table(beta_out$marker.ID, 
            paste0('Lasso/RESULTS/beta_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp'),
            row.names=F, quote=F, col.names=F)
write.table(beta_out, 
            paste0('Lasso/RESULTS/beta_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
            row.names=F, quote=F)

#############
## results ##
#############
print('results')

val_data$PRS_add_res = residuals(lm(PRS_add ~ ., data=val_data %>% select(PRS_add, paste0('PC', 1:20))))

r2_results = val_data %>%
  mutate(ancestry = substr(FID, 1, 3)) %>%
  bind_rows(., val_data %>% mutate(ancestry = 'ALL')) %>%
  bind_rows(., val_data %>% mutate(ancestry = substr(FID, 1, 3)) %>%
              filter(ancestry != 'EUR') %>% mutate(ancestry = 'nonEUR')) %>%
  mutate(ancestry = factor(ancestry, 
                           levels=c('ALL', 'EUR', 'nonEUR', 'AFR', 'EAS', 'SAS', 'AMR'))) %>%
  group_by(ancestry) %>%
  summarize(n=n(),
            R2_add = cor(PRS_add, res)^2,
            R2_add_res = cor(PRS_add_res, res)^2) %>%
  ungroup() %>%
  mutate(train=train,
         rho=rho,
         gr=gr,
         GA=GA,
         rep=rep,
         model=model,
         supp_add = grid %>% slice_max(tun_r2) %>% pull(nb_active) %>% min)

r2_results

saveRDS(val_data, paste0('Lasso/RESULTS/Lasso_val_train', train, '_rho', rho,
                         '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

saveRDS(r2_results, paste0('Lasso/RESULTS/val_train', train, '_rho', rho, 
                           '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

file.remove(paste0('Lasso/PRS/beta_add_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
file.remove(paste0('Lasso/PRS/snps_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
file.remove(paste0('Lasso/PRS/prs_add_train', train, '_rho', rho, 
                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'))
