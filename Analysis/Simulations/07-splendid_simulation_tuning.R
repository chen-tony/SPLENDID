suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(Matrix)
  library(glmnet)
  library(SuperLearner)
})

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = as.numeric(params[[2]])
gr = as.numeric(params[[3]])
GA = as.numeric(params[[4]])
model = params[[5]]

map = readRDS(paste0('PLINK/sim_train', train, '_map.rds'))

##########
## data ##
##########
# data
# train_id = fread(paste0('ID/train', train, '.id'), header=F, col.names=c('FID', 'IID'))
tuning_id = fread(paste0('ID/subtuning.id'), header=F, col.names=c('FID', 'IID'))
validation_id = fread(paste0('ID/validation.id'), header=F, col.names=c('FID', 'IID'))
tunval_id = fread(paste0('ID/subtunval.id'), header=F, col.names=c('FID', 'IID'))

pc = fread(paste0('PCA/all_sim_pca.sscore'),
           col.names=c('FID', 'IID', paste0('PC', 1:20)))


pheno_full = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) 

for (rep in 1:10) {
  if (file.exists(paste0('L0L1/RESULTS/val_train', train, '_rho', rho, 
                         '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) next
  
  if (!file.exists(paste0('L0L1/OUTPUT/model_train', train, '_rho', rho, 
                          '_gr', gr, '_GA', GA, '_rep', rep, '_', model, 
                          '_pval0_lambda1.RDS'))) next
  
  print(rep)
  
  cat('data \n')
  
  pheno = pheno_full %>%
    select(FID, IID, Y=paste0('PHENO', rep)) 
  
  data = pheno %>%
    left_join(pc) %>%
    filter(FID %in% c(tuning_id$FID, validation_id$FID))
  data$res = residuals(lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), data))
  rm(pheno); gc()
  
  fit_tun = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% tuning_id$FID))
  fit_val = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% validation_id$FID))
  
  data$res = 0
  data$res[data$FID %in% tuning_id$FID] = residuals(fit_tun)
  data$res[data$FID %in% validation_id$FID] = residuals(fit_val)
  
  tun_data = data %>%
    filter(FID %in% tuning_id$FID)
  
  val_data = data %>%
    filter(FID %in% validation_id$FID)
  
  rm(data); gc()
  
  #############################
  ## clean beta / grid files ##
  #############################
  cat('beta \n')
  beta = NULL
  grid = NULL
  for (pval in c('5e-6', '5e-8', '0')) {
    for (lambda_ix in 1:5) {
      if (file.exists(paste0('L0L1/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))) {
        # read in
        model0 = readRDS(paste0('L0L1/OUTPUT/model_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))
        beta0 = model0$beta
        rownames(beta0) = model0$snps
        
        grid0 = readRDS(paste0('L0L1/OUTPUT/results_train', train, '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '_pval', pval, '_lambda', lambda_ix, '.RDS'))
        
        # aggregate
        beta = cbind(beta, beta0)
        
        grid = rbind(grid, data.frame(grid0, pval=pval, lambda_ix=lambda_ix))
        
        rm(beta0, grid0, model0); gc()
      }
    }
  }
  
  ix_keep = which(grid$nb_active > 0)
  ix_int = which(grid$nb_interact > 0)
  
  beta_main_col = sapply(ix_keep, FUN=function(i) (i-1)*6 + 1)
  beta_main_row = which(rowSums(beta[,beta_main_col] != 0) > 0)
  beta_main = beta[beta_main_row, beta_main_col]
  colnames(beta_main) = paste0('ADD_', which(grid$nb_active > 0))
  
  beta_main_map = map %>%
    filter(marker.ID %in% rownames(beta_main)) %>%
    select(chromosome, marker.ID, allele1, allele2) %>%
    cbind(data.matrix(beta_main))
  write.table(beta_main_map,
              paste0('L0L1/PRS/beta_add_train', train, '_rho', rho,
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
              row.names=F, quote=F)
  write.table(beta_main_map$marker.ID,
              paste0('L0L1/PRS/snps_train', train, '_rho', rho,
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
              row.names=F, quote=F, col.names=F)
  
  if (length(ix_int) > 0) {
    beta_int_col = sort(c(sapply(ix_int, FUN=function(i) (i-1)*6 + 2:6)))
    beta_int_row = which(rowSums(beta[,beta_int_col] != 0) > 0)
    
    beta_int = beta[beta_int_row, beta_int_col]
    
    if (length(beta_int_row) == 1) {
      beta_int_map = map %>%
        filter(marker.ID %in% rownames(beta)[beta_int_row]) %>%
        select(chromosome, marker.ID, allele1, allele2) %>%
        cbind(matrix(beta_int, nrow=1))
    } else {
      beta_int_map = map %>%
        filter(marker.ID %in% rownames(beta_int)) %>%
        select(chromosome, marker.ID, allele1, allele2) %>%
        cbind(data.matrix(beta_int))
    }
    
    write.table(beta_int_map,
                paste0('L0L1/PRS/beta_int_train', train, '_rho', rho,
                       '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
                row.names=F, quote=F)
    write.table(beta_int_map$marker.ID,
                paste0('L0L1/PRS/ints_train', train, '_rho', rho,
                       '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
                row.names=F, quote=F)
  }
  
  #################
  ## compute PRS ##
  #################
  cat('plink \n')
  # compute PRS
  system(paste0('~/software/plink2 ',
                '--silent ',
                '--bfile all_chr ',
                '--keep ID/subtunval.id ',
                '--extract L0L1/PRS/snps_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt ',
                '--score L0L1/PRS/beta_add_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 2 3 header ',
                'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                '--score-col-nums 5-', ncol(beta_main_map), ' ',
                '--out L0L1/PRS/prs_add_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model))
  
  if (length(ix_int) > 0) {
    system(paste0('~/software/plink2 ',
                  '--silent ',
                  '--bfile all_chr ',
                  '--keep ID/subtunval.id ',
                  '--extract L0L1/PRS/ints_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt ',
                  '--score L0L1/PRS/beta_int_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 2 3 header ',
                  'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                  '--score-col-nums 5-', ncol(beta_int_map), ' ',
                  '--out L0L1/PRS/prs_int_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model))
  }
  
  # read PRS
  cat('prs \n')
  prs = fread(paste0('L0L1/PRS/prs_add_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
              col.names=c('FID', 'IID', paste0('ADD_', which(grid$nb_active > 0)))) %>%
    data.frame()
  
  if (length(ix_int) > 0) {
    prs_int = fread(paste0('L0L1/PRS/prs_int_train', train, '_rho', rho, 
                           '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                    col.names = c('FID', 'IID', c(outer(1:5,which(grid$nb_interact > 0),  
                                                        FUN=function(k,i) paste0('INT', k, '_', i))))) %>%
      data.frame()
    
    # integrate PCs
    Esum = readRDS(paste0('SUMMARY/Esum_train', train, '.RDS'))
    pc_tunval = pc %>%
      filter(FID %in% tunval_id$FID) %>%
      select(paste0('PC', 1:5)) %>%
      data.frame()
    for (k in 1:5) pc_tunval[,k] = (pc_tunval[,k] - Esum$mean[k]) / Esum$sd[k]
    
    for (i in which(grid$nb_interact > 0)) {
      prs[,paste0('ADD_', i)] = prs[,paste0('ADD_', i)] + prs_int[,paste0('INT1_', i)] * pc_tunval[,1] +
        prs_int[,paste0('INT2_', i)] * pc_tunval[,2] + prs_int[,paste0('INT3_', i)] * pc_tunval[,3] +
        prs_int[,paste0('INT4_', i)] * pc_tunval[,4] + prs_int[,paste0('INT5_', i)] * pc_tunval[,5]
    }
  }
  
  prs_tun = prs %>% filter(FID %in% tuning_id$FID)
  prs_val = prs %>% filter(FID %in% validation_id$FID)
  
  #################
  ## grid search ##
  #################
  cat('grid \n')
  # R2
  grid$tun_r2 = 0
  grid$tun_r2[which(grid$nb_active > 0)] = apply(prs_tun %>% select(-FID, -IID), 2, FUN=function(x) 
    ifelse(sd(x)>0, cor(tun_data$res, x)^2, 0))
  
  grid$val_r2 = 0
  grid$val_r2[which(grid$nb_active > 0)] = apply(prs_val %>% select(-FID, -IID), 2, FUN=function(x) 
    ifelse(sd(x)>0, cor(val_data$res, x)^2, 0))
  
  grid %>%
    group_by(pval) %>%
    slice_max(tun_r2) %>%
    slice_head(n=1) %>%
    data.frame() %>% 
    print()
  
  ix_add = grid %>% 
    mutate(ix=row_number()) %>% 
    filter(pval == '0') %>%
    slice_max(tun_r2) %>% 
    slice_min(ix) %>%
    pull(ix)
  
  ix_int = grid %>% 
    mutate(ix=row_number()) %>% 
    filter(pval != '0') %>%
    slice_max(tun_r2) %>% 
    slice_min(ix) %>%
    pull(ix)
  
  ix_grid = grid %>% 
    mutate(ix=row_number()) %>% 
    slice_max(tun_r2) %>% 
    slice_min(ix) %>%
    pull(ix)
  
  val_data$PRS_add = prs_val[,paste0('ADD_', ix_add)]
  
  if (length(ix_int) > 0) {
    val_data$PRS_int = prs_val[,paste0('ADD_', ix_int)]
  } else {
    val_data$PRS_int = val_data$PRS_add
  }
  
  val_data$PRS_grid = prs_val[,paste0('ADD_', ix_grid)]
  
  saveRDS(grid, paste0('L0L1/RESULTS/grid_train', train, '_rho', rho,
                       '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  ################
  ## ensembling ##
  ################
  cat('ensemble \n')
  # pooled-ancestry tuning
  grid0 = grid %>%
    mutate(tun_r2 = ifelse(nb_active > min(100000, 2*grid$nb_active[which.max(grid$tun_r2)]),
                           tun_r2, NA))
  ix_grid_order = intersect(order(-grid$tun_r2), which(grid$nb_active > 0))
  beta_grid = beta[beta_main_row,(ix_grid_order[1]-1)*6 + 1]
  ix_grid_keep = ix_grid_order[1]
  for (i in 2:length(ix_grid_order)) {
    if (is.na(grid0$tun_r2[ix_grid_order[i]] > 0)) {
      beta_grid0 = cbind(beta_grid, beta[beta_main_row,(ix_grid_order[i]-1)*6 + 1])
      if (sum(rowSums(beta_grid0 != 0) > 0) < min(100000, 2*grid$nb_active[which.max(grid$tun_r2)])) {
        beta_grid = beta_grid0
        ix_grid_keep = c(ix_grid_keep, ix_grid_order[i])
      }
    }
  }
  ix_grid_keep = sort(ix_grid_keep)
  
  cat(length(ix_grid_keep), 'PRS kept \n')
  # print(ix_grid %in% ix_grid_keep)
  
  prs_tun_ens = prs_tun %>% select(paste0('ADD_', ix_grid_keep)) %>% data.matrix()
  prs_val_ens = prs_val %>% select(paste0('ADD_', ix_grid_keep)) %>% data.matrix()
  
  cat(' cv \n')
  set.seed(123)
  cv = cv.glmnet(prs_tun_ens, tun_data$res, alpha=0)
  val_data$PRS_ens = as.numeric(predict(cv, prs_val_ens, s='lambda.min'))
  
  # ancestry-specific tuning
  set.seed(123)
  val_data$PRS_ancens = 0
  for (pop in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    ix_pop_tun = which(substr(tun_data$FID, 1, 3) == pop)
    ix_pop_val = which(substr(val_data$FID, 1, 3) == pop)
    cv3 = cv.glmnet(prs_tun[ix_pop_tun,] %>% select(contains('ADD_')) %>% data.matrix(), tun_data$res[ix_pop_tun], alpha=0)
    
    val_data$PRS_ancens[ix_pop_val] = as.numeric(predict(cv3, prs_val[ix_pop_val,] %>% select(contains('ADD_')) %>% data.matrix(), s='lambda.min'))
  }
  
  cat(' pred \n')
  ### print SNPs
  ix_grid_all_add = (ix_grid_keep-1)*6 + 1
  ix_grid_all_int = (ix_grid_keep-1)*6 + 2
  
  ### beta
  cat(' beta \n')
  # grid
  beta_grid = beta[,(ix_grid-1)*6 + 1:6]
  colnames(beta_grid) = c('ADD', paste0('INT', 1:5))
  ix_beta_grid = which(rowSums(beta_grid != 0) > 0)
  beta_grid_map = map %>%
    filter(marker.ID %in% rownames(beta_grid)[ix_beta_grid]) %>%
    select(chromosome, marker.ID, allele1, allele2) %>%
    cbind(data.matrix(beta_grid[ix_beta_grid,]))
  # print(all.equal(beta_grid_map$marker.ID, rownames(beta_grid_map)))
  rownames(beta_grid_map) = NULL
  
  # print(dim(beta_grid_map))
  # print(summary(beta_grid_map))
  # print(head(beta_grid_map))
  
  # sparse ensemble
  ix_ens_row = which(rowSums(beta[,ix_grid_all_add] != 0) > 0)
  ix_ens_col= sort(c(sapply(ix_grid_keep, FUN=function(i) (i-1)*6 + 1:6)))
  
  coef_cv = coef(cv, s='lambda.min')[-1]
  
  beta_ens = matrix(nrow=length(ix_ens_row), ncol=6)
  rownames(beta_ens) = rownames(beta)[ix_ens_row]
  colnames(beta_ens) = c('ADD', paste0('INT', 1:5))
  for (k in 1:6)  {
    beta_ens[,k] = data.matrix(beta[ix_ens_row,(ix_grid_keep-1)*6 + k]) %*% coef_cv
  }
  
  beta_ens_map = map %>%
    filter(marker.ID %in% rownames(beta_ens)) %>%
    select(chromosome, marker.ID, allele1, allele2) %>%
    cbind(beta_ens)
  # print(all.equal(beta_ens_map$marker.ID, rownames(beta_ens_map)))
  rownames(beta_ens_map) = NULL
  
  # print(dim(beta_ens_map))
  # print(summary(beta_ens_map))
  # print(head(beta_ens_map))
  
  supp_ens = sum(beta_ens[,1] != 0)
  het_ens = max(colSums(beta_ens[,2:6] != 0))
  
  saveRDS(list(grid=beta_grid_map, ens=beta_ens_map),
          paste0('L0L1/RESULTS/L0L1_beta_train', train, '_rho', rho,
                 '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  #############
  ## results ##
  #############
  cat('results \n')
  supp_add = grid %>% filter(pval == '0') %>% slice_max(tun_r2) %>% pull(nb_active) %>% head(1)
  if (length(unique(grid$pval)) > 1) {
    supp_int = grid %>% filter(pval != '0') %>% slice_max(tun_r2) %>% pull(nb_active) %>% head(1)
    het_int = grid %>% filter(pval != '0') %>% slice_max(tun_r2) %>% pull(nb_interact) %>% head(1)
  } else {
    supp_int = supp_add
    het_int = 0
  }
  supp_grid = grid %>% slice_max(tun_r2) %>% pull(nb_active) %>% head(1)
  het_grid = grid %>% slice_max(tun_r2) %>% pull(nb_interact) %>% head(1)
  
  val_data$PRS_add_res = residuals(lm(PRS_add ~ ., data=val_data %>% select(PRS_add, paste0('PC', 1:20))))
  val_data$PRS_int_res = residuals(lm(PRS_int ~ ., data=val_data %>% select(PRS_int, paste0('PC', 1:20))))
  val_data$PRS_grid_res = residuals(lm(PRS_grid ~ ., data=val_data %>% select(PRS_grid, paste0('PC', 1:20))))
  val_data$PRS_ens_res = residuals(lm(PRS_ens ~ ., data=val_data %>% select(PRS_ens, paste0('PC', 1:20))))
  val_data$PRS_ancens_res = residuals(lm(PRS_ancens ~ ., data=val_data %>% select(PRS_ancens, paste0('PC', 1:20))))
  
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
              R2_int = cor(PRS_int, res)^2,
              R2_grid = cor(PRS_grid, res)^2,
              R2_ens = cor(PRS_ens, res)^2,
              R2_ancens = cor(PRS_ancens, res)^2,
              R2_add_res = cor(PRS_add_res, res)^2,
              R2_int_res = cor(PRS_int_res, res)^2,
              R2_grid_res = cor(PRS_grid_res, res)^2,
              R2_ens_res = cor(PRS_ens_res, res)^2,
              R2_ancens_res = cor(PRS_ancens_res, res)^2) %>%
    ungroup() %>%
    mutate(train=train,
           rho=rho,
           gr=gr,
           GA=GA,
           rep=rep,
           model=model,
           supp_add = supp_add,
           supp_int = supp_int,
           supp_grid = supp_grid,
           supp_ens = supp_ens,
           het_int = het_int,
           het_grid = het_grid,
           het_ens = het_ens)
  
  # print(data.frame(r2_results))
  print(data.frame(r2_results) %>% select(-contains('_add'), -contains('_int'), -contains('_grid'), -contains('_anc')))
  
  saveRDS(r2_results, paste0('L0L1/RESULTS/val_train', train, '_rho', rho, 
                             '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  saveRDS(val_data, paste0('L0L1/RESULTS/L0L1_val_train', train, '_rho', rho,
                           '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  
  rm(beta, beta_main, beta_int, beta_grid, beta_ens,
     pheno, data, val_data, tun_data,
     prs, prs_int, prs_tun_ens, prs_val_ens, pc_tunval, cv, r2_results)
  gc()
  
  file.remove(paste0('L0L1/PRS/beta_add_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
  file.remove(paste0('L0L1/PRS/beta_int_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
  file.remove(paste0('L0L1/PRS/snps_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
  file.remove(paste0('L0L1/PRS/ints_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'))
  
  file.remove(paste0('L0L1/PRS/prs_add_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'))
  file.remove(paste0('L0L1/PRS/prs_add_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.log'))
  
  file.remove(paste0('L0L1/PRS/prs_int_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'))
  file.remove(paste0('L0L1/PRS/prs_int_train', train, '_rho', rho, 
                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.log'))
}

print('done!')





