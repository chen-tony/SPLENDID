##############
### SERVER ###
##############
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(Matrix)
  library(glmnet)
})

r2_results = data.frame()
snp_results = data.frame()

for (train in c('2a', '4a')) {
  for (model in c('aou')) {
    print(c(train, model))
    for (rho in c(1, 1.5, 2, 3, 4)) {
      for (GA in 1:2) {
        for (gr in c(0.4, 0.6, 0.8)) {
          # print(c(train, model, GA, gr))
          for (rep in 1:10) {
            # STELLAR
            if (file.exists(paste0('L0L1/RESULTS/val_train', train, '_rho', rho, 
                                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              results0 = readRDS(paste0('L0L1/RESULTS/val_train', train, '_rho', rho, 
                                        '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              
              r2_result0 = results0 %>% 
                select(train, rho, gr, GA, rep, ancestry,
                       R2_grid, R2_ens, R2_ancens,
                       R2_grid_res, R2_ens_res, R2_ancens_res) %>%
                gather(key=Model, value=R2, 
                       R2_grid, R2_ens, R2_ancens,
                       R2_grid_res, R2_ens_res, R2_ancens_res) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = case_when(
                  grepl('R2_grid', Model) ~ 'STELLAR Grid',
                  grepl('R2_ens', Model) ~ 'STELLAR Ensemble', 
                  grepl('R2_ancens', Model) ~ 'STELLAR Ancestry Ensemble', 
                  grepl('R2_fullens', Model) ~ 'STELLAR Full Ensemble'))
              
              r2_results = rbind(r2_results, r2_result0)
              
              snp_result0 = results0 %>%
                filter(ancestry == 'ALL') %>%
                select(-contains('R2'), -ends_with('add'), -ends_with('int')) %>%
                gather(key=Model, value=Supp, 
                       supp_grid, supp_ens,
                       # het_int, 
                       het_grid, het_ens) %>%
                mutate(Component = case_when(
                  grepl('supp', Model) ~ 'Main',
                  grepl('het', Model) ~ 'Interaction'
                )) %>%
                mutate(Model = case_when(
                  grepl('_grid', Model) ~ 'STELLAR Grid',
                  grepl('_ens', Model) ~ 'STELLAR Ensemble'
                )) %>%
                select(train, rho, gr, GA, rep, model, Model, Supp, Component)
              snp_results = rbind(snp_results, snp_result0)
            }
            
            # iPGS
            if (file.exists(paste0('Lasso/RESULTS/val_train', train, '_rho', rho,
                                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              results0 = readRDS(paste0('Lasso/RESULTS/val_train', train, '_rho', rho,
                                        '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              
              r2_result0 = results0 %>%
                select(train, rho, gr, GA, rep, ancestry, R2_add, R2_add_res) %>%
                gather(key=Model, value=R2, R2_add, R2_add_res) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = 'Lasso')
              
              r2_results = rbind(r2_results, r2_result0)
              
              snp_result0 = results0 %>%
                filter(ancestry == 'ALL') %>%
                mutate(Model = 'Lasso',
                       Component = 'Main') %>%
                mutate(model=model) %>%
                select(train, rho, gr, GA, rep, model, Model, Supp=supp_add, Component)
              
              snp_results = rbind(snp_results, snp_result0)
            }
            
            # iPGS + Refit
            if (file.exists(paste0('Lasso/RESULTS/valREFIT_train', train, '_rho', rho,
                                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              results0 = readRDS(paste0('Lasso/RESULTS/valREFIT_train', train, '_rho', rho,
                                        '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              
              r2_result0 = results0 %>%
                select(train, rho, gr, GA, rep, ancestry, 
                       R2_refit, R2_refit2, R2_refit_res, R2_refit2_res) %>%
                gather(key=Model, value=R2, contains('R2')) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = case_when(
                  grepl('R2_refit2', Model) ~ 'iPGS + refit (all)',
                  grepl('R2_refit', Model) ~ 'iPGS + refit'))
              r2_results = rbind(r2_results, r2_result0)
              
              snp_result0 = results0 %>%
                filter(ancestry == 'ALL') %>%
                mutate(model=model) %>%
                select(-contains('iPGS')) %>%
                gather(key=Model, value=Supp, supp_add, supp_int, supp_add2, supp_int2) %>%
                mutate(Component = case_when(
                  grepl('add', Model) ~ 'Main',
                  grepl('int', Model) ~ 'Interaction'
                )) %>%
                mutate(Model = ifelse(grepl(2, Model), 'iPGS+refit (all)', 'iPGS+refit')) %>%
                select(-ancestry, -n) %>%
                select(train, rho, gr, GA, rep, model, Model, Supp, Component)
              
              snp_results = rbind(snp_results, snp_result0)
            } else {
              if (file.exists(paste0('Lasso/RESULTS/val_train', train, '_rho', rho,
                                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
                results0 = readRDS(paste0('Lasso/RESULTS/val_train', train, '_rho', rho,
                                          '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
                
                r2_result0 = results0 %>%
                  select(train, rho, gr, GA, rep, ancestry, R2_add, R2_add_res) %>%
                  gather(key=Model, value=R2, R2_add, R2_add_res) %>%
                  mutate(Res = grepl('_res', Model)) %>%
                  mutate(Model = 'iPGS+refit')
                
                r2_results = rbind(r2_results, 
                                   r2_result0, 
                                   r2_result0 %>% mutate(Model = 'iPGS+refit (all)'))
                
                snp_result0 = results0 %>%
                  filter(ancestry == 'ALL') %>%
                  mutate(Model = 'iPGS+refit',
                         Component = 'Main') %>%
                  mutate(model=model) %>%
                  select(train, rho, gr, GA, rep, model, Model, Supp=supp_add, Component)
                
                snp_results = rbind(snp_results, snp_result0, snp_result0 %>% mutate(Model = 'iPGS+refit (all)'))
              }
            }
            
            # PROSPER
            if (file.exists(paste0('PROSPER/Results/results_train_', train, '_rho', rho,
                                   '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              
              results0 = readRDS(paste0('PROSPER/Results/results_train_', train, '_rho', rho, 
                                        '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              
              if (file.exists(paste0('PROSPER/Results/snps_train_', train, '_rho', rho, 
                                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
                snp = readRDS(paste0('PROSPER/Results/snps_train_', train, '_rho', rho, 
                                     '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              } else {
                beta = readRDS(paste0('PROSPER/Results/beta_train_', train, '_rho', rho, 
                                      '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
                
                snp = unique(unlist(lapply(beta, FUN=function(x) x$rsid)))
                
                if (length(snp) == 0) snp = rep(0, nrow(beta[[1]]))
                
                saveRDS(snp, paste0('PROSPER/Results/snps_train_', train, '_rho', rho, 
                                    '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              }
              
              supp = length(snp)
              
              r2_result0 = results0 %>%
                mutate(train=train, rho=rho, gr=gr, GA=GA, rep=rep) %>%
                select(train, rho, gr, GA, rep, ancestry=Ancestry, 
                       R2_PROSPER, R2_PROSPER2, R2_PROSPER_res, R2_PROSPER2_res) %>%
                gather(key=Model, value=R2, contains('R2')) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = case_when(
                  grepl('PROSPER2', Model) ~ 'PROSPER (all)',
                  grepl('PROSPER', Model) ~ 'PROSPER'
                ))
              r2_results = rbind(r2_results, r2_result0)
              
              snp_result0 = data.frame(train=train, 
                                       rho=rho, 
                                       gr=gr,
                                       GA=GA,
                                       rep=rep, 
                                       model=model, 
                                       Model = 'PROSPER',
                                       Supp = supp,
                                       Component = 'Main') 
              snp_results = rbind(snp_results, snp_result0)
              
              rm(results0, beta, snp, supp); gc()
            }
            
            # PRSCSX
            if (file.exists(paste0('PRSCSX/Results/results_train_', train, 
                                   '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              results0 = readRDS(paste0('PRSCSX/Results/results_train_', train, 
                                        '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
              
              r2_result0 = results0 %>%
                mutate(train=train, rho=rho, gr=gr, GA=GA, rep=rep) %>%
                select(train, rho, gr, GA, rep, ancestry=Ancestry, 
                       R2_PRSCSX, R2_PRSCSX2, R2_PRSCSX_res, R2_PRSCSX2_res) %>%
                gather(key=Model, value=R2, contains('R2')) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = case_when(
                  grepl('PRSCSX2', Model) ~ 'PRSCSX (all)',
                  grepl('PRSCSX', Model) ~ 'PRSCSX'
                ))
              r2_results = rbind(r2_results, r2_result0)
              
              {
                snps = unique(c(
                  fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                               '_rep', rep, '_GA', GA, '_', model, '_EUR_snps.txt'), header=F)$V1,
                  fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                               '_rep', rep, '_GA', GA, '_', model, '_AFR_snps.txt'), header=F)$V1,
                  fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                               '_rep', rep, '_GA', GA, '_', model, '_EAS_snps.txt'), header=F)$V1,
                  fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                               '_rep', rep, '_GA', GA, '_', model, '_SAS_snps.txt'), header=F)$V1,
                  fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                               '_rep', rep, '_GA', GA, '_', model, '_AMR_snps.txt'), header=F)$V1
                ))
              }
              
              snp_result0 = data.frame(train=train, 
                                       rho=rho, 
                                       gr=gr,
                                       GA=GA,
                                       rep=rep, 
                                       model=model, 
                                       Model = 'PRSCSX',
                                       Supp = length(snps),
                                       Component = 'Main') 
              snp_results = rbind(snp_results, snp_result0)
              
              rm(results0, snps); gc()
            }
            
            # CTSLEB
            if (file.exists(paste0('CTSLEB/Results/results_train_', train, 
                                   '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) {
              
              results0 = readRDS(paste0('CTSLEB/Results/results_train_', train, 
                                        '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))

              
              r2_result0 = results0 %>%
                mutate(train=train, rho=rho, gr=gr, GA=GA, rep=rep) %>%
                select(train, rho, gr, GA, rep, ancestry=Ancestry, 
                       R2_CTSLEB, R2_CTSLEB2, R2_CTSLEB_res, R2_CTSLEB2_res) %>%
                gather(key=Model, value=R2, contains('R2')) %>%
                mutate(Res = grepl('_res', Model)) %>%
                mutate(Model = case_when(
                  grepl('CTSLEB2', Model) ~ 'CTSLEB (all)',
                  grepl('CTSLEB', Model) ~ 'CTSLEB'
                ))              
              r2_results = rbind(r2_results, r2_result0)
              
              
              {
                snp = unique(c(
                  readRDS(paste0('CTSLEB/PRS/snpEB_AFR_train', train, '_rho', rho, 
                                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')),
                  readRDS(paste0('CTSLEB/PRS/snpEB_EAS_train', train, '_rho', rho, 
                                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')),
                  readRDS(paste0('CTSLEB/PRS/snpEB_SAS_train', train, '_rho', rho, 
                                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS')),
                  readRDS(paste0('CTSLEB/PRS/snpEB_AMR_train', train, '_rho', rho, 
                                 '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '.RDS'))
                ))
              }
              
              supp = length(snp)
              
              snp_result0 = data.frame(train=train, 
                                       rho=rho, 
                                       gr=gr,
                                       GA=GA,
                                       rep=rep, 
                                       model=model, 
                                       Model = 'CTSLEB',
                                       Supp = supp,
                                       Component = 'Main') 
              snp_results = rbind(snp_results, snp_result0)
              
              rm(results0, beta, snp, supp); gc()
            }
            
          }
        }
      }
    }
  }
}


saveRDS(r2_results, paste0('multi_simulation_r2_results_revisions.RDS'))
saveRDS(snp_results, paste0('multi_simulation_snp_results_revisions.RDS'))
print('DONE')

q('no')
