# use METAL to identify heterogeneous SNP effects
# use glmnet for refit analysis

####################
### extract SNPs ###
####################
# ipgs_extract.R
library(data.table)
library(dplyr)

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
# train='1a'; rho=1

sim = '/n/holystore01/LABS/xlin/Everyone/JD_HZ/'
multi = '/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/'
data = '/n/holystore01/LABS/xlin/Lab/tonychen/Data/1000GP_Phase3/'

home = '/n/holystore01/LABS/xlin/Lab/tonychen/2.Multi-Ancestry/Simulation/'
scratch = '/n/holyscratch01/xlin/tonychen/MultiAncestry/Simulation/'

gr_GA = data.frame(
  gr=c(0.6, 0.8, 0.4, 0.6),
  GA=c(1,1,2,2)
)

for (rho in c(1, 1.5, 2, 3)) {
  for (train in c('1a', '2a', '3a', '4a')) {
    snp = c()
    for (i in 1:4) {
      gr = gr_GA$gr[i]
      GA = gr_GA$GA[i]
      
      for (rep in 1:10) {
        snp0 = try(read.table(paste0('METAL/metal_train', train, 
                                     '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_hm3_meta.snp'), 
                              header=F)$V1)
        if (!'try-error' %in% class(snp0)) snp = c(snp, snp0)
      }
    }
    
    snp = unique(snp)
    snp_map = fread(paste0('sim_mapping.txt')) %>%
      filter(rsid %in% snp) %>%
      pull(snp_code)
    write.table(snp_map, paste0('METAL/metal_train', train, '_rho', rho, '_hm3_meta.snp'),
                row.names=F, quote=F, col.names=F)
    
    for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
      system(paste0('~/software/plink2 ',
                    '--bfile ', sim, 'ALL/all_chr ',
                    '--keep ', 'ID/train', train, '_', anc, '.id ',
                    '--extract METAL/metal_train', train, '_rho', rho, '_hm3_meta.snp ',
                    '--make-bed ',
                    '--out PLINK/META/meta_train', train, '_rho', rho, '_', anc))
    }
  }
}

quit('no')

for (train in c('2a', '4a')) {
  train_id = fread(paste0('ID/train', train, '.id'), header=F, col.names=c('FID', 'IID'))
  tunval_id = fread('ID/subtunval.id', header=F, col.names=c('FID', 'IID'))
  all_id = rbind(train_id, tunval_id)
  write.table(all_id, paste0('ID/alltrain', train, '.id'),
              row.names=F, col.names=F, quote=F)
}

######################
### refit analysis ###
######################
suppressPackageStartupMessages({
  library(glmnet)
  library(glmnetUtils)
  library(dplyr)
  library(data.table)
  library(bigsnpr)
  library(tidyr)
  library(foreach)
  library(doParallel)
})

registerDoParallel(6)

silent = function(x) suppressMessages(suppressWarnings(x))

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = as.numeric(params[[3]])
GA = as.numeric(params[[4]])
model = params[[5]]
if (length(params) > 5) {
  rep0 = as.numeric(params[[6]])
} else {
  rep0 = 1:10
}
print(rep0)

# main data
print('data')
sim_map = fread('sim_mapping.txt')
plink = snp_attach(paste0('PLINK/sim_train', train, '_tun.rds'))

train_id = fread(paste0('ID/train', train, '.id'), header=F, col.names=c('FID', 'IID'))
tuning_id = fread('ID/subtuning.id', header=F, col.names=c('FID', 'IID'))
validation_id = fread(paste0('ID/validation.id'), header=F, col.names=c('FID', 'IID'))
Esum = readRDS(paste0('SUMMARY/Esum_train', train, '.RDS'))

pc = fread(paste0('PCA/all_sim_pca.sscore'),
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt')) 

for (rep in rep0) {
  if (file.exists(paste0('Lasso/RESULTS/valREFIT_train', train, '_rho', rho, 
                         '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) next
  
  print(rep)
  
  data = silent(pheno %>%
    select(FID, IID, Y=paste0('PHENO', rep)) %>%
    left_join(pc) %>%
    filter(FID %in% c(train_id$FID)) %>%
    mutate(ANC = substr(FID, 1, 3)))
  data$PC1 = (data$PC1 - Esum$mean[1]) / Esum$sd[1]
  data$PC2 = (data$PC2 - Esum$mean[2]) / Esum$sd[2]
  
  tunval_data = silent(pheno %>%
    select(FID, IID, Y=paste0('PHENO', rep)) %>%
    left_join(pc) %>%
    filter(FID %in% c(tuning_id$FID, validation_id$FID)) %>%
    mutate(ANC = substr(FID, 1, 3)))
  
  fit_tun = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% tuning_id$FID))
  fit_val = lm(paste0('Y~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% validation_id$FID))
  
  tunval_data$res = 0
  tunval_data$res[tunval_data$FID %in% tuning_id$FID] = residuals(fit_tun)
  tunval_data$res[tunval_data$FID %in% validation_id$FID] = residuals(fit_val)
  
  val_data = tunval_data %>% filter(FID %in% validation_id$FID)
  val_data$PC1 = (val_data$PC1 - Esum$mean[1]) / Esum$sd[1]
  val_data$PC2 = (val_data$PC2 - Esum$mean[2]) / Esum$sd[2]
  
  rm(tunval_data); gc()
  
  # compute main effects iPGS
  print('iPGS')
  system(paste0('~/software/plink2 ',
                '--silent ',
                '--bfile all_chr ',
                '--keep ', 'ID/train', train, '.id ',
                '--extract Lasso/RESULTS/beta_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp ',
                '--score Lasso/RESULTS/beta_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 2 3 5 header ',
                'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                '--out Lasso/PRS/prsTRAIN_add_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model))

  system(paste0('~/software/plink2 ',
                '--silent ',
                '--bfile all_chr ',
                '--keep ', 'ID/validation.id ',
                '--extract Lasso/RESULTS/beta_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp ',
                '--score Lasso/RESULTS/beta_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 2 3 5 header ',
                'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                '--out Lasso/PRS/prsVAL_add_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model))
  
  # cv.glmnet with interactions (in training)
  print('CV')
  meta_snps = try(read.table(paste0('METAL/metal_train', train, 
                                '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, '_meta.snp'),
                         header=F)$V1)
  if ('try-error' %in% class(meta_snps)) {
    meta_snps = character(0)
  }
  
  meta_snps_map = sim_map %>% filter(rsid %in% meta_snps)
  
  prs_train = fread(paste0('Lasso/PRS/prsTRAIN_add_train', train, '_rho', rho, 
                           '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                    col.names=c('FID', 'IID', 'iPGS'))
  
  ix_snp = which(plink$map$marker.ID %in% meta_snps_map$snp_code)
  cat(length(ix_snp), ' het. SNPs \n')
  
  if (length(meta_snps_map$rsid) > 0) {
    formula = as.formula(paste0('Y~', paste(paste0('PC', 1:20), collapse='+'), '+',
                                'iPGS + ',
                                paste(meta_snps_map$rsid, collapse='+'), '+',
                                paste(paste0(meta_snps_map$rsid,':PC1'), collapse='+'), '+',
                                paste(paste0(meta_snps_map$rsid,':PC2'), collapse='+')))
  } else {
    formula = as.formula(paste0('Y~', paste(paste0('PC', 1:20), collapse='+'), '+',
                                'iPGS'))
  }
  
  iPGS_coef = foreach(anc = c('EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'ALL'), .combine='c') %dopar% {
    message = paste0(anc, ' started \n')
    cat(message)
    
    if (anc == 'ALL') {
      ix_train = which(plink$fam$fam %in% train_id$FID)
    } else {
      ix_train = which(plink$fam$fam %in% train_id$FID & substr(plink$fam$fam, 1, 3) == anc)
    }

    if (length(ix_snp) > 1) {
      G = plink$genotypes[,ix_snp][ix_train,]
      colnames(G) = meta_snps_map$rsid
      G[is.na(G)] = 0
    } else if (length(ix_snp) == 1) {
      G = matrix(plink$genotypes[,ix_snp][ix_train], ncol=1)
      colnames(G) = meta_snps_map$rsid
      G[is.na(G)] = 0
    } else {
      G = NULL
    }

    if (anc == 'ALL') {
      reg_data = silent(data %>%
        filter(FID %in% train_id$FID) %>%
        left_join(prs_train) %>%
        select(-FID, -IID, -ANC) %>%
        cbind(G))
    } else {
      reg_data = silent(data %>%
        filter(FID %in% train_id$FID & ANC == anc) %>%
        left_join(prs_train) %>%
        select(-FID, -IID, -ANC) %>%
        cbind(G))
    }
    
    rm(G); gc()

    set.seed(123)
    cv = cv.glmnet(formula, data=reg_data,
                   standardize=F,
                   weights=rep(1, nrow(reg_data)),
                   penalty.factor=c(rep(0, 20), # covariates
                                    1, # iPGS,
                                    rep(1.1, length(meta_snps_map$rsid)), # main effects
                                    rep(1.2, length(meta_snps_map$rsid)*2) # interactions
                   ))

    coef_cv = coef(cv, s='lambda.min')

    if (length(ix_snp) > 0) {
      coef_out = data.frame(Var = rownames(coef_cv),
                            Beta = as.numeric(coef_cv))[-c(1:21),] %>%
        separate(Var, into=c('rsid', 'Model'), sep=':') %>%
        mutate(Model = case_when(
          Model == 'PC1' ~ 'INT1',
          Model == 'PC2' ~ 'INT2',
          T ~ 'ADD'
        )) %>%
        spread(key=Model, value=Beta) %>%
        filter(ADD != 0 | INT1 != 0 | INT2 != 0)
      coef_out[is.na(coef_out)] = 0
    } else {
      coef_out = data.frame(Var = rownames(coef_cv),
                            ADD = as.numeric(coef_cv))[-c(1:21),] %>%
        mutate(INT1 = 0, INT2 = 0)
    }

    if (length(ix_snp) > 0) {
      map_out = silent(coef_out %>%
        left_join(data.frame(plink$map[ix_snp,], rsid=meta_snps_map$rsid)) %>%
        select(chromosome, physical.pos, rsid, marker.ID, allele1, allele2, ADD, INT1, INT2))
    } else {
      map_out = silent(coef_out %>%
        mutate(chromosome=NA, physical.pos=NA, rsid=NA, marker.ID=NA,
               allele1='A', allele2='T') %>%
        select(chromosome, physical.pos, rsid, marker.ID, allele1, allele2, ADD, INT1, INT2))
    }

    write.table(map_out, paste0('Lasso/RESULTS/betaREFIT_', anc, '_train', train, '_rho', rho,
                                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'),
                row.names=F, quote=F)
    write.table(map_out$marker.ID, paste0('Lasso/RESULTS/betaREFIT_', anc, '_train', train, '_rho', rho,
                                          '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp'),
                row.names=F, col.names=F, quote=F)

    message = paste0(anc, ' finished \n')
    cat(message)

    x = coef_out$ADD[1]
    names(x) = anc
    x
  }
  iPGS_coef[is.na(iPGS_coef)] = 0

  saveRDS(iPGS_coef, paste0('Lasso/RESULTS/iPGS_coef_train', train, '_rho', rho,
                            '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  iPGS_coef = readRDS(paste0('Lasso/RESULTS/iPGS_coef_train', train, '_rho', rho, 
                             '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  
  print(iPGS_coef)
  
  # compute validation
  print('validation')
  for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    system(paste0('~/software/plink2 ',
                  '--silent ',
                  '--bfile ', sim, 'ALL/all_chr ',
                  '--keep ', 'ID/validation_', anc, '.id ',
                  '--extract Lasso/RESULTS/betaREFIT_', anc, '_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp ',
                  '--score Lasso/RESULTS/betaREFIT_', anc, '_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 4 5 header-read ',
                  'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                  '--score-col-nums 7-9 ',
                  '--out Lasso/PRS/prsVAL_int_', anc, '_train', train, '_rho', rho,
                  '_gr', gr, '_GA', GA, '_rep', rep, '_', model))
  }

  system(paste0('~/software/plink2 ',
                '--silent ',
                '--bfile ', sim, 'ALL/all_chr ',
                '--keep ', 'ID/validation.id ',
                '--extract Lasso/RESULTS/betaREFIT_ALL_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp ',
                '--score Lasso/RESULTS/betaREFIT_ALL_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt 4 5 header-read ',
                'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
                '--score-col-nums 7-9 ',
                '--out Lasso/PRS/prsVAL_int_ALL_train', train, '_rho', rho,
                '_gr', gr, '_GA', GA, '_rep', rep, '_', model))
  
  val_data = val_data %>%
    select(FID, IID, res, paste0('PC', 1:2))
  
  prsVAL_add = fread(paste0('Lasso/PRS/prsVAL_add_train', train, '_rho', rho, 
                            '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                     col.names=c('FID', 'IID', 'iPGS'))
  prsVAL_add = prsVAL_add %>%
    mutate(ancestry = substr(FID, 1, 3)) %>%
    mutate(iPGS_wtd = iPGS * iPGS_coef[ancestry]) %>%
    select(-ancestry)
  
  prsVAL_int = NULL
  for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    if (file.exists(paste0('Lasso/PRS/prsVAL_int_', anc, '_train', train, '_rho', rho, 
                           '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'))) {
      prsVAL_int0 = fread(paste0('Lasso/PRS/prsVAL_int_', anc, '_train', train, '_rho', rho, 
                                 '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                          col.names=c('FID', 'IID', 'ADD', 'INT1', 'INT2'))
      prsVAL_int = rbind(prsVAL_int, prsVAL_int0)
    } 
  }
  
  if (file.exists(paste0('Lasso/PRS/prsVAL_int_ALL_train', train, '_rho', rho, 
                         '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'))) {
    prsVAL_all = fread(paste0('Lasso/PRS/prsVAL_int_ALL_train', train, '_rho', rho, 
                              '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.sscore'),
                       col.names=c('FID', 'IID', 'ADD', 'INT1', 'INT2'))
  } else {
    prsVAL_all = NULL
  }
  
  # ancestry-specific refit
  if (!is.null(prsVAL_int)) {
    val_data_anc = val_data %>%
      left_join(prsVAL_add) %>%
      left_join(prsVAL_int) %>%
      mutate(add_term = ifelse(is.na(ADD), 0, ADD + INT1*PC1 + INT2*PC2)) %>%
      mutate(iPGS_refit = iPGS_wtd + add_term)
  } else {
    val_data_anc = val_data %>%
      left_join(prsVAL_add) %>%
      mutate(iPGS_refit = iPGS_wtd)
  }
  
  # all-ancestry refit
  if (!is.null(prsVAL_all)) {
    val_data_all = val_data %>%
      left_join(prsVAL_add) %>%
      left_join(prsVAL_all) %>%
      mutate(add_term = ifelse(is.na(ADD), 0, ADD + INT1*PC1 + INT2*PC2)) %>%
      mutate(iPGS_refit2 = iPGS_wtd + add_term)
  } else {
    val_data_all = val_data %>%
      left_join(prsVAL_add) %>%
      mutate(iPGS_refit2 = iPGS_wtd)
  }
  
  val_data = val_data_anc %>%
    left_join(val_data_all %>% select(FID, IID, iPGS_refit2))
  
  # extract SNPs
  # print('SNPs')
  iPGS_snp = read.table(paste0('Lasso/RESULTS/beta_train', train, '_rho', rho,
                               '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.snp'), header=F)$V1

  iPGS_refit_beta = NULL
  for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    iPGS_refit_beta = rbind(iPGS_refit_beta,
                            read.table(paste0('Lasso/RESULTS/betaREFIT_', anc, '_train', train, '_rho', rho,
                                              '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'), header=T))
  }
  iPGS_snp_refit_add = iPGS_refit_beta %>% filter(ADD != 0) %>% pull(marker.ID) %>% na.omit() %>% unique()
  iPGS_snp_refit_int = iPGS_refit_beta %>% filter(INT1 != 0 | INT2 != 0) %>% pull(marker.ID) %>% unique()
  
  supp_add = length(unique(c(iPGS_snp, iPGS_snp_refit_add)))
  supp_int = length(unique(iPGS_snp_refit_int))

  
  iPGS_refit_beta2 = read.table(paste0('Lasso/RESULTS/betaREFIT_ALL_train', train, '_rho', rho,
                                       '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.txt'), header=T)

  iPGS_snp_refit_add2 = iPGS_refit_beta2 %>% filter(ADD != 0) %>% pull(marker.ID) %>% na.omit() %>% unique()
  iPGS_snp_refit_int2 = iPGS_refit_beta2 %>% filter(INT1 != 0 | INT2 != 0) %>% pull(marker.ID) %>% unique()
  
  supp_add2 = length(unique(c(iPGS_snp, iPGS_snp_refit_add2)))
  supp_int2 = length(unique(iPGS_snp_refit_int2))
  
  # compile results
  print('results')
  
  val_data = val_data %>% select(-PC1, -PC2) %>% left_join(pc)
  
  val_data$iPGS_res = residuals(lm(iPGS ~ ., data=val_data %>% select(iPGS, paste0('PC', 1:20))))
  val_data$iPGS_refit_res = residuals(lm(iPGS_refit ~ ., data=val_data %>% select(iPGS_refit, paste0('PC', 1:20))))
  val_data$iPGS_refit2_res = residuals(lm(iPGS_refit2 ~ ., data=val_data %>% select(iPGS_refit2, paste0('PC', 1:20))))
  
  r2_results = val_data %>%
    mutate(ancestry = substr(FID, 1, 3)) %>%
    bind_rows(., val_data %>% mutate(ancestry = 'ALL')) %>%
    bind_rows(., val_data %>% mutate(ancestry = substr(FID, 1, 3)) %>%
                filter(ancestry != 'EUR') %>% mutate(ancestry = 'nonEUR')) %>%
    mutate(ancestry = factor(ancestry, 
                             levels=c('ALL', 'nonEUR', 'EUR', 'AFR', 'EAS', 'SAS', 'AMR'))) %>%
    group_by(ancestry) %>%
    summarize(n=n(),
              R2_iPGS = cor(iPGS, res)^2,
              R2_refit = cor(iPGS_refit, res)^2,
              R2_refit2 = cor(iPGS_refit2, res)^2,
              R2_iPGS_res = cor(iPGS_res, res)^2,
              R2_refit_res = cor(iPGS_refit_res, res)^2,
              R2_refit2_res = cor(iPGS_refit2_res, res)^2) %>%
    ungroup() %>%
    mutate(train=train,
           rho=rho,
           gr=gr,
           GA=GA,
           rep=rep,
           model=model,
           supp_iPGS = length(iPGS_snp),
           supp_add = supp_add,
           supp_int = supp_int,
           supp_add2 = supp_add2,
           supp_int2 = supp_int2)
  
  print(data.frame(r2_results))
  
  saveRDS(r2_results, paste0('Lasso/RESULTS/valREFIT_train', train, '_rho', rho, 
                             '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  rm(data, tunval_data, val_data, reg_data, prsVAL_add, prsVAL_int); gc()
}

print('DONE')

quit('no')
