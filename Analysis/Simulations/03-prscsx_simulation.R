# PRS-CSx: https://github.com/getian107/PRScsx
# we ran PRS-CSx in groups of chromosomes (batch job script below)

#################
## RUN PRS-CSX ##
#################

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(foreach)
})

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]]
GA = params[[4]]
rep = as.numeric(params[[5]])
model = params[[6]]
chr_id = as.numeric(params[[7]])

registerDoParallel(cores=4)

chr_vec = c(
  '1,2', '3,4', '5,6,7',
  '8,9,10,11', '12,13,14,15,16', '17,18,19,20,21,22'
)

chr_string = chr_vec[chr_id]
chr_num = as.numeric(strsplit(chr_string, split=',')[[1]])

chr_num_out = NULL
for (chrom in chr_num) {
  files = sapply(c('1e+00', '1e-02', '1e-04', '1e-06'), FUN=function(phi) {
    !file.exists(paste0('PRSCSX/train', train, '_rho', rho, '_gr', gr,
                        '_rep', rep, '_GA', GA, '_', model,
                        '_SAS_pst_eff_a1_b0.5_phi', phi, '_chr', chrom, '.txt'))
  })
  if (sum(files) > 0) {
    chr_num_out = c(chr_num_out, chrom)
  }
}
chr_string = paste(chr_num_out, collapse=',')

print(chr_string)

if (train == '2a') {
  n_gwas = paste(c(7500,7500,7500,70000,7500), collapse=',')
} else if (train == '4a') {
  n_gwas = paste(c(20000,20000,20000,20000,20000), collapse=',')
}
print(n_gwas)

# run PRSCSX
analysis = foreach(phi=c('1e+00', '1e-02', '1e-04', '1e-06')) %dopar% {
  print(phi)
  system(paste0(
    'python ', multi, '/PRSCSX/Package/PRScsx.py',
    ' --ref_dir=', package,
    ' --bim_prefix=PLINK/subtunval ',
    ' --sst_file=', paste(paste0('PRSCSX/GWAS/', c('AFR', 'AMR', 'EAS', 'EUR', 'SAS'), 
                                 '_gwas_train', train, '_rho', rho, '_gr', gr, 
                                 '_rep', rep, '_GA', GA, '_', model, '.txt'), 
                          collapse=','),
    ' --n_gwas=', n_gwas,
    ' --pop=AFR,AMR,EAS,EUR,SAS',
    ' --out_dir=', out,
    ' --out_name=train', train, '_rho', rho, '_gr', gr, '_rep', rep, '_GA', GA, '_', model, 
    ' --phi=', phi,
    ' --seed 123',
    ' --chrom=', chr_string
  ), intern=T)
}

#############
## SCORING ##
#############

params = commandArgs(trailingOnly=TRUE)
train = params[[1]]
rho = params[[2]]
gr = params[[3]] 
GA = params[[4]]
model = params[[5]]

tun_id = fread('ID/subtuning.id', header=F,
               col.names=c('FID', 'IID'))
val_id = fread('ID/validation.id', header=F,
               col.names=c('FID', 'IID'))

pc = fread('PCA/all_sim_pca.sscore',
           col.names=c('FID', 'IID', paste0('PC', 1:20)))

pheno = fread(paste0('PHENO/phenotypes_rho', rho, '_gr', gr, '_GA', GA, '_', model, '.txt'))

phi_vec = c('1e+00', '1e-02', '1e-04', '1e-06')
phi_vec2 = c('1e.00', '1e.02', '1e.04', '1e.06')

for (rep in 1:10) {
  print(rep)
  
  if (file.exists(paste0('PRSCSX/Results/results_train_', train, 
                         '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))) next
  
  n_files = length(list.files(path = 'PRSCSX/Analysis/', 
                              pattern=glob2rx(paste0('train', train, '_rho', rho, '_gr', gr, 
                                                     '_rep', rep, '_GA', GA, '_', model, '*'))))

  prs_tunval = fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                            '_rep', rep, '_GA', GA, '_', model, '_EUR.sscore'),
                     col.names=c('FID', 'IID', paste0('PRS_EUR_', phi_vec))) %>%
    full_join(fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                           '_rep', rep, '_GA', GA, '_', model, '_AFR.sscore'),
                    col.names=c('FID', 'IID', paste0('PRS_AFR_', phi_vec)))) %>%
    full_join(fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                           '_rep', rep, '_GA', GA, '_', model, '_EAS.sscore'),
                    col.names=c('FID', 'IID', paste0('PRS_EAS_', phi_vec)))) %>%
    full_join(fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                           '_rep', rep, '_GA', GA, '_', model, '_SAS.sscore'),
                    col.names=c('FID', 'IID', paste0('PRS_SAS_', phi_vec)))) %>%
    full_join(fread(paste0('PRSCSX/PRS/train', train, '_rho', rho, '_gr', gr, 
                           '_rep', rep, '_GA', GA, '_', model, '_AMR.sscore'),
                    col.names=c('FID', 'IID', paste0('PRS_AMR_', phi_vec))))
  
  tunval_data = pheno %>%
    select(FID, IID, Pheno=paste0('PHENO', rep)) %>%
    filter(FID %in% c(tun_id$FID, val_id$FID)) %>%
    left_join(pc) %>%
    mutate(Ancestry = substr(FID, 1, 3)) %>%
    left_join(prs_tunval) %>%
    data.frame()
  
  # compute residuals for tuning
  fit_tun = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% tun_id$FID))
  fit_val = lm(paste0('Pheno~', paste(paste0('PC', 1:20), collapse='+')), tunval_data %>% filter(FID %in% val_id$FID))
  
  tunval_data$Y = 0
  tunval_data$Y[tunval_data$FID %in% tun_id$FID] = residuals(fit_tun)
  tunval_data$Y[tunval_data$FID %in% val_id$FID] = residuals(fit_val)
  
  tunval_data$PRSCSX = 0
  
  tun_data = tunval_data %>%
    filter(FID %in% tun_id$FID) %>%
    select(FID, IID, Ancestry, Y, contains('PRS'))
  
  val_data = tunval_data %>%
    filter(FID %in% val_id$FID) %>%
    select(FID, IID, Ancestry, Y, contains('PRS'))
  
  # ancestry-specific tuning
  for (pop in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
    r2 = rep(0, 4)
    fit_lm = list()
    for (i in 1:4) {
      if (train == '0') {
        fit_lm[[i]] = lm(Y~., data = tun_data %>% filter(Ancestry==pop) %>% 
                           select(Y, paste0('PRS_', phi_vec2[i])))
        
      } else {
        fit_lm[[i]] = lm(Y~., data = tun_data %>% filter(Ancestry==pop) %>% 
                           select(Y, paste0('PRS_', c('EUR', 'AFR', 'EAS', 'SAS', 'AMR'), '_', phi_vec2[i])))
        
      }
      
      r2[i] = summary(fit_lm[[i]])$r.squared
    }
    # print(r2)
    
    val_data$PRSCSX[val_data$Ancestry==pop] = predict(fit_lm[[which.max(r2)]], val_data %>% filter(Ancestry==pop))
  }
  
  # all-ancestry tuning
  r2 = rep(0, 4)
  fit_lm = list()
  for (i in 1:4) {
    if (train == '0') {
      fit_lm[[i]] = lm(Y~., data = tun_data %>% 
                         select(Y, paste0('PRS_', phi_vec2[i])))
    } else {
      fit_lm[[i]] = lm(Y~., data = tun_data %>% 
                         select(Y, paste0('PRS_', c('EUR', 'AFR', 'EAS', 'SAS', 'AMR'), '_', phi_vec2[i])))
    }
    
    r2[i] = summary(fit_lm[[i]])$r.squared
  }
  val_data$PRSCSX2 = predict(fit_lm[[which.max(r2)]], val_data)
  
  prscsx_result = val_data %>%
    bind_rows(., val_data %>% mutate(Ancestry = 'ALL')) %>%
    mutate(Ancestry = factor(Ancestry, 
                             levels=c('ALL', 'EUR', 'AFR', 'EAS', 'SAS', 'AMR'))) %>%
    group_by(Ancestry) %>%
    summarize(n=n(),
              R2_PRSCSX = cor(PRSCSX, Y)^2,
              R2_PRSCSX2 = cor(PRSCSX2, Y)^2)
  print(prscsx_result)
  
  saveRDS(val_data, paste0('PRSCSX/Results/val_train_', train, 
                           '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
  saveRDS(prscsx_result, paste0('PRSCSX/Results/results_train_', train, 
                                '_rho', rho, '_gr', gr, '_GA', GA, '_rep', rep, '_', model, '.RDS'))
  
}


##################
## BATCH SCRIPT ##
##################
## PRS-CSx
# model=aou
# 
# chr_max=( 0 2 4 7 11 16 22 )
# 
# for rho in {1..4}; do
# for train in 2a 4a; do
# for chr_id in {1..6}; do
# for rep in {1..10}; do
# 
# GA=1
# for gr in 0.6 0.8; do
# chr_max0=${chr_max[$chr_id]}
# files=$(ls PRSCSX/Analysis/train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}_EUR_pst_eff_a1_b0.5_phi*_chr${chr_max0}.txt | wc -l)
# if ! [[ $files -eq 4 ]]; then
# job=$(sbatch --parsable -o jobs/prscsx_${train}_${rho}_${gr}_${GA}_${rep}_${model}_${chr_id}.out prscsx_sim.sh $train $rho $gr $GA $rep $model $chr_id)
# sleep 1
# fi
# done
# 
# # heterogeneous effects
# GA=2
# for gr in 0.4 0.6; do
# chr_max0=${chr_max[$chr_id]}
# files=$(ls PRSCSX/Analysis/train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}_EUR_pst_eff_a1_b0.5_phi*_chr${chr_max0}.txt | wc -l)
# if ! [[ $files -eq 4 ]]; then
# job=$(sbatch --parsable -o jobs/prscsx_${train}_${rho}_${gr}_${GA}_${rep}_${model}_${chr_id}.out prscsx_sim.sh $train $rho $gr $GA $rep $model $chr_id)
# sleep 1
# fi
# done
# 
# done
# done
# done
# done

## PRS-CSx Scoring
# for train in 2a 4a; do
# for rho in 1 2 3 4; do
# GA=1
# for gr in 0.6 0.8; do
# sbatch -o jobs/prscsx_score_${rho}_${gr}_${GA}.out prscsx_score.sh $train $rho $gr $GA $model
# sleep 1
# done
# 
# GA=2
# for gr in 0.6 0.4; do
# sbatch -o jobs/prscsx_score_${rho}_${gr}_${GA}.out prscsx_score.sh $train $rho $gr $GA $model
# sleep 1
# done
# done
# done