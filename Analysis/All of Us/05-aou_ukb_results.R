library(dplyr)
library(data.table)
library(glmnet)
library(tidyr)

#############
### SETUP ###
#############
# convert allele pairs for possible strand flip
alleles = c('A', 'C', 'G', 'T')
allele_pair = NULL
for (a in 1:4) {
  for (b in setdiff(1:4, a)) {
    allele_pair = rbind(allele_pair, c(
      paste0(alleles[a], alleles[b]), 
      paste0(alleles[b], alleles[a]), 
      paste0(alleles[5-a], alleles[5-b]), 
      paste0(alleles[5-b], alleles[5-a])
    ))
  }
}
rownames(allele_pair)=allele_pair[,1]

strand_match = Vectorize(function(a, b) {
  if (!(a %in% rownames(allele_pair))) {
    return(F)
  } else {
    b %in% allele_pair[a,]
  }
})

############
### DATA ###
############
print('data')

split1_id = fread(paste0('split1.id'))
split2_id = fread(paste0('split2.id'))
all_id = fread(paste0('all.id'))

args = commandArgs(trailingOnly = T)
trait = args[[1]]
# trait = 'LDL'

pheno = fread(paste0('ukb_multi.pheno')) %>%
  select(FID, IID, Y=ends_with(trait))
cov = fread(paste0('ukb_multi.cov')) %>%
  select(FID, IID, age, female)
anc = fread(paste0('ukb_multi.anc')) %>%
  select(FID, IID, ancestry=predicted, contains('P.'))
aou_pc = fread('/n/holystore01/LABS/xlin/Lab/tonychen/Data/1000GP_Phase3/ukb_aou_pca.sscore') %>%
  select(FID=`#FID`, IID, PC=contains('SCORE'))

data = pheno %>%
  filter(FID %in% all_id$FID) %>%
  left_join(anc) %>%
  left_join(cov) %>%
  left_join(aou_pc) %>%
  na.omit() %>% 
  data.frame()

fit_tun = lm(paste0('Y~age+female+', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% split1_id$FID))
fit_val = lm(paste0('Y~age+female+', paste(paste0('PC', 1:20), collapse='+')), data %>% filter(FID %in% split2_id$FID))

data$res = 0
data$res[data$FID %in% split1_id$FID] = residuals(fit_tun)
data$res[data$FID %in% split2_id$FID] = residuals(fit_val)

split1 = data %>%
  filter(FID %in% split1_id$FID)

# downsampled EUR tuning
split1a = data %>%
  filter(FID %in% split1_id$FID) %>%
  group_by(ancestry) %>%
  slice_sample(n=5000) %>%
  ungroup() %>%
  arrange(FID)

table(split1$ancestry)
table(split1a$ancestry)

table(split1$ancestry) / nrow(split1)
table(split1a$ancestry) / nrow(split1a)

split2 = data %>%
  filter(FID %in% split2_id$FID) %>%
  select(FID, IID, ancestry, Y, res)

rm(pheno, cov, anc, aou_pc); gc()

snps = matrix(0, nrow=2, ncol=14)
rownames(snps) = c('main', 'int')
colnames(snps) = c('PROSPER', 'CTSLEB', 'PRSCSX', 'iPGS', 'iPGS + Refit', 'SPLENDID')

##############
### PRSCSX ###
##############
print('PRSCSX')
phi = c('1e-06', '1e-04', '1e-02', '1e+00')
prscsx = list()
prscsx_snps = list()
for (i in 1:4) {
  prscsx[[i]] = fread(paste0('PRSCSX/Results/prs_', trait, '_', phi[i], '_EUR.sscore'),
                      col.names=c('FID', 'IID', 'NALLELE', 'DOSAGESUMS', 'EUR')) %>%
    select(-NALLELE, -DOSAGESUMS) %>%
    left_join(
      fread(paste0('PRSCSX/Results/prs_', trait, '_', phi[i], '_AFR.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGESUMS', 'AFR')) %>%
        select(-NALLELE, -DOSAGESUMS)
    ) %>%
    left_join(
      fread(paste0('PRSCSX/Results/prs_', trait, '_', phi[i], '_EAS.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGESUMS', 'EAS')) %>%
        select(-NALLELE, -DOSAGESUMS)
    ) %>%
    left_join(
      fread(paste0('PRSCSX/Results/prs_', trait, '_', phi[i], '_SAS.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGESUMS', 'SAS')) %>%
        select(-NALLELE, -DOSAGESUMS)
    ) %>%
    left_join(
      fread(paste0('PRSCSX/Results/prs_', trait, '_', phi[i], '_AMR.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGESUMS', 'AMR')) %>%
        select(-NALLELE, -DOSAGESUMS)
    ) %>%
    data.frame() %>%
    left_join(data %>% select(FID, IID, res, ancestry))
  
  prscsx_snps[[i]] = unique(c(
    fread(paste0('PRSCSX/Results/', trait, '_EUR_pst_eff_a1_b0.5_phi', phi[i], '_beta.snp'),header=F)$V1,
    fread(paste0('PRSCSX/Results/', trait, '_AFR_pst_eff_a1_b0.5_phi', phi[i], '_beta.snp'),header=F)$V1,
    fread(paste0('PRSCSX/Results/', trait, '_EAS_pst_eff_a1_b0.5_phi', phi[i], '_beta.snp'),header=F)$V1,
    fread(paste0('PRSCSX/Results/', trait, '_AMR_pst_eff_a1_b0.5_phi', phi[i], '_beta.snp'),header=F)$V1,
    fread(paste0('PRSCSX/Results/', trait, '_SAS_pst_eff_a1_b0.5_phi', phi[i], '_beta.snp'),header=F)$V1
  ))
}

# grid search: linear-combination of ancestry-specific models
split2$PRSCSX = 0
for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
  prscsx_anc = list()
  r2_anc = rep(0, 4)
  
  for (i in 1:4) {
    prscsx_anc[[i]] = lm(res ~ EUR + AFR + EAS + SAS + AMR, prscsx[[i]] %>% 
                           filter(FID %in% split1$FID) %>% filter(ancestry==anc))
    r2_anc[i] = summary(prscsx_anc[[i]])$r.squared
  }
  
  split2$PRSCSX[split2$ancestry==anc] = predict(prscsx_anc[[which.max(r2_anc)]], prscsx[[which.max(r2_anc)]] %>% 
                                                  filter(FID %in% split2$FID) %>% filter(ancestry==anc))
}

snps[1,'PRSCSX'] = length(unique(unlist(prscsx_snps)))

rm(prscsx, prscsx_anc); gc()

###############
### CT-SLEB ###
###############
print('CT-SLEB')
ctsleb_split2 = rbind(
  readRDS(paste0('CTSLEB/aou_', trait, '_EUR_prs.RDS')) %>%
    select(FID, IID, CTSLEB),
  readRDS(paste0('CTSLEB/aou_', trait, '_AFR_prs.RDS')) %>%
    select(FID, IID, CTSLEB),
  readRDS(paste0('CTSLEB/aou_', trait, '_EAS_prs.RDS')) %>%
    select(FID, IID, CTSLEB),
  readRDS(paste0('CTSLEB/aou_', trait, '_SAS_prs.RDS')) %>%
    select(FID, IID, CTSLEB),
  readRDS(paste0('CTSLEB/aou_', trait, '_AMR_prs.RDS')) %>%
    select(FID, IID, CTSLEB)
) %>% 
  arrange(FID, IID)

snp_ctsleb = read.table(paste0('CTSLEB/', trait, '_eb_snps.txt'), header=F)$V1
snps[1,'CTSLEB'] = length(unique(snp_ctsleb))

split2 = split2 %>%
  left_join(ctsleb_split2)

rm(ctsleb_split2); gc()

###############
### PROSPER ###
###############
print('PROSPER')
# PROSPER PRS for split2
prosper_split2 = rbind(
  fread(paste0('Prosper/AOU/', trait, '/PROSPER/tmp/sample_scores_EUR/after_ensemble_testing.txt'),
        col.names=c('FID', 'IID', 'PROSPER')),
  fread(paste0('Prosper/AOU/', trait, '/PROSPER/tmp/sample_scores_AFR/after_ensemble_testing.txt'),
        col.names=c('FID', 'IID', 'PROSPER')),
  fread(paste0('Prosper/AOU/', trait, '/PROSPER/tmp/sample_scores_EAS/after_ensemble_testing.txt'),
        col.names=c('FID', 'IID', 'PROSPER')),
  fread(paste0('Prosper/AOU/', trait, '/PROSPER/tmp/sample_scores_SAS/after_ensemble_testing.txt'),
        col.names=c('FID', 'IID', 'PROSPER')),
  fread(paste0('Prosper/AOU/', trait, '/PROSPER/tmp/sample_scores_AMR/after_ensemble_testing.txt'),
        col.names=c('FID', 'IID', 'PROSPER'))
)

prosper_snps = list(
  eur=fread(paste0('Prosper/AOU/', trait, '/PROSPER/after_ensemble_EUR/PROSPER_prs_file.txt'))$rsid,
  afr=fread(paste0('Prosper/AOU/', trait, '/PROSPER/after_ensemble_AFR/PROSPER_prs_file.txt'))$rsid,
  eas=fread(paste0('Prosper/AOU/', trait, '/PROSPER/after_ensemble_EAS/PROSPER_prs_file.txt'))$rsid,
  sas=fread(paste0('Prosper/AOU/', trait, '/PROSPER/after_ensemble_SAS/PROSPER_prs_file.txt'))$rsid,
  amr=fread(paste0('Prosper/AOU/', trait, '/PROSPER/after_ensemble_AMR/PROSPER_prs_file.txt'))$rsid
)
print(length(unique(unlist(prosper_snps))))
snps[1,'PROSPER'] = length(unique(unlist(prosper_snps))) 

split2 = split2 %>%
  left_join(prosper_split2)

rm(prosper_split2); gc()

################
### SPLENDID ###
################
print('L0L1')
aou_map = readRDS(paste0('arrays_map.rds'))
ukb_map = fread(paste0('aou_ukb_map_final.txt'))

grid = readRDS(paste0('Beta/', trait, '_L0L1_results.RDS'))

## compute PRS
beta = readRDS(paste0('Beta/', trait, '_L0L1_beta.RDS'))

beta_map = aou_map %>%
  filter(marker.ID %in% rownames(beta)) %>%
  left_join(ukb_map, by=c('marker.ID'='aou_code')) %>%
  mutate(allele_aou = paste0(allele1, allele2)) %>%
  mutate(allele_ukb = paste0(alt_ukb, ref_ukb)) %>%
  mutate(matched = strand_match(allele_allele_ukb)) %>%
  filter(matched) %>%
  select(chromosome, marker.ID, rsid, allele1, allele2)

ix_keep = which(rownames(beta) %in% beta_map$marker.ID &
                  (rowSums(beta != 0) > 0))

all.equal(beta_map$marker.ID, rownames(beta))

### write files
ix_main = sapply(1:nrow(grid), FUN=function(i) (i-1)*6 + 1)
beta_main = cbind(beta_map[ix_keep,], data.matrix(beta[ix_keep,ix_main]))
write.table(beta_main, paste0('Beta/', trait, '_beta_main.txt'),
            row.names=F, quote=F)
write.table(beta_main$rsid, paste0('Beta/', trait, '_beta_main.snp'),
            col.names=F, row.names=F, quote=F)
n_main = ncol(beta_main)
rm(beta_main); gc()

system(paste0('~/software/plink2 ',
              '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
              '--keep /n/holyscratch01/xlin/tonychen/UKB/PROSPER/ID/all.id ',
              '--extract Beta/', trait, '_beta_main.snp ',
              '--score Beta/', trait, '_beta_main.txt 3 4 header list-variants cols=+scoresums,-scoreavgs ',
              '--score-col-nums 6-', n_main, ' ',
              '--out Beta/', trait, '_full_main'))

ix_int_col = sort(c(sapply(which(grid$nb_interact != 0), FUN=function(i) (i-1)*6 + 2:6)))
ix_int_row = which(Matrix::rowSums(beta[,ix_int_col] != 0) > 0)
beta_int = cbind(beta_map[ix_int_row,], matrix(data.matrix(beta[ix_int_row,ix_int_col]), nrow=length(ix_int_row)))
write.table(beta_int, paste0('Beta/', trait, '_beta_int.txt'),
            row.names=F, quote=F)
write.table(beta_int$rsid, paste0('Beta/', trait, '_beta_int.snp'),
            col.names=F, row.names=F, quote=F)
n_int = ncol(beta_int)
rm(beta_int); gc()

system(paste0('~/software/plink2 ',
              '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
              '--keep /n/holyscratch01/xlin/tonychen/UKB/PROSPER/ID/all.id ',
              '--extract Beta/', trait, '_beta_int.snp ',
              '--score Beta/', trait, '_beta_int.txt 3 4 header list-variants cols=+scoresums,-scoreavgs ',
              '--score-col-nums 6-', n_int, ' ',
              '--out Beta/', trait, '_full_int'))

### additive PRS
prs = fread(paste0('Beta/', trait, '_full_main.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                        paste0('ADD_', 1:nrow(grid)))) %>%
  data.frame() %>%
  filter(FID %in% data$FID)

### interaction PRS
Esum = readRDS(paste0('Beta/', trait, '_Esum.RDS'))

PC_std = sapply(1:5, FUN=function(k)
  (data[,paste0('PC', k)] - Esum$mean[k]) / Esum$sd[k])
rm(Esum); gc()

# integrate PCs
prs_int = fread(paste0('Beta/', trait, '_full_int.sscore'),
                col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                            c(outer(1:5,which(grid$nb_interact != 0),
                                    FUN=function(k,i) paste0('INT', k, '_', i))))) %>%
  data.frame() %>%
  filter(FID %in% data$FID)

for (i in which(grid$nb_interact != 0)) {
  prs[,paste0('ADD_', i)] = prs[,paste0('ADD_', i)] + prs_int[,paste0('INT1_', i)] * PC_std[,1] +
    prs_int[,paste0('INT2_', i)] * PC_std[,2] + prs_int[,paste0('INT3_', i)] * PC_std[,3] +
    prs_int[,paste0('INT4_', i)] * PC_std[,4] + prs_int[,paste0('INT5_', i)] * PC_std[,5]
}
rm(prs_int); gc()

saveRDS(prs, paste0('Beta/', trait, '_L0L1_prs.RDS'))
prs = readRDS(paste0('Beta/', trait, '_L0L1_prs.RDS'))

prs1 = prs %>% filter(FID %in% split1$FID) %>% select(contains('ADD')) %>% data.matrix()
prs1a = prs %>% filter(FID %in% split1a$FID) %>% select(contains('ADD')) %>% data.matrix()
prs2 = prs %>% filter(FID %in% split2$FID) %>% select(contains('ADD')) %>% data.matrix()

## tuning / validation
grid$tun_r2 = apply(prs1, 2, FUN=function(x) cor(x, split1$res)^2)
grid$tun0_r2 = apply(prs1a, 2, FUN=function(x) cor(x, split1a$res)^2)
grid$val_r2 = apply(prs2, 2, FUN=function(x) cor(x, split2$res)^2)
saveRDS(grid, paste0('AOU/Results/', trait, '_L0L1_grid.RDS'))
grid = readRDS(paste0('AOU/Results/', trait, '_L0L1_grid.RDS'))

print('Ensemble')
ix_grid_order = grid %>%
  select(lambda0_G, lambda1_G, pval, nb_active, nb_interact, tun_r2) %>%
  mutate(ix = row_number()) %>%
  group_by(pval) %>%
  arrange(desc(tun_r2)) %>%
  mutate(pval_order = row_number()) %>%
  ungroup() %>%
  arrange(pval_order, desc(tun_r2), nb_active)

max_supp = grid %>% 
  slice_max(tun_r2) %>%
  pull(nb_active)

beta_grid = beta[,(ix_grid_order$ix[1]-1)*6 + 1]
ix_grid_keep = ix_grid_order$ix[1]
for (i in 2:nrow(ix_grid_order)) {
  beta_grid0 = cbind(beta_grid, beta[,(ix_grid_order$ix[i]-1)*6 + 1])
  if (sum(rowSums(beta_grid0 != 0) > 0) < min(10000, 2*max_supp)) {
    beta_grid = beta_grid0
    ix_grid_keep = c(ix_grid_keep, ix_grid_order$ix[i])
  }
}
ix_grid_keep = sort(ix_grid_keep)

set.seed(123)
cv = cv.glmnet(prs1[,ix_grid_keep], split1$res, alpha=0)
split2$SPLENDID = as.numeric(predict(cv, prs2[,ix_grid_keep], s='lambda.min'))

coef_ens = coef(cv, s='lambda.min')[-1]
beta_grid_ens = cbind(
  beta[,(ix_grid_keep-1)*6 + 1] %*% coef_ens,
  beta[,(ix_grid_keep-1)*6 + 2] %*% coef_ens,
  beta[,(ix_grid_keep-1)*6 + 3] %*% coef_ens,
  beta[,(ix_grid_keep-1)*6 + 4] %*% coef_ens,
  beta[,(ix_grid_keep-1)*6 + 5] %*% coef_ens,
  beta[,(ix_grid_keep-1)*6 + 6] %*% coef_ens
)
beta_grid_ens = beta_grid_ens[which(beta_grid_ens[,1] != 0),]
colSums(beta_grid_ens != 0)

saveRDS(beta_grid_ens, paste0('Results/', trait, '_beta_ridge.RDS'))

snps[1,'SPLENDID'] = sum(beta_grid_ens[,1] != 0)
snps[2,'SPLENDID'] = max(colSums(beta_grid_ens != 0)[-1])

############################
### Downsampled Ensemble ###
############################
print('Downsampled EUR Ensemble')

ix_grid_order = grid %>%
  select(lambda0_G, lambda1_G, pval, nb_active, nb_interact, tun0_r2) %>%
  mutate(ix = row_number()) %>%
  group_by(pval) %>%
  arrange(desc(tun0_r2)) %>%
  mutate(pval_order = row_number()) %>%
  ungroup() %>%
  arrange(pval_order, desc(tun0_r2), nb_active)

max_supp = grid %>% 
  slice_max(tun0_r2) %>%
  pull(nb_active)

beta_grid = beta[,(ix_grid_order$ix[1]-1)*6 + 1]
ix_grid_keep = ix_grid_order$ix[1]
for (i in 2:nrow(ix_grid_order)) {
  beta_grid0 = cbind(beta_grid, beta[,(ix_grid_order$ix[i]-1)*6 + 1])
  if (sum(rowSums(beta_grid0 != 0) > 0) < min(10000, 2*max_supp)) {
    beta_grid = beta_grid0
    ix_grid_keep = c(ix_grid_keep, ix_grid_order$ix[i])
  }
}
ix_grid_keep = sort(ix_grid_keep)

set.seed(123)
cv = cv.glmnet(prs1a[,ix_grid_keep], split1a$res, alpha=0)
split2$SPLENDID_sub = as.numeric(predict(cv, prs2[,ix_grid_keep], s='lambda.min'))

rm(cv); gc()

############################
### L0L1 (External GWAS) ###
############################
print('L0L1 Ext')
aou_map = readRDS(paste0('arrays_map.rds'))
ukb_map = fread(paste0('aou_ukb_map_final.txt'))

grid = readRDS(paste0('Beta/', trait, '_L0L1_ext_results.RDS'))

## compute PRS
beta = readRDS(paste0('Beta/', trait, '_L0L1_ext_beta.RDS'))

beta_map = aou_map %>%
  filter(marker.ID %in% rownames(beta)) %>%
  left_join(ukb_map, by=c('marker.ID'='aou_code')) %>%
  mutate(allele_aou = paste0(allele1, allele2)) %>%
  mutate(allele_ukb = paste0(alt_ukb, ref_ukb)) %>%
  mutate(matched = strand_match(allele_allele_ukb)) %>%
  filter(matched) %>%
  select(chromosome, marker.ID, rsid, allele1, allele2)

ix_keep = which(rownames(beta) %in% beta_map$marker.ID &
                  (rowSums(beta != 0) > 0))

all.equal(beta_map$marker.ID, rownames(beta))

### write files
ix_main = sapply(1:nrow(grid), FUN=function(i) (i-1)*6 + 1)
beta_main = cbind(beta_map[ix_keep,], data.matrix(beta[ix_keep,ix_main]))
write.table(beta_main, paste0('Beta/', trait, '_ext_beta_main.txt'),
            row.names=F, quote=F)
write.table(beta_main$rsid, paste0('Beta/', trait, '_ext_beta_main.snp'),
            col.names=F, row.names=F, quote=F)
n_main = ncol(beta_main)
rm(beta_main); gc()

system(paste0('~/software/plink2 ',
              '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
              '--keep /n/holyscratch01/xlin/tonychen/UKB/PROSPER/ID/ALL.id ',
              '--extract Beta/', trait, '_ext_beta_main.snp ',
              # '--chr 22 ',
              '--score Beta/', trait, '_ext_beta_main.txt 3 4 header list-variants cols=+scoresums,-scoreavgs ',
              '--score-col-nums 6-', n_main, ' ',
              '--out Beta/', trait, '_ext_full_main'))

ix_int_col = sort(c(sapply(which(grid$nb_interact != 0), FUN=function(i) (i-1)*6 + 2:6)))
ix_int_row = which(Matrix::rowSums(beta[,ix_int_col] != 0) > 0)
beta_int = cbind(beta_map[ix_int_row,], matrix(data.matrix(beta[ix_int_row,ix_int_col]), nrow=length(ix_int_row)))
write.table(beta_int, paste0('Beta/', trait, '_ext_beta_int.txt'),
            row.names=F, quote=F)
write.table(beta_int$rsid, paste0('Beta/', trait, '_ext_beta_int.snp'),
            col.names=F, row.names=F, quote=F)
n_int = ncol(beta_int)
rm(beta_int); gc()

system(paste0('~/software/plink2 ',
              '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
              '--keep /n/holyscratch01/xlin/tonychen/UKB/PROSPER/ID/ALL.id ',
              '--extract Beta/', trait, '_ext_beta_int.snp ',
              # '--chr 22 ',
              '--score Beta/', trait, '_ext_beta_int.txt 3 4 header list-variants cols=+scoresums,-scoreavgs ',
              '--score-col-nums 6-', n_int, ' ',
              '--out Beta/', trait, '_ext_full_int'))

### additive PRS
prs = fread(paste0('Beta/', trait, '_ext_full_main.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                        paste0('ADD_', 1:nrow(grid)))) %>%
  data.frame() %>%
  filter(FID %in% data$FID)

### interaction PRS
Esum = readRDS(paste0('Beta/', trait, '_Esum.RDS'))

PC_std = sapply(1:5, FUN=function(k)
  (data[,paste0('PC', k)] - Esum$mean[k]) / Esum$sd[k])
rm(Esum); gc()

# integrate PCs
prs_int = fread(paste0('Beta/', trait, '_ext_full_int.sscore'),
                col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                            c(outer(1:5,which(grid$nb_interact != 0),
                                    FUN=function(k,i) paste0('INT', k, '_', i))))) %>%
  data.frame() %>%
  filter(FID %in% data$FID)

for (i in which(grid$nb_interact != 0)) {
  prs[,paste0('ADD_', i)] = prs[,paste0('ADD_', i)] + prs_int[,paste0('INT1_', i)] * PC_std[,1] +
    prs_int[,paste0('INT2_', i)] * PC_std[,2] + prs_int[,paste0('INT3_', i)] * PC_std[,3] +
    prs_int[,paste0('INT4_', i)] * PC_std[,4] + prs_int[,paste0('INT5_', i)] * PC_std[,5]
}
rm(prs_int); gc()

saveRDS(prs, paste0('Beta/', trait, '_L0L1_ext_prs.RDS'))
prs = readRDS(paste0('Beta/', trait, '_L0L1_ext_prs.RDS'))

prs1 = prs %>% filter(FID %in% split1$FID) %>% select(contains('ADD')) %>% data.matrix()
prs2 = prs %>% filter(FID %in% split2$FID) %>% select(contains('ADD')) %>% data.matrix()

## tuning / validation
grid$tun_r2 = apply(prs1, 2, FUN=function(x) cor(x, split1$res)^2)
grid$val_r2 = apply(prs2, 2, FUN=function(x) cor(x, split2$res)^2)
saveRDS(grid, paste0('AOU/Results/', trait, '_L0L1_ext_grid.RDS'))

grid = readRDS(paste0('AOU/Results/', trait, '_L0L1_ext_grid.RDS'))

print('Ensemble Ext')

# full ensemble (ridge / lasso) w/ limited support
ix_grid_order = grid %>%
  select(lambda0_G, lambda1_G, pval, nb_active, nb_interact, tun_r2) %>%
  mutate(ix = row_number()) %>%
  group_by(pval) %>%
  arrange(desc(tun_r2)) %>%
  mutate(pval_order = row_number()) %>%
  ungroup() %>%
  arrange(pval_order, desc(tun_r2), nb_active)

max_supp = grid %>% 
  slice_max(tun_r2) %>%
  pull(nb_active)

beta_grid = beta[,(ix_grid_order$ix[1]-1)*6 + 1]
ix_grid_keep = ix_grid_order$ix[1]
for (i in 2:nrow(ix_grid_order)) {
  beta_grid0 = cbind(beta_grid, beta[,(ix_grid_order$ix[i]-1)*6 + 1])
  if (sum(rowSums(beta_grid0 != 0) > 0) < min(10000, 2*max_supp)) {
    beta_grid = beta_grid0
    ix_grid_keep = c(ix_grid_keep, ix_grid_order$ix[i])
  }
}
ix_grid_keep = sort(ix_grid_keep)

set.seed(123)
cv = cv.glmnet(prs1[,ix_grid_keep], split1$res, alpha=0)
split2$SPLENDID_ext = as.numeric(predict(cv, prs2[,ix_grid_keep], s='lambda.min'))

############
### iPGS ###
############
print('Lasso')
## compute PRS
beta = readRDS(paste0('Beta/', trait, '_lasso_beta.RDS'))
grid = readRDS(paste0('Beta/', trait, '_lasso_results.RDS'))

beta_map = aou_map %>%
  filter(marker.ID %in% rownames(beta)) %>%
  left_join(ukb_map, by=c('marker.ID'='aou_code')) %>%
  mutate(allele_aou = paste0(allele1, allele2)) %>%
  mutate(allele_ukb = paste0(alt_ukb, ref_ukb)) %>%
  mutate(matched = strand_match(allele_allele_ukb)) %>%
  filter(matched) %>%
  select(chromosome, marker.ID, rsid, allele1, allele2)

ix_keep = which(rownames(beta) %in% beta_map$marker.ID)
beta = beta[ix_keep,]

all.equal(beta_map$marker.ID, rownames(beta))

### write files
beta_main = cbind(beta_map, data.matrix(beta))
write.table(beta_main, paste0('Beta/', trait, '_lasso_beta_main.txt'),
            row.names=F, quote=F)
write.table(beta_main$rsid, paste0('Beta/', trait, '_lasso_beta_main.snp'),
            col.names=F, row.names=F, quote=F)

n_main = ncol(beta_main)
rm(beta_main); gc()

system(paste0('~/software/plink2 --silent ',
              '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
              '--keep /n/holyscratch01/xlin/tonychen/UKB/PROSPER/ID/all.id ',
              '--extract Beta/', trait, '_lasso_beta_main.snp ',
              '--score Beta/', trait, '_lasso_beta_main.txt 3 4 header list-variants cols=+scoresums,-scoreavgs ',
              '--score-col-nums 6-', n_main, ' ',
              '--out Beta/', trait, '_lasso_full_main'))

### additive PRS
prs = fread(paste0('Beta/', trait, '_lasso_full_main.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                        paste0('ADD_', 1:nrow(grid)))) %>%
  data.frame() %>% 
  filter(FID %in% data$FID)

prs1 = prs %>% filter(FID %in% split1$FID) %>% select(contains('ADD')) %>% data.matrix()
prs1a = prs %>% filter(FID %in% split1a$FID) %>% select(contains('ADD')) %>% data.matrix()
prs2 = prs %>% filter(FID %in% split2$FID) %>% select(contains('ADD')) %>% data.matrix()

## tuning / validation
grid$tun_r2 = apply(prs1, 2, FUN=function(x) cor(x, split1$res)^2)
grid$tun0_r2 = apply(prs1a, 2, FUN=function(x) cor(x, split1a$res)^2)
grid$val_r2 = apply(prs2, 2, FUN=function(x) cor(x, split2$res)^2)

saveRDS(grid, paste0('AOU/Results/', trait, '_Lasso_grid.RDS'))
grid= readRDS(paste0('AOU/Results/', trait, '_Lasso_grid.RDS'))

split2$PRS_lasso = prs2[,grid %>% mutate(ix=row_number()) %>% slice_max(tun_r2) %>% head(1) %>% pull(ix)]
split2$PRS_lasso0 = prs2[,grid %>% mutate(ix=row_number()) %>% slice_max(tun0_r2) %>% head(1) %>% pull(ix)]

grid %>%
  slice_max(tun_r2) %>%
  head(1)

snps[1,'iPGS'] = grid %>% slice_max(tun_r2) %>% head(1) %>% pull(nb_active)

rm(beta); gc()

####################
### iPGS + Refit ###
####################
print('iPGS w/ Refit')

grid = readRDS(paste0('Results/', trait, '_Lasso_grid.RDS'))

prs = fread(paste0('Beta/', trait, '_lasso_full_main.sscore'),
            col.names=c('FID', 'IID', 'NALLELE', 'DOSAGE',
                        paste0('ADD_', 1:nrow(grid)))) %>%
  data.frame() %>% 
  filter(FID %in% data$FID)

supp_iPGS = grid %>% slice_max(tun_r2) %>% slice_min(nb_active) %>% pull(nb_active)
iPGS_val = prs[,c(1, 2, which.max(grid$tun_r2)+4)]; names(iPGS_val)[3] = 'iPGS'

## compute PRS
main_snps = c()
het_snps = c()
iPGS_coef = rep(0, 5)
names(iPGS_coef) = c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')
for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
  beta_refit = readRDS(paste0('Beta/iPGS_refit_coef_', trait, '_', tolower(anc), '.RDS')) %>%
    data.matrix() %>%
    data.frame() %>%
    tibble::rownames_to_column('Var') %>%
    filter(grepl('rs', Var) | Var == 'iPGS') %>%
    separate(Var, into=c('V1', 'V2'), sep=':') %>%
    mutate(SNP = ifelse(V1 %in% c('PC1', 'PC2'), V2, V1),
           PC = ifelse(V1 %in% c('PC1', 'PC2'), V1, 'ADD')) %>%
    select(SNP, PC, s1) %>%
    spread(key=PC, value=s1) %>%
    filter(ADD != 0 | PC1 != 0 | PC2 != 0)
  beta_refit[is.na(beta_refit)] = 0
  
  beta_map = aou_map %>%
    left_join(ukb_map, by=c('marker.ID'='aou_code')) %>%
    filter(rsid %in% beta_refit$SNP) %>%
    select(chromosome, marker.ID, rsid, allele1, allele2) %>%
    left_join(beta_refit, by=c('rsid'='SNP'))
  
  if ('iPGS' %in% beta_refit$SNP) iPGS_coef[anc] = beta_refit$ADD[beta_refit$SNP == 'iPGS']
  
  write.table(beta_map, paste0('Beta/iPGS_refit_coef_', trait, '_', anc, '.txt'),
              row.names=F, quote=F)
  
  main_snps = c(main_snps, beta_map$rsid)
  het_snps = c(het_snps, beta_map$rsid[beta_map$PC1 != 0 | beta_map$PC2 != 0])
}
main_snps = unique(main_snps)
het_snps = unique(het_snps)
write.table(main_snps, paste0('Beta/iPGS_refit_coef_', trait, '.snp'),
            row.names=F, quote=F, col.names=F)

for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
  system(paste0('~/software/plink2 ',
                '--silent ',
                '--bfile /n/holyscratch01/xlin/tonychen/Lasso/aou_common ',
                '--keep ', 'ID/split2.id ',
                '--extract Beta/iPGS_refit_coef_', trait, '.snp ',
                '--score Beta/iPGS_refit_coef_', trait, '_', anc, '.txt 3 4 header-read list-variants cols=+scoresums,-scoreavgs ',
                '--score-col-nums 6-8 ',
                '--out Beta/', trait, '_iPGS_refit_', anc))
}

### interaction PRS
Esum = readRDS(paste0('Beta/', trait, '_Esum.RDS'))

# integrate PCs
split2$iPGS_refit = 0
for (anc in c('EUR', 'AFR', 'EAS', 'SAS', 'AMR')) {
  ix_IID = split2 %>% filter(ancestry == anc) %>% pull(IID)
  
  PC_std = sapply(1:2, FUN=function(k)
    (data[(data$IID %in% ix_IID),paste0('PC', k)] - Esum$mean[k]) / Esum$sd[k])
  
  if (file.exists(paste0('Beta/', trait, '_iPGS_refit_', anc, '.sscore'))) {
    prs_anc = fread(paste0('Beta/', trait, '_iPGS_refit_', anc, '.sscore')) %>%
      filter(IID %in% ix_IID)
    
    split2$PRS_refit[split2$IID %in% ix_IID] = split2$PRS_lasso[split2$IID %in% ix_IID]*iPGS_coef[anc] + 
      prs_anc$ADD_SUM +
      prs_anc$PC1_SUM*PC_std[,1] + prs_anc$PC2_SUM*PC_std[,2]
  } else {
    split2$iPGS_refit[split2$IID %in% ix_IID] = split2$PRS_lasso[split2$IID %in% ix_IID]
  }
  
}

rm(prs_anc, Esum, PC_std); gc()

iPGS_beta0 = fread(paste0('Beta/', trait, '_iPGS.txt'))

snps[1,'iPGS + Refit'] = length(unique(c(iPGS_beta0$rsid, main_snps)))
snps[2,'iPGS + Refit'] = length(het_snps)

rm(beta); gc()

###############
### RESULTS ###
###############

low_prob90 = ancestry %>%
  filter(pmax(P.AFR, P.EUR, P.EAS, P.AMR, P.SAS) < 0.9) %>%
  pull(IID)

low_prob70 = ancestry %>%
  filter(pmax(P.AFR, P.EUR, P.EAS, P.AMR, P.SAS) < 0.7) %>%
  pull(IID)

cat(trait, ' START\n')
split2 = silent(readRDS(paste0(multi, 'AOU/Results/', trait, '_validation.RDS')) %>%
                  left_join(cov) %>% 
                  left_join(aou_pc) %>%
                  data.frame())

# residualize PRS
for (prs in c('PRSCSX', 'CTSLEB', 'PROSPER', 'iPGS_refit', 'SPLENDID', 'SPLENDID_sub', 'SPLENDID_ext')) {
  split2[,paste0(prs, '_res')] = residuals(lm(formula(paste0(prs, '~age+female+', paste(paste0('PC', 1:20), collapse='+'))), data=split2))
}

# bootstrap
val_data = split2 %>%
  select(IID, ancestry, res, 
         PRSCSX, PRSCSX_res,
         CTSLEB, CTSLEB_res, 
         PROSPER, PROSPER_res,
         contains('PRS_'))

main_r2 = silent(val_data %>% 
                   bind_rows(., val_data %>% mutate(ancestry = 'ALL')) %>%
                   bind_rows(., val_data %>% filter(ancestry != 'EUR') %>% mutate(ancestry='nonEUR')) %>%
                   bind_rows(., val_data %>% filter(IID %in% low_prob90) %>% mutate(ancestry='Other90')) %>%
                   bind_rows(., val_data %>% filter(IID %in% low_prob70) %>% mutate(ancestry='Other70')) %>%
                   mutate(ancestry = factor(ancestry, 
                                            levels=c('ALL', 'EUR', 'nonEUR', 'Other90', 'Other70', 'AFR', 'AMR', 'EAS', 'SAS'))) %>%
                   group_by(ancestry) %>%
                   summarize(n=n(),
                             ## original PRS
                             # GWAS
                             R2_prscsx = cor(PRSCSX, res)^2,
                             R2_ctsleb = cor(CTSLEB, res)^2,
                             R2_prosper = cor(PROSPER, res)^2,
                             # individual-level
                             R2_iPGS = cor(iPGS_refit, res)^2,
                             R2_SPLENDID=cor(SPLENDID, res)^2,
                             # downsample
                             R2_SPLENDID.sub=cor(SPLENDID_sub, res)^2,
                             # external
                             R2_SPLENDID.ext=cor(SPLENDID_ext, res)^2,
                             ## residualized PRS (partial correlation)
                             # GWAS
                             R2_prscsx = cor(PRSCSX_res, res)^2,
                             R2_ctsleb = cor(CTSLEB_res, res)^2,
                             R2_prosper = cor(PROSPER_res, res)^2,
                             # individual-level
                             R2_iPGS = cor(iPGS_refit_res, res)^2,
                             R2_SPLENDID=cor(SPLENDID_res, res)^2,
                             # downsample
                             R2_SPLENDID.sub=cor(SPLENDID_sub_res, res)^2,
                             # external
                             R2_SPLENDID.ext=cor(SPLENDID_ext_res, res)^2,
                   data.frame())

set.seed(123)
B = 10000
boot_mat = matrix(sample(1:nrow(val_data), nrow(val_data)*B, replace=T), ncol=B)
boot_r2 = NULL
# boot_r2 = foreach(b=1:B, .combine = 'rbind') %dopar% {
for (b in 1:B) {
  if (b %% 100 == 0) cat(paste0(trait, b), '.')
  if (b %% 1000 == 0) cat('\n')
  # ix_boot = sample(1:nrow(val_data), nrow(val_data), replace=T)
  ix_boot = boot_mat[,b]
  
  boot_r2_b = silent(val_data[ix_boot,] %>% 
                       bind_rows(., val_data[ix_boot,] %>% mutate(ancestry = 'ALL')) %>%
                       bind_rows(., val_data[ix_boot,] %>% filter(ancestry != 'EUR') %>% mutate(ancestry='nonEUR')) %>%
                       bind_rows(., val_data[ix_boot,] %>% filter(IID %in% low_prob90) %>% mutate(ancestry='Other90')) %>%
                       bind_rows(., val_data[ix_boot,] %>% filter(IID %in% low_prob70) %>% mutate(ancestry='Other70')) %>%
                       mutate(ancestry = factor(ancestry, 
                                                levels=c('ALL', 'EUR', 'nonEUR', 'Other90', 'Other70', 'AFR', 'AMR', 'EAS', 'SAS'))) %>%
                       group_by(ancestry) %>%
                       summarize(n=n(),
                                 ## original PRS
                                 # GWAS
                                 R2_prscsx = cor(PRSCSX, res)^2,
                                 R2_ctsleb = cor(CTSLEB, res)^2,
                                 R2_prosper = cor(PROSPER, res)^2,
                                 # individual-level
                                 R2_iPGS = cor(iPGS_refit, res)^2,
                                 R2_SPLENDID=cor(SPLENDID, res)^2,
                                 # downsample
                                 R2_SPLENDID.sub=cor(SPLENDID_sub, res)^2,
                                 # external
                                 R2_SPLENDID.ext=cor(SPLENDID_ext, res)^2,
                                 ## residualized PRS (partial correlation)
                                 # GWAS
                                 R2_prscsx = cor(PRSCSX_res, res)^2,
                                 R2_ctsleb = cor(CTSLEB_res, res)^2,
                                 R2_prosper = cor(PROSPER_res, res)^2,
                                 # individual-level
                                 R2_iPGS = cor(iPGS_refit_res, res)^2,
                                 R2_SPLENDID=cor(SPLENDID_res, res)^2,
                                 # downsample
                                 R2_SPLENDID.sub=cor(SPLENDID_sub_res, res)^2,
                                 # external
                                 R2_SPLENDID.ext=cor(SPLENDID_ext_res, res)^2,
                       data.frame() %>%
                       mutate(b = b))
  
  boot_r2 = silent(rbind(boot_r2, boot_r2_b))
  rm(boot_r2_b); gc()
}
cat('\n')

saveRDS(list(main = main_r2,
             boot = boot_r2), 
        paste0('AOU/Results/', trait, '_boot_results_revision.RDS'))


cat(trait, ' DONE\n')
rm(main_r2, boot_r2, val_data, val, split2); gc()
