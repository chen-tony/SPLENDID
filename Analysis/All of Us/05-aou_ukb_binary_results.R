library(dplyr)
library(data.table)
library(glmnet)
library(tidyr)
library(ROCnReg)

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

AUC = function(Var, Pheno, Cov, Data) {
  if (Data %>% pull(Var) %>% sd() > 0) {
    partial_auc = AROC.sp(formula.h = formula(paste0(Var, '~', paste(Cov, collapse='+'))), 
                          group=Pheno, tag.h=0, data=Data, B=1)
    return(c(data.matrix(partial_auc$AUC[1])))
  } else {
    return(0.5)
  }
}

silent = function(x) suppressMessages(suppressWarnings(x))

############
### DATA ###
############
print('data')

split2_id = fread(paste0('split2.id'))
all_id = fread(paste0('all.id'))

args = commandArgs(trailingOnly = T)
trait = args[[1]]

pheno = fread(paste0('ukb_multi.pheno')) %>%
  select(FID, IID, Y=ends_with(trait))
cov = fread(paste0('ukb_multi.cov')) %>%
  select(FID, IID, age, female)
anc = fread(paste0('ukb_multi.anc')) %>%
  select(FID, IID, ancestry=predicted, contains('P.'))
aou_pc = fread('ukb_aou_pca.sscore') %>%
  select(FID=`#FID`, IID, PC=contains('SCORE'))

data = pheno %>%
  filter(FID %in% all_id$FID) %>%
  left_join(anc) %>%
  left_join(cov) %>%
  left_join(aou_pc) %>%
  na.omit() %>%
  data.frame()

if (trait %in% c('Breast', 'Prostate')) {
  split2 = data %>%
    # filter(FID %in% split2_id$FID) %>%
    select(FID, IID, ancestry, Y, age, paste0('PC', 1:20), contains('P.'))
} else {
  split2 = data %>%
    # filter(FID %in% split2_id$FID) %>%
    select(FID, IID, ancestry, Y, age, female, paste0('PC', 1:20), contains('P.'))
}
if (sd(split2$female) == 0) split2$female = NULL

snps = matrix(0, nrow=2, ncol=3)
rownames(snps) = c('main', 'int')
colnames(snps) = c('PRSCSX', 'SPLENDID (Linear)', 'SPENDID (Logistic)')

###################
### PRSCSX-auto ###
###################
print('PRSCSX-auto')

prscsx_beta = data.frame()
prscsx_snps = c()
for (chr in 1:22) {
  prscsx_chr = fread(paste0('PRSCSX/Results_AUTO/', trait,
                            '_META_pst_eff_a1_b0.5_phiauto_chr', chr, '.txt'),
                     col.names=c('chr', 'rsid', 'pos', 'a1', 'a0', 'beta'))
  prscsx_snps = c(prscsx_snps, prscsx_chr$rsid)
  prscsx_beta = rbind(prscsx_beta, prscsx_chr)
}

fwrite(prscsx_beta, paste0('PRSCSX/prscsx_auto_', trait, '_beta.txt'), sep='\t')
write.table(prscsx_snps, paste0('PRSCSX/prscsx_auto_', trait, '_snps.txt'),
            row.names=F, col.names=F, quote=F)
rm(prscsx_beta, prscsx_chr); gc()

prscsx_snps = read.table(paste0('PRSCSX/prscsx_auto_', trait, '_snps.txt'))$V1
snps[1,'PRSCSX'] = length(prscsx_snps)

system(paste0('~/software/plink2 ',
              '--bfile aou_common ',
              '--extract PRSCSX/prscsx_auto_', trait, '_snps.txt ',
              # '--keep split2.id ',
              '--keep all.id ',
              '--score ',  'PRSCSX/prscsx_auto_', trait, '_beta.txt 2 4 6 header ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--out PRSCSX/prscsx_auto_', trait))

prscsx_prs = fread(paste0('PRSCSX/prscsx_auto_', trait, '.sscore'),
                   col.names=c('FID', 'IID', 'PRSCSX')) %>%
  data.frame() %>%
  filter(FID %in% split2$FID)

split2 = left_join(split2, prscsx_prs)

#########################
### SPLENDID (Linear) ###
#########################
# print('L0L1')
aou_map = readRDS(paste0('arrays_map.rds'))
ukb_map = fread(paste0('aou_ukb_map_final.txt'))

## compute PRS
beta = readRDS(paste0('Beta/Linear/', trait, '_beta_cv.RDS'))

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

ix_match = match(beta_map$marker.ID, rownames(beta))
beta = beta[ix_match,]

all.equal(beta_map$marker.ID, rownames(beta))

### write files
colnames(beta) = c('ADD', paste0('INT', 1:5))
beta_main = cbind(beta_map, data.matrix(beta))
write.table(beta_main, paste0('Beta/Linear/', trait, '_beta_auto.txt'),
            row.names=F, quote=F)
write.table(beta_main$rsid, paste0('Beta/Linear/', trait, '_beta_auto.snp'),
            col.names=F, row.names=F, quote=F)
n_main = ncol(beta_main)
beta_main = read.table(paste0('Beta/Linear/', trait, '_beta_auto.txt'), header=T)

snps[1,'SPLENDID (Linear)'] = sum(abs(beta_main$ADD) > 1e-20)
snps[2,'SPLENDID (Linear)'] = sum(abs(beta_main$INT1) > 1e-20  | abs(beta_main$INT2) > 1e-20 |
                       abs(beta_main$INT3) > 1e-20 | abs(beta_main$INT4) > 1e-20 |
                       abs(beta_main$INT5) > 1e-20 )

system(paste0('~/software/plink2 ',
              '--bfile aou_common ',
              '--keep all.id ',
              '--extract Beta/Linear/', trait, '_beta_auto.snp ',
              '--score Beta/Linear/', trait, '_beta_auto.txt 3 4 header list-variants ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 6-11 ',
              '--out Beta/Linear/', trait, '_full_auto'))

# PRS
prs = fread(paste0('Beta/Linear/', trait, '_full_auto.sscore'),
            col.names=c('FID', 'IID', 'ADD', paste0('INT', 1:5)))%>%
  data.frame() %>%
  filter(FID %in% split2$FID)

Esum = readRDS(paste0('Beta/', trait, '_Esum.RDS'))
PC_std = sapply(1:5, FUN=function(k)
  (data[which(data$FID %in% split2$FID), paste0('PC', k)] - Esum$mean[k]) / Esum$sd[k])

prs$SPLENDID_Linear = prs$ADD + prs$INT1 * PC_std[,1] +
  prs$INT2 * PC_std[,2] + prs$INT3 * PC_std[,3] +
  prs$INT4 * PC_std[,4] + prs$INT5 * PC_std[,5]


###########################
### SPLENDID (Logistic) ###
###########################
# print('L0L1')
## compute PRS
beta = readRDS(paste0('Beta/Logistic/', trait, '_beta_cv.RDS'))

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

ix_match = match(beta_map$marker.ID, rownames(beta))
beta = beta[ix_match,]

all.equal(beta_map$marker.ID, rownames(beta))

### write files
colnames(beta) = c('ADD', paste0('INT', 1:5))
beta_main = cbind(beta_map, data.matrix(beta))
write.table(beta_main, paste0('Beta/Logistic/', trait, '_beta_auto.txt'),
            row.names=F, quote=F)
write.table(beta_main$rsid, paste0('Beta/', trait, '_beta_auto.snp'),
            col.names=F, row.names=F, quote=F)
n_main = ncol(beta_main)
beta_main = read.table(paste0('Beta/Logistic/', trait, '_beta_auto.txt'), header=T)

snps[1,'SPLENDID (Logistic)'] = sum(abs(beta_main$ADD) > 1e-20)
snps[2,'SPLENDID (Logistic)'] = sum(abs(beta_main$INT1) > 1e-20  | abs(beta_main$INT2) > 1e-20 |
                       abs(beta_main$INT3) > 1e-20 | abs(beta_main$INT4) > 1e-20 |
                       abs(beta_main$INT5) > 1e-20 )

system(paste0('~/software/plink2 ',
              '--bfile aou_common ',
              '--keep all.id ',
              '--extract Beta/Logistic/', trait, '_beta_auto.snp ',
              '--score Beta/Logistic/', trait, '_beta_auto.txt 3 4 header list-variants ',
              'cols=+scoresums,-scoreavgs,-dosagesum,-nallele ',
              '--score-col-nums 6-11 ',
              '--out Beta/Logistic/', trait, '_full_auto'))

# PRS
prs = fread(paste0('Beta/Logistic/', trait, '_full_auto.sscore'),
            col.names=c('FID', 'IID', 'ADD', paste0('INT', 1:5)))%>%
  data.frame() %>%
  filter(FID %in% split2$FID)

Esum = readRDS(paste0('Beta/', trait, '_Esum.RDS'))
PC_std = sapply(1:5, FUN=function(k)
  (data[which(data$FID %in% split2$FID), paste0('PC', k)] - Esum$mean[k]) / Esum$sd[k])

prs$SPLENDID_Logistic = prs$ADD + prs$INT1 * PC_std[,1] +
  prs$INT2 * PC_std[,2] + prs$INT3 * PC_std[,3] +
  prs$INT4 * PC_std[,4] + prs$INT5 * PC_std[,5]

saveRDS(prs, paste0('Results/', trait, '_auto_validation.RDS'))

############
### DATA ###
############
print('data')

split2_id = fread(paste0('split2.id'))
all_id = fread(paste0('all.id'))

args = commandArgs(trailingOnly = T)
trait = args[[1]]
# trait = 'T2D'

cat(trait, 'START \n')
val_data = silent(readRDS(paste0('Results/', trait, '_auto_validation.RDS')) %>%
                    select(FID, IID, ancestry, Y, PRSCSX_AUTO=PRSCSX, SPLENDID_Linear, SPLENDID_Logistic) %>%
                    left_join(cov) %>%
                    left_join(aou_pc))

for (prs in c('PRSCSX_AUTO', 'SPLENDID_Linear', 'SPLENDID_Logistic')) {
  val_data[,paste0(prs, '_res')] = residuals(lm(formula(paste0(prs, '~age+female+', paste(paste0('PC', 1:20), collapse='+'))), data=val_data))
}

low_prob90 = ancestry %>%
  filter(pmax(P.AFR, P.EUR, P.EAS, P.AMR, P.SAS) < 0.9) %>%
  pull(IID)

low_prob70 = ancestry %>%
  filter(pmax(P.AFR, P.EUR, P.EAS, P.AMR, P.SAS) < 0.7) %>%
  pull(IID)

main_auc = silent(val_data %>%
                    bind_rows(., val_data %>% mutate(ancestry = 'ALL')) %>%
                    bind_rows(., val_data %>% filter(ancestry != 'EUR') %>% mutate(ancestry='nonEUR')) %>%
                    bind_rows(., val_data %>% filter(IID %in% low_prob90) %>% mutate(ancestry='Other90')) %>%
                    bind_rows(., val_data %>% filter(IID %in% low_prob70) %>% mutate(ancestry='Other70')) %>%
                    mutate(ancestry = factor(ancestry, 
                                             levels=c('ALL', 'EUR', 'nonEUR','Other90', 'Other70',  'AFR', 'AMR', 'EAS', 'SAS'))) %>%
                    group_by(ancestry) %>%
                    summarize(n=n(), case = sum(Y), control = sum(1-Y)) %>%
                    data.frame() %>% 
                    ungroup())
main_auc[,c('PRSCSX_AUTO', 'SPLENDID_Linear', 'SPLENDID_Logistic', 'PRSCSX_AUTO_res', 'SPLENDID_Linear_res', 'SPLENDID_Logistic_res')] = NA

for (pop in c('ALL', 'nonEUR', 'Other90', 'Other70', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS')) {
  if (pop == 'ALL') {
    val = val_data
  } else if (pop == 'nonEUR') {
    val = val_data %>% filter(ancestry != 'EUR')
  } else if (pop == 'Other90') {
    val = val_data %>% filter(IID %in% low_prob90)
  } else if (pop == 'Other70') {
    val = val_data %>% filter(IID %in% low_prob70)
  } else {
    val = val_data %>% filter(ancestry == pop)
  }
  
  for (prs in c('PRSCSX_AUTO', 'L0L1_AUTO', 'L0L1_AUTOLOG', 'PRSCSX_AUTO_res', 'L0L1_AUTO_res', 'L0L1_AUTOLOG_res')) {
    auc = NA
    if (trait %in% c('Breast', 'Prostate')) {
      auc = AUC(Var=prs, Pheno='Y', Cov=c('age', paste0('PC', 1:20)), Data=val)
    } else {
      auc = AUC(Var=prs, Pheno='Y', Cov=c('age', 'female', paste0('PC', 1:20)), Data=val)
    }
    
    if (!'try-error' %in% class(auc)) {
      main_auc[main_auc$ancestry==pop, prs] = auc
    }
  }
}

set.seed(123)
B = 1000

boot_mat = matrix(sample(1:nrow(val_data), nrow(val_data)*B, replace=T), ncol=B)

# boot_auc = NULL
# for (b in 1:B) {
boot_auc = foreach(b=1:B, .combine = 'rbind') %dopar% {
  if (b %% 10 == 0) cat(paste0(trait, b), '.')
  if (b %% 100 == 0) cat('\n')
  # ix_boot = sample(1:nrow(val_data), nrow(val_data), replace=T)
  ix_boot = boot_mat[,b]
  
  boot = silent(val_data[ix_boot,] %>%
                  bind_rows(., val_data[ix_boot,] %>% mutate(ancestry = 'ALL')) %>%
                  bind_rows(., val_data[ix_boot,] %>% filter(ancestry != 'EUR') %>% mutate(ancestry='nonEUR')) %>%
                  bind_rows(., val_data[ix_boot,] %>% filter(IID %in% low_prob90) %>% mutate(ancestry='Other90')) %>%
                  bind_rows(., val_data[ix_boot,] %>% filter(IID %in% low_prob70) %>% mutate(ancestry='Other70')) %>%
                  mutate(ancestry = factor(ancestry, 
                                           levels=c('ALL', 'EUR', 'nonEUR','Other90', 'Other70',  'AFR', 'AMR', 'EAS', 'SAS'))) %>%
                  group_by(ancestry) %>%
                  summarize(n=n(), case = sum(Y), control = sum(1-Y)) %>%
                  data.frame() %>% 
                  ungroup() %>%
                  mutate(b=b))
  boot[,c('PRSCSX_AUTO', 'L0L1_AUTO', 'L0L1_AUTOLOG', 'PRSCSX_AUTO_res', 'L0L1_AUTO_res', 'L0L1_AUTOLOG_res')] = NA
  
  for (pop in c('ALL', 'nonEUR', 'Other90', 'Other70', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS')) {
    if (pop == 'ALL') {
      val = val_data[ix_boot,]
    } else if (pop == 'nonEUR') {
      val = val_data[ix_boot,] %>% filter(ancestry != 'EUR')
    } else if (pop == 'Other90') {
      val = val_data[ix_boot,] %>% filter(IID %in% low_prob90)
    } else if (pop == 'Other70') {
      val = val_data[ix_boot,] %>% filter(IID %in% low_prob70)
    } else {
      val = val_data[ix_boot,] %>% filter(ancestry == pop)
    }
    
    if (length(unique(val$Y)) == 2) {
      for (prs in c('PRSCSX_AUTO', 'L0L1_AUTO', 'L0L1_AUTOLOG', 'PRSCSX_AUTO_res', 'L0L1_AUTO_res', 'L0L1_AUTOLOG_res')) {
        auc = NA
        if (trait %in% c('Breast', 'Prostate')) {
          auc = AUC(Var=prs, Pheno='Y', Cov=c('age', paste0('PC', 1:20)), Data=val)
        } else {
          auc = AUC(Var=prs, Pheno='Y', Cov=c('age', 'female', paste0('PC', 1:20)), Data=val)
        }
        
        if (!'try-error' %in% class(auc)) {
          boot[boot$ancestry==pop, prs] = auc
        }
      }
    }
  }
  # gc()
  
  # boot_auc = rbind(boot_auc, boot)
  
  boot
}
cat('\n')

cat(trait, 'DONE \n')

saveRDS(list(main=main_auc, boot=boot_auc), 
        paste0('AOU/Results/', trait, '_auto_boot_results_revisions.RDS'))


