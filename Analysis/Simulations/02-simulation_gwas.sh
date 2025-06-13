# rho = percentage of causal variants
# gr = genetic correlation
# GA = heterogeneity pattern
# train = training size
# model = SNP set
# anc = ancestry


plink2 \
--bfile all_chr \
--keep ID/train${train}_${anc}.id \
--extract Generate/aou_snp.txt \
--linear hide-covar cols=+a1freq \
--covar-variance-standardize \
--pheno PHENO/phenotypes_rho${rho}_gr${gr}_GA${GA}_${model}.txt \
--pheno-col-nums 3-12 \
--covar PCA/all_sim_pca.sscore \
--out GWAS/gwas_${anc}${train}_rho${rho}_gr${gr}_GA${GA}_${model}