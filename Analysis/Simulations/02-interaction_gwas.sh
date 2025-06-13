# rho = percentage of causal variants
# gr = genetic correlation
# GA = heterogeneity pattern
# train = training size
# model = SNP set

plink2 \
--bfile all_chr \
--keep ID/train${train}.id \
--extract Generate/aou_snp.txt \
--linear interaction hide-covar \
--covar-variance-standardize \
--pheno PHENO/phenotypes_rho${rho}_gr${gr}_GA${GA}_${model}.txt \
--pheno-col-nums 3-12 \
--covar PCA/all_sim_pca.sscore \
--parameters 1-26 \
--tests 22-26 \
--out GWAS/gwas_INT${train}_rho${rho}_gr${gr}_GA${GA}_${model}