train=${1}
rho=${2}
gr=${3}
GA=${4}
rep=${5}
model=${6}

n=5
v=1
package=Prosper/Package

mkdir PROSPER/Analysis/train_${train}_rho${rho}_gr${gr}_GA${GA}_rep${rep}_${model}

outdir=PROSPER/Analysis/train_${train}_rho${rho}_gr${gr}_GA${GA}_rep${rep}_${model}

# LASSOSUM2
Rscript ${package}/scripts/lassosum2.R \
--PATH_package ${package} \
--PATH_out ${outdir}/lassosum2 \
--PATH_plink plink2 \
--FILE_sst PROSPER/GWAS/EUR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/AFR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/EAS_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/SAS_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/AMR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt \
--pop EUR,AFR,EAS,SAS,AMR \
--bfile_tuning PROSPER/PLINK/EUR_tuning,PROSPER/PLINK/AFR_tuning,PROSPER/PLINK/EAS_tuning,PROSPER/PLINK/SAS_tuning,PROSPER/PLINK/AMR_tuning \
--pheno_tuning PROSPER/Pheno/EUR_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/AFR_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/EAS_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/SAS_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/AMR_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--bfile_testing PROSPER/PLINK/EUR_validation,PROSPER/PLINK/AFR_validation,PROSPER/PLINK/EAS_validation,PROSPER/PLINK/SAS_validation,PROSPER/PLINK/AMR_validation \
--pheno_testing PROSPER/Pheno/EUR_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/AFR_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/EAS_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/SAS_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt,PROSPER/Pheno/AMR_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--testing TRUE \
--NCORES $n \
--verbose $v

# PROSPER
Rscript ${package}/scripts/PROSPER.R \
--PATH_package ${package} \
--PATH_out ${outdir}/PROSPER \
--FILE_sst PROSPER/GWAS/EUR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/AFR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/EAS_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/SAS_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt,PROSPER/GWAS/AMR_gwas_train${train}_rho${rho}_gr${gr}_rep${rep}_GA${GA}_${model}.txt \
--pop EUR,AFR,EAS,SAS,AMR \
--lassosum_param ${outdir}/lassosum2/EUR/optimal_param.txt,${outdir}/lassosum2/AFR/optimal_param.txt,${outdir}/lassosum2/EAS/optimal_param.txt,${outdir}/lassosum2/SAS/optimal_param.txt,${outdir}/lassosum2/AMR/optimal_param.txt \
--NCORES $n \
--verbose $v

# Ancestry-Specific Ensemble
for anc in EUR AFR EAS SAS AMR; do
Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink plink2 \
--PATH_out ${outdir}/PROSPER \
--prefix ${anc} \
--testing TRUE \
--SL_library SL.glmnet,SL.ridge \
--bfile_tuning PROSPER/PLINK/${anc}_tuning \
--pheno_tuning PROSPER/Pheno/${anc}_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--bfile_testing PROSPER/PLINK/${anc}_validation \
--pheno_testing PROSPER/Pheno/${anc}_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--cleanup F \
--NCORES $n \
--verbose $v
done


# Pooled Ancestry Ensemble
for anc in ALL; do
Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink plink2 \
--PATH_out ${outdir}/PROSPER \
--prefix ${anc} \
--testing TRUE \
--SL_library SL.glmnet,SL.ridge \
--bfile_tuning PROSPER/PLINK/${anc}_tuning \
--pheno_tuning PROSPER/Pheno/${anc}_tun_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--bfile_testing PROSPER/PLINK/${anc}_validation \
--pheno_testing PROSPER/Pheno/${anc}_val_pheno_rho${rho}_gr${gr}_GA${GA}_${model}_rep${rep}.txt \
--cleanup F \
--NCORES $n \
--verbose $v
done