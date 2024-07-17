# SPLENDID

SPLENDID is a biobank-scale penalized regression framework to model shared and heterogeneity genetic effects, such as across diverse ancestries. The associated preprint can be found on [bioRxiv](url). The software (detailed below) requires R, C++, and plink2, and entirely run through command line arguments. Much of the package borrows code from the [bigstatsr](https://github.com/privefl/bigstatsr) package for efficient analysis of large genetic data . Feedback and suggestions are always welcome to improve code functionality and usability!

# Tutorial
## Required R packages
```
Rscript -e 'install.packages(c('optparse', 'Rcpp','RcppArmadillo', 'dplyr', 'data.table','tidyr', 'bigsnpr', 'Matrix', 'bigparallelr', 'foreach', 'glmnet'))'
```

## Create file-backed matrix for genotype data
Inputs require the following formats:
- BFILE: plink file without .bed/.bim/.fam extensions
- OUT: header for output file name (keep consistent throughout)
- KEEP.txt: text file with two columns (FID and IID) for all samples to use in analysis
- EXTRACT.txt: text file with SNPs to use in analysis
- NCORES: number of cores for parallel computing
- VERBOSE: 0 (no messaging), 1 (some messaging), 2 (more messsaging)
```
Rscript splendid_bigsnpr \
--bfile BFILE \
--out OUT \
--keep KEEP.txt \
--extract EXTRACT.txt \
--ncores NCORES \
--verbose VERBOSE
```
This will create an FBM for efficient analysis of individual-level data. If there are missing values, we can either ignore NAs (impute with 0), or impute with various options within the bigsnpr package. 

## Run initial summaries of data
Heterogeneity can be measured by meta-analysis of external GWAS or an interaction GWAS within the training data. For interaction GWAS, inputs require the following formats:
- BFILE: plink file without .bed/.bim/.fam extensions
- OUT: header for output file name (keep consistent throughout)
- PLINK2: path to plink2 software
- KEEP.txt: text file with two columns (FID and IID) for training samples
- EXTRACT.txt: text file with SNPs to use in analysis
- COVAR.txt: text file with FID, IID, and covariates
- COVAR_NAMES: comma-separated list of covariates of interest
- INTERACT_NAMES: command-separated list of covariates of interest for interactions
- PHENO.txt: text file with FID, IID, and phenotype(s)
- PHENO_NAME: name of phenotype of interest
- NCORES: number of cores for multi-threading in plink2
- VERBOSE: 0 (no messaging), 1 (some messaging), 2 (more messsaging)
```
Rscript splendid_gwas \
--bfile BFILE \
--out OUT \
--plink2 PLINK2 \
--keep KEEP.txt \
--extract EXTRACT.txt \
--covar COVAR.txt \
--covar-names COVAR_NAMES \
--interact-names INTERACT_NAMES \
--pheno PHENO.txt \
--pheno-name PHENO_NAME \
--ncores NCORES \
--verbose VERBOSE
```
This gives a text file with p-values for heterogeneity used for thresholding in regression. 

Then, we compute summaries within the training data to initialize in regression, with the following inputs:
Inputs require the following formats:
- BIGSNPR: FBM file for genotype data (with .rds extension)
- OUT: header for output file name (keep consistent throughout)
- PATH: path with SPLENDID package code
- KEEP.txt: text file with two columns (FID and IID) for all samples to use in analysis
- EXTRACT.txt: text file with SNPs to use in analysis
- COVAR.txt: text file with FID, IID, and covariates
- COVAR_NAMES: comma-separated list of covariates of interest
- INTERACT_NAMES: command-separated list of covariates of interest for interactions
- PHENO.txt: text file with FID, IID, and phenotype(s)
- PHENO_NAME: name of phenotype of interest
- NCORES: number of cores for multi-threading in plink2
- VERBOSE: 0 (no messaging), 1 (some messaging), 2 (more messsaging)
- IMPUTED: 0 to replace NA with 0, 1 if the FBM is already imputed (using bigsnpr package)
```
Rscript splendid_summaries \
--bigsnpr BIGSNPR \
--out OUT \
--plink2 PLINK2 \
--path PATH \
--keep KEEP.txt \
--extract EXTRACT.txt \
--covar COVAR.txt \
--covar-names COVAR_NAMES \
--interact-names INTERACT_NAMES \
--pheno PHENO.txt \
--pheno-name PHENO_NAME \
--ncores NCORES \
--verbose VERBOSE \
--imputed IMPUTED
```

## Run group-L0L1 regression
Inputs follow the same format as splendid_summaries.R, in addition to the following:
- KEEP.text: text file with FID and IID of training samples
- KEEP_TUN.text: text file with FID and IID of tuning samples (optional)
- SUMMARIES: summary file produced by splendid_summaries.R
- HET_GWAS: heterogeneity GWAS file with three columns:
  - marker.ID: SNP names
  - P_ADD: p-values for main effects
  - P_HET: p-values for heterogeneity
- LAMBDA: comma-separated list of lambda1 values
  - default to 5 log-scale values between 1e-2 and 1e-4
- PVAL: command-separate list of p-values for thresholding
- NLAMBDA0: number of lambda0 values (for each lambda1 and p-value)
- NLAMBDA0_MIN: minimum number of lambda0 values to test (with tuning data)
- NCHECK: number of active set checks
- NABORT: number of extract lambda0 to test (with tuning data)
- MAXSUPP: maximum number of SNPs to include
```
Rscript splendid_linear \
--bigsnpr BIGSNPR \
--out OUT \
--plink2 PLINK2 \
--path PATH \
--keep KEEP.txt \
--keep-tun KEEP_TUN.txt \
--extract EXTRACT.txt \
--covar COVAR.txt \
--covar-names COVAR_NAMES \
--interact-names INTERACT_NAMES \
--pheno PHENO.txt \
--pheno-name PHENO_NAME \
--lambda LAMBDA \
--pval PVAL \
--nlambda0 NLAMBDA0 \
--nlambda0-min NLAMBDA0_MIN \
--ncheck NCHECK \
--nabort NABORT \
--maxsupp MAXSUPP \
--ncores NCORES \
--verbose VERBOSE \
--imputed IMPUTED
```

## Run grid search and ensemble learning
Inputs follow the same format as before, in addition to the following:
- KEEP_TUN.text: text file with FID and IID of tuning samples
- KEEP_VAL.text: text file with FID and IID of validation samples (optional)
- ESUM: mean/sd for interaction covariates (from splendid_summaries.R)
- MAXSUPP: maximum number of SNPs to include in ensemble
- MAXFACTOR: maximum number of SNPs (relative to grid search PRS) to include in ensemble
- CLEANUP: 1 to remove intermediate files
```
Rscript splendid_tuning \
--bigsnpr BIGSNPR \
--bfile BFILE \
--out OUT \
--plink2 PLINK2 \
--path PATH \
--keep-tun KEEP_TUN.txt \
--extract EXTRACT.txt \
--covar COVAR.txt \
--covar-names COVAR_NAMES \
--interact-names INTERACT_NAMES \
--pheno PHENO.txt \
--pheno-name PHENO_NAME \
--Esum ESUM \
--maxsupp MAXSUPP \
--maxfactor MAXFACTOR \
--cleanup CLEANUP \
--verbose VERBOSE
```

## Apply to external validation data
Inputs follow the same format as before, in addition to the following:
- KEEP_VAL.text: text file with FID and IID of validation samples (optional)
```
Rscript splendid_validation \
--bfile BFILE \
--out OUT \
--plink2 PLINK2 \
--path PATH \
--keep-val KEEP_VAL.txt \
--covar COVAR.txt \
--covar-names COVAR_NAMES \
--interact-names INTERACT_NAMES \
--pheno PHENO.txt \
--pheno-name PHENO_NAME \
--Esum ESUM \
--cleanup CLEANUP \
--verbose VERBOSE
```

# Notes
- SPLENDID was primarily built for multi-ancestry PRS analysis, where we use interactions with ancestry PCs to model heterogeneity. However, ancestry PCs can be replaced with other variables, such as social risk factors or environmental exposures, depending on interest. The current framework models all interactions simultaneously, but future work can refine the penalty to sparse model interactions within groups. 
- SPLENDID is currently built assuming a training - tuning - validation data split structure. However, the same code can be used for cross-validation-type analysis by running it several times on different folds of the training data. Then, one can use the multiple folds to choose specific tuning parameters, or combine models across folds similar to the CMSA strategy in bigstatsr. 
