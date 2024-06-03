# Slingshot

Slingshot is a biobank-scale penalized regression framework to model shared and heterogeneity genetic effects, such as across diverse ancestries. The associated preprint can be found on [bioRxiv](url). The software (detailed below) requires R, C++, and plink2, and entirely run through command line arguments. Much of the package borrows code from the bigstatsr package for efficient analysis of large genetic data. Feedback and suggestions are always welcome to improve code functionality and usability!

# Tutorial
## Create file-backed matrix for genotype data

## Run initial summaries of data

## Run group-L0L1 regression

## Run grid search and ensemble learning

## Apply PRS to external validation data

# Notes
- Slingshot was primarily built for multi-ancestry PRS analysis, where we use interactions with ancestry PCs to model heterogeneity. However, ancestry PCs can be replaced with other variables, such as social risk factors or environmental exposures, depending on interest. The current framework models all interactions simultaneously, but future work can refine the penalty to sparse model interactions within groups. 
- Slingshot is currently built assuming a training - tuning - validation data split structure. However, the same code can be used for cross-validation-type analysis by running it several times on different folds of the training data. Then, one can use the multiple folds to choose specific tuning parameters, or combine models across folds similar to the CMSA strategy in bigstatsr. 
