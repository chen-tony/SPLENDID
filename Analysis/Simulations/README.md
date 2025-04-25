The order of files are
1) Generate data
   1) Subset genotypes into training, tuning, and validation
   2) Create bigsnpr files
   3) Conduct PCA on genotypes
   4) Generate effect sizes
   5) Generate phenotypes
2) Conduct GWAS (ancestry-specific and GxPC
   1) Re-format for PRS-CSx, CT-SLEB, and PROSPER
3) Run PRS-CSx
   1) Run PRS-CSx
   2) Tuning / obtain results
4) Run CT-SLEB
   1) 2-way Clumping+threshold and empirical Bayes
   2) Super learning / obtain results
5) Run PROSPER
   1) Run PROSPER
   2) Obtains results
6) Run iPGS+refit
   1) iPGS (Lasso)
   2) Tuning
   3) iPGS + refit / obtain results
7) Run SPLENDID
   1) Compute summaries
   2) Run penalized regression
   3) Tuning
   4) Ensembling / obtain results
8) Compile results
