We used simulated multi-ancestry genotypes accessible at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP

The order of files are
1) Generate data
   1) Subset genotypes into training, tuning, and validation
   2) Conduct PCA on genotypes
   3) Generate effect sizes
   4) Generate phenotypes
   5) Create bigsnpr files
6) Conduct GWAS (ancestry-specific and GxPC
   1) Re-format for PRS-CSx, CT-SLEB, and PROSPER
7) Run PRS methods
   1) PRS-CSx, CT-SLEB, and PROSPER w/ tuning
   2) iPGS w/ tuning + refitting
   3) SPLENDID w/ tuning
8) Apply PRS to validation data
