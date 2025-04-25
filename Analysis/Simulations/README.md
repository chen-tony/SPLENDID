We used simulated multi-ancestry genotypes accessible at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/COXHAP

The order of files are
1) Subset genotypes into training, tuning, and validation
2) Conduct PCA on genotypes
3) Generate effect sizes
4) Generate phenotypes
5) Conduct GWAS (ancestry-specific and GxPC
   a) Re-format for PRS-CSx, CT-SLEB, and PROSPER
6) Run PRS methods
   a) PRS-CSx, CT-SLEB, and PROSPER w/ tuning
   b) iPGS w/ tuning + refitting
   c) SPLENDID w/ tuning
7) Apply PRS to validation data
