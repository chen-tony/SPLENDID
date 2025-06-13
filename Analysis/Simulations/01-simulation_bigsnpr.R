## CONVERTING PLINK TO BIGSNPR FORMAT

library(bigsnpr)
library(dplyr)
library(data.table)
library(tidyr)

train = params[[1]] # '2a' or '4a'

### TRAINING-TUNING SPLIT
train.id = fread(paste0('ID/train', train, '.id'), header=F)$V1
tunval.id = fread(paste0('ID/tunval.id'), header=F)$V1

all_id = unique(c(train.id, tunval.id))

ind.row = which(fam$V1 %in% all_id)
ind.col = 1:nrow(bim)

table(substr(fam$V1[ind.row], 1, 3))

print(length(ind.row))
print(length(ind.col))

snp_readBed2(
  paste0('all_chr.bed'),
  backingfile = paste0('PLINK/sim_train', train, '_tun'),
  ind.row = ind.row,
  ind.col = ind.col,
  ncores = 5
)

plink = snp_attach(paste0('PLINK/sim_train', train, '_tun.rds'))

map = plink$map
saveRDS(map, paste0('PLINK/sim_train', train, '_tun_map.rds'))