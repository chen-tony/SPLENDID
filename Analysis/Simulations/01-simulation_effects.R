library(dplyr)
library(data.table)
library(tidyr)
library(mvtnorm)

GenSigma = function(sigma,n1,n2,n3,n4,n5,gr){
  vecsigma = sigma*c(1/n1, gr/sqrt(n1*n2), gr/sqrt(n1*n3), gr/sqrt(n1*n4), gr/sqrt(n1*n5), 
                     gr/sqrt(n2*n1), 1/n2, gr/sqrt(n2*n3), gr/sqrt(n2*n4), gr/sqrt(n2*n5),
                     gr/sqrt(n3*n1), gr/sqrt(n3*n2), 1/n3, gr/sqrt(n3*n4), gr/sqrt(n3*n5), 
                     gr/sqrt(n4*n1), gr/sqrt(n4*n2), gr/sqrt(n4*n3), 1/n4, gr/sqrt(n4*n5),
                     gr/sqrt(n5*n1), gr/sqrt(n5*n2), gr/sqrt(n5*n3), gr/sqrt(n5*n4), 1/n5)
  Sigma = matrix(vecsigma,5,5)
  return(Sigma)
}

aou = fread('AOU/aou_ukb_map_final.txt')

bim_map = fread('sim_mapping.txt') %>%
  filter(rsid %in% aou$rsid)

set.seed(123)
causal_snps = list()
causal_snps[[4]] = sample(bim_map$rsid, length(bim_map$rsid)*0.02) # 12096
causal_snps[[3]] = sample(causal_snps[[4]], length(bim_map$rsid)*0.01) # 6048
causal_snps[[2]] = sample(causal_snps[[3]], length(bim_map$rsid)*0.005) # 3024
causal_snps[[1]] = sample(causal_snps[[2]], length(bim_map$rsid)*0.001) # 604
saveRDS(causal_snps, 'Generate/causal_snps_aou.RDS')

causal_snps = list(`1`=604, `2`=3024, `3`=6048, `4`=12096)
causal_proportions = list(`1`='0.1%', `2`='0.5%', `3`='1%', `4`='2%')

eth = c('EUR', 'AFR', 'EAS', 'AMR', 'SAS')

causal_snps = readRDS(paste0(home, 'Generate/causal_snps_aou.RDS'))

mix_prop = c(0, 0.5) # proportion of mixed effects
# GA=1: all homogeneous or all heterogeneous
# GA=2: 50% homogeneous / 50% heterogeneous

grGA = data.frame(
  GA = c(1, 1, 2, 2),
  gr = c(0.6, 0.8, 0.4, 0.6)
)

for (rho in 1:4) {
  print(rho)
  freq_cau = bim_map %>%
    filter(rsid %in% causal_snps[[rho]])
  write.table(freq_cau$snp_code,
              file = paste0('Generate/select.cau_rho',rho,'_aou.snp'),
              row.names = F,col.names = F,quote=F)
  
  n.total.snp = length(causal_snps[[rho]])
  
  sigma = 0.4
  
  GA = 1
  for (gr in c(0.6, 0.8)) {
    Sigma = GenSigma(sigma,n.total.snp,n.total.snp,n.total.snp,
                     n.total.snp,n.total.snp,gr)
    
    beta = rmvnorm(n.total.snp, rep(0, 5), Sigma)
    
    for (i in 1:length(eth)) {
      select.cau = data.frame(SNP=freq_cau$snp_code, BETA=beta[,i])
      write.table(select.cau,file = paste0('Generate/Beta/select.cau_', eth[i], '_rho',rho,
                                           '_gr',gr, '_GA', GA, '_aou'),
                  # '_gr',gr, '_GA', GA, '_hm3'),
                  row.names = F,col.names = F,quote=F)
    }
    print(colSums(beta^2))
  }
  
  GA = 2
  for (gr in c(0.4, 0.6)) {
    Sigma = GenSigma(sigma,n.total.snp,n.total.snp,n.total.snp,
                     n.total.snp,n.total.snp,gr)
    
    n.mix.snps = ceiling(n.total.snp * 0.5)
    n.shared.snps = n.total.snp-n.mix.snps
    beta_mix = rmvnorm(n.mix.snps, rep(0, 5), Sigma)
    beta_shared = rnorm(n.shared.snps, 0, sqrt(sigma / n.total.snp)) %>%
      replicate(5,.)
    
    beta = rbind(beta_mix, beta_shared)
    
    for (i in 1:length(eth)) {
      select.cau = data.frame(SNP=freq_cau$snp_code, BETA=beta[,i])
      write.table(select.cau,file = paste0('Generate/Beta/select.cau_', eth[i], '_rho',rho,
                                           '_gr',gr, '_GA', GA, '_aou'),
                  # '_gr',gr, '_GA', GA, '_hm3'),
                  row.names = F,col.names = F,quote=F)
    }
    print(colSums(beta^2))
  }
}

