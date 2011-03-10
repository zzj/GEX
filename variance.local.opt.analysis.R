##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(emma)
source('lib/imageplot.R')
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  a=1
  b=19
  print("No arguments supplied.")
  genotypefolder="result/precc_liver_gene_expression/genotype/"
  kinshipfolder="result/precc_liver_gene_expression/kinship/"
  chrid=3
  genestart=105762278
  geneend=105772678
  rangestart=102910198
  rangeend=112949064
  phenotypename="ENSMUSG00000000001"
  phenotypefile="result/precc_liver_gene_expression/gene_expression/0/ENSMUSG00000000001"
  datafolder="result/precc_liver_gene_expression/emma/0/"
  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


Y=data.matrix(read.table(phenotypefile))
Y
# 10000000 249
step=10000000
size=249
markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
x=data.matrix(read.table(paste(genotypefolder,1,'.genotype',sep=""),na.string='N'))
nk=ncol(x)

# return density of sigma
density.sigma.all <- function(sigma){
#  print(Omega)
  for (k in 1:idx){
    if(sigma[k]!=lastsigma[k])
      Omega<-Omega + (sigma[k]-lastsigma[k])*K[k,,]
  }
  lastsigma <- sigma
  if (sum(sigma)==0) return (Inf)
  t <- -(-nk/2*log(2*pi)-(log(det(Omega)))/2-t(currentY) %*% solve(Omega) %*% currentY/2 )
  if (is.nan(t)){
    print(det(Omega))
    return (100000000000)
  }
  print(sigma)
  print(t)
  t
}

density.sigma.all.grad <- function(sigma){
  for (k in 1:idx){
    if(sigma[k]!=lastsigma[k])
      Omega<-Omega + (sigma[k]-lastsigma[k])*K[k,,]
  }
  lastsigma <- sigma
  invOmega <- solve(Omega)
  detOmega <- det(Omega)
  ret <- array(0,dim=c(idx))
  for (k in 1:(idx)){
    ret[k] <- sum(rowSums((invOmega)* (t(K[k,,]))))/2-t(currentY) %*% invOmega %*% K[k,,] %*% invOmega %*% currentY/2
  }
  ret
}

cauchy.p <- 10
load('variance.10M.kinship.Rdata')
sigma <- array(0.00,dim=c(idx))
sigma[idx] <- var(Y)
# calculate total variance covariance matrix
TotalVariance<-array(0,dim=c(nk,nk))
for (k in 1:idx){
  if(k!=idx)
    K[k,,] <- K[k,,]-diag(nk)
  TotalVariance<-TotalVariance + sigma[k]*K[k,,]
}
Omega <- TotalVariance
lastsigma <- sigma
current<-chrid
niter<-10
currentY <- Y-mean(Y)
lower <- array(0,c(idx))
upper <- array(var(Y)*3,c(idx))
optresult <- optim(sigma, density.sigma.all, density.sigma.all.grad, method="L-BFGS-B", lower=lower)
result <- optresult$par
save(chrid, genestart,result, optresult,file=paste(datafolder,phenotypename,"_local_bfgs.Rdata",sep=""))
