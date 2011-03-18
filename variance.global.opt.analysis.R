##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
options(scipen=4)
library(emma)
library('DEoptim')
source('lib/imageplot.R')
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  a=1
  b=19
  print("No arguments supplied.")
  a=1
  b=19
  a=1
  b=19
  genotypefolder="result/precc_liver_gene_expression/genotype/"
  kinshipfolder="result/precc_liver_gene_expression/kinship/"
  chrid=3
  genestart=105762278
  geneend=105772678
  rangestart=102910198
  rangeend=112949064
  phenotypename="ENSMUSG00000000001"
  phenotypefile="result/precc_liver_gene_expression/gene_expression/0/ENSMUSG00000000001"
  step=10000000
  datafolder="result/precc_liver_gene_expression/emma/0/"
  kinshipfolder="result/gse13870_female/kinship/local_10000000/local.kinships.Rdata"
  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


density.sigma <- function(x){
  currentcomp=K[current,,]
  Omega=x*currentcomp+R
#  print(Omega)
  t <- -(-size/2*log(2*pi)-(det(Omega,T))/2-t(currentY) %*% solve(Omega) %*% currentY/2)
  if (is.infinite(t)){
    print(Omega)
    print(t)
    print(current)
  }
  t
}

                                        # return density of sigma
density.sigma.all <- function(sigma){
#  print(Omega)
  for (k in 1:idx){
    if(sigma[k]!=lastsigma[k])
      Omega<-Omega + (sigma[k]-lastsigma[k])*K[k,,]
  }
  lastsigma <- sigma
  if (sum(sigma)==0) return (Inf)
  t <- -(-size/2*log(2*pi)-(log(det(Omega)))/2-t(currentY) %*% solve(Omega) %*% currentY/2 )
  if (is.infinite(t)){
    return (Inf)
  }
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


density.sigma.batch <- function(x){
  exp(sapply(x,density.sigma))
}

Y=data.matrix(read.table(phenotypefile))
is.known=(!is.na(Y))
Y <- Y[is.known]
Y <- Y*100
size=length(Y)
markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
cauchy.p <- 10
load(paste(kinshipfolder,'/local_',step,'/local.kinships.Rdata',sep=""))
K <- K[,is.known, is.known]
idx <- dim(K)[1]
sigma <- array(0.00,dim=c(idx))
sigma[idx] <- var(Y)
# calculate total variance covariance matrix
TotalVariance<-array(0,dim=c(size,size))
for (k in 1:idx){
  TotalVariance<-TotalVariance + sigma[k]*K[k,,]
}
Omega <- TotalVariance
lastsigma <- sigma
current<-chrid
niter<-10
currentY <- Y-mean(Y)
lower <- array(0,c(idx))
upper <- array(var(Y)*3,c(idx))
optresult <- nlminb(sigma, density.sigma.all, density.sigma.all.grad, lower=lower)
print(optresult)
result <- optresult$par
save(chrid, genestart,result, optresult,file=paste(datafolder,phenotypename,"_global_nlminb.Rdata",sep=""))

#r <- DEoptim(density.sigma.all,lower,upper,control = DEoptim.control(NP = 50))
# r$member$bestmemit[200,]
