source('gex.init.R')
source('lib/multi_emma.R')
current <- 0
load('variance.kinship.Rdata')
load('variance.10M.kinship.Rdata')

current<-chrid
sigma <- array(0,dim=c(b+1))
sigma[b+1] <- 1
niter<-200
currentY <- Y-mean(Y)

temp.fun <- function(K,Y,X){
  print(dim(K))
  emma.REML.t(Y,X,K)
}
result <- apply(K,1,temp.fun,Y,(rep(1,length(Y))))
save(chrid, genestart,result,file=paste(datafolder,phenotypename,"_emma_10M.Rdata",sep=""))
