##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(emma)
source('lib/imageplot.R')
source('lib/slice.R')
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
  phenotypename="ENSMUSG00000000561"
  phenotypefile="result/precc_liver_gene_expression/gene_expression/0/ENSMUSG00000000561"
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

temp <- function(){
  K<-array(0,dim=c(size,ncol(x), ncol(x)))
  chr <- array(0,dim=c(size))
  startpos <- array(0,dim=c(size))
  
  idx <- 1
  for ( i in a:b){
    x=data.matrix(read.table(paste(genotypefolder,i,'.genotype',sep=""),na.string='N'))
    x=x/2
    
    minpos=min(markers[markers[,2]==i,3])
    maxpos=max(markers[markers[,2]==i,3])
    while(minpos<maxpos){
      selected=which((markers[,3]>=minpos) & (markers[,3]<minpos+step) & (markers[,2]==i))
      XX=X[selected,]
      K[idx,,] <- emma.kinship(XX)
      chr[idx] <- i
      startpos[idx] <- minpos
      idx <- idx+1
      minpos <- minpos+step
    }
  }
  K[idx,,]=diag(nk)
  stop('a')
}
# return density of sigma
density.sigma <- function(x){
  currentcomp=K[current,,]
  
  Omega=x*currentcomp+R
  -nk/2*log(2*pi)-log(det(Omega))/2-t(currentY) %*% solve(Omega) %*% currentY/2 +(dcauchy(x,location=0,scale=cauchy.p,log=T))
}

density.sigma.batch <- function(x){
  exp(sapply(x,density.sigma))
}
current <- 0
load('variance.10M.kinship.Rdata')
current<-chrid
sigma <- array(0,dim=c(idx))
sigma[idx] <- var(Y)
niter<-50
currentY <- Y-mean(Y)
K[idx,,]=diag(nk)

#warm up
list.cauchy.p<- c(0.05,0.005,0.0005,0.00005)
result <- array(0,c(length(list.cauchy.p),niter,length(sigma)))
TotalVariance<-array(0,dim=c(nk,nk))

for (i in 1:idx){
  TotalVariance<-TotalVariance + sigma[i]*K[i,,]
}
ip <- 1
for (cauchy.p in list.cauchy.p){
  for (i in 1:niter){
    for (current in 1:idx){
      TotalVariance <- TotalVariance-sigma[current]*K[current,,]
      R <- TotalVariance
      sigma[current] <- uni.slice(sigma[current],density.sigma,w=0.1,lower=0)
      TotalVariance <- TotalVariance+sigma[current]*K[current,,]
    }
    result[ip,i,]=sigma
  }
  ip <- ip+1
}
result

save(chrid, genestart,result,file=paste(datafolder,phenotypename,"_10M_other.Rdata",sep=""))
