## First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library('abind')

library(emma)
library(matlab)
source('lib/imageplot.R')
## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  a=1
  b=19
  genotypefolder="result/precc_liver_gene_expression/genotype/"
  datafolder="result/precc_liver_gene_expression/kinship/"
  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

gex.kinship.image <- function(matrix, datafolder, filename,zlim){
  pdf(paste(datafolder,'/',filename,".pdf",sep=""))
  myImagePlot(as.matrix(matrix),zlim=zlim)
  dev.off()
}

options(scipen=4)

gex.build.kinship.by.chr <- function(genotypefolder, datafolder, a,b){
  chrkinships <- vector("list",b-a+1)
  # build local kinship matrix
  for ( i in a:b){
    x=data.matrix(read.table(paste(genotypefolder,i,'.genotype',sep=""),na.string='N'))
    x=x/2
    if (i==a)
      X=x
    else
      X=rbind(X,x)
    chrkinships[[i]]=emma.kinship(x)
    gex.kinship.image(chrkinships[[i]], datafolder, i)
  }
  # build global kinship matrix
  globalkinship <- emma.kinship(X)
  gex.kinship.image(globalkinship,datafolder,'all')
  save(chrkinships, globalkinship,file=paste(datafolder,"kinships.Rdata",sep=""));
}

gex.build.kinship.by.region <- function(genotypefolder, datafolder,a,b,step){
  idx <- 1
  size <- 1
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  K<-array(0,dim=c(size,ncol(x), ncol(x)))
  datafolder <- file.path(datafolder, "local_1M")
  print(datafolder)
  dir.create(datafolder,showWarnings=F)
  zlim=c(0.3,1)
  for ( i in a:b){
    print(i)
    x=data.matrix(read.table(paste(genotypefolder,i,'.genotype',sep=""),na.string='N'))
    x=x/2
    if (i==1)
      X=x
    else X=rbind(X,x)
    minpos=min(markers[markers[,2]==i,3])
    maxpos=max(markers[markers[,2]==i,3])
    
    while(minpos<maxpos){
      selected=which((markers[,3]>=minpos) & (markers[,3]<minpos+step) & (markers[,2]==i))
      if (length(selected)>0){
        XX=X[selected,]
        tempK <- array(0,dim=c(1,ncol(x), ncol(x)))
        tempK[1,,]=emma.kinship(XX)
        if (idx==1){
          K=tempK
          K.chr=c(i)
          K.startpos=c(minpos)
          K.num=c(length(selected))
        }
        else{
          K=abind(K,tempK,along=1)
          K.chr <- append(K.chr, i)
          K.startpos <- append(K.startpos, minpos)
          K.num <- append(K.num, length(selected))
        }
        gex.kinship.image(K[idx,,],datafolder,paste("chr_",i,"_local_",idx,"markers_",K.num[idx],sep=""),zlim)
        idx <- idx+1
      }
      minpos <- minpos+step
    }
  }
  cbind(K,diag(nk))
  K.chr <- append(K.chr, -1)
  K.startpo <- append(K.startpos, -1)
  K.num <- append(K.num, -1)
  print(K.num)
  save(K,file=paste(datafolder,"local.kinships.Rdata",sep=""))
}
  
#gex.build.kinship(genotypefolder,datafolder,a,b);
gex.build.kinship.by.region(genotypefolder, datafolder,a,b,1000000)

