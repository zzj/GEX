##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(emma)
library(matlab)
library(gap)
library(lars)
source('lib/imageplot.R')
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  a=1
  b=19
  a=1
  b=19
  genotypefolder="result/precc_liver_gene_expression/genotype/"
  kinshipfolder="result/precc_liver_gene_expression/kinship/"
  chrid=3
  genestart=107910198
  geneend=107949064
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

chrgenotypes <- vector("list",b-a+1)
Y=data.matrix(read.table(phenotypefile))
markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
for ( i in a:b){
  x=data.matrix(read.table(paste(genotypefolder,i,'.genotype',sep=""),na.string='N'))
  x=x/2
  if (i==a){
    X=x
  } else{
    X=rbind(X,x)
  }
  chrgenotypes[[i-a+1]]=x
  if (chrid==i){
    chrg <- x
  }
}


gex.mhtplot <- function(markers,chrid,genestart,geneend,ps,phenotypename,logscale=TRUE){
  mydata=(as.matrix(cbind(as.integer(markers[,2]),as.integer(markers[,3]),as.numeric(ps))))
  selected=which((markers[,3]>=genestart) & (markers[,3]<=geneend) & (markers[,2]==chrid))
  if (length(selected)==1)
    hdata=as.data.frame(cbind(t(mydata[selected,]),(rep(phenotypename,length(selected)))))
  else 
    hdata=as.data.frame(cbind((mydata[selected,]),(rep(phenotypename,length(selected)))))
  print(dim(hdata))
  if (dim(hdata)[2]==0 || dim(hdata)[1]==0){
    hdata <- NULL
  }
  print(hdata)
  color <- rep(c("lightgray","gray"),11)
  glen <- length(selected)
  hcolor <- rep("red",glen)
  print(hcolor)
  par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
  ops <- mht.control(colors=color,yline=1.5,xline=3,cex=10)
  hops <- hmht.control(data=hdata,colors=hcolor)
  mhtplot(mydata,ops,hops,pch=19)
  abline(v=max(which(markers[,2]==chrid & markers[,3]<genestart)),col='red')
  axis(2,pos=2,at=1:16)
  title(paste("Manhattan plot with gene highlighted","chr=",chrid,"genepos=",genestart),cex.main=1.8)
}

