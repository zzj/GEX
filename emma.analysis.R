##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(emma)
library(matlab)
library(gap)
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
  phenotypename="ENSMUSG00000000561"
  phenotypefile="result/precc_liver_gene_expression/gene_expression/0/ENSMUSG00000000561"
  datafolder="result/precc_liver_gene_expression/emma/0/"

  ##supply default values
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

gex.emma.mapping <- function(phenotypename, phenotypefile, chrid, genotypefolder, kinshipfolder, datafolder, a, b,rangestart, rangeend, genestart,geneend){
  load(paste(kinshipfolder,"kinships.Rdata",sep=""))
  Y=data.matrix(read.table(phenotypefile))
  Y
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  for ( i in a:b){
    x=data.matrix(read.table(paste(genotypefolder,i,'.genotype',sep=""),na.string='N'))
    x=x/2
    if (i==a){
      X=x
    } else{
      X=rbind(X,x)
    }
    if (chrid==i){
      chrg <- x
    }
  }
  selected=which((markers[,3]>=rangestart) & (markers[,3]<=rangeend) & (markers[,2]==chrid))
  XX=X[selected,]
#  positions=markers[selected,3]
#  ret=emma.REML.t(Y,X,globalkinship)
#  retchr=emma.REML.t(Y,X,chrkinships[[chrid]])
  retren=emma.REML.t(Y,X,emma.kinship(XX))
#  save(ret,retchr,retren,file=paste(datafolder,phenotypename,"_append.Rdata",sep=""))
  save(retren,file=paste(datafolder,phenotypename,"_append.Rdata",sep=""))
}


gex.emma.mhtplot <- function(markers,chrid,genestart,geneend,ps,phenotypename){
  mydata=(as.matrix(cbind(as.integer(markers[,2]),as.integer(markers[,3]),as.numeric(ps))))
  selected=which((markers[,3]>=genestart) & (markers[,3]<=geneend) & (markers[,2]==chrid))
  if (length(selected)==1)
    hdata=as.data.frame(cbind(t(mydata[selected,]),(rep(phenotypename,length(selected)))))
  else 
    hdata=as.data.frame(cbind((mydata[selected,]),(rep(phenotypename,length(selected)))))
  print(dim(hdata))
  if (dim(hdata)[2]==0){
    hdata <- NULL
  }
  print(hdata)
  color <- rep(c("lightgray","gray"),11)
  glen <- length(selected)
  hcolor <- rep("red",glen)
  print(hcolor)
  par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
  ops <- mht.control(colors=color,yline=1.5,xline=3)
  hops <- hmht.control(data=hdata,colors=hcolor)
  mhtplot(mydata,ops,hops,pch=19)
  abline(v=max(which(markers[,2]==chrid & markers[,3]<genestart)),col='red')
  axis(2,pos=2,at=1:16)
  title(paste("Manhattan plot with gene highlighted","chr=",chrid,"genepos=",genestart,"max=",max((mydata[,3])),", min=",min(mydata[,3])),cex.main=1.8)
}

gex.emma.mapping.plot <- function(phenotypename, phenotypefile, chrid, genotypefolder, kinshipfolder, datafolder, a, b,rangestart, rangeend, genestart,geneend){
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  load(file=paste(datafolder, phenotypename, ".Rdata",sep=""))
  load(file=paste(datafolder, phenotypename, "_append.Rdata",sep=""))
  png(paste(datafolder,phenotypename,"_global.png",sep="") ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,ret$ps,phenotypename)
  dev.off();
  png(paste(datafolder,phenotypename,"_local.png",sep="")  ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,retren$ps,phenotypename)
  dev.off();
  png(paste(datafolder,phenotypename,"_global_ves.png",sep="")  ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,ret$ves,phenotypename)
  dev.off();
  png(paste(datafolder,phenotypename,"_global_vgs.png",sep="")  ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,ret$vgs,phenotypename)
  dev.off();
  png(paste(datafolder,phenotypename,"_local_ves.png",sep="")  ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,retren$ves,phenotypename)
  dev.off();
  png(paste(datafolder,phenotypename,"_local_vgs.png",sep="")  ,width=1440)
  gex.emma.mhtplot(markers,chrid,genestart,geneend,retren$vgs,phenotypename)
  dev.off();
}
#}
comment <- function(){

  plot( positions,-log10(ret$ps), type="o", col="blue")
  lines(positions,-log10(retchr$ps), type="o", pch=22, lty=2, col="red")
  lines(positions,-log10(retren$ps), type="o", pch=23, lty=3, col="black")
  abline(v=genestart)
  abline(v=geneend)
  dev.off();
  pdf(paste(datafolder,phenotypename,"global.pdf",sep="") )
  plot( positions,-log10(ret$ps), type="o", col="blue")
  abline(v=genestart)
  abline(v=geneend)
  dev.off();
  pdf(paste(datafolder,phenotypename,"chr.pdf",sep="") )
plot(positions,-log10(retchr$ps), type="o", pch=22, lty=2, col="red")
  abline(v=genestart)
  abline(v=geneend)
  dev.off();
  pdf(paste(datafolder,phenotypename,"local.pdf",sep="") )
  plot(positions,-log10(retren$ps), type="o", pch=23, lty=3, col="black")
  abline(v=genestart)
  abline(v=geneend)

}

gex.emma.mapping.plot(phenotypename, phenotypefile, chrid, genotypefolder, kinshipfolder, datafolder, a,b,rangestart,rangeend, genestart,geneend)


