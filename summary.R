##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library('multicore')
library(emma)
source('lib/imageplot.R')
source('lib/slice.R')
source('methods.R')
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
library(matlab)
if(length(args)==0){
  print("No arguments supplied.")
  genotypefolder="result/precc_liver_gene_expression/genotype/"
  kinshipfolder="result/precc_liver_gene_expression/kinship/"
  geneexpfolder="result/precc_liver_gene_expression/gene_expression/"
  resultfolder="result/precc_liver_gene_expression/summary/"
  a=1
  b=19
}else{
   for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


gex.prepare.plot.eqtl <- function(filelist,qtlresultfun,totalinterval,groups){
  Img <- list()
  for ( i in 1:length(filelist)){
    result <- qtlresultfun(filelist[i])
    if (is.null(groups)) t <- result[]
    else t <- ((as.numeric(tapply(result,groups,max))))
    if (max(t)!=0) t=t/max(t)
    Img[[i]] <- t
  }
  Img
}


# geneexpfolder : gene expression data
# resultfolder : gene result folder

markers <- data.frame(read.table(paste(genotypefolder,'marker_list',sep="")),stringsAsFactors =FALSE)
idx <- 1
groups <- rep(0,dim(markers)[1])
for ( i in 1:19){
  minpos=min(markers[markers[,2]==i,3])
  maxpos=max(markers[markers[,2]==i,3])
  while(minpos<maxpos){
    selected=which((markers[,3]>=minpos) & (markers[,3]<minpos+step) & (markers[,2]==i))
    if (length(selected)>10){
      groups[selected] <- idx
      idx <- idx+1
    }
    minpos <- minpos+step
  }
}

stdfilelist <- c()
emmafilelist <- c()
variancefilelist <- c()
lassofilelist <- c()
emmafilelist <- c()
#names[i,1] name
#names[i,2] folder
#names[i,3] chr
#names[i,4] gene strat position
#names[i,5] gene end position
names <- (read.table(paste(geneexpfolder,"gene_list",sep=""), stringsAsFactors =FALSE))
for (i in 1:dim(names)[1]){
  filename=paste(root,'/std/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    stdfilelist <- c(stdfilelist,filename)
  }
  filename=paste(root,'/lasso/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    lassofilelist <- c(lassofilelist,filename)
  }    
  filename=paste(root,'/emma/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    emmafilelist <- c(emmafilelist,filename)
  }    
  filename=paste(root,'/variance/',names[i,2],'/',names[i,1],'_global_nlminb.Rdata',sep='')
  if (file.exists(filename)){
    variancefilelist <- c(variancefilelist,filename)
  }
}


if (fun=='variance') {
  resultfilename <- 'eqtl_variance'
  Img <- pvec((variancefilelist),gex.prepare.plot.eqtl,variance.qtl.result,idx,NULL)
#  Img <- gex.prepare.plot.eqtl(variancefilelist, variance.qtl.result, idx,NULL)
  Img <- do.call(rbind, Img)
  save(Img, file=paste(resultfolder,'/',resultfilename,'.Rdata',sep=""))
  mycolors <- jet.colors(16)
  pdf(paste(resultfolder,'/',resultfilename,'.pdf',sep=""))
  imagesc((Img)*16, xlab = "pos", ylab = "gene", col = mycolors )
  dev.off()
}else  if (fun=='lasso') {
  resultfilename <- 'eqtl_lasso'
  Img <- pvec(lassofilelist,gex.prepare.plot.eqtl,lasso.qtl.result,idx,groups)
  Img <- do.call(rbind, Img)
  save(Img, file=paste(resultfolder,'/',resultfilename,'.Rdata',sep=""))
  mycolors <- jet.colors(16)
  pdf(paste(resultfolder,'/',resultfilename,'.pdf',sep=""))
  imagesc((Img/max(Img))*16, xlab = "pos", ylab = "gene", col = mycolors )
  dev.off()
} else  if (fun=='std') {
  resultfilename <- 'eqtl_std'
  Img <- pvec(stdfilelist,gex.prepare.plot.eqtl,std.qtl.result,idx,groups)
  Img <- do.call(rbind, Img)
  save(Img, file=paste(resultfolder,'/',resultfilename,'.Rdata',sep=""))
  mycolors <- jet.colors(16)
  pdf(paste(resultfolder,'/',resultfilename,'.pdf',sep=""))
  imagesc((Img/max(Img))*16, xlab = "pos", ylab = "gene", col = mycolors )
  dev.off()
} else  if (fun=='emma') {
  resultfilename <- 'eqtl_emma'
  Img <- pvec(emmafilelist,gex.prepare.plot.eqtl,emma.qtl.result,idx,groups)
  Img <- do.call(rbind, Img)
  save(Img, file=paste(resultfolder,'/',resultfilename,'.Rdata',sep=""))
  mycolors <- jet.colors(16)
  pdf(paste(resultfolder,'/',resultfilename,'.pdf',sep=""))
  imagesc((Img/max(Img))*16, xlab = "pos", ylab = "gene", col = mycolors )
  dev.off()
} else if (fun=='true'){
  resultfilename <- 'eqtl_true'
  selected <- seq(1,dim(names)[1],by=rough.scan)
  true <- names[selected,]
  print((as.numeric(true[,3])))
  true <- true[which(!is.na(as.numeric(true[,3]))),]
  Img<-mat.or.vec(dim(true)[1],idx-1)
  for ( i in 1:dim(true)[1]){
    t <- mat.or.vec(1,idx-1)
      right=min(which((markers[,3]>=true[i,4]) & (markers[,3]<true[i,4]+step) & (markers[,2]==true[i,3])))
      if (is.infinite(right)){
        print(true[i,])
      }
      t[groups[right]]=1
      Img[i,]=t[]
  }
  pdf(paste(resultfolder,'/',resultfilename,'.pdf',sep=""))
  imagesc((Img/max(Img))*16, xlab = "pos", ylab = "gene", col = jet.colors(16) )
  dev.off()
}
