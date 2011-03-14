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

# geneexpfolder : gene expression data
# resultfolder : gene result folder

markers <- data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
idx <- 1
step=10000000
groups <- rep(0,dim(markers)[1])
for ( i in 1:19){
  minpos=min(markers[markers[,2]==i,3])
  maxpos=max(markers[markers[,2]==i,3])
  while(minpos<maxpos){
    selected=which((markers[,3]>=minpos) & (markers[,3]<minpos+step) & (markers[,2]==i))
    if (length(selected)!=0){
      groups[selected] <- idx
      idx <- idx+1
    }
    minpos <- minpos+step
  }
}

#names[i,1] name
#names[i,2] folder
#names[i,3] chr
stdfilelist <- c()
emmafilelist <- c()
variancefilelist <- c()
lassofilelist <- c()

names <- data.matrix(read.table(paste(geneexpfolder,"gene_list",sep="")))
for (i in 1:dim(names)[1]){
  filename=paste(root,'/std/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    stdfilelist <- c(stdfilelist,filename)
  }
  filename=paste(root,'/lasso/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    lassofilelist <- c(lassofilelist,filename)
  }    
}


gex.prepare.plot.eqtl <- function(filelist,qtlresultfun,totalinterval,groups){
  mycolors <- jet.colors(16)
  Img <- mat.or.vec(length(filelist),totalinterval-1)
  for ( i in 1:length(filelist)){
    result <- qtlresultfun(filelist[i])
    t <- ((as.numeric(tapply(result,groups,max))))
    t <- t/max(t)
    Img[i,]=t[]
  }
  Img
}

resultfilename <- 'eqtl_lasso'
Img <- pvec(lassofilelist, gex.prepare.plot.eqtl,lasso.qtl.result,idx,groups)
save(Img, file=paste(resultfilename,'.Rdata',sep=""))
pdf(paste(resultfilename,'.pdf',sep=""))
imagesc(Img*16, xlab = "pos", ylab = "gene", col = mycolors )
dev.off()


