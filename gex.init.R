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
markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
x=data.matrix(read.table(paste(genotypefolder,1,'.genotype',sep=""),na.string='N'))
nk=ncol(x)
