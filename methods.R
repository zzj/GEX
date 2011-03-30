
lasso.analysis <- function(X,Y,datafolder,phenotypename){
  choose=complete.cases(X)
  X=X[choose,]
  object3 <- lars(t(X),Y,type="lar",use.Gram=FALSE) # Can use abbreviations
  save(object3,choose,file=paste(datafolder,phenotypename,".Rdata",sep=""))

}

lasso.plot <- function(datafolder, phenotypename, markers,chrid,genestart,geneend){
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  load(file=paste(datafolder, phenotypename, ".Rdata",sep=""))
  png(paste(datafolder,phenotypename,"_lasso.png",sep="") ,width=1440)
  #gex.mhtplot(markers[choose,],chrid,genestart,geneend,abs(object3$beta[which.min(object3$Cp),]),phenotypename,logscale=FALSE)
  gex.mhtplot(markers[choose,],chrid,genestart,geneend,abs(object3$beta[10,]),phenotypename,logscale=FALSE)
  dev.off();
}

lasso.qtl.result <- function(filename){
  load(file=filename)
  a=abs(object3$beta[dim(object3$beta[,])[1]-1,])
  
  return(a)
}

variance.qtl.result <- function(filename){
  load(file=filename)
  return (optresult$par[-length(optresult$par)])
}


std.analysis <- function(datafolder, phenotypename,X,Y){
  num.snps=dim(X)[1]
  p=mat.or.vec(num.snps,1)
  for (j in 1:(num.snps)){
    snp=as.matrix(X[j,])
    s=(anova(glm(Y~ snp)))
    p[j]=pchisq(s$Deviance[2],1, lower.tail=F)
  }
  save(p,file=paste(datafolder,phenotypename,".Rdata",sep=""))
}

std.plot <- function(datafolder, phenotypename, markers,chrid,genestart,geneend){
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  load(file=paste(datafolder, phenotypename, ".Rdata",sep=""))
  png(paste(datafolder,phenotypename,"_std.png",sep="") ,width=1440)
  gex.mhtplot(markers,chrid,genestart,geneend,p,phenotypename)
  dev.off();
}

std.pvalue.hist.plot <- function(datafolder, phenotypename, markers, chrid, genestart, geneend){
  markers=data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
  load(file=paste(datafolder, phenotypename, ".Rdata",sep=""))
  png(paste(datafolder,phenotypename,"_std_hist.png",sep=""))
  gex.pvalue.hist(p,phenotypename)
  dev.off();
}

std.qtl.result <- function(filename, log=T){
  load(file=filename)
  if (log)
    return (-log(p))
  else
    return (p)
}
