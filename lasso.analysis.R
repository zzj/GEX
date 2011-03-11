source('gex.header.R')
chrgenotype=0
selected=which((markers[,3]>=rangestart) & (markers[,3]<=rangeend) & (markers[,2]==chrid))
#%XX=X[selected,]
#XX=XX[complete.cases(XX),]

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


if (fun=='analysis'){
  choose=complete.cases(X)
  X=X[choose,]
  object3 <- lars(t(X),Y,type="lar",use.Gram=FALSE) # Can use abbreviations
  save(object3,choose,file=paste(datafolder,phenotypename,".Rdata",sep=""))
  #Can not call function, because R does not support pass by reference
  #lasso.analysis(X,Y,datafolder,phenotypename)
}else if (fun=='plot'){
  lasso.plot(datafolder,phenotypename,markers,chrid,genestart,geneend)
}
#which(object3$beta[which.min(object3$Cp),]>0)
