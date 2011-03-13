

source('gex.header.R')
source('methods.R')
chrgenotype=0
selected=which((markers[,3]>=rangestart) & (markers[,3]<=rangeend) & (markers[,2]==chrid))
#%XX=X[selected,]
#XX=XX[complete.cases(XX),]

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
