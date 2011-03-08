source('gex.header.R')
selected=which((markers[,3]>=rangestart) & (markers[,3]<=rangeend) & (markers[,2]==chrid))
XX=X[selected,]

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

if (fun=='analysis'){
  std.analysis(datafolder, phenotypename,X,Y)
}else if (fun=='plot'){
  std.plot(datafolder, phenotypename,markers,chrid,genestart,geneend)
}

