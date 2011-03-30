source('gex.header.R')
source('methods.R')
selected=which((markers[,3]>=rangestart) & (markers[,3]<=rangeend) & (markers[,2]==chrid))
XX=X[selected,]

if (fun=='analysis'){
  std.analysis(datafolder, phenotypename,X,Y)
}else if (fun=='plot'){
  std.plot(datafolder, phenotypename,markers,chrid,genestart,geneend)
  std.pvalue.hist.plot(datafolder, phenotypename,markers,chrid,genestart,geneend)
}


