library(matlab)

genotypefolder="result/precc_liver_gene_expression/genotype/"
markers <- data.matrix(read.table(paste(genotypefolder,'marker_list',sep="")))
names <- read.table('result/precc_liver_gene_expression/gene_expression/gene_list')
info <- read.csv('result/precc_liver_gene_expression/gene_expression/gene_info')
gene.index <- data.frame()
gene.index <- info$Chromosome.Name
names(gene.index) <- info$Ensembl.Gene.ID
gene.position <- data.frame()
gene.position <- info$Gene.Start..bp.
names(gene.position) <- info$Ensembl.Gene.ID
right <- 0
all <- 0
filelist <- c()
namelist <- c()
vfilelist <- c()
for (i in 1:dim(names)[1]){
  filename=paste('result/precc_liver_gene_expression/std/',names[i,2],'/',names[i,1],'.Rdata',sep='')
  if (file.exists(filename)){
    filelist <- c(filelist,filename)
    namelist <- c(namelist,paste(names[i,1]))
    vfilelist <- c(vfilelist,paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_global_nlminb.Rdata',sep=''))
  }    
}

gex.variance.summary.chr <- function(){
  for (i in 1:dim(names)[1]){
    filename=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_chr.Rdata',sep='')
    if (file.exists(filename)){
      load(filename)
      m=apply(result[60:(60+test*10),1:19],2,mean);
      m=result[1,]
      order=sort.int(m,decreasing=T,index.return=T)
      if(chrid %in% order$ix[1:3]){
        right <- right+1
      }
      all <- all+1
    }
  }
  list(right,all)
}
gex.variance.summary.emma.chr <- function(){
  for (i in 1:dim(names)[1]){
    filename=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_emma_chr.Rdata',sep='')
    if (file.exists(filename)){
      load(filename)
      m <- rep(0,19)
      for ( i in 1:19){
        m[i]=result[[i]]$vgs
      }
      order=sort.int(m,decreasing=T,index.return=T)
      if(chrid %in% order$ix[1:1]){
        right <- right+1
      }
      all <- all+1
    }
  }
   list(right,all)
}
gex.variance.summary.10M <- function(){
  load('variance.10M.kinship.Rdata')
  for (k in 1:length(namelist)){
    #filename=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_10M_other.Rdata',sep='')
    filename <- vfilelist[k]
    if (file.exists(filename)){
      load(filename)
      #print(dim(result[4,,]))
      m=apply(result[2,,1:(dim(result)[3]-1)],2,mean);
      #print(dim(m))
      order=sort.int(m,decreasing=T,index.return=T)
      chrid=gene.index[paste(namelist[k])]
      pos=gene.position[paste(namelist[k])]
      order=sort.int(m,decreasing=T,index.return=T)
      
      for (s in 1:1){
        i <- order$ix[s]
        if (chrid==chr[i] && startpos[i]-10000000 <= genestart && startpos[i]+20000000>=genestart){
          right <- right+1
        }
      }
      all <- all+1
    }
  }
   list(right,all)
}

gex.variance.summary.local.bfgs <- function(){
  load('variance.10M.kinship.Rdata')
  t <- NA
  for (k in 1:length(namelist)){
    #filename=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_10M_other.Rdata',sep='')
    filename <- vfilelist[k]
    if (file.exists(filename)){
      load(filename)
      #print(dim(result[4,,]))
      m=(result[15,1:(dim(result)[2]-1)])
      #m=result[-length(result)]
      #print(dim(m))
      order=sort.int(m,decreasing=T,index.return=T)
      chrid=gene.index[paste(namelist[k])]
      pos=gene.position[paste(namelist[k])]
      order=sort.int(m,decreasing=T,index.return=T)
      print(c(order$x[1],order$ix[1],order$ix[2]))
      for (s in 1:1){
        i <- order$ix[s]
        if (chrid==chr[i] && startpos[i]-10000000 <= genestart && startpos[i]+20000000>=genestart){
          right <- right+1
          if (is.nan(t)){
           p <- c(genestart+10000000-startpos[i])
          }else {
            p <- c(p,genestart+10000000-startpos[i])
          }
          break
        }
      }
      all <- all+1
    }
  }
   list(right,all)
}

gex.variance.summary.emma.10M <- function(){
  load('variance.10M.kinship.Rdata')
  for (i in 1:dim(names)[1]){
    filename=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_emma_10M.Rdata',sep='')
    if (file.exists(filename)){
      load(filename)
      nk=dim(K)[1]
      m <- rep(0,nk)
      for ( i in 1:nk){
        m[i]=result[[i]]$vgs
      }
      order=sort.int(m,decreasing=T,index.return=T)
      i <- which.max(m)
      if (chrid==chr[i] && startpos[i] <= genestart && startpos[i]+10000000>=genestart){
        right <- right+1
      }
      all <- all+1
    }
  }
  list(right,all)
}

gex.std.summary.chr <- function(){
  for (i in 1:dim(names)[1]){
    filename=paste('result/precc_liver_gene_expression/std/',names[i,2],'/',names[i,1],'.Rdata',sep='')
    if (file.exists(filename)){
      load(filename)
      r <- -log(p)
      t <- which.max(-log(p))
      if (r[t]>15){
        chrid=gene.index[paste(names[i,1])]
        pos=gene.position[paste(names[i,1])]
        if (chrid==markers[t,2] && pos-15000000<=markers[t,3] && pos+15000000>=markers[t,3]){  
          right <- right+1
        }
        all <- all+1
      }
    }
  }
  list(right,all)
}

gex.plot.eqtl <- function(){
# generate gene list
  filelist <- c()
  namelist <- c()
  chrlist <- c()
  poslist <- c()
  vfilelist <- c()
  idx <- 1
  step=10000000
  groups=rep(0,dim(markers)[1])
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
  for (i in 1:dim(names)[1]){
    filename=paste('result/precc_liver_gene_expression/std/',names[i,2],'/',names[i,1],'.Rdata',sep='')
    filename1=paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_global_nlminb.Rdata',sep='')
    if (file.exists(filename) && file.exists(filename1)){
      filelist <- c(filelist,filename)
      vfilelist <- c(vfilelist,paste('result/precc_liver_gene_expression/variance/',names[i,2],'/',names[i,1],'_global_nlminb.Rdata',sep=''))
      namelist <- c(namelist,paste(names[i,1]))
      chrlist <- c(chrlist,gene.index[paste(names[i,1])])
      poslist <- c(poslist,gene.position[paste(names[i,1])])
    }    
  }

  orders=order(chrlist[namelist],poslist[namelist])
  filelist <- filelist[orders]
  chrlist <- chrlist[orders]
  poslist <- poslist[orders]
  vfilelist <- vfilelist[orders]
  mycolors <- jet.colors(16)
  #mycolors[-16] <- mycolors[1]
  if (!file.exists("eqtl_result/std_plot_max_2.pdf")){
    Img <- mat.or.vec(length(filelist),idx-1)
    for ( i in 1:length(filelist)){
      load(filelist[i])
      data <- -log(p)
      t <- ((as.numeric(tapply(data,groups,max))))
      t <- t/max(t)
      Img[i,]=t[]
    }
    save(Img, file="temp.Rdata")
    pdf("eqtl_result/std_plot_max_2.pdf")
    imagesc(Img*16, xlab = "pos", ylab = "gene", col = mycolors )
    dev.off()
  }
  if (!file.exists("eqtl_result/true_plot.pdf")){
    Img <- mat.or.vec(length(filelist),idx-1)
    for ( i in 1:length(filelist)){
      t <- mat.or.vec(1,idx-1)
      right=min(which((markers[,3]>=poslist[i]) & (markers[,3]<poslist[i]+step) & (markers[,2]==chrlist[i])))
      t[groups[right]]=1
      print(dim(Img))
      print(length(t))
      Img[i,]=t[]
    }
    save(Img, file="eqtl_result/true_plot.Rdata")
    pdf("eqtl_result/true_plot.pdf")
    imagesc(Img*2, xlab = "pos", ylab = "gene", col = jet.colors(2))
    dev.off()
  }

  if (!file.exists("eqtl_result/variance_10M_plot_global_nlminb_2.pdf")){
    Img <- mat.or.vec(length(vfilelist),idx)
    for ( i in 1:length(vfilelist)){
      load(vfilelist[i])
#      t=apply(result[4,10:50,-dim(result)[3]],2,mean);
      #t=(result[15,-dim(result)[2]])
      t=(result[])
      if (max(t) != 0)
        t <- t/max(t)
      else
        t[] <- 0
      Img[i,]=t[]
    }
    print(Img)
    save(Img, file="eqtl_result/variance_10M_plot_global_nlminb_2.Rdata")
    pdf("eqtl_result/variance_10M_plot_global_nlminb_2.pdf")
    imagesc(Img*16, xlab = "pos", ylab = "gene", col =mycolors)
    dev.off()
  }
}

#print("The standard t-test")
#print(gex.std.summary.chr())
# 762/1119
#print("The variance 10M test")
#print(gex.variance.summary.emma.chr())
#print(gex.variance.summary.10M())
print(gex.variance.summary.local.bfgs())
(gex.plot.eqtl())
