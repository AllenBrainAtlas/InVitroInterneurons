library(limma)
DE.genes <- function(dat,cl){  
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  de.df=list()
  for(x in tmp.cl){
    ctr <<- paste(x, "- (", paste(setdiff(tmp.cl,x), collapse="+"),")/",(length(tmp.cl)-1))
    contrasts.matrix <- makeContrasts(ctr,  levels=design)
    fit2 = contrasts.fit(fit, contrasts.matrix)
    fit2 = eBayes(fit2)
    padj = apply(fit2$p.value, 2, p.adjust)
    lfc = coef(fit2)
    mA = rowMeans(dat[, cl==x,drop=F])
    mB = rowMeans(dat[, cl!=x,drop=F])
    de.df[[x]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1],meanA=mA, meanB=mB)
    row.names(de.df[[x]])= row.names(dat)
  }
  return(de.df)
}


DE.genes.pw <- function(dat,cl){  
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  cl.means = tapply(1:ncol(dat), cl, function(x){
    rowMeans(dat[,x,drop=F])
  })
  cl.means = cl.means[tmp.cl]
  de.df=list()
  for(i in 1:(length(tmp.cl)-1)){
    for(j in (i+1):length(tmp.cl)){
      x=tmp.cl[i]
      y=tmp.cl[j]
      ctr <<- paste(x, "- ", y)
      contrasts.matrix <- makeContrasts(ctr,  levels=design)
      fit2 = contrasts.fit(fit, contrasts.matrix)
      fit2 = eBayes(fit2)
      padj = apply(fit2$p.value, 2, p.adjust)
      lfc = coef(fit2)
      pair=paste(x,y,sep="_")
      de.df[[pair]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1],meanA=cl.means[[i]], meanB=cl.means[[j]])
      row.names(de.df[[pair]])= row.names(dat)
    }
  }
  return(de.df)
}

DE.genes.one.vs.other <- function(dat,cl, x.cl){  
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  cl.means = tapply(1:ncol(dat), cl, function(x){
    rowMeans(dat[,x,drop=F])
  })
  cl.means = cl.means[tmp.cl]
  de.df=list()
  for(y.cl in setdiff(tmp.cl, x.cl)){
    x=min(x.cl, y.cl)
    y=max(x.cl, y.cl)
    ctr <<- paste(x, "- ", y)
    contrasts.matrix <- makeContrasts(ctr,  levels=design)
    fit2 = contrasts.fit(fit, contrasts.matrix)
    fit2 = eBayes(fit2)
    padj = apply(fit2$p.value, 2, p.adjust)
    lfc = coef(fit2)
    pair=paste(x,y,sep="_")
    de.df[[pair]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1],meanA=cl.means[[x]], meanB=cl.means[[y]])
    row.names(de.df[[pair]])= row.names(dat)
  }
  return(de.df)
}



DE.genes.groups <- function(dat,cls1, cls2){
  cl = c(cls1, cls2)
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  dat = dat[,names(cl)]
  fit = lmFit(dat , design=design)
  de.df=list()
  tmp.cl1 = sort(unique(cls1))
  tmp.cl2 = sort(unique(cls2))
  comb.cl1 = paste("(", paste(tmp.cl1, collapse="+"), ")/",length(tmp.cl1))
  comb.cl2 = paste("(", paste(tmp.cl2, collapse="+"), ")/",length(tmp.cl2))
  for(x in c(tmp.cl1,tmp.cl2)){
    if(x %in% tmp.cl1){
      ctr <<- paste(x, "- ", comb.cl2)
    }
    else{
      ctr <<- paste(x, "- ", comb.cl1)
    }
    contrasts.matrix <- makeContrasts(ctr,  levels=design)
    fit2 = contrasts.fit(fit, contrasts.matrix)
    fit2 = eBayes(fit2)
    padj = apply(fit2$p.value, 2, p.adjust)
    lfc = coef(fit2)
    pair=paste(x,y,sep="_")
    de.df[[x]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1])
    row.names(de.df[[x]])= row.names(dat)
  }
  return(de.df)
}



DESeq.genes<- function(dat,cl,...)
  {
    require("DESeq2")
    require("parallel")
    df= data.frame(cl=cl)
    dds = DESeqDataSetFromMatrix(dat, colData=df, design=~ cl)
    dds=DESeq(dds,...)
    df = results(dds)
    colnames(df)[2] = "lfc"
    return(df)
  }

DESeq.genes.pw <- function(dat,cl, dds.file="dds.rda",mc.cores=4){
  require("DESeq2")
  require("parallel")
  if(is.null(dds.file) || !file.exists(dds.file)){
    df= data.frame(cl=cl)
    dds = DESeqDataSetFromMatrix(dat, colData=df, design=~ cl)
    dds=DESeq(dds)
    if(!is.null(dds.file)){
      save(dds, file=dds.file)
    }
  }
  tmp.cl = sort(unique(cl))
  pairs = data.frame(X=as.character(rep(tmp.cl,length(tmp.cl))), Y=as.character(rep(tmp.cl, rep(length(tmp.cl), length(tmp.cl)))), stringsAsFactors=F)
  pairs = pairs[pairs[,1]<pairs[,2],]
  de.df=sapply(1:nrow(pairs), function(i){
    print(pairs[i,])
    x=pairs[i,1]
    y=pairs[i,2]
    res=results(dds, contrast = c("cl", x,y))
    colnames(res)[2]="lfc"
    res
  },simplify=F)
  names(de.df) = paste(pairs[,1],pairs[,2],sep="_")
  return(de.df)
}
