source("de.genes.R")
library(WGCNA)
library(limma)
library(mclust)
library(flashClust)
library(matrixStats)

findCor <- function(gene, dat, cells, n=10)
{
  tmp=cor(t(dat[,cells]), dat[gene,cells])
  tmp = tmp[order(tmp[,1],decreasing=T),]
  head(tmp, n)
}

MyTOM <- function(adj,sparse.cutoff=10^-4)
{
  require(Matrix)
  adj[adj < sparse.cutoff] = 0
  adj = Matrix(adj, sparse=T)
  print("matrix multiplication")
  olap= crossprod(adj)
  print("normalize")
  adj.counts = rowSums(adj)-1
  nsize = nrow(adj)
  TOM= .C("tomSimilarity", as.double(as.matrix(olap)), as.double(as.matrix(adj)), as.integer(nsize))
  tmp1 = matrix(rep(adj.counts, ncol(olap)), nrow=length(adj.counts))
  N = pmin(tmp1, t(tmp1))
  TOM = (olap - adj)/(N + 1 - adj)
  diag(TOM)=1
  TOM
}


getGeneTree <- function(norm.dat, ercc_counts, counts, select.cells, min.cells=10, padj.th=0.05, prefix="cl", softPower=4,select.genes=NULL,type="unsigned",method="average",maxGenes=6000)
  {
    if(is.null(select.genes)){
      tmp.cells = select.cells[colMedians(ercc_counts[,select.cells]) > 0]
      if(length(tmp.cells)>20){
        vg = run_brennecke(ercc_counts[,tmp.cells],counts[,tmp.cells])#,plot.fig=paste0("DESeq.",prefix,".pdf"))
      }
      else{
        vg = run_brennecke(counts[,select.cells],counts[,select.cells])#,plot.fig=paste0("DESeq.",prefix,".pdf"))
      }
      select.genes = row.names(norm.dat)[which(rowSums(norm.dat[,select.cells] > 0) >= min.cells & vg$padj < padj.th)]
      select.genes = head(select.genes[order(vg[select.genes, "padj"],-vg[select.genes, "cv2"])],maxGenes)
      cat("Select genes", length(select.genes),"\n")
    }
    if(length(select.genes) < 5){return(NULL)}
    dat = norm.dat[select.genes,select.cells]
    adj=adjacency(t(dat), power = softPower,type=type)
    adj[is.na(adj)]=0
    TOM = TOMsimilarity(adj,TOMType=type)
    dissTOM = as.matrix(1-TOM)
    row.names(dissTOM)= colnames(dissTOM) = row.names(dat)
    geneTree = flashClust(as.dist(dissTOM), method = method)
    #save(geneTree, dissTOM, file=paste(prefix, "geneTree.rda",sep="."))
    return(list(geneTree, dissTOM))
  }

scoreGeneMod <-  function(norm.dat, select.cells, gene.mod, min.cells=3, padj.th=0.01, lfc.th=1,method="average"){
    if(length(gene.mod)==0){
      return(NULL)
    }
    gene.mod.val=sapply(names(gene.mod), function(y){
      x= gene.mod[[y]]
      tmp.dat = norm.dat[x,select.cells]
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      if(method=="average"){
        tmp.cl = cutree(hclust(dist(t(tmp.dat)),method="average"),2)
      }
      else if(method=="ward"){
        tmp.cl = cutree(hclust(dist(t(tmp.dat)),method="ward"),2)
      }
      else if(method=="kmeans"){
        tmp.cl = kmeans(t(tmp.dat), 2)$cluster
      }
      else if(method=="mclust"){
        tmp = Mclust(t(tmp.dat), G=2)
        if(!is.na(tmp[[1]])){
          tmp.cl = tmp$classification
        }
        else{
          return(list(c(0,0,0),NULL))
        }
      }
      tb=table(tmp.cl)
      if(min(tb) < min.cells){
        return(list(c(min(tb), 0,0),NULL))
      }
      de.df = DE.genes.pw(tmp.dat, paste("cl",tmp.cl, sep=""))[[1]]
      de.df$padj[which(de.df$padj < 10^-20)]=10^-20
      select = with(de.df, which(padj < padj.th & abs(lfc)> lfc.th))
      sc = -sum(log10(de.df$padj)[select])
      return(list(c(min(tb), sc, length(select)),de.df))         
    },simplify=F)
    return(gene.mod.val)
  }


tomHC <- function(dat,power=4,method="ward")
  {
    adj = adjacency(dat, power = power,type="signed")
    tom = TOMsimilarity(adj, TOMType="signed")
    rm(adj)
    gc()
    row.names(tom)=colnames(tom)= colnames(dat)
    return(tom)
  }



cutGeneTree <- function(norm.dat, select.cells, geneTree,  dissTOM,  minModuleSize=20, cutHeight=0.99,plotMod=T,...)
  {
                                        # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight=cutHeight,
      deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
    gene.mod = split(row.names(dissTOM), dynamicMods)
    gene.mod = gene.mod[setdiff(names(gene.mod),"0")]
    if(is.null(gene.mod)){return(NULL)}
    filterGeneMod(norm.dat, select.cells, gene.mod, minModuleSize=minModuleSize, plotMod=plotMod,...)
  }

filterGeneMod <- function(norm.dat, select.cells, gene.mod, minModuleSize=10,min.cells=3, padj.th= 0.01, lfc.th = 1, min.padj=40, maxSize=200, bykME=FALSE,prefix="cl", max.mod=NULL,plotMod=T,...)
  {
    if(length(select.cells) > 1000){
      method="kmeans"
    }
    else{
      method=c("average","ward","kmeans")
    }
    gene.mod.val = sapply(method, function(m){
      tmp=scoreGeneMod(norm.dat, select.cells, gene.mod, min.cells, padj.th, lfc.th,method=m)
    },simplify=F)
    save(gene.mod.val, gene.mod, file=paste0(prefix,".gene.mod.rda"))
    mod.score = rep(0, length(gene.mod))
    select.mod.val = list()
    for(val in gene.mod.val){
      x = do.call("cbind", sapply(val, function(x)x[[1]],simplify=F))
      tmp=which(x[1,] >= min.cells & x[2,] > min.padj & x[3,] > minModuleSize & x[2,] > mod.score)
      mod.score[tmp] = x[2, tmp]
      select.mod.val[tmp] = val[tmp]
    }
    select.mod = mod.score > 0
    if(sum(select.mod)>0){
      gene.mod = gene.mod[select.mod]
      gene.mod.val = select.mod.val[select.mod]
      gene.mod.sc = do.call("cbind", sapply(gene.mod.val, function(x)x[[1]],simplify=F))
      ord = order(gene.mod.sc[2,],decreasing=T)
      if(!is.null(max.mod)){
        ord = head(ord, max.mod)
      }
      gene.mod = gene.mod[ord]
      gene.mod.val = gene.mod.val[ord]
      gene.mod.sc = gene.mod.sc[,ord]
      
      tmp.dat = norm.dat[unlist(gene.mod),select.cells]
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      #colnames(tmp.dat) = sample.df[colnames(tmp.dat),"short_name"]
    
      for(x in 1:length(gene.mod)){
        g = gene.mod[[x]]
        if(length(g) > maxSize){
          if(bykME){
            eigen = prcomp(t(norm.dat[g,]))$x[,1]
            kME=cor(t(norm.dat[g,]), eigen)
            g = head(g[order(abs(kME),decreasing=T)], maxSize)
            gene.mod[[x]] <- g
          }
          else{
            de.df = gene.mod.val[[x]][[2]]
            g = head(g[order(de.df[g, "padj"])], maxSize)
            gene.mod[[x]] <- g
          }
        }
      }
      if(plotMod){
        plotGeneMod(norm.dat, select.cells,gene.mod,prefix=prefix,...)
      }
      return(list(gene.mod, gene.mod.sc))
    }
    return(NULL)
  }


plotGeneMod <- function(norm.dat, select.cells,gene.mod,prefix="cl",method="ward",...)
  { 
    pdf(paste(prefix,"gene.mod.pdf",sep="."), height=9, width=9)
    for(x in 1:length(gene.mod)){
      g = gene.mod[[x]]
      tmp.dat = norm.dat[g,select.cells]
      tmp.dat = tmp.dat - rowMeans(tmp.dat)
      hc = hclust(dist(t(tmp.dat[g,])), method=method)
      gene.hc = hclust(dist(tmp.dat[g,]), method=method)
      cexCol = min(70/ncol(tmp.dat),1)
      cexRow = min(60/length(g),1)
      heatmap.3(tmp.dat[g,],  Colv=as.dendrogram(hc), Rowv=as.dendrogram(gene.hc),col=blue.red(100), trace="none",  dendrogram="column", cexCol=cexCol,cexRow=cexRow,...)
    }
    dev.off()
  }



dePair <- function(df,dat1, dat2, padj.th=0.01, lfc.th=1, wilcox.th=NULL, q=NULL, q1=NULL, q2=NULL)
  {
    df = df[order(df$pval),]
    df$padj[df$padj < 10^-20]=10^-20
    select=with(df, padj < padj.th & abs(lfc)>lfc.th)
    select=row.names(df)[which(select)]
    select = intersect(select, row.names(dat1))
    if(!is.null(wilcox.th)){
      pairs = unlist(strsplit(s, "_"))
      pairs = gsub("cl", "",pairs)
      tmp.wilcox = sapply(select, function(i){
        tmp=wilcox.test(dat1[i,], dat2[i,])
        tmp$p.value
      })
      select=select[tmp.wilcox < wilcox.th]
    }
    if(is.null(select) | length(select)==0){
      return(list(score=0, num=0, genes= NULL))
    }
    if (!is.null(q)){
      if(is.null(q1)){
        q1=rowQuantiles(dat1[select,,drop=F], probs=c(1-q,q))
      }
      if(is.null(q2)){
        q2=rowQuantiles(dat2[select,,drop=F], probs=c(1-q,q))
      }
      if(length(select)>1){
        select1= q1[select,2] +0.001 < q2[select,1]| q2[select,2] + 0.001 < q1[select,1]
      }
      else{
        select1= q1[2] +0.001 < q2[1]| q2[2] + 0.001 < q1[1]      
      }
      select = select[select1]
    }
    if(length(select)==0){
      list(score=0, num=0, genes= NULL)
    }
    else{
      list(score=sum(-log10(df[select,"padj"])), num=length(select), genes= select)
    }
  }


deScore <- function(norm.dat, cl, select.cells=names(cl), de.df=NULL, min.cells=4, ...)
  {
    cl.size = table(cl)
    cl.n = names(cl.size)
    cl.small = cl.n[cl.size < min.cells]
    cl.big =  setdiff(cl.n,cl.small)
    select.cells =names(cl)[!cl %in% cl.small]
    cl = cl[select.cells]
    de.genes=list()
    select.genes= row.names(norm.dat)[rowSums(norm.dat > 0) > min.cells]
    if(length(cl.big)>1){
      if(is.null(de.df)){
        de.df = DE.genes.pw(norm.dat[select.genes,select.cells], paste0("cl",cl[select.cells]))
      }
      ##merge clusters are have too few de.genes
      de.genes = sapply(names(de.df), function(x){
        if(is.null(de.df[[x]])){
          return(list(score=0, num=0, genes=NULL))          
        }
        pair = gsub("cl", "", unlist(strsplit(x, "_")))
        i = as.integer(pair[1])
        j = as.integer(pair[2])
        dat1 = norm.dat[select.genes,select.cells[cl==i]]
        dat2 = norm.dat[select.genes,select.cells[cl==j]]
        tmp=dePair(de.df[[x]],dat1, dat2, ...)
      },simplify=F)
    }
    for(i in cl.small){
      for(j in cl.n[cl.n!=i]){
        x = paste0("cl",min(i,j))
        y = paste0("cl",max(i,j))
        pair=paste(x,y,sep="_")
        de.genes[[pair]]=list(score=0, num=0, genes=NULL)
        de.df[[pair]]=list()
      }
    }
    list(de.df, de.genes)
  }


deScore.one.vs.other <- function(norm.dat, cl, x.cl, select.cells=names(cl), de.df=NULL, min.cells=4, ...)
  {
    cl.size = table(cl)
    cl.n = names(cl.size)
    cl.small = cl.n[cl.size < min.cells]
    cl.big =  setdiff(cl.n,cl.small)
    select.cells =names(cl)[cl %in% cl.big]
    de.genes=list()
    if(x.cl %in% cl.big & length(select.cells)>0 & length(cl.big)>1){
      cl = cl[select.cells]
      if(is.null(de.df)){
        de.df = DE.genes.one.vs.other(norm.dat[,select.cells], paste0("cl",cl[select.cells]), paste0("cl",x.cl))
      }
      de.genes = sapply(names(de.df), function(s){
        if(is.null(de.df[[s]])){
          return(list(score=0, num=0, genes=NULL))
        }
        pair = gsub("cl", "", unlist(strsplit(s, "_")))
        i = as.integer(pair[1])
        j = as.integer(pair[2])
        dat1 = norm.dat[,select.cells[cl==i]]
        dat2 = norm.dat[,select.cells[cl==j]]
        dePair(de.df[[s]],dat1, dat2, ...)
      },simplify=F)
    }
    if(x.cl %in% cl.small){
      for(i in cl.big){
        x = paste0("cl",min(i,x.cl))
        y = paste0("cl",max(i,x.cl))
        pair=paste(x,y,sep="_")
        de.genes[[pair]]=list(score=0, num=0, genes=NULL)
        de.df[[pair]]=list()
      }
    }
    for(i in setdiff(cl.small,x.cl)){
      x = paste0("cl",min(i,x.cl))
      y = paste0("cl",max(i,x.cl))
      pair=paste(x,y,sep="_")
      de.genes[[pair]]=list(score=0, num=0, genes=NULL)
      de.df[[pair]]=list()
    }
    list(de.df, de.genes)
  }


mergeCl <- function(norm.dat, cl, hc=NULL, min.padj=40, min.genes=5, min.cells=3, eigen=NULL, rm.eigen=NULL,rm.th=0.3, ...)
  {
    select.cells = names(cl)
    if(!is.null(rm.eigen)){
      rm.kME=cor(t(norm.dat[,select.cells]), rm.eigen)
      rm.kME[is.na(rm.kME)] =0
      select = rowMaxs(rm.kME) < rm.th
      if (!is.null(eigen)){
        kME=cor(t(norm.dat[,select.cells]), eigen)
        kME[is.na(kME)] =0          
        select = select | rowMaxs(abs(kME)) > rowMaxs(rm.kME)
      }
      #tmp = row.names(rm.kME)[!select]
      #tmp = tmp[order(rowMaxs(rm.kME[tmp,]),decreasing=T)]
      norm.dat=norm.dat[select,]
    }
    de.score = deScore(norm.dat, cl, select.cells, min.cells=min.cells,...)
    gc()
    de.df= de.score[[1]]
    de.genes=de.score[[2]]
    if(is.null(de.genes)|length(de.genes)==0){
      return(NULL) 
    }
    else{
      sc = sapply(de.genes, function(x)x$score)
      num = sapply(de.genes, function(x)x$num)
      to.merge = sc < min.padj | num < min.genes    
      pairs = names(to.merge)[to.merge]
      print(sort(sc))
    }
    print(table(cl))
    print(sc[pairs])
    pairs=gsub("cl", "",pairs)
    pairs = do.call(rbind,strsplit(pairs, "_"))
    n.cl = length(unique(cl))
    while(length(pairs)>0 & n.cl > 1){
      if(!is.null(hc)){
        ###Determine which clusters can be merged at the current level
        tmp.cl = cutree(hc, n.cl-1)
        tb= table(cl, tmp.cl)
        tmp.df = as.data.frame(tb)
        tmp.df=tmp.df[tmp.df$Freq > 0,]
        merge.map = setNames(tmp.df$tmp.cl, tmp.df$cl)
        merge.pairs = pairs[merge.map[pairs[,1]]==merge.map[pairs[,2]],,drop=F]
        ##Only merge the most smilar pair one time.
      }
      else{
        merge.pairs = pairs
      }
      if(nrow(merge.pairs)>0){
        to.merge= merge.pairs[1,]
        cat("Merge ",to.merge[1], to.merge[2],"\n")
        tmp=do.call(rbind,strsplit(gsub("cl", "",names(de.df)), "_"))
        select = tmp[,1] %in%  to.merge | tmp[,2] %in% to.merge
        de.df = de.df[!select]
        de.genes = de.genes[!select]
        cl[cl==to.merge[2]] = to.merge[1]
        ###Update de.score
        if(length(unique(cl)) < 2){
          return(NULL)
        }
        tmp = deScore.one.vs.other(norm.dat, cl, x=to.merge[1], min.cells=min.cells,...)
        de.df = c(de.df, tmp[[1]])
        de.genes= c(de.genes, tmp[[2]])
        print(names(de.df))
        sc = sapply(de.genes, function(x)x$score)
        num = sapply(de.genes, function(x)x$num)
        to.merge = sc < min.padj | num < min.genes
        print(table(cl))
        print("To merge")
        print(sc[to.merge])
        ord = order(sc[to.merge])
        pairs = names(to.merge)[to.merge][ord]        
        pairs=gsub("cl", "",pairs)
        pairs = do.call(rbind,strsplit(pairs, "_"))
        if(is.null(hc)){
          n.cl = length(unique(cl))
        }
      }
      else{
        if(!is.null(hc)){
          n.cl = n.cl-1
        }
      }
    }
    return(list(cl=cl, de.df=de.df,sc=sc))
  }


selectMarkers <- function(norm.dat, cl, de.df, n.markers=20, ...)
  {
    de.score=deScore(norm.dat, cl, de.df=de.df, ...)[[2]]
    select.genes = unique(unlist(sapply(names(de.score), function(s){
      genes = de.score[[s]]$genes
    },simplify=F)))
    de.genes = sapply(names(de.score), function(s){
      genes = de.score[[s]]$genes
      genes =genes[genes %in% select.genes]
      up = genes[de.df[[s]][genes, "lfc"] > 0]
      down = genes[de.df[[s]][genes, "lfc"] < 0]
      c(head(up,n.markers), head(down,n.markers))
    },simplify=F)
    markers = unique(unlist(de.genes))
  }

displayCl <- function(norm.dat, cl, de.df, prefix=NULL,hc=NULL, gene.hc=NULL,markers=NULL, centered=T,labels=names(cl),rm.eigen=NULL,eigen=NULL,kME.th=0.3,rm.kME.th=0.7, by.cl=FALSE,ColSideColors=NULL,maxValue=5,min.sep=10,main="",default.markers=NULL,...)
  {
    select.cells = names(cl)
    select = rep(TRUE, nrow(norm.dat))
    de.df=de.df[sapply(de.df, length)>0]
    if (is.null(markers)){
      norm.dat= norm.dat[row.names(de.df[[1]]),]
      if(!is.null(eigen)){
        kME=cor(t(norm.dat[,select.cells]), eigen)
        kME[is.na(kME)] =0
        select = rowMaxs(abs(kME)) > kME.th
      }
      if(!is.null(rm.eigen)){
        rm.kME=cor(t(norm.dat[,select.cells]), rm.eigen)
        rm.kME[is.na(rm.kME)] =0
        select =  select & rowMaxs(abs(rm.kME))< rm.kME.th
        if(!is.null(eigen)){
          select = select & rowMaxs(abs(kME)) > rowMaxs(abs(rm.kME))
        }
      }
      if(sum(!select)>0){
        de.df = sapply(de.df, function(x)x[select,],simplify=F)
        norm.dat = norm.dat[select,]
      }
      markers <- selectMarkers(norm.dat, cl, de.df, ...)
    }
    if(!is.null(markers) & !is.null(prefix)){
      markers=union(markers, default.markers)
      tmp.dat = norm.dat[markers,select.cells]
      if(centered){
        tmp.dat = tmp.dat - rowMeans(tmp.dat)
        breaks=c(min(min(tmp.dat)-0.1,-maxValue),  seq(-maxValue,maxValue, length.out=99), max(max(tmp.dat)+1))
      }
      else{
        tmp.dat = tmp.dat/rowMaxs(tmp.dat)
        breaks=c(0, seq(0.05, 1, length.out=100))
      }
      colnames(tmp.dat)=labels
      dendro = "column"
      cexCol = min(70/ncol(tmp.dat),1)
      cexRow = min(60/nrow(tmp.dat),1)
      if(is.null(gene.hc)){
        gene.hc = flashClust(dist(tmp.dat), method="ward")
      }
      else{
        dendro="both"
      }
      if(is.null(hc)){
        hc = flashClust(dist(t(tmp.dat)), method="ward")
      }

      pdf(paste(prefix,"pdf",sep="."), height=13, width=9)
      if(by.cl){
        ord = order(cl, order(hc$order))
        sep = cl[ord]
        sep=which(sep[-1]!=sep[-length(sep)])
        sep = c(sep[1], sep[which(sep[-1] - sep[-length(sep)] >=min.sep)+1])
        heatmap.3(tmp.dat[,ord],Rowv=as.dendrogram(gene.hc), Colv=NULL, col=blue.red(100), trace="none", dendrogram="none", cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors[,ord],breaks=breaks,colsep=sep, sepcolor="black",main=main)
      }
      else{
        heatmap.3(tmp.dat,Rowv=as.dendrogram(gene.hc), Colv=as.dendrogram(hc), col=blue.red(100), trace="none", dendrogram=dendro, cexCol=cexCol,cexRow=cexRow,ColSideColors=ColSideColors,breaks=breaks,main=main)
      }
      dev.off()        
    }
    return(markers)
  }
    
    
getEigen <- function(gene.mod, norm.dat, select.cells, prefix=NULL,method="ward",hc=NULL,...)
  {
    gene.vector = setNames(rep(names(gene.mod), sapply(gene.mod, length)), unlist(gene.mod))
    eigen = moduleEigengenes(t(norm.dat[names(gene.vector),select.cells]), gene.vector)[[1]]
    eigen=eigen[,paste0("ME",names(gene.mod)),drop=F]
    row.names(eigen)= select.cells
    kME=cor(t(norm.dat[names(gene.vector),select.cells]), eigen)
    hub=sapply(names(gene.mod), function(i){
      x = gene.mod[[i]]
      x[which.max(kME[x, paste0("ME",i)])]
    })
    hub = setNames(hub, paste0("ME",names(gene.mod)))
    colnames(eigen)=paste(colnames(eigen),hub[colnames(eigen)])
    row.names(eigen)=select.cells
    if(is.null(hc)){
      hc = hclust(dist(eigen), method=method)
    }
    if(!is.null(prefix) & ncol(eigen)>1){
      colv = as.dendrogram(hc)
      pdf(paste0(prefix, ".eigen.pdf"),height=6,width=7)
      heatmap.3(t(eigen), Colv= colv, col=blue.red(100), trace="none", dendrogram="column",...)
      dev.off()
    }
    return(list(eigen,hc))
  }


select.markers  <-  function(norm.dat, cl, de.df=NULL, n=2, q=0.8,padj.th = 0.01,wilcox.th=0.01)
  {
    select.cells = names(cl)
    if(is.null(de.df)){
      de.df = DE.genes.pw(norm.dat[,select.cells], paste("cl",cl,sep=""))
    }
    names(de.df) = gsub("cl", "", names(de.df))
    de.genes = sapply(names(de.df), function(s){
      pairs = unlist(strsplit(s, "_"))
      x = de.df[[s]]
      select=with(x, padj < padj.th & abs(lfc)>1)
      tmp=row.names(x)[which(select)]
      tmp.wilcox = sapply(tmp, function(i){
        tmp=wilcox.test(norm.dat[i,names(cl)[cl==pairs[1]]], norm.dat[i,names(cl)[cl==pairs[2]]])
        tmp$p.value
      })
      tmp[tmp.wilcox < wilcox.th]
    },simplify=F)
    select.genes = sort(unique(unlist(de.genes)))
    cl.cells = split(select.cells, as.character(cl[select.cells]))
    tmp.dat = norm.dat[select.genes,]
    select.genes.q = sapply(cl.cells, function(x){
      print(length(x))
      rowQuantiles(tmp.dat[,x,drop=F], q)
    })
    select.genes.lq = sapply(cl.cells, function(x){
      print(length(x))
      rowQuantiles(tmp.dat[,x,drop=F], 1-q)
    })
    row.names(select.genes.q) = row.names(select.genes.lq) = select.genes
    sig = sapply(names(de.df), function(s){
      x = de.df[[s]][select.genes,]
      pair = unlist(strsplit(s, "_"))
      print(pair)
      q1= select.genes.q[select.genes, pair[[1]]]
      q2= select.genes.q[select.genes, pair[[2]]]
      lq1= select.genes.lq[select.genes, pair[[1]]]
      lq2= select.genes.lq[select.genes, pair[[2]]]
      select = (q1 < lq2 | q2 < lq1)
      select = select & with(x, abs(lfc) > 1 & padj < padj.th)
      y = setNames(rep(0, length(select.genes)), select.genes)
      y[which(select & x$lfc>0)]=1
      y[which(select & x$lfc<0)]=1
      print(s)
      print(table(y))
      y
    })
    row.names(sig)= select.genes
    sc = sapply(names(de.df), function(s) {
      x = de.df[[s]][select.genes,]
      abs(-log10(x$padj) * sig[, s])
    })    
    de.genes = sapply(names(de.df), function(s){
      x= de.df[[s]][select.genes,]
      x = x[sig[,s]!=0,]
      if(nrow(x)==0){
        return(NULL)
      }
      x = x[order(rowSums(abs(sig[row.names(x),,drop=F])),rowSums(sc[row.names(x),,drop=F]),decreasing=T),]
      c(head(row.names(x)[x$lfc>0],n),head(row.names(x)[x$lfc<0],n))
    },simplify=F)
    markers = sort(unique(unlist(de.genes)))
    ###find markers for each cell type
    cl.markers=sapply(unique(cl), function(i){
      print(i)
      pairs = sort(c(grep(paste("^", i, "_",sep=""), colnames(sig)),grep(paste("_", i, "$",sep=""), colnames(sig))))
      select = c()
      while(length(pairs) > 0){
        tmp=order(rowSums(abs(sig[,pairs,drop=F])),rowSums(sc[,pairs,drop=F]),decreasing=T)
        tmp = tmp[rowSums(sig[tmp, pairs,drop=F]!=0) > 0]
        if(length(tmp) > 0){
          select = c(select, head(row.names(sig)[tmp],n))
          pairs = pairs[colSums(abs(sig[select,pairs,drop=F]))==0]
        }
        else{
          select = c(select, paste("Fail:", colnames(sig)[pairs],collapse=" "))
          break
        }
      }
      select
    },simplify=F)
    return(list(markers, cl.markers))
  }



test.cv <- function(norm.dat, cl, n.bin=5,padj.th = 0.01, n.markers=10){
  require(randomForest)
  bins = sample(1:n.bin, length(cl),replace=T)
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  rf.pred.cl = setNames(rep(NA, length(cl)), names(cl))
  rf.pred.prob = matrix(0, nrow=length(cl), ncol=length(unique(cl)))
  row.names(rf.pred.prob) = names(cl)

  for(i in 1:n.bin){
    print(i)
    select.cells = names(cl)[bins!=i]
    de.df = DE.genes.pw(norm.dat[,select.cells], paste("cl",cl[select.cells],sep=""))
    de.genes = sapply(de.df, function(x){
      x = x[order(x$padj),]
      up = x$padj < padj.th & x$lfc > 1
      down = x$padj < padj.th & x$lfc < -1
      c(head(row.names(x)[up],n.markers), head(row.names(x)[down],n.markers))},simplify=F)
    markers = unique(unlist(de.genes))
    rf.result<-randomForest(t(norm.dat[markers, select.cells]),as.factor(cl[select.cells]),ntree=2000)
    tmp = predict(rf.result, t(norm.dat[markers, names(cl)[bins==i]]),type="prob")
    rf.pred.prob[bins==i,] = tmp
    rf.pred.cl[bins ==i]= apply(tmp, 1, which.max)
  }
  list(class=rf.pred.cl, prob = rf.pred.prob)
}
    
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param counts 
##' @param ercc.counts 
##' @param select.cells 
##' @param prefix 
##' @param remove.gene.mod 
##' @param pos.gene.list 
##' @param col 
##' @param DESeq.padj.th 
##' @param min.cells 
##' @param softPower 
##' @param minModuleSize 
##' @param cutHeight 
##' @param mod.min.padj 
##' @param cl.min.padj 
##' @param de.padj.th 
##' @param kME.th 
##' @return 
##' @author 
WGCNA.cluster <- function(counts, ercc.counts, select.cells, prefix, norm.dat=log2(counts[,select.cells]+1), rm.gene.mod=NULL, pos.gene.list=NULL, col=NULL, DESeq.padj.th=0.05,min.cells=4, softPower=4,minModuleSize=6,cutHeight=0.995,mod.min.padj=80,cl.min.padj=40,de.padj.th=0.05, lfc.th=1, kME.th=0.25,n.markers=30,max.mod=NULL,gene.mod=NULL,dissTOM=NULL, geneTree=NULL,maxGenes=6000,rm.th=0.8,sampleSize=1500)
  {
    if(is.null(gene.mod)){
      if(length(select.cells)>sampleSize){
        tmp.cells = sample(select.cells, pmin(length(select.cells),sampleSize))
      }
      else{
        tmp.cells = select.cells
      }
      if(is.null(dissTOM)| is.null(geneTree)){
        print("Get Gene Tree")
        tmp = getGeneTree(norm.dat, ercc.counts, counts, select.cells=tmp.cells,min.cells=min.cells, padj.th=DESeq.padj.th, prefix=prefix, softPower=4,select.genes=NULL,type="unsigned",maxGenes=maxGenes)
        if(is.null(tmp)){return(NULL)}
        geneTree = tmp[[1]]
        dissTOM= tmp[[2]]
      }
      print("Get Gene Module")
      tmp=cutGeneTree(norm.dat, tmp.cells, geneTree, dissTOM, minModuleSize=6, min.cells=min.cells, cutHeight=cutHeight,min.padj=mod.min.padj, prefix=prefix, ColSideColors=col[,tmp.cells],padj.th=de.padj.th,lfc.th=lfc.th, max.mod=max.mod)
      gene.mod= tmp[[1]]
      gene.mod.val=tmp[[2]]
      #save(gene.mod, gene.mod.val, file=paste0(prefix, ".gene.mod.rda"))
      if(is.null(gene.mod)){return(NULL)}
      if(!is.null(pos.gene.list)){
        tmp=sapply(gene.mod, function(x)length(intersect(x, pos.gene.list)))
        gene.mod = gene.mod[tmp > 0] 
      }
      if(is.null(gene.mod) | length(gene.mod)==0){return(NULL)}
    }
    names(gene.mod)=1:length(gene.mod)    
    print("Cluster cells")
    eigen = getEigen(gene.mod, norm.dat,select.cells,prefix=prefix,cexRow=1)[[1]]
    if(!is.null(rm.gene.mod)){
      rm.eigen = getEigen(rm.gene.mod, norm.dat,select.cells)[[1]]    
      rm.cor=cor(eigen, rm.eigen)
      select1 = rowMaxs(rm.cor) < rm.th
      select2 = apply(sapply(gene.mod, function(x){
        sapply(rm.gene.mod, function(y){
          m=length(intersect(x,y))
          max(m/length(x),m/length(y)) < 0.5
        })
      }), 2, all)     
      select = select1 & select2
      eigen= eigen[,select,drop=F]
      gene.mod = gene.mod[select]
    }
    else{
      rm.eigen=NULL
    }
    if(is.null(gene.mod) | length(gene.mod)==0){return(NULL)}
    if(length(gene.mod) <=3){
      markers=unlist(gene.mod)
      hc = hclust(dist(t(norm.dat[markers,select.cells])),method="ward")
    }
    else{    
      hc = hclust(dist(eigen),method="ward")
      pdf(paste0(prefix, ".eigen.pdf"),height=6,width=7)
      heatmap.3(t(eigen), Colv= as.dendrogram(hc), col=blue.red(100), trace="none", dendrogram="column",ColSideColors=col[,select.cells])
      dev.off()
    }
    cl = cutree(hc, length(gene.mod)+4)
    merge.result = mergeCl(norm.dat, hc, cl = cl, min.padj=cl.min.padj,padj.th=de.padj.th,lfc.th=lfc.th, eigen=eigen, rm.eigen=rm.eigen,min.cells=min.cells)
    if(is.null(merge.result))return(NULL)
    #save(merge.result, file=paste0(prefix, ".merge.rda"))
    sc = merge.result$sc    
    cl = merge.result$cl
    cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
    tmp.col =rbind(col[,select.cells],cl.col)
    de.df = merge.result$de.df
    if(ncol(eigen)>1){
      eigen.de = DE.genes.pw(t(eigen),paste0("cl",cl))
      select = rowMins(sapply(eigen.de, function(x)x$padj)) < de.padj.th
      eigen = eigen[,select,drop=F]
      gene.mod= gene.mod[select]
    }
    print(sort(merge.result$sc))
    if(length(gene.mod) >= 2){
      z = matrix(0, nrow=length(select.cells),ncol=length(unique(cl)))
      row.names(z)=select.cells
      colnames(z)= sort(unique(cl))
      for(i in unique(cl)){
        z[select.cells[cl==i],i] = 1
      }
      cl.me=me("VVI", eigen,z)
      if(attr(cl.me, "returnCode") >= 0){
        row.names(cl.me$z) = row.names(z)
        colnames(cl.me$z) = colnames(z)
        z = cl.me$z
        new.cl = setNames(colnames(z)[apply(z, 1, which.max)], select.cells)
        merge.result = mergeCl(norm.dat, cl = new.cl, min.padj=cl.min.padj,padj.th=de.padj.th,lfc.th=lfc.th, eigen=eigen, rm.eigen=rm.eigen,min.cells=min.cells)
        new.sc = merge.result$sc
        cat("Score ", sum(sc), "mclust", sum(new.sc), "\n")
        if(sum(new.sc) > sum(sc)){
          cl=merge.result$cl
          de.df = merge.result$de.df
        }
      }
    }
    markers = displayCl(norm.dat, cl, de.df=de.df, n.markers=n.markers,eigen=eigen, rm.eigen=rm.eigen, kME.th=kME.th, min.cells=min.cells, padj.th= de.padj.th,lfc.th=lfc.th,ColSideColors=col[,select.cells])
    cl = setNames(factor(cl), names(cl))
    levels(cl) = 1:length(levels(cl))
    if(length(gene.mod)>=2){
      cl.means = do.call("cbind",tapply(names(cl), cl, function(x){
        rowMeans(norm.dat[markers,x])
      }))
      cl.hc = hclust(dist(t(cl.means)),method="single")
      cl = setNames(factor(as.integer(cl), levels= colnames(cl.means)[cl.hc$order]), names(cl))
    }
    cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
    tmp.col =rbind(col[,select.cells],cl.col)
    markers = displayCl(norm.dat, cl, de.df=de.df, ColSideColors=tmp.col, prefix=prefix, markers=markers, by.cl=TRUE)
    result=list(cl=cl, markers=markers,gene.mod= gene.mod, eigen=eigen)
    return(result)
  }


recursive_cluster<-function(split.size=10,counts, select.cells,prefix,norm.dat=log2(counts[,select.cells]+1), col,result=NULL,display=TRUE,...)
  {
    print(prefix)
    if(is.null(result)){
      result=WGCNA.cluster(counts=counts, select.cells=select.cells, prefix=prefix,norm.dat=norm.dat, col=col,...)
      gc()
      save(result, file=paste0(prefix,".rda"))
    }
    if(!is.null(result)){
      cl = result$cl
      gene.mod = result$gene.mod
      markers=result$markers
      cl = setNames(as.integer(cl),names(cl))
      cl.size = table(cl)
      print(cl.size)
      to.split = names(cl.size)[cl.size >=split.size]
      if(length(to.split)>0){
        n.cl = max(as.integer(cl))
        for(x in to.split){
          tmp.cells = names(cl)[cl==x]
          tmp.prefix = paste(prefix, x, sep=".")
          tmp.result=recursive_cluster(split.size=split.size, counts=counts, select.cells=tmp.cells, prefix=tmp.prefix,col=col,norm.dat=norm.dat,...)
          if(!is.null(tmp.result)){            
            tmp.cl = tmp.result$cl
            if(length(unique(tmp.cl)>1)){
              cat("Expand",tmp.prefix, "\n")
              print(table(tmp.cl))
              cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
              gene.mod = c(gene.mod, tmp.result$gene.mod)
              markers = union(markers,tmp.result$markers)
            }
          }
          n.cl = max(as.integer(cl))
        }
        cl.means = do.call("cbind",tapply(names(cl), cl, function(x){
          rowMeans(norm.dat[markers,x])
        }))
        cl.hc = hclust(dist(t(cl.means)))
        cl = setNames(factor(as.integer(cl), levels= colnames(cl.means)[cl.hc$order]), names(cl))
        levels(cl) = 1:length(levels(cl))
      }
      if(display){
        cl.col = jet.colors(length(unique(cl)))[as.factor(cl)]
        tmp.col =rbind(col[,select.cells],cl.col)      
        tmp.markers = displayCl(norm.dat[,select.cells], cl, ColSideColors=tmp.col, prefix=prefix, markers=markers, de.df=NULL,by.cl=TRUE)
      }
      return(list(cl=cl, gene.mod=gene.mod, markers=markers))
    }
    return(NULL)
  }


