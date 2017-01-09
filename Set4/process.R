###This is a simplified the script for clustering that produces similar but not identical results as the published clustering result. The published clustering result were produced using eariler version of the script involving more manual intervention, and iteractive refinement for cleaner idenfication of clusters boundaries. 

library(limma)
source("de.genes.R")
source("DESeq.var.R")
source("heatmap.R")
source("iWGNCA.R")

load("sample.df.rda")
load("counts.rda")
load("ercc.counts.rda")
load("tpm.rda")
load("col.rda")
load("rm.gene.mod.rda")
load("all.cells.rda")

norm.dat = log2(tpm+1)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))

cl = c()

for(day in levels(sample.df$day)){
  prefix= paste(day, sep=".")
  select.cells = all.cells[as.character(sample.df[all.cells, "day"])==day]
  print(prefix)
  print(length(select.cells))
  min.cells=8
  if(day=="D125"){
    min.cells=4
  }
  ###split progenitors and neurons
  result1=WGCNA.cluster(counts=counts, ercc.counts=ercc.counts, select.cells, prefix, norm.dat=norm.dat,rm.gene.mod=rm.gene.mod, col=col[,select.cells],DESeq.padj.th=0.5,de.padj.th=0.01,min.cells=min.cells, mod.min.padj=1000, cl.min.padj=5000,maxGenes=6000, rm.th=0.8)
  ###recursive clustering
  result=recursive_cluster(split.size=20,counts=counts, ercc.counts=ercc.counts, select.cells, prefix, norm.dat=norm.dat,rm.gene.mod=rm.gene.mod, col=col[,select.cells],DESeq.padj.th=0.5,de.padj.th=0.01,min.cells=min.cells, mod.min.padj=40, cl.min.padj=80,maxGenes=6000, rm.th=0.7,result=result1)
  tmp.cl = setNames(as.integer(result$cl), names(result$cl))
  rm.eigen = getEigen(rm.gene.mod, norm.dat,select.cells)[[1]]
  ###Merge similar clusters
  result2=mergeCl(norm.dat[, select.cells], cl=tmp.cl, rm.eigen=rm.eigen, min.cells=min.cells, padj.th=0.01,min.padj=80,q=0.5,rm.th=0.5)
  cl[select.cells] = paste0(day, ".",result2$cl[select.cells])
}
save(cl, file="cl.rda")













  
