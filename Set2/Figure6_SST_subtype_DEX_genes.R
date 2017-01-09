# Load functions
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)  # https://bioconductor.org/packages/release/bioc/html/limma.html
require(pheatmap)  # https://cran.r-project.org/web/packages/pheatmap/index.html
require(RColorBrewer)  # https://cran.r-project.org/web/packages/RColorBrewer/index.html

DexVsAll <-  function(f.levels) {
  n <- length(f.levels)
  design <- matrix(-1 / (n - 1), n, n)
  rownames(design) <- f.levels
  colnames(design) <- f.levels
  diag(design) <- 1
  return(design)
}

# Load single cell data
load(file = "interneuron_ss.RData", verbose = TRUE)
expr.dat <- as.data.frame(log2(datExprR + 1))

# Select clusters
diff.cl <- c("D54.6", "D54.7", "D54.8", "D54.9",  # MGE (SST+)
             "D100.8", "D100.9", "D100.10", "D125.7", "D125.8",  # MGE (SST+)
             # "D26.5", "D26.6", 
             "D54.11", "D125.6",  # Early MGE (SST-)
             # "D26.7", "D26.8", 
             "D54.10", "D100.12")  # LGE, striatal
keep.cells <- which(cl %in% diff.cl[1:4])
expr.subset <- expr.dat[, keep.cells]
samp.subset <- droplevels(samp.dat[keep.cells, ])
cl.subset <- droplevels(cl[keep.cells])

# DEX genes
design <- model.matrix(~ 0 + cl.subset)
colnames(design) <- gsub("cl.subset", "", colnames(design), fixed = TRUE)

fit <- lmFit(expr.subset, design)
cont.matrix <- DexVsAll(colnames(design))

fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
top1 <- topTable(fit2, number = Inf, p.value = 0.05, adjust = "BH", 
                 sort.by = "F")

dex.df <- sapply(colnames(top1)[1:4], function(x) {  # colnames(design)
  rownames(top1)[order(abs(top1[, x]), decreasing = TRUE)]
})

write.csv(dex.df, file = "dex.csv", row.names = FALSE)

dex.genes <- rownames(top1)[1:50]


# Plot heatmap of gene sets
hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
div.colors <- colorRampPalette(c("#0C8CFF", "white", "#FF3300"))(100)
cl.label.df <- as.data.frame(matrix(NA, length(diff.cl), 1))
rownames(cl.label.df) <- diff.cl
colnames(cl.label.df) <- "group"
cl.label.df$group <- c(rep("MGE_SST", 9), rep("early_MGE", 2), 
                       rep("LGE_striatal", 2))

plot.genes <- dex.genes[dex.genes %in% rownames(top1)]
expr.bycl2 <- expr.bycl[plot.genes, diff.cl[1:4]]
colnames(expr.bycl2) <- c("nCTX.54.SST", "nCTX.54.GRIA4", "nCTX.54.CORT", 
                          "nCTX.54.CORT.NEFL")
pheatmap(expr.bycl2, 
         scale = "row", border = NA, fontsize_row = 8,
         # clustering_method = "ward.D2",
         # annotation_col = cl.label.df,  
         color = div.colors,
         main = "(cluster z-score)")
