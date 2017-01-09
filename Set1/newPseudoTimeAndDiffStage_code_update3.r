########################################################################################
print("Read in and prepare the data set.")

zFolder    = "\\\\aibsdata/0285/Personal_folders/Zizhen/projects/Inh_Shared/2016-01-22/"
cFolder    = "\\\\aibsdata/0285/Personal_folders/Zizhen/projects/Inh_Shared/2016_02_25/"
mainFolder = "\\\\aibsdata/humancelltypes/JeremyM/singleCellAnalyses/interneuronAnalysis/"
moduleCall = "\\\\AIBSFILEPRINT/Public/JennieC/Interneuron paper_6-1-2016/Misc/sc_cell_type/cl.2016-03-07.rda"
clInfoFn   = "\\\\AIBSFILEPRINT/Public/JennieC/Interneuron paper_6-1-2016/Misc/sc_cell_type/cl.df.2016-03-07.csv"

setwd(mainFolder)

source("normalizationScript_additionalFunctions.r")
anac = function(x) as.numeric(as.character(x))
library(WGCNA)
library(gplots)

load(paste(zFolder,"counts.rda",sep=""))
genes = rownames(counts)

load(moduleCall)
cl = all.cl
cellAsgn  = cl
cellNames = names(cl)
cellCols  = cl.col
names(cellCols) = names(cl) 

clInfo    = read.csv(clInfoFn,row.names=1)
v1Marks   = read.csv("adultV1InterneuronMarkers.csv")
geneCalls = read.csv("cellTypeMarkers.csv")

load(paste(zFolder,"sample.df.rda",sep=""))
sampleInfo = sample.df[cellNames,]

load(paste(zFolder,"rpkm.rda",sep=""))
colnames(rpkm) = colnames(counts)

########################################################################################
print("Subset the data to only include genes with reasonable expression in at least a small subset of cells.")

qFrac     = 0.95  
minThresh = 1    # Require RPKM=1 in at least 5% of the in vitro cells
kpGene    = apply(rpkm[,cellNames],1,quantile,qFrac)>=minThresh
datExprIn = log2(rpkm[kpGene,cellNames]+1)


########################################################################################
print("Read in and format subpopulation data.")

subpopFn   = "\\\\aibsdata2/ivhct/RNAseq/Nextera/Inh_subpop/subpopulations.csv"
datSubpop0 = read.csv(subpopFn,row.names=1)
datSubpop0 = datSubpop0[,1:40]
cn = scan(subpopFn,nlines=1,sep=",",what="character")[2:41]
cn = gsub("100PS_","",cn)
cn = gsub("50PS_","",cn)
cn = paste(cn,1:40,sep="")
cn[c(15,17)] = cn[c(17,15)] 
cn[c(12,16)] = cn[c(16,12)] 
colnames(datSubpop0) = cn
kpSamp = !is.element(colnames(datSubpop0),c("-_D54_13","+_D54_33","+_D54_y37"))
MDSplot(datSubpop0[,kpSamp],main="tmp",col=labels2colors(substr(cn[kpSamp],3,6))) # This looks good

datSubpop = datSubpop0[,kpSamp]
# Require RPKM=1 in at least 5% of the populations
kpGeneP   = apply(datSubpop,1,quantile,qFrac)>=minThresh
datSubpop = log2(datSubpop[kpGeneP,]+1)
popStage  = substr(colnames(datSubpop),3,6)
popStage  = gsub("_","",popStage)
popHasDCX = paste("DCX",substr(colnames(datSubpop),1,1),sep="")
MDSplot(datSubpop,main="tmp",col=labels2colors(popStage),log2dat=FALSE)

########################################################################################
print("Read in and format the in vivo data for comparison with in vitro")

datElliot   = read.csv("\\\\Aibsdata/0285/Single_cell_rna_seq/ERT/SmartSeq2/3Brn_100C/Analysis/tpm_rename_with4023.csv",row.names=1)
	# DOUBLE CHECK THAT THIS IS CORRECT DATA SET.  THE NAME CHANGED
rownames(datElliot) = rownames(rpkm)  # Correct the MARCH / SEPT error that excel introduces.
gateElliot  = substr(colnames(datElliot),1,3)
brainElliot = substr(colnames(datElliot),5,8)
agesElliot  = c(96,115,101,113)
regsElliot  = c("anterior/frontal","unknown","anterior/ventral","unknown")
names(agesElliot) <- names(regsElliot) <- brainsElliot <- c("4023","4013","4015","4005")
ageElliot   = agesElliot[brainElliot]
regElliot   = regsElliot[brainElliot]
colEll      = labels2colors(gateElliot)
colEll[colEll=="yellow"] = "orange"
textEll     = substr(gateElliot,3,3)
textEll2    = substr(brainsElliot,3,4)
textEll3    = ifelse(textEll=="9",19,NA)
datCompB    = log2(datElliot[kpGeneP,]+1)
datCompB    = as.matrix(datCompB)
meanElliot  = findFromGroups(t(datElliot),gateElliot)


########################################################################################
print("Generate a bar plot showing some key genes across all of the gating types to prove P9 is interneurons.")

tmpPlot = meanElliot[,order(meanElliot["SST",])]
gn3     = intersect(rownames(meanElliot),c("SOX2","PAX6","VIM","MKI67","NRN1","NRGN","GAD1","SLC17A6",
  "SLC32A1","GJA1","GJB6","LHX6","EOMES","SATB2","TBR1","C1QB","PLLP","OLIG1","SST","CDH5","PDGFRA",
  "HOPX","DCX","AQP4","DLX1","DLX2","GAD2","GABRA6","STARD8","HBA1","SOX6","NEUROD6",
  "SP8","HTR3A"))
zPlot  = apply(tmpPlot[gn3,],1,function(x) return((x-mean(x))/sd(x)))
cn = c("P04","P11","P10","P07","P06","P08","P09")
pdf("primaryCellHeatmap_update.pdf",width=4,height=8)
heatmap.2(t(pmax(pmin(zPlot,2),-1))[,cn],col=blueWhiteRed(20),Colv=NA,trace="none",scale="none")
dev.off()

########################################################################################
print("Perform principal component analysis to cluster populations by stage and differentiation status, then map in in vivo data.")

percents  = c(0.5,1,2,3,5,10,25,50,100,1)  # We are going with 1%, so put this at the end of the list so the variables end up with these values
varCounts = round(percents*dim(datSubpop)[1]/100)
gnPlotCt  = 25
colModsB  = labels2colors(popStage,colorSeq =c("turquoise","brown","orange"))
colModsB2 = labels2colors(popHasDCX)
pchModsB  = 9*(popHasDCX=="DCX+")+8
maxPopFC  = findFromGroups(t(datSubpop),paste(colModsB,colModsB2))
maxPopFC  = apply(maxPopFC,1,max) - apply(maxPopFC,1,min)
maxPopFC  = -sort(-maxPopFC)

pdf("bb_pcaMappingTestPlots_subpop.pdf",height=6,width=11.5)
for (i in 1:length(varCounts)){
 gnUse     = names(maxPopFC)[1:varCounts[i]]
 datUse    = datSubpop[gnUse,]
 pcaUse    = prcomp(t(datUse))
 datUseNew = datCompB[gnUse,]
 datUseNew[is.na(datUseNew)] = 0
 pcaUseNew = scale(t(datUseNew), pcaUse$center, pcaUse$scale) %*% pcaUse$rotation 
 varExp    = round(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10
 mn = paste("Top ",varCounts[i]," FC genes (",percents[i],"% of expressed genes)",sep="")
 px = pcaUse$rotation[,1]
 py = pcaUse$rotation[,2]
 kp = unique(c(names(sort(px))[1:gnPlotCt],names(sort(-px))[1:gnPlotCt],names(sort(py))[1:gnPlotCt],names(sort(-py))[1:gnPlotCt]))
 xl = paste("PC 1: ",varExp[1],"% var explained")
 yl = paste("-PC 2: ",varExp[2],"% var explained")
 
 par(mfrow=c(1,2))
 plot(pcaUse$x[,1],-pcaUse$x[,2],pch=pchModsB,col=colModsB,cex=2.5,xlab=xl,ylab=yl,
     xlim=range(c(pcaUse$x[,1],pcaUseNew[,1]))*1.2,ylim=range(-c(pcaUse$x[,2],pcaUseNew[,2]))*1.2,main=mn)
 points(pcaUse$x[,1],-pcaUse$x[,2],pch=19,cex=0.3)
 points(pcaUseNew[,1],-pcaUseNew[,2],pch=textEll3,col="green",cex=1.5,xlim=range(pcaUseNew[,1])*1.2,ylim=range(pcaUseNew[,2])*1.2)
 if(i==length(varCounts)) abline(v=-10.5,lty="dotted")
 plot(px,py,col="grey",xlab="PC1",ylab="PC2",main="PC loading",pch=19,cex=0.3)
 text(px[kp],py[kp],kp,col="black",cex=0.8)
}
dev.off()
loadingGenes = kp

ke = textEll3!=""
ke[is.na(ke)] = FALSE
xv = c(paste(popHasDCX,popStage),textEll3[ke])
op =  c("DCX- D26","DCX- D54","DCX- D100","N","DCX+ D26","DCX+ D54","DCX+ D100","19") 
x  = apply(cbind(xv,xv),1, function(x,op) return(which(op==x[1])[1]),op)
x[x>7.5] = x[x>7.5]+1
x[x>4.5] = x[x>4.5]+1
x[x>3.5] = x[x>3.5]+1
set.seed(10)
xj = jitter(x,1)
cl = c(colModsB,rep("green",sum(ke)))

pdf("bb_topLoadingGenes.pdf",height=5,width=8)
par(mfrow=c(1,2))
for (gn in loadingGenes){
 y  = c(as.numeric(datSubpop[gn,]),log2(as.numeric(datElliot[gn,ke])+1))
 y  = c(as.numeric(datSubpop[gn,]),log2(as.numeric(datElliot[gn,ke])+1))
 zy = findFromGroupsVector(y,x)
 zx = as.numeric(names(zy))
 plot(xj,y,pch=19,col=cl,xlim=c(0,12),xlab="",ylab="log2(RPKM+1)",main=gn)
 rect(zx-0.5,zy-0.05,zx+0.5,zy+0.05,border=NA,col="black")
 abline(v=6)
 abline(v=c(4,10),lty="dotted")
 plot(xj,2^(y-1),pch=19,col=cl,xlim=c(0,12),xlab="",ylab="RPKM",main=gn)
 rect(zx-0.5,2^(zy-1)-0.05,zx+0.5,2^(zy-1)+0.05,border="black",col="black")
 abline(v=6)
 abline(v=c(4,10),lty="dotted")
}
dev.off()



########################################################################################
print("Perform principal component analysis on single cells to cluster samples by stage and differentiation status.")

percents  = c(0.5,1,2,3,5,10,100,1)  # We are going with 1%, so put this at the end of the list so the variables end up with these values
varCounts = round(percents*dim(datExprIn)[1]/100)
gnPlotCt  = 50
colMods   = labels2colors(sampleInfo$day)
colMods[colMods=="yellow"] = "orange"
colMods2  = labels2colors(sampleInfo$line)
pchMods   = 9*(sampleInfo$line=="Dp")+8
cexMods   = 0*(sampleInfo$line!="Dp")+1
maxCelFC  = findFromGroups(t(datExprIn),paste(colMods,colMods2))
maxCelFC  = apply(maxCelFC,1,max) - apply(maxCelFC,1,min)
maxCelFC  = -sort(-maxCelFC)

pdf("bb_pcaCellTypePlots.pdf",height=8.5,width=8)
for (i in 1:length(varCounts)){
 gnUse  = names(maxCelFC)[1:varCounts[i]]
 datUse    = datExprIn[gnUse,]
 pcaUse    = prcomp(t(datUse))
 varExp    = round(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10
 mn = paste("Top ",varCounts[i]," FC genes (",percents[i],"% of expressed genes)",sep="")
 px = pcaUse$rotation[,1]
 py = pcaUse$rotation[,2]
 kp = unique(c(names(sort(px))[1:gnPlotCt],names(sort(-px))[1:gnPlotCt],names(sort(py))[1:gnPlotCt],names(sort(-py))[1:gnPlotCt]))
 xl = paste("PC 1: ",varExp[1],"% var explained")
 yl = paste("PC 2: ",varExp[2],"% var explained")
 
 plot(pcaUse$x[,1],pcaUse$x[,2],pch=pchMods,col=colMods,cex=cexMods,xlab=xl,ylab=yl,
     xlim=range(c(pcaUse$x[,1],pcaUseNew[,1])),ylim=range(c(pcaUse$x[,2],pcaUseNew[,2])),main=mn)
 points(pcaUse$x[,1],pcaUse$x[,2],pch=19,col="black",cex=0.1)
}
dev.off()

diffStage  = pcaUse$x[,1];
pseudoTime = pcaUse$x[,2];
ldct       = 4

pdf("bb_pcaLoadingGenes_eachAxis.pdf",height=3.2,width=10)
linearPlotValues(datExprIn[c(names(sort(py))[1:ldct],names(sort(-py))[1:ldct]),],pseudoTime,main="Top genes loading by pseudotime")
linearPlotValues(datExprIn[c(names(sort(px))[1:ldct],names(sort(-px))[1:ldct]),],diffStage,main="Top genes loading by diff stage")
linearPlotValues(datExprIn[c("PTN","BCAN","C1orf61","FAM60A","CRABP2","DLK1"),],pseudoTime,main="Select genes loading by pseudotime")
linearPlotValues(datExprIn[c("GAD1","STMN2","SST","S1PR1","PDPN","VIM"),],diffStage,main="Select genes loading by diff stage")
dev.off()

testGenes = read.csv("testGenesForPlotting.csv")
testGenes = testGenes[is.element(testGenes[,1],rownames(datExprIn)),]

pdf("bb_testGenes_eachAxis.pdf",height=5,width=10)
for (r in unique(testGenes$Reason)){
 gnTmp = as.character(testGenes$Gene)[testGenes$Reason==r]
 tmpI  = which(testGenes$Reason==r)[1]
 if(testGenes$variable[tmpI]=="diffStage") { value   = diffStage; kpValue = pseudoTime; }
 if(testGenes$variable[tmpI]!="diffStage") { kpValue = diffStage; value   = pseudoTime; }
 kpSam = (kpValue>testGenes$minVal[tmpI])&(kpValue<testGenes$maxVal[tmpI])
 linearPlotValues(datExprIn[gnTmp,kpSam],value[kpSam],main=r,ylim=c(0,15))
}
dev.off()


########################################################################################
print("Determine the significance of correlations between each gene and diffStage and pseudoTime.")

corPval <- function(x,...){
  cr = cor.test(x, ...)
  return(cr$p.value)
}

dsCor  = apply(datExprIn,1,cor,diffStage)
ptCor  = apply(datExprIn,1,cor,pseudoTime)
dsPval = apply(datExprIn,1,corPval,diffStage)
ptPval = apply(datExprIn,1,corPval,pseudoTime)
out    = data.frame(Gene_Symbol = names(dsCor), Diff_Status_Correlation = dsCor, Diff_Status_Pvalue = dsPval,
         Pseudotime_Correlation = ptCor, Pseudotime_Pvalue = ptPval) 
write.csv(out,"DSandPT_correlations.csv",row.names=FALSE)

# dsPermMax <- dsPermMin <- ptPermMax <- ptPermMin <- NULL
# for (i in 1:9){
 # set.seed(i)
 # dsTmp = apply(datExprIn,1,cor,sample(diffStage,length(diffStage)))
 # ptTmp = apply(datExprIn,1,cor,sample(pseudoTime,length(pseudoTime))) 
 # print(i)
# }


########################################################################################
print("Compare the maximum FC for the population and single cell data.")

percent     = 1
maxPopFC2   = maxPopFC
maxCelFC2   = maxCelFC
varCountPop = round(percent*dim(datSubpop)[1]/100)
varCountCel = round(percent*dim(datExprIn)[1]/100)
valCountPop = maxPopFC2[varCountPop]
valCountCel = maxCelFC2[varCountCel]
maxPopFC2   = maxPopFC2[intersect(names(maxPopFC2),names(maxCelFC2))]
maxCelFC2   = maxCelFC2[intersect(names(maxPopFC2),names(maxCelFC2))]
popAndCel   = intersect(names(maxPopFC[1:varCountPop]),names(maxCelFC[1:varCountCel]))
print(length(popAndCel))
colPC       = rep("black",length(maxCelFC2))
colPC[is.element(names(maxCelFC2),popAndCel)] = "green"
pdf("PopulationVsSingleCell_ComparisonPlotFC.pdf")
verboseScatterplot(maxPopFC2,maxCelFC2,xlab="Max log2(FC) in population",
  ylab="Max log2(FC) in single cells",pch=19,col=colPC,cex=0.5)
abline(v=valCountPop,lty="dashed",col="grey")
abline(h=valCountCel,lty="dashed",col="grey")
verboseScatterplot(rank(maxPopFC2),rank(maxCelFC2),xlab="Rank Max log2(FC) in population",
  ylab="Rank Max log2(FC) in single cells",pch=19,col=colPC,cex=0.5)
dev.off()


########################################################################################
print("Plot all genes expressed in at least 1% of cells in MDS space")

plotDir   ="\\\\aibsdata/hct/HCT_RNAseq/inVitroInterneuronGenes/"
dir.create(plotDir)
cexMods2  = 0.6*(sampleInfo$line!="Dp")+1.3
datExpr   = log2(rpkm[,cellNames]+1)
onCount   = rowSums(datExpr>0)
plotGenes = sort(names(onCount)[onCount>=(dim(datExpr)[2]*0.01)])
plotGenes = plotGenes[!is.element(1:length(plotGenes),grep("/",plotGenes))]
for(gn in plotGenes){
 fn = paste(plotDir,gn,".pdf",sep="")
 pdf(fn,height=8.5,width=8)
 x    = datExpr[gn,]
 isG0 = x>0
 gnCol= numbers2colors(x[isG0])
 mn   = paste(gn,"- log2(FPKM+1) range =",signif(min(x[isG0]),2),"to",signif(max(x),2))
 plot(diffStage[!isG0],pseudoTime[!isG0],xlab="PC 1 (differentiation stage)",ylab="PC 2 (pseudotime)",
   main=mn,col="grey",pch=19,cex=0.3,xlim=range(diffStage),ylim=range(pseudoTime))
 points(diffStage[isG0],pseudoTime[isG0],pch=pchMods[isG0],col=gnCol,cex=cexMods2[isG0])
 dev.off()
}


########################################################################################
print("Redefine progenitors and differentiation stage based on the PCA, and then find top correlates to both metrics.")

dVal    = 4
pVal    = -16
aVal    = -5
newLine = rep("intermediate",length=length(diffStage))
newLine[diffStage<pVal] = "progenitor"
newLine[diffStage>dVal] = "differentiated"
newStage = rep("early",length=length(diffStage))
newStage[pseudoTime>aVal] = "late"
lines  = c("progenitor","intermediate","differentiated")
stages = c("early","late")

pdf("bb_pcaCellTypePlot_Final.pdf",height=8.5,width=8)
for (i in length(varCounts)){
 gnUse     = names(maxCelFC)[1:varCounts[i]]
 datUse    = datExprIn[gnUse,]
 pcaUse    = prcomp(t(datUse))
 varExp    = round(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10
 mn = paste("Top ",varCounts[i]," FC genes (",percents[i],"% of expressed genes)",sep="")
 px = pcaUse$rotation[,1]
 py = pcaUse$rotation[,2]
 kp = unique(c(names(sort(px))[1:gnPlotCt],names(sort(-px))[1:gnPlotCt],names(sort(py))[1:gnPlotCt],names(sort(-py))[1:gnPlotCt]))
 xl = paste("PC 1: ",varExp[1],"% var explained")
 yl = paste("PC 2: ",varExp[2],"% var explained")
 
 plot(pcaUse$x[,1],pcaUse$x[,2],pch=pchMods,col=colMods,cex=cexMods,xlab=xl,ylab=yl,
     xlim=range(c(pcaUse$x[,1],pcaUseNew[,1])),ylim=range(c(pcaUse$x[,2],pcaUseNew[,2])),main=mn)
 points(pcaUse$x[,1],pcaUse$x[,2],pch=19,col="black",cex=0.1)
 abline(h=aVal,lty="dotted")
 abline(v=c(dVal,pVal),lty="dotted")
}
dev.off()

allCorrs = data.frame(
 pseudoTimeAll = apply(datExprIn,1,cor,pseudoTime),
 pseudoTimeProgenitors = apply(datExprIn[,newLine=="progenitor"],1,cor,pseudoTime[newLine=="progenitor"]),
 pseudoTimeNeurons = apply(datExprIn[,newLine=="differentiated"],1,cor,pseudoTime[newLine=="differentiated"]),
 diffStageAll   = apply(datExprIn,1,cor,diffStage),
 diffStageEarly = apply(datExprIn[,newStage=="early"],1,cor,diffStage[newStage=="early"]),
 diffStageLate  = apply(datExprIn[,newStage=="late"],1,cor,diffStage[newStage=="late"]))

write.csv(allCorrs,"CorrelationsWithPseudotimeAndDiffStage.csv")
write.csv(cbind(pseudoTime,diffStage),"pseudotimeAndDiffStageValues.csv")


########################################################################################
print("Compare pseudotime with expression changes in MGE of prenatal NHP and human.")

## Read in and subset the NHP data
nhpFolder = paste(mainFolder,"NHPdata/",sep="")  # Data downloaded from http://blueprintnhpatlas.org/static/download into this folder
nhpDat    = read.csv(paste(nhpFolder,"expression_matrix.csv",sep=""),row.names=1,header=FALSE)
nhpGene   = read.csv(paste(nhpFolder,"rows_metadata.csv",sep=""),row.names=3)
# Updated gene/probe list is Table S2 from http://www.nature.com/nature/journal/v535/n7612/full/nature18637.html, renamed as BakkenMiller_TableS2.csv
nhpGene2  = read.csv(paste(nhpFolder,"BakkenMiller_TableS2.csv",sep=""),row.names=1)
nhpGene2  = nhpGene2[rownames(nhpGene),]
nhpSample = read.csv(paste(nhpFolder,"columns_metadata.csv",sep=""),row.names=1)
kpGene    = nhpGene2$keep_for_analysis
nhpDat3   = nhpDat[kpGene,]
rownames(nhpDat3) <- nhpGene2[kpGene,"macaque_genesymbol"]
nhpAge3   = factor(as.character(nhpSample$age),levels=unique(nhpSample$age))
isGerm    = is.element(nhpSample$structure_acronym,c("V1sz","V1szi","V1szo","V1vz","V1vzi","V1vzo"))&(substr(nhpAge3,1,1)=="E")#
cpRegs    = c("MGE","CGE","LGE","LGEcx")#c("V1-5","V1-6","V1cp","V1cpi")#,"V1sp","V1-2","V1-2/3","V1-3","V1-4","V1-4A","V1-4B","V14Ca","V14Cb","V1cpo"
isCP      = is.element(nhpSample$structure_acronym,cpRegs)&(substr(nhpAge3,1,1)=="E")

nhpMeanGerm = findFromGroups(t(nhpDat3[,isGerm]),paste(nhpAge3,nhpSample$donor_name,sep="_")[isGerm],function(x) return(mean(x,na.rm=TRUE)))
nhpMeanCP   = findFromGroups(t(nhpDat3[,isCP]),paste(nhpAge3,nhpSample$donor_name,sep="_")[isCP],function(x) return(mean(x,na.rm=TRUE)))
nhpAgeGerm  = as.character(lapply(colnames(nhpMeanGerm),function(x) return(strsplit(x,"_")[[1]][1])))
nhpAgeCP    = as.character(lapply(colnames(nhpMeanCP),function(x) return(strsplit(x,"_")[[1]][1])))


## Read in and subset the BS data
bsFolder = paste(mainFolder,"BrainSpanData/exonArray/",sep="")  # Data downloaded from http://blueprintbsatlas.org/static/download into this folder
bsDat    = read.csv(paste(bsFolder,"expression_matrix.csv",sep=""),row.names=1,header=FALSE)
#bsGene   = read.csv(paste(bsFolder,"rows_metadata.csv",sep=""),row.names=3)
bsGene   = read.csv(paste(bsFolder,"rows_metadata.csv",sep=""))
bsSample = read.csv(paste(bsFolder,"columns_metadata.csv",sep=""),row.names=1)
bsAge3   = as.character(bsSample$age)
isV1     = is.element(bsSample$structure_acronym,c("V1C","Ocx"))&(substr(bsAge3,nchar(bsAge3)-2,nchar(bsAge3))=="pcw") #is.element(1:length(bsAge3),grep("cortex",bsSample$structure_name))
isV1     = isV1&(!is.element(bsAge3,c("35 pcw","37 pcw")))
kpSample = !is.element(bsGene$gene_symbol,names(table(bsGene$gene_symbol))[table(bsGene$gene_symbol)>1])
bsMeanV1 = bsDat[kpSample,isV1]
rownames(bsMeanV1) = bsGene$gene_symbol[kpSample]
colnames(bsMeanV1) = paste(bsAge3,bsSample$donor_name,sep="_")[isV1]
bsAgeV1  = as.character(lapply(colnames(bsMeanV1),function(x) return(strsplit(x,"_")[[1]][1])))

# Function to convert macaque and human ages to human pcd
convertAge <- function(ageIn,speciesIn="human",speciesOut="macaque"){
 # speciesOut can also be "eventScore"
 if(speciesIn=="human")     { I1 = 3.167; I2 = 3.72;  }
 if(speciesIn=="macaque")   { I1 = 3.27;  I2 = 2.413; }
 if(speciesIn=="mouse")     { I1 = 2.145; I2 = 1.894; }
 eventScore = (log(2^ageIn)-I1)/I2
 if(speciesOut=="eventScore") return(eventScore)
 if(speciesOut=="human")    { R1 = 3.167; R2 = 3.72;  }
 if(speciesOut=="macaque")  { R1 = 3.27;  R2 = 2.413; }
 if(speciesOut=="mouse")    { R1 = 2.145; R2 = 1.894; }
 ageOut = log2(exp(eventScore*R2+R1))
 return(ageOut)
}
cnvAge <- function(x){
 if(substr(x[1],1,1)=="E") return(convertAge(as.numeric(substr(x,2,nchar(x))),"macaque","human"))
 return(as.numeric(gsub(" pcw","",x))*7)
}
nhpAgeG = cnvAge(nhpAgeGerm)
nhpAgeC = cnvAge(nhpAgeCP)
bsAgeV  = cnvAge(bsAgeV1)


## Add PCs to plot
n     = 500
pt    = pseudoTimeAll
nhpPT = rep("grey",dim(nhpMeanGerm)[1])
nhpPT[is.element(rownames(nhpMeanGerm),names(sort(pt))[1:n])]  = "Down"
nhpPT[is.element(rownames(nhpMeanGerm),names(sort(-pt))[1:n])] = "Up"
bsPT  = rep("grey",dim(bsMeanV1)[1])
bsPT[is.element(rownames(bsMeanV1),names(sort(pt))[1:n])]  = "Down"
bsPT[is.element(rownames(bsMeanV1),names(sort(-pt))[1:n])] = "Up"
mesG  = t(moduleEigengenes(t(nhpMeanGerm),nhpPT,excludeGrey=TRUE)$eigengene)
mesC  = t(moduleEigengenes(t(nhpMeanCP),nhpPT,excludeGrey=TRUE)$eigengene)
mesH  = t(moduleEigengenes(t(bsMeanV1),bsPT,excludeGrey=TRUE)$eigengene)


## Plot selected genes as line plots
someGenes = c("BCAN","PTN","FAM60A","CRABP2","DLK1")  # "C1orf61",
nmg = nhpMeanGerm[is.element(rownames(nhpMeanGerm),someGenes),]
nmc = nhpMeanCP[is.element(rownames(nhpMeanCP),someGenes),]
bsv = bsMeanV1[is.element(rownames(bsMeanV1),someGenes),]
gt  = c(rownames(nmg),"C1orf61"); nmg = rbind(nmg,0); rownames(nmg) = gt
gt  = c(rownames(nmc),"C1orf61"); nmc = rbind(nmc,0); rownames(nmc) = gt
nmg = nmg[someGenes,]
nmc = nmc[someGenes,]
bsv = bsv[someGenes,]
bsv = log2(bsv+1)
nmg = rbind(nmg,mesG)
nmc = rbind(nmc,mesC)
colnames(mesH) = colnames(bsv)
bsv = rbind(bsv,mesH)
someGenes = rownames(nmc)

library("ggplot2")
source("multiplot.R")
geneplots <- vector("list", 0)
j=0
for (i in 1:length(someGenes)){
  j=j+1	 
  modsAll2 = paste("Human V1:",someGenes[i])
  region   = factor(rep("V1",length(bsAgeV1)),levels="V1")
  datAndInfo <- data.frame(age=bsAgeV1, age.numeric=bsAgeV, region=region,subregion=region, 
     module=rep(modsAll2,length(region)),value=as.numeric(bsv[i,]))
  geneplots[[j]] <- qplot(data=datAndInfo, x=age.numeric, y=value, #stat="summary", #fun.y="mean", 
     facets=region~., geom=c("point"), col=subregion, ylab="log2(expression)", 
     shape=subregion, xlab="Age (log2 PSD)", main=modsAll2) + theme_bw() + ylim(range(bsv[i,])) + geom_smooth(se=FALSE)+
     #theme(legend.position="right", text = element_text(size=10)) + 
     xlim(c(50,200))# range(bsAgeV))
  j=j+1
  modsAll2 = paste("NHP VZ/SZ:",someGenes[i])
  region   = factor(rep("VZ/SZ",length(nhpAgeGerm)),levels="VZ/SZ")
  datAndInfo <- data.frame(age=nhpAgeGerm, age.numeric=nhpAgeG, region=region,subregion=region, 
     module=rep(modsAll2,length(region)),value=as.numeric(nmg[i,]))
  geneplots[[j]] <- qplot(data=datAndInfo, x=age.numeric, y=value, #stat="summary", #fun.y="mean", 
     facets=region~., geom=c("point"), col=subregion, ylab="log2(expression)", 
     shape=subregion, xlab="Age (log2 PSD)", main=modsAll2) + theme_bw() + ylim(range(nmg[i,])) + geom_smooth(se=FALSE)+
     #theme(legend.position="right", text = element_text(size=10)) +  
     xlim(c(50,200))#range(nhpAgeG))
  j=j+1
  modsAll2 = paste("NHP GEs:",someGenes[i])
  region   = factor(rep("GEs",length(nhpAgeCP)),levels="GEs")
  datAndInfo <- data.frame(age=factor(nhpAgeC,levels=sort(unique(nhpAgeC))), age.numeric=nhpAgeC, region=region,subregion=region, 
     module=rep(modsAll2,length(region)),value=as.numeric(nmc[i,]))
  geneplots[[j]] <- qplot(data=datAndInfo, x=age.numeric, y=value, #stat="summary", #fun.y="mean", 
     facets=region~., geom=c("point"), col=subregion, ylab="log2(expression)", 
     shape=subregion, xlab="Age (log2 PSD)", main=modsAll2) + theme_bw() + ylim(range(nmc[i,])) + geom_smooth(se=FALSE)+
     #theme(legend.position="right", text = element_text(size=10)) + 
     xlim(c(50,200))#range(nhpAgeC))
}

pdf("NHPandBS_comparisonPlots.pdf",height=20,width=16);   r=length(someGenes);  c=3;
multiplot(plotlist=geneplots,cols=c,layout=matrix(1:(r*c),nrow=r, byrow=TRUE))
dev.off()




