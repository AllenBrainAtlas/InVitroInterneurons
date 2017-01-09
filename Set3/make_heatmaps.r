require(gplots);
fulldata<-read.csv("Population_rnaseq_data_tpm.csv",as.is=T,row.names=1,check.names=F)
rowcols<-c(rep("purple",2),rep("grey",2),rep("cyan",2),rep("red",2),rep("navy",2),rep("tan",2),rep("green",4),rep("orange"),rep("blue",4),rep("black",4));
subinds<-grep("D0|D6|D9|D12|D19|D26|D54",colnames(fulldata))
subgenedat<-fulldata[genelist,subinds];
subcolvec<-do.call('rbind',strsplit(colnames(subgenedat),"-"))[,3]
genelist<-c("ACTB","GAPDH","NANOG","POU5F1","OTX2","DACH1","SOX10","FOXD3","GATA1","SOX17","RUNX1","GATA2","PAX6","EMX2","LHX2","TBR1","FOXG1","DLX1","ASCL1","NKX2-1","GAD1","EN2","PAX7","TFAP2B","TH");

####collapse by date###
newdat<-subgenedat[,1:length(unique(subcolvec))];
for (ii in 1:length(unique(subcolvec))) {
	keepinds<-which(subcolvec==unique(subcolvec)[ii]);
	if (length(keepinds)>1) {
		newdat[,ii]<-apply(subgenedat[,keepinds],1,mean);
	} else {
		newdat[,ii]<-subgenedat[,keepinds];
	}
}
colvec2<-unique(subcolvec);
colnames(newdat)<-unique(subcolvec)
rownames(newdat)<-paste(rownames(newdat),apply(newdat,1,mean),sep="_")
pdf("gene_exp_heatmaps.pdf");
heatmap.2(as.matrix(log10(newdat+1)),scale='row',col=bluered,trace='none',Rowv=FALSE,Colv=FALSE,RowSideColors=rowcols);
heatmap.2(as.matrix(log10(newdat+1)),scale='none',trace='none',Rowv=FALSE,Colv=FALSE,RowSideColors=rowcols);
newmat<-as.matrix(log10(newdat+1))
newmat<-sweep(newmat,1,apply(newmat,1,max),'/')
heatmap.2(newmat,scale='none',trace='none',Rowv=FALSE,Colv=FALSE,RowSideColors=rowcols);
dev.off();

