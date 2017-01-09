## This file contains additional functions required for normalization of RNA-Seq data and outlier detection

linearPlotValues <- function(exprMat,val,numBins=40,cdxRange=c(0.5,0.5),firstPlot=FALSE,
         nm=rownames(exprMat),orderGenes=TRUE,main="",ylim=NA){
 if(!is.matrix(exprMat)) exprMat = t(exprMat)
 exprMat = exprMat[,order(val)]
 if(orderGenes){
  ord = apply(exprMat,1,function(x) return(min(which(cumsum(x)/sum(x)>0.5))))
  #ord = apply(exprMat,1,cor,val)
  exprMat = exprMat[order(ord),]
 }
 val   = sort(val)
 len   = dim(exprMat)[1]
 if(is.na(ylim[1])) ylim = c(0,len+1)
 yVals = 1:len
 plot(0,0,xlim=c(-numBins/4,numBins),ylim=ylim,axes=FALSE,xlab="Expression",ylab="Genes",col="white",main=main)
 for (i in 1:len){
  expr = exprMat[i,]
  minV   = min(val)
  rangeV = diff(range(val))
  xN     = 0:numBins
  valN   = round(numBins*(val-minV)/rangeV+0.5)
  valN   = factor(valN,levels=xN)
  cntN   = table(valN)
  exprN  = findFromGroupsVector(expr,valN)
  colN   = numbers2colors(exprN)
  cdxN   = (cntN-min(cntN))/diff(range(cntN))
  cdxN   = (cdxN*diff(cdxRange))+min(cdxRange)
  rect(xN-0.5,yVals[i]-cdxN,xN+0.5,yVals[i]+cdxN,col=colN,border="lightgrey")
  text(-numBins/8,yVals[i],nm[i])
 }
}



anac = function(x) as.numeric(as.character(x))

findFromGroupsVector <- function(x,y,...) { 
 x = as.numeric(x)
 x = cbind(x,x)
 f = findFromGroups(x,y,...)
 return(colMeans(f))
}


mouse2human2 <- function (mouse, m2h){
  # Convert mouse to human symbols
  rownames(m2h) = m2h$Mou
  noHumn = which(!(mouse%in%m2h$Mouse_Symbol))
  humn = which((mouse%in%m2h$Mouse_Symbol))
  mouse[humn] = as.character(m2h[mouse[humn],1])
  mouse[noHumn] = toupper(mouse[noHumn])
  return(mouse)
}

findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}

getAnovaPvalforApply <- function(x,varLabels,varWeights=NULL){
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))
  return(aov.out[[1]]$'Pr(>F)'[1])
}



isFirst <- function(x){
 out = rep(TRUE,length(x));
 for (i in 2:length(x)) if(is.element(x[i],x[1:(i-1)])) out[i]=FALSE
 return(out)
}

MDSplot <- function(dat,main,log2dat=TRUE,col=NA,useText=TRUE,matchCol=NA,
    useCol="black",yellowToOrange=TRUE,distmat=NA, pch=19, cexText=0.6, ...) {
 nm = colnames(dat);
 if(is.na(col)[1]) col = rep(useCol,length(nm))
 if(yellowToOrange) col[col=="yellow"] = "orange"
 if(is.na(distmat[1])) {
  if(log2dat)    dat = log2(dat+1)
  distmat = sqrt(1-cor(dat,use="na.or.complete")^2);   
 }
 loc     = cmdscale(distmat,k=2, eig=TRUE);   x = -loc[[1]][,1];   y = loc[[1]][,2]
 loc1    = cmdscale(distmat,k=1, eig=TRUE);   # For percent variance explained calculation
 if(!is.na(matchCol[1])){
  names(matchCol) = col
  matchCol = matchCol[isFirst(matchCol)]
  plot(rep(0,length(matchCol)+2),0:(length(matchCol)+1),col="white",main=main,axes=FALSE,xlab="",ylab="",...)
  text(rep(0,length(matchCol)),1:length(matchCol),matchCol,col=names(matchCol),cex=cexText)
 }
 if(useText) {
  plot(x, y, type="n", main=main,xlim=range(x)*1.2,ylim=range(y)*1.2,
   ylab = paste("% Var explained, Dim. 2:  ",100*signif(loc$GOF[1]-loc1$GOF[1],3)),
   xlab = paste("% Var explained, Dim. 1:  ",100*signif(loc1$GOF[1],3)),...)
  text(x, y, nm, col = col,cex=cexText)
 } else {
  plot(x, y, pch=pch, main=main,xlim=range(x)*1.2,ylim=range(y)*1.2,
   ylab = paste("% Var explained, Dim. 2:  ",100*signif(loc$GOF[1]-loc1$GOF[1],3)),
   xlab = paste("% Var explained, Dim. 1:  ",100*signif(loc1$GOF[1],3)),col = col,...)
 }
}

HCplot <- function(dat,main,plotText=NULL,isNumeric=NULL,log2dat=TRUE) {
 if(log2dat)    dat = log2(dat+1)
 corFits   = cor(dat)
 clustFits = hclust(as.dist(1-(1+corFits)/2))
 colFits   = data.frame(t(t(rep("white",dim(corFits)[1]))))
 if (length(plotText)>0){
  for (i in 1:length(isNumeric)) { if(isNumeric[i]){
    colFits[,i] = numbers2colors(as.numeric(as.character(plotText[,i])))
   } else { colFits[,i] = labels2colors(as.character(plotText[,i])); }
  }
  colnames(colFits) = colnames(plotText)
 }
 plotDendroAndColors(clustFits,colFits,main=main)
}

IACplot <- function(dat,main,thirdVar,ROI,log2dat=TRUE,thirdVarName="Batch",ylim=NA) {
 nm = colnames(dat);  rois = sort(unique(ROI))
 if(log2dat)    dat = log2(dat+1)
 IACall <- IACreg <- IACrt <- rep(0,length(nm))
 names(IACall) <- names(IACreg) <- names(IACrt) <- nm
 numRois = 1:length(rois);            names(numRois) = rois
 IAC     = cor(dat);                  IACall  = rowMeans(IAC)

 for (i in 1:length(IACall)){
  kpReg = (ROI==ROI[i])&((1:length(IACall))!=i)
  if(sum(kpReg)==0) kpReg = (1:length(IACall))==i
  kpTum = kpReg&(thirdVar==thirdVar[i])
  if(sum(kpTum)==0) kpTum = (1:length(IACall))==i
  IACreg[i] = mean(as.numeric(IAC[i,kpReg]))
  IACrt[i]  = mean(as.numeric(IAC[i,kpTum]))
 }
 
 l = length(rois)+1
 if(is.na(ylim[1])) ylim = c(min(c(IACall,IACreg,IACrt)),1)
 plot(numRois[ROI],IACall,col="white",ylab="IAC",main="Overall IAC",xlab="",xlim=c(0,l),ylim=ylim)
 text(numRois[ROI],IACall,nm,cex=0.6)
 plot(numRois[ROI],IACreg,col="white",ylab="IAC",main=paste(main,"- In Region IAC"),xlab="",xlim=c(0,l),ylim=ylim)
 text(numRois[ROI],IACreg,nm,cex=0.6); abline(h=1,col="red",lwd=2)
 main3 = paste("In Reg./",thirdVarName," IAC",sep="")
 plot(numRois[ROI],IACrt,col="white",ylab="IAC",xlim=c(0,l),ylim=ylim,xlab="",main=main3)
 text(numRois[ROI],IACrt,nm,cex=0.6); abline(h=1,col="red",lwd=2)
}




######################################################################################################################

normalizeByTbT <- function(datInput, dexVector){
  # normalizeByTbT - This function performs TbT normalization on RNA-Seq data.  It is a version of the R scripts 
  #  generated by Robinson and Oshlack (Genome Biol., 11:R25, 2010).  Robinson and Oshlack get full credit for 
  #  the TbT normalization algorithm used here.  Note that it requires several R libraries to be installed.


  ## Load libraries and get appropriate input data
  # library(baySeq);  library(edgeR);  library(DESeq);  library(NBPSeq);  library(ROC)  # Libraries loaded in main script
  # source("http://bioinf.wehi.edu.au/folders/tmm_rnaseq/functions.R")
    # R scripts provided by Robinson and Oshlack, Genome Biol., 11:R25, 2010, and copied below

  data    <- datInput
  data.cl <- dexVector
  RPM     <- sweep(data, 2, 1000000/colSums(data), "*") # Convert to RPM
  groups  <- list(NDE=rep(1, length(data.cl)), DE=data.cl)

  ##  Generation of initial TMM-normalized data
  d          <- DGEList(counts=data, group=data.cl)
  d          <- calcNormFactors(d)
  norm_f_TMM <- d$samples$norm.factors
  RPM_TMM    <- sweep(RPM, 2, 1/norm_f_TMM, "*")

  ##  TbT normalization
  data <- round(RPM)
  hoge <- new("countData", data=as.matrix(data), replicates=data.cl, libsizes=colSums(data)*norm_f_TMM, groups=groups)
  hoge.NB <- getPriors.NB(hoge, samplesize=2000, estimation="QL", cl=NULL)   # THIS IS A SLOW STEP!
  out <- getLikelihoods.NB(hoge.NB, pET="BIC", cl=NULL)                      # THIS IS A VERY SLOW STEP!
  PDEG <- out@estProps[2]
  rank_bayseq <- rank(-out@posteriors[,2])
  DEGnum <- (nrow(data) * PDEG)

  data <- RPM_TMM                                         # For calculating PA value
  meanA <- log2(apply(data[,data.cl==1], 1, mean))        # For calculating PA value
  meanB <- log2(apply(data[,data.cl==2], 1, mean))        # For calculating PA value
  y_axis <- meanB - meanA                                 # For calculating PA value
  PA <- sum(y_axis[rank_bayseq < DEGnum] < 0)/DEGnum      # For calculating PA value

  obj_DEGn <- (rank_bayseq >= DEGnum)
  data <- datInput[obj_DEGn,]
  d <- DGEList(counts=data, group=data.cl)
  d <- calcNormFactors(d)
  norm_f_TbT_RAW <- 1000000/(colSums(data)*d$samples$norm.factors)
  norm_f_TbT <- d$samples$norm.factors*colSums(data)/colSums(datInput)
  RPM_TbT <- sweep(RPM, 2, 1/norm_f_TbT, "*")
  return(RPM_TbT)
}

## The scripts below are additional functions required for TbT normlization and are copied 
##   directly from this web link: http://bioinf.wehi.edu.au/folders/tmm_rnaseq/functions.R


plotReverseCumDist <- function(x, xlab="Tag count X", ylab="# tags >= X", add=FALSE, ...) {
  v <- ecdf(x)
  matplot( knots(v), (1-v(knots(v)))*sum(D[,1]), log="xy", xlab=xlab, ylab=ylab, add=add, ... )
}

generateDataset2 <- function(commonTags=15000, uniqueTags=c(1000,3000), group=c(1,2), libLimits=c(.9,1.1)*1e6, 
                            empiricalDist=NULL, lengthDist=NULL, pDifferential=.05, pUp=.5, foldDifference=2, nreps=c(2,2)) {
                            
  # some checks
  stopifnot( length(group) == length(uniqueTags) )
  stopifnot( length(group) == length(nreps) )
  stopifnot( length(empiricalDist) == length(lengthDist) )
  group <- as.factor(rep(group,nreps))
  stopifnot( nlevels(group) == 2 ) 
  
  print(group)

  #exampleCounts <- empiricalDist/lengthDist
  exampleCounts <- empiricalDist
  exampleLambda <- exampleCounts/sum(exampleCounts)
  exampleIds <- seq_len( length(empiricalDist) )
 
  # set up libraries
  nLibraries <- sum( nreps )
  libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )

  # vector of starts/stops for the unique Tags
  en <- commonTags + cumsum(uniqueTags)
  st <- c(commonTags+1,en[-nLibraries]+1)

  # create matrix of LAMBDA(=relative expression levels)
  LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
  
  ID <- rep(0, max(en))
  ID[1:commonTags] <- sample(exampleIds, commonTags, replace=TRUE)
  LAMBDA[1:commonTags,] <- exampleLambda[ ID[1:commonTags] ]

  # set unique tag totals
  for(i in 1:length(nreps))
    if(uniqueTags[i] > 0) {
      ID[st[i]:en[i]] <- sample(exampleIds, uniqueTags[i], replace=TRUE)
      LAMBDA[st[i]:en[i],group==levels(group)[i]] <- exampleLambda[ ID[st[i]:en[i]] ]
    }
    
  g <- group == levels(group)[1]
  ind <- seq_len(floor(pDifferential*commonTags))
  if(length(ind)>0) {
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    LAMBDA[ind,g] <- LAMBDA[ind,g]*exp(log(foldDifference)/2*fcDir)
    LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
  }

  sampFactors <- colSums(LAMBDA)

  sampFactorsM <- outer( rep(1,max(en)), sampFactors )
  libSizesM <- outer(  rep(1,max(en)), libSizes )

  # create observed means
  MEAN <- LAMBDA / sampFactorsM * libSizesM  # to get the totals to sum to 1

  # sample observed data (column sums will be *close* to set library sizes)
  DATA <- matrix(0, nr=nrow(LAMBDA), ncol=nLibraries)
  DATA <- matrix(rpois(length(MEAN), lambda=MEAN),ncol=nLibraries)

  trueFactors <- colSums(MEAN[1:commonTags,])
  trueFactors <- trueFactors/trueFactors[1]
  
  colnames(DATA) <- paste(paste("group",group,sep=""),1:ncol(DATA),sep=".")
  
  list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, trueFactors=trueFactors, group=group, libSizes=libSizes,  
       differentialInd=c(ind,(commonTags+1):nrow(DATA)), commonInd=1:commonTags, ID=ID, length=lengthDist[ID])
}

takeSubset <- function(obj, subsetInd) {
  allInd <- 1:nrow(obj$DATA)
  commonInd <- allInd %in% obj$commonInd
  differentialInd <- allInd %in% obj$differentialInd
  list(DATA=obj$DATA[subsetInd,], LAMBDA=obj$LAMBDA[subsetInd,], trueFactors=obj$trueFactors, group=obj$group, 
       libSizes=obj$libSizes, differentialInd=which(differentialInd[subsetInd]), commonInd=which(commonInd[subsetInd]),
	   ID=obj$ID[subsetInd], length=obj$length[subsetInd])
}

generateDataset <- function(commonTags=15000, uniqueTags=c(1000,3000), group=c(1,2), libLimits=c(.9,1.1)*1e6, 
                            empiricalDist=NULL, randomRate=1/100, pDifferential=.05, pUp=.5, foldDifference=2) {
                            
  # some checks
  group <- as.factor(group)
  stopifnot( length(group) == length(uniqueTags) )
  #stopifnot( length(group) == 2 ) # code below only works for 2 samples
  stopifnot( nlevels(group) == 2 ) 

  # define where to take random sample from (empirical distribution OR random exponential)
  if(is.null(empiricalDist))
    exampleCounts <- ceiling(rexp(commonTags,rate=randomRate))
  else
    exampleCounts <- empiricalDist
  
  exampleLambda <- exampleCounts/sum(exampleCounts)
 
  # set up libraries
  nLibraries <- length(uniqueTags)
  libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )

  # vector of starts/stops for the unique Tags
  en <- commonTags + cumsum(uniqueTags)
  st <- c(commonTags+1,en[-nLibraries]+1)

  # create matrix of LAMBDA(=relative expression levels)
  LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
  LAMBDA[1:commonTags,] <- sample(exampleLambda, commonTags, replace=TRUE)

  # set unique tag totals
  for(i in 1:nLibraries)
    if(uniqueTags[i] > 0)
      LAMBDA[st[i]:en[i],i] <- sample(exampleLambda, uniqueTags[i])    
    
  ind <- seq_len(floor(pDifferential*commonTags))
  g <- group == levels(group)[1]
  if(length(ind)>0) {
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    LAMBDA[ind,g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*fcDir)
    LAMBDA[ind,!g] <- LAMBDA[ind,!g]*exp(log(foldDifference)/2*(-fcDir))
  }
  
  sampFactors <- colSums(LAMBDA)

  sampFactorsM <- outer( rep(1,max(en)), sampFactors )
  libSizesM <- outer(  rep(1,max(en)), libSizes )

  # create observed means
  MEAN <- LAMBDA / sampFactorsM * libSizesM

  # sample observed data (column sums will be *close* to set library sizes)
  DATA <- matrix(rpois(length(MEAN), lambda=MEAN),ncol=nLibraries)

  trueFactors <- colSums(MEAN[1:commonTags,])
  trueFactors <- trueFactors/trueFactors[1]
  list(DATA=DATA, LAMBDA=LAMBDA, MEAN=MEAN, trueFactors=trueFactors, group=group, libSizes=libSizes,  differentialInd=c(ind,(commonTags+1):nrow(DATA)), 
  commonInd=1:commonTags)
}

calcFactor <- function(obs, ref, trim=.45) {
  logR <- log2(obs/ref)
  fin <- is.finite(logR)
  2^mean(logR[fin],trim=trim)
}

Poisson.model <- function(MA,group1,group2){

  require(limma)
  Poisson.glm.pval <- vector()
  Fold.changes <- vector()
  
  CS <- colSums(MA$M[,c(group1,group2)])

  for (i in 1:(nrow(MA))){
    S1 <- MA$M[i,group1] 
    S2 <- MA$M[i,group2] 
    In <- c(S1,S2)
    sample.f <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
    In <- as.vector(unlist(In))
    GLM.Poisson <- glm(In ~ 1 + sample.f + offset(log(CS)),family=poisson)
    Poisson.glm.pval[i] <- anova(GLM.Poisson,test="Chisq")[5][2,1]
    Fold.changes[i] <- exp(GLM.Poisson$coefficients[1])/(exp(GLM.Poisson$coefficients[1]+GLM.Poisson$coefficients[2]))
  }
  
  output <- matrix(ncol=2,nrow=nrow(MA$M))
  output[,1] <- Poisson.glm.pval
  output[,2] <- Fold.changes
  output <- as.data.frame(output)
  names(output) <- c("pval","FC")
  output
}

Poisson.model.new <- function(countMatrix,group1,group2, ref=1, calcFactor=TRUE){

  Poisson.glm.pval <- vector()
  Fold.changes <- vector()
  
  props <- countMatrix / outer( rep(1,nrow(countMatrix)), colSums(countMatrix) )
  
  refS <- colSums(countMatrix[,c(group1,group2)])
  
  if( calcFactor ) {
    require(edgeR)
	  CS <- calcNormFactors(countMatrix[,c(group1,group2)])
  } else {
    CS <- rep(1,length(group1)+length(group2))
  }
  
  offsets <- log(CS)+log(refS)

  sample.f <- factor(c(rep(1,length(group1)),rep(2,length(group2))))
  
  for (i in 1:(nrow(countMatrix))){
    S1 <- countMatrix[i,group1] 
    S2 <- countMatrix[i,group2] 
    In <- c(S1,S2)
    In <- as.vector(unlist(In))
    GLM.Poisson <- glm(In ~ 1 + sample.f + offset(offsets),family=poisson)
    Poisson.glm.pval[i] <- anova(GLM.Poisson,test="Chisq")[5][2,1]
    Fold.changes[i] <- exp(GLM.Poisson$coefficients[1])/(exp(GLM.Poisson$coefficients[1]+GLM.Poisson$coefficients[2]))
    if(i %% 100==0) cat(".")
  }
  cat("\n")
  
  #output <- matrix(ncol=2,nrow=nrow(countMatrix))
  #output[,1] <- Poisson.glm.pval
  #output[,2] <- Fold.changes
  #output <- as.data.frame(output)
  #names(output) <- c("pval","FC")
  
  list(stats=data.frame(pval=Poisson.glm.pval, FC=Fold.changes),offsets=offsets,factors=CS)
}


exactTestPoisson <- function(dataMatrix, meanMatrix, group1Ind, group2Ind, verbose=TRUE) {
  
  y1 <- rowSums(dataMatrix[,group1Ind])
  y2 <- rowSums(dataMatrix[,group2Ind])
  m1 <- rowSums(meanMatrix[,group1Ind])
  m2 <- rowSums(meanMatrix[,group2Ind])
  
  N <- rowSums( dataMatrix[,c(group1Ind,group2Ind)] )
  
  pvals <- rep(NA, nrow(dataMatrix))
  
  for (i in 1:length(pvals)) {
    v <- 0:N[i]
    p.top <- dpois(v, lambda=m1[i]) * dpois(N[i]-v, lambda=m2[i])
    p.obs <- dpois(y1[i], lambda=m1[i]) * dpois(y2[i], lambda=m2[i])
    p.bot <- dpois(N[i], lambda=m1[i]+m2[i])
    keep <- p.top <= p.obs
    pvals[i] <- sum(p.top[keep]/p.bot)
    if (verbose)
        if (i%%1000 == 0)
          cat(".")
  }
  if (verbose)
    cat("\n")
    
  pvals

}

calcFactorRLM <- function(obs, ref, logratioTrim=.20, sumTrim=0.01) {

  if( all(obs==ref) )
    return(1)

  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))         # log ratio of expression, accounting for library size
  p0 <- obs/nO
  pR <- ref/nR
  
  x <- log2(p0)-log2(pR)
  x <- x[ !is.na(x) & is.finite(x) ]
  
  r <- rlm(x~1, method="MM")
  2^r$coef
  
}

calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05) {

  if( all(obs==ref) )
    return(1)

  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))         # log ratio of expression, accounting for library size
  absE <- log2(obs/nO) + log2(ref/nR)     # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref  # estimated asymptotic variance
  
  fin <- is.finite(logR) & is.finite(absE)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  
}

calcFactor2 <- function(obs, ref) {
  logR <- log2(obs/ref)
  fin <- is.finite(logR)
  d<-density(logR,na.rm=TRUE)
  2^d$x[which.max(d$y)]
}


fdPlot <- function( score, indDiff, add=FALSE, xlab="Number of Genes Selected", 
                    ylab="Number of False Discoveries", lwd=4, type="l", ... ) {
  o <- order(score)
  w <- o %in% indDiff
  x <- 1:length(indDiff)
  y <- cumsum(!w[indDiff])
  matplot(x, y, xlab=xlab, ylab=ylab, lwd=lwd, type=type, add=add, ... )
}

