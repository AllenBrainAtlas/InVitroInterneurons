####Identify high variability genes using Brennecke's method#####
run_brennecke<-function(refv,sampv, plot.fig=NULL) {
	require("DESeq");
	require("statmod") 
	require("matrixStats")
	#refv<-alldat$e;
	#sampv<-alldat$g;
	remrow<-grep("tdTomato",rownames(refv));
	if (length(remrow)>0) {refv<-refv[-remrow,]}
	sfref <- estimateSizeFactorsForMatrix(refv)
	sfsamp <- estimateSizeFactorsForMatrix(sampv)
	nCountsref <- t( t(refv) / sfref )
	nCountssamp <- t( t(sampv) / sfsamp )
	meansref <- rowMeans( nCountsref )
	varsref <- rowVars(nCountsref)
	cv2ref <- varsref / meansref^2

	minMeanForFit <- unname( quantile( meansref[ which( cv2ref > .3 ) ], .80 ) )

	useForFit <- meansref >= minMeanForFit
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansref[useForFit] ),cv2ref[useForFit] )

	xi <- median( 1 / sfref )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"] - xi )
	c( a0, a1 )

	###plot fit###
        if(!is.null(plot.fig)){
          pdf(plot.fig)
          plot( NULL, xaxt="n", yaxt="n",log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
          axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
          axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10","100" ), las=2 )
          abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
          
          # Add the data points
          points( meansref, cv2ref, pch=20, cex=1, col="blue" )
	  # Plot the fitted curve
          xg <- 10^seq( -2, 6, length.out=1000 )
          lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
	  # Plot quantile lines around the fit
          df <- ncol(refv) - 1
          lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df, col="#FF000080", lwd=2, lty="dashed" )
          lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df, col="#FF000080", lwd=2, lty="dashed" )
        }

	#####test samples####
	meanssamp <- rowMeans( nCountssamp )
	varssamp <- apply(nCountssamp,1,var)
	varssamp <- rowVars(nCountssamp)
	cv2samp <- varssamp / meanssamp^2
	psia1theta <- mean( 1 / sfsamp ) + a1 * mean( sfref / sfsamp )
	minBiolDisp <- 0.25^2
	m <- ncol(sampv)
	cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
	testDenom <- ( meanssamp * psia1theta + meanssamp^2 * cv2th ) / ( 1 + cv2th/m )
	p <- 1 - pchisq( varssamp * (m-1) / testDenom, m-1 )
	padj <- p.adjust( p, "BH" )
	sig <- padj < .05
	sig[is.na(sig)] <- FALSE
	print(table( sig ))

	#####plot samples####
        if(!is.null(plot.fig)){
          #Plot the plant genes, use a different color if they are highly variable
          points( meanssamp, cv2samp, pch=20, cex=.2,col = ifelse( padj < .05, "#C0007090", "darkgreen" ))
	  # Add the technical noise fit, as before
          xg <- 10^seq( -2, 6, length.out=1000 )
          lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
	  # Add a curve showing the expectation for the chosen biological CV^2 thershold
          lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df + minBiolDisp, col="#FF000080", lwd=2, lty="dashed" )
          lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )
          dev.off()
        }

	####Table of genes####
	fullVarTable <- data.frame(geneSymbol = rownames(sampv),padj = padj,meanNormCount = meanssamp, cv2 = cv2samp)
        return(fullVarTable)
      }

