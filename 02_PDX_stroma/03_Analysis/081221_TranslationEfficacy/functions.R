library(anota2seq)
library(RColorBrewer)
################################################################################
################################################################################
#' Plot PCA components in relation to a given factor
#'@description
#' this function does a parplot of the desired n of components and color 
#' them according to a given factor.
#'
#'prints a list of spearman correlations of the desired metadata factor
#'with the IC that best correlate with it, from the ICAs with the different number of components 
#'
#'@param pca_toplot data frame containing the PCs and the factors 
#'you want to use to correlate
#'@param feat name of the feature to color
#'@param ncomp number of PCs to plot
#'@param dotsize to adjust the dotsize manually
#'
#'TBD 
#'continuous coloring
#'

plot_PCs <- function(pca_toplot, feat, ncomp, dotsize){
  col_factor <- as.factor(pca_toplot[[feat]])
  col_n <- nlevels(col_factor)
  cols <- brewer.pal(col_n, "Spectral")
  cols <- colorRampPalette(cols)(col_n)
  pairs(pca_toplot[, paste("PC",1:ncomp,sep ="")],
        col = cols[as.numeric(col_factor)],
        pch = 19,
        cex = dotsize,
        lower.panel = NULL)
  
  par(xpd = TRUE)
  x11()  
  plot.new()
  legend(x = "center",fill = cols, legend = levels(col_factor), horiz = F)
}

################################################################################

################################################################################
# Translation efficacy analysis
## To do:
## - More efficient lm test, 
## - retreieve a matrix with details about more than one outlier (give back more than one matrix)
calculaTE <- function(x)
{
  dplyr::filter(normHost_cyt, EnsemblID == x)[-1] %>% unlist() -> cyt
  dplyr::filter(normHost_pol, EnsemblID == x)[-1] %>% unlist() -> pol
  fit <- lm(pol~cyt)
  
  slope <- fit$coefficients[["cyt"]]
  residuals <- fit$residuals
  homoscedasticity <- bptest(fit,studentize = TRUE) # Koenker–Bassett test (homoscedasticity is H0)
  normality <- shapiro.test(residuals)
  #cooksd <- cooks.distance(fit)
  #outl <- any(cooksd > (3 * mean(cooksd, na.rm = TRUE)))
  #outl <- any(cooksd > 4/21)
  outl <- car::outlierTest(fit)
  return(c(x, residuals, slope, homoscedasticity$p.value[[1]], normality$p.value, outl$bonf.p[[1]], names(outl$bonf.p)[[1]]))
}
################################################################################
################################################################################
## Anota2seq-like set of functions to
## perform gene-wise outlier test
## see https://rdrr.io/bioc/anota2seq/src/R/anota2seqResidOutlierTest.R

anota2seqResidOutlierPlotAll <- function(all=NULL, xsAll=xsAll,  env=env, obtained, expected, obtRelExpected, confInt){
  allMax=max(all)
  allMin=min(all)
  xsMin=min(xsAll)
  xsMax=max(xsAll)
  ##Get number of outliers per rankposition
  allLog <- all<env[,1] | all>env[,2]
  allLogSum <- apply(allLog, 1, sum)
  allLogSumP <- 100*(allLogSum / dim(allLog)[2])
  allLogSumP <- round(allLogSumP, digits=3)
  ##plot
  plot(x=c(xsMin,xsMax), y=c(allMin-0.2, allMax+0.4), pch="", axes=TRUE, xlab="Quantiles of standard normal", ylab="R", main="Summary of all residuals")
  for(i in 1:dim(all)[2]){
    points(x=xsAll[,i], y=all[,i], pch=16, cex=0.2)
  }
  if(is.null(obtained)==FALSE){
    text(x=xsMin+0.4, y=allMax-0.4-(0.05*allMax), labels=paste("Expected outliers: ", confInt*100, "%", sep=""))
    text(x=xsMin+0.4, y=allMax-0.4-(0.1*allMax), labels=paste("Obtained outliers: ", round(obtRelExpected*confInt*100, digits=3), "%", sep=""))
  }
  segments(xsAll[,i]-0.05, env[,1], xsAll[,i]+0.05, env[,1], col=2, lwd=1.5)
  segments(xsAll[,i]-0.05, env[,2], xsAll[,i]+0.05, env[,2], col=2, lwd=1.5)
  text(x=xsAll[,1], y=env[,2]+0.2, labels=paste(allLogSumP, "%", sep=""))
}

OutlierTest <- function(residualMatrix, confInt=0.01, iter=5, nGraphs=200, #need testing
                        generateSummaryPlot=TRUE, residFitPlot=TRUE){
  
  # Get data
  nData <- dim(residualMatrix)[1]
  geneNames <- rownames(residualMatrix)
  
  # Data structures
  residualMatrixOutlierSum <- matrix(ncol=iter,                                  # to save how many samples are outlier per gene
                                     nrow=nData,
                                     dimnames=list("rownames"=rownames(residualMatrix)))
  
  residualMatrixOutlier <- matrix(ncol=dim(residualMatrix)[2],                   # to save for each gene wether a sample was an outlier or not
                                  nrow=nData,
                                  dimnames=list("rownames"=rownames(residualMatrix), 
                                                "colnames"=colnames(residualMatrix)))
  allSortScaleTrue <- allXs <- matrix(nrow=dim(residualMatrix)[2], ncol=nData)
  residualMatrixOutlierSumP <- rep(NA, iter)
  for(j in 1:iter){                                                               # Now iter times this:
    ##generate the random normally distributed data.                              ## sample 99(when confint == 0.01)*your n samples from a normal distribution.
    rnormMat <- matrix(data=rnorm(dim(residualMatrix)[2]*((1/confInt)-1)),        ## samples as rows, 99 iter as cols
                       nrow=dim(residualMatrix)[2], ncol=((1/confInt)-1)) 
    
    ##Calculate the upper and lower limits of the data and scale
    rnormMat <- apply(scale(rnormMat),2,sort)                                     ## this scales each "randomized" sample to mean 0 std 1
    env <- t(apply(rnormMat, 1, range))                                           ## This gets the max and min of the rNorm once scaled (they call it envelope)
    for(i in 1:nData){                                                            ## per gene i
      ##Scale true data per gene sort and cbind to rnorm data set
      trueVec <- sort(scale(residualMatrix[i,]))                                  ### scale and sort the gene i (order is kept as it is a named vector)
      sampMat <- cbind(trueVec, rnormMat)                                         ### create a new mat with this vector as 1st column (genes as columns)
      ##get real data set qq
      rs <- sampMat[,1]
      xs <- qqnorm(rs, plot=FALSE)$x                                              ### get the theoretical quantiles of your gene (x axis of a qqplot)
      ##get range of the sampled distribution sort position
      ##Calculate if obtained residuals falls outside expected from rnorm
      rsLog <- rs<env[,1] | rs>env[,2]                                            ### Is the real sample between the envelopes of the simulated data?
      residualMatrixOutlierSum[i,j] <- sum(rsLog)                                 ### retrieve the sum of samples where this is unfulfilled in this j iteration
      residualMatrixOutlier[i,] <- rsLog                                          ### This saves in the matrix which samples did not fullfilled in this iteration
      ##save true data
      allSortScaleTrue[,i] <- trueVec                                             ### This saves the sorted scaled vector (loose sample and gene names)
      allXs[,i] <- xs                                                             ### This saves the theoretical quantiles (loose sample and gene names)
    }
    ##Collect data for single
    residualMatrixOutlierSumP[j] <- sum(residualMatrixOutlier>0)                  ## Save how many sample/genes are outliers in the iteration FROM DOWN HERE IT JUST DESCRIBES and compare with expected.
  }                                                                               ##not usefull for filtering purposes
  ################################################
  ##calcaulte obtained expected
  ##only create full summary for last iteration
  residualMatrixOutlierLog <- residualMatrixOutlier>0                             # why this line if this is already boolean!
  residualMatrixOutlierSumP <- residualMatrixOutlierSumP/(nData*dim(residualMatrixOutlier)[2]) # proportion of outliers of the total per iteration
  obtVsExpected <- residualMatrixOutlierSumP[j]/confInt                           # ratio compared with the 1% expected
  expected <- nData *dim(residualMatrixOutlier)[2] *confInt                       
  obtained <- sum(residualMatrixOutlierLog)
  ################################################
  # outputList <- new("Anota2seqResidOutlierTest",
  #                   confInt = confInt,
  #                   inputResiduals = residualMatrix,
  #                   rnormIter = iter,
  #                   outlierMatrixLog = residualMatrixOutlierLog,
  #                   meanOutlierPerIteration = residualMatrixOutlierSumP,
  #                   obtainedComparedToExpected = obtVsExpected,
  #                   nExpected = expected,
  #                   nObtained = obtained)
  #################################################
  #Plotting summary
  if(generateSummaryPlot==TRUE){
    # jpeg("ANOTA2SEQ_residual_distribution_summary.jpeg", width=800, height=800, quality=100)
    anota2seqResidOutlierPlotAll(all=allSortScaleTrue, xsAll=allXs, env=env, obtained=obtained, expected=expected, obtRelExpected=obtVsExpected, confInt=confInt)
    # dev.off()
  }
  ##plot fitted vs residuals
  #   if(residFitPlot==TRUE){
  #     jpeg("ANOTA2SEQ_residual_vs_fitted.jpeg", width=900, height=900, quality=100)
  #     par(mfrow=c(2,1))
  #     plot(x=as.vector(Anota2seqDataSet@qualityControl@fittedValues), y=as.vector(Anota2seqDataSet@qualityControl@residuals), ylab="residuals", xlab="Fitted values", main="Residual vs fitted values")
  #     dev.off()
  #   }
  # # 
  #   Anota2seqDataSet@residOutlierTest <- outputList
  #   return(Anota2seqDataSet)
}
