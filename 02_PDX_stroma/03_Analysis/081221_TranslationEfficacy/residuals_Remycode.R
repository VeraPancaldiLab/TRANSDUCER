# Remy Original code
library(tidyverse)
library(MASS)
library(car)
### KNOWN DIFFERENCES
#  - keep genes with some 0s (they would get discarded afterwards with the lm)

setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")


# dataT = cpm(calcNormFactors(DGEList(counts = Hrawallrna[,subsamps]), method = "TMM"), # normalize by TMM & log2 transform
#             normalized.lib.sizes = T, log = T)
# dataP = cpm(calcNormFactors(DGEList(counts = Hrawriborna[,subsamps]), method = "TMM"),
#             normalized.lib.sizes = T, log = T)[rownames(dataT),]

# here I load mine to understand the rest of the code
read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv") %>%
  column_to_rownames("EnsemblID") %>% data.matrix() -> dataT

read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv") %>%
  column_to_rownames("EnsemblID") %>% data.matrix() %>% .[rownames(dataT), colnames(dataT)] -> dataP


#residualCompute=mclapply(I,function(i){                                           # Parallel Versions of lapply and mapply using Forking. This iterate through genes
residualCompute=mclapply(rownames(dataT),function(i){
  print(i)
  geneID=i                                                    # get gene names
  
  lmfit=lm(dataP[i,]~dataT[i,])                                                   # do a lm(y~x)
  
  rlmfit=rlm(dataP[i,]~dataT[i,])                                                 # do a robust lm(y~x)
  
  
  residus=resid(lmfit,type="response")                                            # extract the lm residuals (no rlm)
  
  residu=dataP[i,]-lmfit$fitted.values                                            # extract the residuals in another way (checked by all.equal)
  r2=1-(sum(residu^2)  / sum( (dataP[i,] -mean(dataP[i,]))^2 ))                   # Coeficient of determination (R2) of the model (how well observed outcomes are replicated by the model)
  
  
  shapitest=shapiro.test(residu)                                                  # Normality check
  
  slope=lmfit$coefficients[2]                                                     # slope
  
  lmaov <- anova(lmfit)                                                           # do an anova of the LM (how the total affect Polysome)
  MSE=lmaov[1,3]                                                                  ## get standard error
  Fval=lmaov[1,4]                                                                 ## get F value (difference between groups/inside groups)
  Pval=lmaov[1,5]                                                                 ## significance (of what?) we are not testing the influence of anything (perhaps using dataT as categorical?)
  
  
  ncvfit=car::ncvTest(lmfit)                                                      # Score Test for Non-Constant Error Variance (heterocedasticity) go through residuals relaetd to vars
  nonconstvarP=ncvfit$p                                                           ## save Pval
  nonconstvarStat=ncvfit$ChiSquare                                                ## save stat
  
  
  outtest=car::outlierTest(lmfit)                                                 # Bonferroni p-values for testing each observation in turn to be a mean-shift outlier,
  
  
  residuQC=data.frame(                                                            # Save everything in a dataframe and
    GeneName= geneID,                                       ## take from an outside table somehow the gene name associated to id
    ResiduSD=sd(residu),                                                          ## Res SD
    ResiduNormalityStat=shapitest$statistic,                                      ## Stats of normality
    ResiduNormalityPvalue=shapitest$p.value,                                      ## Pval of normality
    ResiduNormalityLogPvalue=-log10(shapitest$p.value),                           ## Pval of normality (in log)
    Slope=slope,                                                                  ## Slope
    # R2=summary(lmfit)$r.squared,
    R2=r2,                                                                        ## the R2
    MSE=MSE,AOVFval=Fval,AOVPval=Pval,AOVLogPval=-log10(Pval),                    ## all the anota stuff stat pval and logpval
    NonConstStat=nonconstvarStat,NonConstPval=nonconstvarP,NonConstLogPval=-log10(nonconstvarP),  ## all the heterocedasticity stuff stat pval and logpval
    CountOutlier=sum(outtest$bonf.p<0.05),                                        ## 1 if the bonferronu pval is significative 0 if it isnt 95% IC
    row.names=geneID,stringsAsFactors = F)
  outtest$bonf.p[colnames(dataT)]<0.05
  
  outliers=outtest$bonf.p[colnames(dataT)]<0.05
  names(outliers) <- colnames(dataT)
  
  list(residuals=residu,robustresdiuals=rlmfit$residuals,QC=residuQC,outliers=outliers)
},mc.cores=14)                                                                    # This produces a list of list

                                                                                 

outlierMatrix=do.call(rbind,lapply(residualCompute,function(x)x$outliers))        # which here are turn into columns of a df
residualQC=do.call(rbind,lapply(residualCompute,function(x)x$QC))
residuals=do.call(rbind,lapply(residualCompute,function(x)x$residuals))
robustresiduals=do.call(rbind,lapply(residualCompute,function(x)x$robustresdiuals))

rownames(robustresiduals)=rownames(residuals)=rownames(outlierMatrix)=rownames(dataT)


OKresid=rownames(residualQC)[which(residualQC$CountOutlier==0 & residualQC$NonConstPval >0.01 & residualQC$ResiduNormalityPvalue >0.01)]
selectedResiduals=residuals[OKresid,]
