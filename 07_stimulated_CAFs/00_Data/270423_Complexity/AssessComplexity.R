library(tidyverse)

################################################################################
mySampleFun <- function(x, names, sampSize=1e6){
  data <- x
  ##The idea here is to sample a fixed number of reads to enable a more comparable analysis
  ##The input is a matrix so we will use an apply function
  tmpData <- rep(names, data)
  if(length(tmpData)>sampSize){
    tmpSamp <- sample(tmpData, sampSize)
    tmpSum <- table(tmpSamp)
    tmpSum2 <- as.vector(tmpSum)
    names(tmpSum2) <- names(tmpSum)
    tmpOut <- c(rep(0, length(names)))
    names(tmpOut) <- names
    tmpOut[names(tmpSum2)] <- tmpSum2[names(tmpSum2)]
  } else {
    tmpOut <- c(rep(0, length(names)))
    names(tmpOut) <- names
  }
  return(tmpOut)
}

PlotComplexity <- function(){
  
}

################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/270423_Complexity/")


# load expression data
datExpr <-  read_tsv("../rawcounts.tsv") %>% 
  column_to_rownames("Geneid") %>% 
  as.matrix()

# load Metadata
datInfo <- read_tsv("../Data_RNA_sample_Jacobo.tsv")

Pol_pos <- which(str_detect(colnames(datExpr), "F8"))
Tot_pos <- which(str_detect(colnames(datExpr), "Input"))

sampSizes <- c(0.1e5,0.2e5,0.3e5,0.4e5,0.5e5,0.6e5,0.7e5,0.8e5,0.9e5,1e5,1.5e5,2e5,2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5,6e5,6.5e5,7e5,7.5e5,8e5,8.5e5,9e5,9.5e5,1e6)

testOut <- list()
for(t in 1:length(sampSizes)){
  testOut[[t]] <- apply(datExpr,2,FUN=mySampleFun,names=rownames(datExpr),sampSize= sampSizes[t])
}

plotTable <- matrix(nrow=ncol(datExpr),ncol=28)
colnames(plotTable) <- paste("test",format(sampSizes,scientific = T),sep="_")
rownames(plotTable) <- colnames(datExpr)
for(d in 1:length(testOut)){
  for(r in 1:nrow(plotTable)){
    plotTable[r,d] <- nrow(testOut[[d]][which(testOut[[d]][,r] >= 2),])
  }
}

par(mar=c(6,6,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.8,cex.main=0.8,cex.lab=1)
plot(sampSizes[which(plotTable[1,]>0)],as.numeric(plotTable[1,][which(plotTable[1,]>0)]),xlab='', ylab='',pch=20,ylim=c(0,18000),xlim=c(0,1e6),col='#F22F08',lwd=2,bty="n",yaxt="n",xaxt="n")

axis(format(sampSizes[c(seq(1,28))],scientific = T,digits = 2),side=1,at=sampSizes[c(seq(1,28))],font=2,las=2,lwd=2)
axis(seq(0,18000,3000),side=2,at=seq(0,18000,3000),font=2,las=2,lwd=2)
mtext(side=2, line=4, 'Number of genes detected', col="black", font=2, cex=1.7)
mtext(side=1, line=5, 'Millions reads sampled', col="black", font=2, cex=1.7)

lines(sampSizes[which(plotTable[1,]>0)],as.numeric(plotTable[1,][which(plotTable[1,]>0)]),lwd=2,col='#F22F08')
for(i in Tot_pos){
  print(i)
  lines(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),lwd=2,col='#F22F08')
  points(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),pch=20,col='#F22F08')
}

for(i in Pol_pos){
  print(i)
  lines(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),lwd=2,col='#8D2F23')
  points(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),pch=20,col='#8D2F23')
}

legend(1,19000,bty='n', fill=c('#8D2F23','#F22F08'), horiz=TRUE, xpd=T, cex=1.5,legend=c('polyRNA','totalRNA'))


complexity_toplot <- as_tibble(plotTable, rownames ="sample_name") %>%
    pivot_longer(cols = -sample_name, names_to = "reads_sampled", values_to = "detected_genes") %>% 
    mutate(reads_sampled = as.numeric(str_remove(reads_sampled, "test_"))) %>% 
    left_join(datInfo, "sample_name")

## Fraction
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = Fraction, group =  sample_name) + 
  geom_line() +
  theme_bw()

## Batch  
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = Batch, group =  sample_name) + 
  geom_line() +
  theme_bw()

## CAF
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = CAF, group =  sample_name) + 
  geom_line() +
  theme_bw()

## Experimentalist
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = Experimentalist, group =  sample_name) + 
  geom_line() +
  theme_bw()

## Condition
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = Condition, group =  sample_name) + 
  geom_line() +
  theme_bw()

## Sample
ggplot(complexity_toplot) +
  aes(x = reads_sampled, y = detected_genes, color = sample, group =  sample_name) + 
  geom_line() +
  theme_bw()

