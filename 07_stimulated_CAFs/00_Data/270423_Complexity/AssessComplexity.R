library(tidyverse)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/00_Data/270423_Complexity/")

# load expression data
datExpr <-  read_tsv("../rawcounts.tsv") %>% 
  column_to_rownames("Geneid") %>% 
  as.matrix()

Pol_pos <- which(str_detect(colnames(datExpr), "F8"))
Tot_pos <- which(str_detect(colnames(datExpr), "Input"))

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

sampSizes <- read_tsv("MarkDuplicates/MarkDuplicatesStats.tsv")  %>% 
  dplyr::mutate(FILE = str_remove_all(FILE,"_Aligned.sortedByCoord.out.bam.txt|__Aligned.sortedByCoord.out.bam.txt") %>%
                  fct(levels = colnames(datExpr))) %>%
  dplyr::arrange(FILE) %>%
  dplyr::select(ESTIMATED_LIBRARY_SIZE) %>%
  deframe()
  

testOut <- list()
for(t in 1:length(sampSizes)){
  testOut[[t]] <- apply(datExpr,2,FUN=mySampleFun,names=rownames(datExpr),sampSize= sampSizes[t])
}

plotTable <- matrix(nrow=ncol(datExpr),ncol=length(sampSizes))
colnames(plotTable) <- paste("test",format(sampSizes,scientific = T),sep="_")
rownames(plotTable) <- colnames(datExpr)
for(d in 1:length(testOut)){
  for(r in 1:nrow(plotTable)){
    plotTable[r,d] <- nrow(testOut[[d]][which(testOut[[d]][,r] >= 2),])
  }
}

limx <- c(min(sampSizes)-1000, max(sampSizes) +1000)
limy <- c(min(plotTable[which(plotTable>0)])-1000, max(plotTable)+ 1000)

par(mar=c(6,6,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.8,cex.main=0.8,cex.lab=1)
plot(sampSizes[which(plotTable[1,]>0)],as.numeric(plotTable[1,][which(plotTable[1,]>0)]),xlab='', ylab='',pch=20,ylim=limy,xlim=limx,col='#F22F08',lwd=2,bty="n",yaxt="n",xaxt="n")

axis(format(sampSizes[c(seq(1,ncol(datExpr)))],scientific = T,digits = 2),side=1,at=sampSizes[c(seq(1,ncol(datExpr)))],font=2,las=2,lwd=2)
axis(seq(0,limy[2],3000),side=2,at=seq(0,limy[2],3000),font=2,las=2,lwd=2)
mtext(side=2, line=4, 'Number of genes detected', col="black", font=2, cex=1.7)
mtext(side=1, line=5, 'Millions reads sampled', col="black", font=2, cex=1.7)

lines(sampSizes[which(plotTable[1,]>0)],as.numeric(plotTable[1,][which(plotTable[1,]>0)]),lwd=2,col='#F22F08')
for(i in Tot_pos){
  print(colnames(datExpr)[i])
  lines(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),lwd=2,col='#F22F08')
  points(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),pch=20,col='#F22F08')
}

for(i in Pol_pos){
  print(colnames(datExpr)[i])
  lines(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),lwd=2,col='#8D2F23')
  points(sampSizes[which(plotTable[i,]>0)],as.numeric(plotTable[i,][which(plotTable[i,]>0)]),pch=20,col='#8D2F23')
}

legend(10e5,24000,bty='n', fill=c('#8D2F23','#F22F08'), horiz=TRUE, xpd=T, cex=1.5,legend=c('polyRNA','totalRNA'))
