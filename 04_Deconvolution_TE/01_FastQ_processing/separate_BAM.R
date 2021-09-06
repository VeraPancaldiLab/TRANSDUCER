wd <- "/home/jacobo/Documents/02_TRANSDUCER/03_Deconvolution/01_FastQ_processing"
setwd(wd)
install.packages("https://bioconductor.org/packages/3.5/bioc/src/contrib/Rsubread_1.26.1.tar.gz", repos=NULL, type="source") # Try with the version of 2017

#library(devtools)
install_github("cit-bioinfo/SMAP")
library("SMAP")
source("01_Input/SMAP/R/SMAPcount.R") #


# Separation of reads
counts = SMAPcount(SMAPBAM = paste(wd, "01_Input", sep="/"),
                   GTF = paste(wd, "02_Output/combined_reference/combined.gtf", sep = "/"))
countsTumor = counts$FCcountsTumor
countsHost = counts$FCcountsHost



## Example data
d = system.file("extdata","Example_STAR_OUTPUT", package = "SMAP")
dcounts = SMAPcount(SMAPBAM=d, GTF = "/home/jacobo/Documents/02_TRANSDUCER/03_Deconvolution/01_FastQ_processing/02_Output/combined_reference/combined.gtf")
countsTumor = counts$FCcountsTumor
countsHost = counts$FCcountsHost

# Normalizations
## TMM
library(edgeR)
countsTumor.edgeR <- DGEList(countsTumor)
countsTumor.edgeR <- calcNormFactors(countsTumor.edgeR,
                                     method = "TMM")
tail(cpm(countsTumor.edgeR))
tail(cpm(countsTumor))
tail(countsTumor/countsTumor.edgeR$samples$norm.factors)

countsTumor.tmm <- cpm(countsTumor.edgeR)
countsTumor.tmm.l <- cpm(countsTumor.edgeR, log = TRUE) # calculate mean dif between these two

## Upper Quartile
library(edgeR)
countsTumor.edgeR <- DGEList(countsTumor)
countsTumor.edgeR <- calcNormFactors(countsTumor.edgeR,
                                     method = "upperquartile")

countsTumor.uq <- cpm(countsTumor.edgeR)
countsTumor.uq.l <- cpm(countsTumor.edgeR, log = TRUE)

## TPM
#https://support.bioconductor.org/p/91218/
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
