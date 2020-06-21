#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("Logolas")

library(Logolas)

# Load the required packages
require(ggplot2)
require(ggseqlogo)

char_data <- read.csv("pssm.csv", stringsAsFactors = F,check.names=FALSE,row.names=1)
char_data

# R format for accesing data 
# char_data[,2:3]

modulus=15
for (i in 1:ncol(char_data)) {
  if ( (i %% modulus) & (i != ncol(char_data)) ){
    next
  }
    str=sprintf("residue_%d.png",i)
    png(filename=str)
    end=i    
    start=( i-modulus)
    logo_pssm(char_data[,start:end], control = list(round_off = 0))
    dev.off()
}
