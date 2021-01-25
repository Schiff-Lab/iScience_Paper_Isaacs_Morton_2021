#deconvolution of RNA
load("eSet.rda")
load("raw_filter_normalize100.rda")

library(yarn)
library(immunedeconv)
library(pheatmap)
library(dplyr)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

es2 <- filterLowGenes(eSet,minSamples=18)
es2 <- es2[-which(is.na(fData(es2)$gene_name)),]
cnts <- exprs(es2)
rownames(cnts) <- fData(es2)$gene_name
width = (fData(es2)$end-fData(es2)$start)
tpm <- tpm3(counts = cnts,len = width)
xcell <- immunedeconv::deconvolute(tpm, 'xcell')
x2 <- xcell
xcell <- x2 %>% as.matrix
rownames(xcell) <- xcell[,1]
xcell = xcell[,-1]
for(i in seq(ncol(xcell))) xcell[,i] = as.numeric(as.character(xcell[,i]))

for(i in 1:nrow(pData(es2))){
  if(pData(es2)$IPM_ID[i] %in% NS_meta$IPM_ID){
    pData(es2)$new_paeni[i] <- NS_meta$Col_OR_PSU_Paeni[NS_meta$IPM_ID==pData(es2)$IPM_ID[i]]    
  }else{
    pData(es2)$new_paeni[i] <- "NA"
  }
}

xcell <- xcell[,colnames(xcell) %in% NS_meta$IPM_ID]

df <- array(NA,dim=c(nrow(xcell),ncol(xcell)))
for(i in seq(ncol(xcell))) df[,i] = as.numeric(as.character(xcell[,i]))
rownames(df) <- rownames(xcell);  colnames(df) <- colnames(xcell);
annodf <- data.frame(Paeni=factor(pData(es2)$new_paeni),Hydrocephalus=pData(es2)$Hydrocephalus,row.names = colnames(df))
pheatmap(df,annotation_col = annodf)

paeni <- NS_meta %>% filter(Col_OR_PSU_Paeni=="Y")
xcell_paeni <- xcell[,which(colnames(xcell) %in% paeni$IPM_ID)]
df <- array(NA,dim=c(nrow(xcell_paeni),ncol(xcell_paeni)))
for(i in seq(ncol(xcell_paeni))) df[,i] = as.numeric(as.character(xcell_paeni[,i]))
rownames(df) <- rownames(xcell_paeni);  colnames(df) <- colnames(xcell_paeni);
annodf <- data.frame(Hydrocephalus=pData(es2)$Hydrocephalus[pData(es2)$IPM_ID %in% paeni$IPM_ID],row.names = colnames(df))
pheatmap(df,annotation_col = annodf)