# Set working directory if necessary
setwd("/Users/albertisaacs/Dropbox/Research Projects/CONSHA/Manuscripts/PIH Proteogenomics/Data Analysis")

##------------------------------LIMMA DIFFERENTIAL EXPRESSION ANALYSIS----------------------------------------------------------------------------

# Install Limma and load libraries
#BiocManager::install("Limma")
library(limma)
library(rafalib)
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

# Read the RDS file
obj4 <- readRDS("pih_proteo_batch_normalized_reduced.rds")

# Perform Limma (by Paeni status)
obj4$CSFCellCount <- as.numeric(numextract(obj4$CSFCellCount))
obj4$CSFCellCountFactor <- factor(ifelse(obj4$CSFCellCount>5,yes="high",no="low"),levels=c("low","high"))
obj4$Paeni_status = factor(obj4$Paeni_status,levels=c("NPIH","NonPaeni_PIH", "Paeni_PIH"))

design <- model.matrix(~ Paeni_status,data=pData(obj4))
fit <- lmFit(exprs(obj4), design)
fit <- eBayes(fit)
res <- topTable(fit,coef=3,number=nrow(obj4),sort.by="logFC")
volcanoplot(fit,coef=3)

#Save output datashe of DE genes
write.csv(res, file = 'differential_abundance_analysis_Paeni_status_coeff3.csv')
