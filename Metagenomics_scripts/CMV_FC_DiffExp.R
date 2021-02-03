##RNASeqHelpSaraAlbert
library(limma)
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

# Read the RDS file
obj4 <- readRDS("pih_proteo_batch_normalized_reduced_071019.rds")

# Limma by Paeni status
obj4$CSFCellCount <- as.numeric(numextract(obj4$CSFCellCount))
obj4$CSFCellCountFactor <- factor(ifelse(obj4$CSFCellCount>5,yes="high",no="low"),levels=c("low","high"))
obj4$Paeni_status = factor(obj4$Paeni_status,levels=c("NPIH","NonPaeni_PIH", "Paeni_PIH"))
obj4$CMV = factor(obj4$M+CMV,levels=c("positive","negative"))

# temporarily remove CSFCellCountFactor - missing sample
# temporarily removed factor(Paeni_status) + factor(Region)+ WBC
design <- model.matrix(~ Paeni_status+CMV,data=pData(obj4))
fit <- lmFit(exprs(obj4), design)
fit <- eBayes(fit)
res <- topTable(fit,coef=3,number=nrow(obj4))
volcanoplot(fit,coef=3)
write.csv(res, file = 'differential_abundance_analysis_Paeni_status_coeff_V_WITHCMV.csv')


