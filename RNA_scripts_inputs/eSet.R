#espression set: one with all data, another with filtered
library(Biobase)
library(dplyr)

#load all data, make raw dataset
# NS_data <- read.delim("data100/Sepsis_100RNA.genes.STAR.results.txt", row.names = 1, header=T, stringsAsFactors = F)
# NS_data <- NS_data[grep("^ENSG", rownames(NS_data)),]
# NS_data <- round(NS_data)
load("raw_filter_normalize100.rda")
NS_matrix <- NS_data %>% as.matrix()
minimalSet <- ExpressionSet(assayData=NS_matrix)

NS_data <- NS_data[,colnames(NS_data) %in% NS_meta$IPM_ID]
for(i in 1:ncol(NS_data)){
  colnames(NS_data)[i] <- NS_meta$IPM_ID[NS_meta$IPM_ID==colnames(NS_data)[i]]
}

### Check that sample names match in both files- this is crucial or it will not work
all(names(NS_data) %in% rownames(NS_meta))
all(names(NS_data) == rownames(NS_meta))

#create datasets, start with metadata
metadata <- colnames(NS_meta) %>% as.data.frame()
metadata[,2] <- c("patient_ID", "study_ID", "hydrocephalus_status", "csf_cell_count", "age", "sex", "blood_WBC", "HIV_status", "home_region", "BCx", "CSFCx", "GA", "fever", "Apgar", "Mort", "meta_ID")

phenoData <- new("AnnotatedDataFrame", data=NS_meta, varMetadata=metadata)
head(pData(phenoData))

experimentData <- new("MIAME", name="Sarah Morton", lab="CONSHA", contact="sarah.morton@childrens.harvard.edu", title="100_PIH", abstract="Initial Dataset",  other=list(notes="Created from counts files"))

#next step is to add features of ENSG genes
library(refGenome)
# ens <- ensemblGenome()
# read.gtf(ens, "../../../../Downloads/Homo_sapiens.GRCh38.94.gtf")
# save(ens, file="ens_genes_hg38.rda")
load("ens_genes_hg38.rda")
genes <- getGenePositions(ens)
NS_genes <- rownames(NS_data) %>% as.data.frame()
NS_genes$short <- sapply(strsplit(as.character(NS_genes$.), "\\."), "[", 1)
NS_ens <- genes %>% filter(gene_id %in% NS_genes$short)

NS_feature <- merge(NS_genes, NS_ens, by.x = "short", by.y = "gene_id", all.x = T)
rownames(NS_feature) <- NS_feature$.

featureData <- new("AnnotatedDataFrame", data=NS_feature)
head(pData(featureData))

head((NS_matrix))
head(pData(phenoData))
all(names(NS_matrix) %in% names(phenoData))

eSet <- ExpressionSet(assayData=NS_matrix, phenoData=phenoData, experimentData=experimentData, featureData = featureData)
#save(eSet, file = "eSet.rda")

#ExpressionSet(NS_data, phenoData=annotatedDataFrameFrom(NS_meta, byrow=FALSE), featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), experimentData=MIAME(), annotation=character(), protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE))

#now redo with filtered data
NS_matrixF <- NS_data1 %>% as.matrix()

NS_featureF <-NS_feature %>% filter(. %in% rownames(NS_data1))
rownames(NS_featureF) <- NS_featureF$.
featureDataF <- new("AnnotatedDataFrame", data=NS_featureF)
eSetF <- ExpressionSet(assayData=NS_matrixF, phenoData=phenoData, experimentData=experimentData, featureData = featureDataF)
#save(eSetF, file = "eSetF.rda")

#now redo with normalized, filtered data
NS_matrixNF <- NS_data1_LOG2 %>% as.matrix()

NS_featureNF <-NS_feature %>% filter(. %in% rownames(NS_data1_LOG2))
rownames(NS_featureNF) <- NS_featureNF$.
featureDataNF <- new("AnnotatedDataFrame", data=NS_featureNF)
eSetNF <- ExpressionSet(assayData=NS_matrixNF, phenoData=phenoData, experimentData=experimentData, featureData = featureDataNF)
#save(eSetNF, file = "eSetNF.rda")

#save(eSet, eSetF, eSetNF, file = "all_eSets.rda")
