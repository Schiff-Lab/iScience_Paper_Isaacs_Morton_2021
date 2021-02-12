# set working directory
setwd("/Users/albertisaacs/Dropbox/Research Projects/CONSHA/Manuscripts/PIH Proteogenomics/Data Analysis")

#install programs
BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")

#Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(rafalib)
library(EnhancedVolcano)
library(pheatmap)

#Read Limma .csv output file of DE genes
d <- read.csv("/Users/albertisaacs/Dropbox/Research Projects/CONSHA/Manuscripts/PIH Proteogenomics/Data Analysis/differential_abundance_analysis_Paeni_status_coeff3.csv")
geneList <- d[,2]
names(geneList) <- as.character(d[,1])

#Create volcano plot
EnhancedVolcano(d,
                lab = as.character(d[,1]),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(-3, 3),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 4.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.9,
                drawConnectors = FALSE,
                widthConnectors = 0.3)


#Use groupGO function for gene classification based on GO distribution
gene.df <- bitr(names(geneList), fromType = "SYMBOL",
                toType = c("ENTREZID","SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
geneEntrez <- gene.df$ENTREZID
names(geneList) = geneEntrez

cc = groupGO(geneEntrez,org.Hs.eg.db,ont='CC')
bp = groupGO(geneEntrez,org.Hs.eg.db,ont='BP')


#Create dot plot for BP
x_pval_only <- enrichGO(names(geneList[which(d$adj.P.Val<0.05)]),OrgDb = org.Hs.eg.db,ont='BP')
x <- enrichGO(names(geneList[which(d$adj.P.Val<0.05 & abs(d$logFC)>1)]),OrgDb = org.Hs.eg.db,ont='BP')

enrichplot::dotplot(x)
p1 <- enrichplot::dotplot(x, showCategory=30) + ggtitle("dotplot for ORA")

#Plot Gene-concept Networks
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)

#Plot Heat plots
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
heatplot(edox, foldChange=geneList)

#Plot enrichment map
emapplot(edox)

