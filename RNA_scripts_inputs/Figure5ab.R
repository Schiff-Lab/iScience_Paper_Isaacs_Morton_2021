#WGCNA: identify DEG based on Paeni status
library(dplyr)
options(stringsAsFactors = F)
library(edgeR)
library(ggplot2)
library(WGCNA)
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)

#load our data, filtered and log normalized
load("raw_filter_normalize100.rda")
#now change to gene names
NS_data1_LOG2$ensembl <- sapply(strsplit(rownames(NS_data1_LOG2), "\\."), "[", 1)
library(org.Hs.eg.db)
ensembl_map <- NS_data1_LOG2$ensembl %in% mappedRkeys(org.Hs.egENSEMBL) #10662 of 11114
map <- select(org.Hs.eg.db, key=NS_data1_LOG2$ensembl,columns=c("ENTREZID", "SYMBOL"),keytype="ENSEMBL")
map <- map[!duplicated(map$ENSEMBL),]
NS_data1_LOG2 <- NS_data1_LOG2 %>% filter(ensembl %in% map$ENSEMBL)

for(i in 1:nrow(NS_data1_LOG2)){
  NS_data1_LOG2$symbol[i] <- map$SYMBOL[map$ENSEMBL==NS_data1_LOG2$ensembl[i]]
}

#subset to genes from edgeR analysis
DEG <- read.delim("Paeni_status_TOTAL_SIGNIFICANT_29June20.txt")

NS_data1_LOG2 <- NS_data1_LOG2 %>% filter(NS_data1_LOG2$symbol %in% rownames(DEG))
############################step 1: data input and cleaning
#set up data and cluster samples
datExpr0 = as.data.frame(t(NS_data1_LOG2[, -c(99:101)]));
names(datExpr0) = NS_data1_LOG2$symbol;
rownames(datExpr0) = names(NS_data1_LOG2)[-c(99:101)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust) #
par(cex = 0.6);
par(mar = c(0,4,2,0))
sampleTree$cluster <- clust
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2, lab=sampleTree$cluster)

#metadata
brief_meta <- NS_meta[,c("IPM_ID", "Hydrocephalus", "Col_OR_PSU_Paeni", "CSF_cells", "CMV_CSF", "Score", "WBC_blood", "CMV_Blood"),]

dim(brief_meta)
names(brief_meta)

# Form a data frame analogous to expression data that will hold the clinical traits.
NSSamples = rownames(brief_meta);
traitRows = match(NSSamples, brief_meta$IPM_ID);
datTraits = brief_meta[traitRows, -1];
rownames(datTraits) = brief_meta[traitRows, 1];
collectGarbage();
datTraits_sub <- datTraits %>% filter(rownames(datTraits) %in% rownames(datExpr))

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
datTraits_sub$group <- datTraits_sub$Hydrocephalus
datTraits_sub$group <- gsub("NPIH", 2, datTraits_sub$group)
datTraits_sub$group <- gsub("PIH", 1, datTraits_sub$group)
datTraits_sub$paeni_status <- datTraits_sub$Col_OR_PSU_Paeni
datTraits_sub$paeni_status <- gsub("Y", 1, datTraits_sub$paeni_status)
datTraits_sub$paeni_status <- gsub("N", 2, datTraits_sub$paeni_status)
datTraits_sub$cmv_status <- datTraits_sub$CMV_CSF
datTraits_sub$cmv_status <- gsub("Y", 1, datTraits_sub$cmv_status)
datTraits_sub$cmv_status <- gsub("N", 2, datTraits_sub$cmv_status)
datTraits_sub$cmv_blood<- datTraits_sub$CMV_Blood
datTraits_sub$cmv_blood <- gsub("Y", 1, datTraits_sub$cmv_blood)
datTraits_sub$cmv_blood <- gsub("N", 2, datTraits_sub$cmv_blood)

#remove uniform values and non-numeric
datTraits_sub <- datTraits_sub[,c("CSF_cells", "Score", "WBC_blood", "paeni_status", "cmv_status", "cmv_blood")]

traitColors = numbers2colors(as.numeric(datTraits_sub$paeni_status), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits)[2], 
                    main = "Sample dendrogram and trait heatmap")

#automatic network assembly
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "cluster2", 
                       verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#relate modules to gene expression
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits_sub, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               #xLabels = names(datTraits_sub),
               xLabels = c("CSF_Cells", "CT_Score", "Blood_WBCs", "Paeni_Status", "CSF_CMV", "Blood_CMV"),
               #yLabels = names(MEs),
               yLabels = c("Module1", "Module2"),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#make corrplot using WGCNA data for Paeni vs non-Paeni from newmeta_oct19.R

library(dplyr)
library(ggcorrplot)

correlations_m1 <- c(0.56, 0.31, -0.37, -0.15, -0.25)
correlations_m2 <- c(-0.11, 0.10, 0.42, -0.07, -0.19)
corr <- rbind(correlations_m1, correlations_m2) %>% as.data.frame()
rownames(corr) <- c("Module 1", "Module 2")
colnames(corr) <- c("CSF Cell Count", "Blood WBC Count", "Paeni Status", "CSF CMV Status", "Blood CMV Status")

pvalues_m1 <- c(0.0007, 0.08, 0.036, 0.42, 0.16)
pvalues_m2 <- c(0.56, 0.58, 0.016, 0.69, 0.29)
p.mat <- rbind(pvalues_m1, pvalues_m2) %>% as.data.frame()
rownames(p.mat) <- c("Module 1", "Module 2")
colnames(p.mat) <- c("CSF Cell Count", "Blood WBC Count", "Paeni Status", "CSF CMV Status", "Blood CMV Status")

ggcorrplot(corr, method = "circle")

library(pheatmap)
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("red", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "blue"))(length(my.breaks)/2))
pheatmap(corr, cluster_rows = F, 
         color = my.colors,
         breaks = my.breaks
)
