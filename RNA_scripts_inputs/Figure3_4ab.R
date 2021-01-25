#DEG by Paeni, CMV status
library(dplyr)
options(stringsAsFactors = F)
library(edgeR)
library(ggplot2)
library(yarn)

##eSet
load("all_eSets.rda")
pData(eSet)$cohort_paeni <- paste(pData(eSet)$Sample_Def, pData(eSet)$Col_OR_PSU_Paeni, sep = ":")
pData(eSet)$legend_data <- pData(eSet)$cohort_paeni
pData(eSet)$legend_data <- gsub("CaseControl:N", "NPIH", pData(eSet)$legend_data)
pData(eSet)$legend_data <- gsub("CaseControl:Y", "NPIH", pData(eSet)$legend_data)
pData(eSet)$legend_data <- gsub("Case:N", "Paeni-negative", pData(eSet)$legend_data)
pData(eSet)$legend_data <- gsub("Case:Y", "Paeni-positive", pData(eSet)$legend_data)

##PCA
myColors <- c("blue","red", "green3")
palette(myColors)
plotCMDS(eSet,pch=21, bg=factor(pData(eSet)$legend_data), cex=2, ylim=c(-40,40), xlim = c(-250,250))
legend("topright",legend=levels(factor(eSet$legend_data)),fill=1:3, box.col=NA, cex = 0.75)
palette(alpha(myColors,0.4))
plotCMDS(eSet,pch=21, bg=factor(pData(eSet)$legend_data), cex=2, ylim=c(-40,40), xlim = c(-250,250))
legend("topright",legend=levels(factor(eSet$legend_data)),fill=1:3, box.col=NA, cex = 0.75)

#load in data
load("raw_filter_normalize100.rda")

for(i in 1:nrow(NS_data1)){
  NS_data1$sum[i] <- sum(as.numeric(NS_data1[1:98,i]))
}
NS_data1$ensembl_brief <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1)

#look at duplicated ensembl
NS_dup <- NS_data1[!duplicated(NS_data1$ensembl_brief),]

#prepare edgeR object
y <- DGEList(counts=NS_data1[,1:98], genes=NS_data1$ensembl_brief)

#add gene names
library(org.Hs.eg.db)
idfound <- y$genes$genes %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
map <- select(org.Hs.eg.db, key=y$genes$genes,columns=c("ENTREZID", "SYMBOL"),keytype="ENSEMBL")
m <- match(y$genes$genes, map$ENSEMBL)
y$genes$Symbol <- map$SYMBOL[m]

#filter and normalize (TMM)
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]

y$samples$lib.size <- colSums(y$counts)
rownames(y$counts) <- rownames(y$genes) <- y$genes$Symbol
y <- edgeR::calcNormFactors(y)

#cpm
logcounts <- cpm(y,log=TRUE)

#assign patient attributes
table(NS_meta$Col_OR_PSU_Paeni)
design <- model.matrix(object = ~NS_meta$Col_OR_PSU_Paeni)
rownames(design) <- rownames(NS_meta)

#estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

#DE genes: Paeni vs Non
fit <- edgeR::glmFit(y, design)
lrt <- edgeR::glmLRT(fit)
edgeR::topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt)) #

#get DEG list
out <- topTags(lrt, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05
sig <- out[keep,] %>% as.data.frame()

padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold_NS <- sig$FDR < padj.cutoff & abs(sig$logFC) > lfc.cutoff
length(which(threshold_NS == TRUE)) #15
sig$threshold <- threshold_NS 
NS_sig<- subset(sig, threshold == TRUE)
sig_genes <- NS_sig %>% as.data.frame()

#heatmap of 500 most DEG
library(gplots)
library(RColorBrewer)
logcpm <- cpm(y, prior.count=2, log=TRUE)
sort_sig <- sig[order(sig$logFC, decreasing = T),]
sig_names <-  sort_sig$Symbol[1:500]
highly_DE_lcpm <- logcpm[sig_names,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("red","green")[as.factor(NS_meta$Col_OR_PSU_Paeni)]

heatmap.2(highly_DE_lcpm,col=rev(morecols(50)),trace="none", main="DE genes across samples",ColSideColors=col.cell,scale="row")

#volcano plot
ggplot(sig) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  ggtitle('NS') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

#clusterProfiler for GO
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
#these were already in ensembl so not needed but nice to have gene names
mart <- useDataset("hsapiens_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'grch37.ensembl.org'))

sig_genes_ensembl <- getBM(filters = "ensembl_gene_id", 
                           values = sig_genes$genes,
                           attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = mart)
input_genes <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1) %>% as.data.frame()
all_genes <- getBM(filters = "ensembl_gene_id", 
                   values = input_genes$.,
                   attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = mart)

ego <- enrichGO(gene=sig_genes_ensembl$ensembl_gene_id, universe=all_genes$ensembl_gene_id, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", pAdjustMethod = "BH", qvalueCutoff =0.05, readable=TRUE)
ego2 <- as.data.frame(ego)
cluster_summary <- summary(ego2)
dotplot(ego, showCategory=25)
dev.copy(png,'dotplot')
dev.off()
enrichMap(ego, n=25, vertex.label.font=3)
cnetplot(ego, categorySize="pvalue", showCategory = 5, vertex.label.font=3)

############################################################
####now repeat but with NPIH/PIH status as co-variate
#assign patient attributes
cl <- ifelse(NS_meta$Hydrocephalus=='NPIH',yes='NPIH',no=ifelse(NS_meta$Col_OR_PSU_Paeni=='Y',yes='PIH_Paeni',no='PIH_noPaeni'))
cl <- factor(cl)

design <- model.matrix(object = ~cl)
rownames(design) <- rownames(NS_meta)

#estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

#DE genes: Paeni vs Non
fit <- edgeR::glmFit(y, design)
lrt <- edgeR::glmLRT(fit, coef = 3)
edgeR::topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt)) 

#get DEG list
out <- topTags(lrt, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.01
sig <- out[keep,] %>% as.data.frame()

padj.cutoff <- 0.01
lfc.cutoff <- 1

threshold_NS <- sig$FDR < padj.cutoff
length(which(threshold_NS == TRUE)) #1905
sig$threshold <- threshold_NS 
NS_sig<- subset(sig, threshold == TRUE)
sig_genes <- NS_sig %>% as.data.frame()

#heatmap of 500 most DEG
library(gplots)
library(RColorBrewer)
logcpm <- cpm(y, prior.count=2, log=TRUE)
sort_sig <- sig[order(sig$logFC, decreasing = T),]
sig_names <-  sort_sig$Symbol[1:500]
highly_DE_lcpm <- logcpm[sig_names,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[as.factor(NS_meta$Col_OR_PSU_Paeni)]

heatmap.2(highly_DE_lcpm,col=rev(morecols(50)),trace="none", main="DE genes across samples",ColSideColors=col.cell,scale="row")

ggplot(sig) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  ggtitle('NS') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

#clusterProfiler for GO
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
#these were already in ensembl so not needed but nice to have gene names
mart <- useDataset("hsapiens_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'grch37.ensembl.org'))

sig_genes_ensembl <- getBM(filters = "ensembl_gene_id", 
                           values = sig_genes$genes,
                           attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = mart)
input_genes <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1) %>% as.data.frame()
all_genes <- getBM(filters = "ensembl_gene_id", 
                   values = input_genes$.,
                   attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = mart)

ego <- enrichGO(gene=sig_genes_ensembl$ensembl_gene_id, universe=all_genes$ensembl_gene_id, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", pAdjustMethod = "BH", qvalueCutoff =0.05, readable=TRUE)
ego2 <- as.data.frame(ego)
cluster_summary <- summary(ego2)
dotplot(ego, showCategory=25)
dev.copy(png,'PIH_Paeni_dotplot')
dev.off()
enrichMap(ego, n=25, vertex.label.font=3)
cnetplot(ego, categorySize="pvalue", showCategory = 5, vertex.label.font=3)

####now repeat but with blood CMV  status as co-variate
#assign patient attributes
cl <- ifelse(NS_meta$Hydrocephalus=='NPIH',yes='NPIH',no=ifelse(NS_meta$CMV_Blood=='N',yes='PIH_noBlood_CMVi',no='PIH_Blood_CMV'))
cl <- factor(cl)

design <- model.matrix(object = ~cl)
rownames(design) <- rownames(NS_meta)

#estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

#DE genes: Paeni vs Non
fit <- edgeR::glmFit(y, design)
lrt <- edgeR::glmLRT(fit, coef = 3)
edgeR::topTags(lrt)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt)) #

#get DEG list
out <- topTags(lrt, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.01
sig <- out[keep,] %>% as.data.frame()

padj.cutoff <- 0.01
lfc.cutoff <- 1

threshold_NS <- sig$FDR < padj.cutoff
length(which(threshold_NS == TRUE)) #4900
sig$threshold <- threshold_NS 
NS_sig<- subset(sig, threshold == TRUE)
sig_genes <- NS_sig %>% as.data.frame()

#heatmap of 500 most DEG
library(gplots)
library(RColorBrewer)
logcpm <- cpm(y, prior.count=2, log=TRUE)
sort_sig <- sig[order(sig$logFC, decreasing = T),]
sig_names <-  sort_sig$Symbol[1:500]
highly_DE_lcpm <- logcpm[sig_names,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[as.factor(NS_meta$CMV_Blood)]

heatmap.2(highly_DE_lcpm,col=rev(morecols(50)),trace="none", main="DE genes across samples",ColSideColors=col.cell,scale="row")

ggplot(sig) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  ggtitle('NS') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

#clusterProfiler for GO
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
#these were already in ensembl so not needed but nice to have gene names
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
sig_genes_ensembl <- getBM(filters = "ensembl_gene_id", 
                           values = sig_genes$genes,
                           attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = ensembl)
input_genes <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1) %>% as.data.frame()
all_genes <- getBM(filters = "ensembl_gene_id", 
                   values = input_genes$.,
                   attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = ensembl)

ego <- enrichGO(gene=sig_genes_ensembl$ensembl_gene_id, universe=all_genes$ensembl_gene_id, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", pAdjustMethod = "BH", qvalueCutoff =0.05, readable=TRUE)
ego2 <- as.data.frame(ego)
cluster_summary <- summary(ego2)
dotplot(ego, showCategory=25)
dev.copy(png,'PIH_Paeni_dotplot')
dev.off()
enrichMap(ego, n=25, vertex.label.font=3)
cnetplot(ego, categorySize="pvalue", showCategory = 5, vertex.label.font=3)

######################################now with CMV status alone
table(NS_meta$CMV_CSF)
design <- model.matrix(object = ~NS_meta$CMV_CSF)
rownames(design) <- rownames(NS_meta)

#estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

#DE genes: Paeni vs Non
fit <- edgeR::glmFit(y, design)
lrt <- edgeR::glmLRT(fit)
edgeR::topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt)) #

#get DEG list
out <- topTags(lrt, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05
sig <- out[keep,] %>% as.data.frame()

padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold_NS <- sig$FDR < padj.cutoff & abs(sig$logFC) > lfc.cutoff
length(which(threshold_NS == TRUE)) #15
sig$threshold <- threshold_NS 
NS_sig<- subset(sig, threshold == TRUE)
sig_genes <- NS_sig %>% as.data.frame()

#heatmap of 500 most DEG
library(gplots)
library(RColorBrewer)
logcpm <- cpm(y, prior.count=2, log=TRUE)
sort_sig <- sig[order(sig$logFC, decreasing = T),]
sig_names <-  sort_sig$Symbol[1:500]
highly_DE_lcpm <- logcpm[rownames(logcpm) %in% sig_names,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[as.factor(NS_meta$Col_OR_PSU_Paeni)]

heatmap.2(highly_DE_lcpm,col=rev(morecols(50)),trace="none", main="DE genes across samples",ColSideColors=col.cell,scale="row")

ggplot(sig) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  ggtitle('NS') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

#clusterProfiler for GO
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
#these were already in ensembl so not needed but nice to have gene names
mart <- useDataset("hsapiens_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'grch37.ensembl.org'))

sig_genes_ensembl <- getBM(filters = "ensembl_gene_id", 
                           values = sig_genes$genes,
                           attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = mart)
input_genes <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1) %>% as.data.frame()
all_genes <- getBM(filters = "ensembl_gene_id", 
                   values = input_genes$.,
                   attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = mart)

ego <- enrichGO(gene=sig_genes_ensembl$ensembl_gene_id, universe=all_genes$ensembl_gene_id, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", pAdjustMethod = "BH", qvalueCutoff =0.05, readable=TRUE)
ego2 <- as.data.frame(ego)
cluster_summary <- summary(ego2)
dotplot(ego, showCategory=25)
dev.copy(png,'dotplot_csfCMV')
dev.off()
enrichMap(ego, n=25, vertex.label.font=3)
cnetplot(ego, categorySize="pvalue", showCategory = 5, vertex.label.font=3)

### CMV positive versus negative in the subset of Paenibacillus spp. positive infants
NS_meta$Col_OR_PSU_Paeni <- gsub("Y", "Paeni", NS_meta$Col_OR_PSU_Paeni)
NS_meta$Col_OR_PSU_Paeni <- gsub("N", "Non_Paeni", NS_meta$Col_OR_PSU_Paeni) 
NS_meta$CMV_CSF <- gsub("Y", "CMV", NS_meta$CMV_CSF)
NS_meta$CMV_CSF <- gsub("N", "Non_CMV", NS_meta$CMV_CSF) 
Group <- factor(paste(NS_meta$Col_OR_PSU_Paeni,NS_meta$CMV_CSF,sep="."))
cbind(NS_meta,Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
rownames(design) <- rownames(NS_meta)
fit <- glmQLFit(y, design)

my.contrasts <- makeContrasts(Paeni.CMVvsNot = Paeni.CMV-Paeni.Non_CMV, CMV.PaenivsNot = Paeni.CMV-Non_Paeni.CMV,PaeniVsCMV = (Paeni.CMV-Paeni.Non_CMV)-(Non_Paeni.CMV-Non_Paeni.Non_CMV), levels=design)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Paeni.CMVvsNot"])

edgeR::topTags(qlf)
colnames(design)
o <- order(qlf$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(qlf)) #
plotMD(qlf)
abline(h=c(-1, 1), col="blue")

#get DEG list
out <- topTags(qlf, n=Inf, adjust.method="BH")
keep <- out$table$FDR <= 0.05
sig <- out[keep,] %>% as.data.frame()

padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold_NS <- sig$FDR < padj.cutoff & abs(sig$logFC) > lfc.cutoff
length(which(threshold_NS == TRUE)) #15
sig$threshold <- threshold_NS 
NS_sig<- subset(sig, threshold == TRUE)
sig_genes <- NS_sig %>% as.data.frame()

#heatmap of 500 most DEG
library(gplots)
library(RColorBrewer)
logcpm <- cpm(y, prior.count=2, log=TRUE)
sort_sig <- sig[order(sig$logFC, decreasing = T),]
sig_names <-  sort_sig$Symbol[1:500]
highly_DE_lcpm <- logcpm[rownames(logcpm) %in% sig_names,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange")[as.factor(NS_meta$Col_OR_PSU_Paeni)]

heatmap.2(highly_DE_lcpm,col=rev(morecols(50)),trace="none", main="DE genes across samples",ColSideColors=col.cell,scale="row")

ggplot(sig) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  xlim(c(-10,10)) +
  ylim(c(0,10)) +
  ggtitle('NS') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

#clusterProfiler for GO
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
#these were already in ensembl so not needed but nice to have gene names
mart <- useDataset("hsapiens_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'grch37.ensembl.org'))

sig_genes_ensembl <- getBM(filters = "ensembl_gene_id", 
                           values = sig_genes$genes,
                           attributes = c("ensembl_gene_id", "external_gene_name"),
                           mart = mart)
input_genes <- sapply(strsplit(rownames(NS_data1), "\\."), "[", 1) %>% as.data.frame()
all_genes <- getBM(filters = "ensembl_gene_id", 
                   values = input_genes$.,
                   attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = mart)

ego <- enrichGO(gene=sig_genes_ensembl$ensembl_gene_id, universe=all_genes$ensembl_gene_id, OrgDb=org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", pAdjustMethod = "BH", qvalueCutoff =0.05, readable=TRUE)
ego2 <- as.data.frame(ego)
cluster_summary <- summary(ego2)
dotplot(ego, showCategory=25)
dev.copy(png,'dotplot_CMVposNeg_PaeniPos')
dev.off()
enrichMap(ego, n=25, vertex.label.font=3)
cnetplot(ego, categorySize="pvalue", showCategory = 5, vertex.label.font=3)