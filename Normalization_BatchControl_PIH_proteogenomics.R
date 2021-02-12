# Set working directory if necessary
setwd("/Users/albertisaacs/Dropbox/Research Projects/CONSHA/Manuscripts/PIH Proteogenomics/Data Analysis")


#Install programs
if (!requireNamespace("BiocManager", quietly = TRUE))    install.packages("BiocManager")
BiocManager::install("yarn")

#load libraries
library(Biobase)
library(readxl)	
library(dplyr)
library(sva)
library(calibrate)
library(rafalib)
library(yarn)

#------------------------------BIOCONUDCTOR WORKING ENVIRONMENT----------------------------------------------------------------------------

# Load files for Bioconductor environment
cnts <- read_excel("measurements.xlsx") %>% data.frame
rownames(cnts) <- cnts[,1]
cnts <- cnts[,-1]
fd   <- read_excel("fdata.xlsx") %>% data.frame
pd   <- read_excel("pdata.xlsx") %>% data.frame

# Annotate the phenotype and the protein features (the Biobase library is required for this step)
rownames(pd) <- pd[,1]
rownames(fd) <- fd[,1]

# Set up the Bioconductor ExpressionSet
obj <- ExpressionSet(assayData=as.matrix(cnts),
    phenoData=AnnotatedDataFrame(pd),
    featureData=AnnotatedDataFrame(fd))

#------------------------------YARN----------------------------------------------------------------------------

# Plot Raw Denisty
yarn::plotDensity(obj,groups="Batch",legendPos="topright",xlab="Abundance",ylab="Density")

# Assess for batch effects with CMDS Plot CMDS
x = yarn::plotCMDS(obj,pch=21,bg=factor(pData(obj)$Batch))
textxy(x[,1],x[,2],pData(obj)$Batch)

# Filter low genes and observe effect
obj2 <- filterLowGenes(obj,minSamples = 18)
yarn::plotDensity(obj2,groups="Batch",legendPos="topright",xlab="Abundance",ylab="Density")

# Normalize the data and observe effect
obj3 <- yarn::normalizeTissueAware(obj2,normalizationMethod="quantile",groups=NULL)
nmat <- assayData(obj3)$normalizedMatrix
yarn::plotDensity(nmat,groups=factor(pData(obj3)$Batch),legendPos="topright")

# Adjust batch effects
combat_dat <- ComBat(dat=nmat, batch=factor(pData(obj)$Batch),
    mod=model.matrix(~Hydrocephalus,data=pData(obj3)),
    par.prior=TRUE, prior.plots=TRUE)

#Visualize normalized data with batch correction
bigpar(1)
plotDensity(combat_dat,normalized=TRUE,groups=factor(pData(obj3)$Batch),
    legendPos="topright",xlab="Abundance",ylab="Density")
    cbmat<-log2(combat_dat+1)

# Save results
storageMode(obj3) <- "environment"
assayData(obj3)$batchnormalized <- nmat
storageMode(obj3) <- "lockedEnvironment"

keepGenes = -which(is.na(rowSums(cbmat)))
keepSamples = -which(ifelse(obj3$Hydrocephalus %in% c("PIH","NPIH"),1,0)==0)
reducedMat = cbmat[keepGenes,keepSamples]

fdr = fData(obj3)[keepGenes,]
pdr = pData(obj3)[keepSamples,]


objr <- ExpressionSet(assayData=as.matrix(reducedMat),
                     phenoData=AnnotatedDataFrame(pdr),
                     featureData=AnnotatedDataFrame(fdr))

#Plot MDS to view sample clustering
bigpar(brewer.name = "Set1")
pdf("~/Dropbox/Research Projects/CONSHA/Manuscripts/PIH Proteogenomics/Data Analysis/MDS_plot.pdf")
    plotCMDS(objr,pch=21,bg=factor(pData(objr)$Paeni_status),cex=2)
legend("topright",legend=levels(factor(objr$Paeni_status)),fill=1:3, box.col=NA)
dev.off()

#Save RDS output
saveRDS(objr,file="pih_proteo_batch_normalized_reduced.rds")
saveRDS(obj3,file="pih_proteo_batch_normalized.rds")

x <- assayData(obj3)$batchnormalized
