setwd("analysis")
library(dplyr)
library(shiny)

#update December 2022: new patient status
new_PCR <- read.csv("../new_data_Sept2022/PIH_CSF_summary_wconfirm_final_09072021.csv", stringsAsFactors = F)
metadata <- read.delim("../data/CURE_CONSHA DATABASE v2.2.0_08162022.txt", stringsAsFactors = F)
table(metadata$Hydrocephalus)
table(metadata$HYDRO_GP)
metadata$Hydrocephalus <- metadata$HYDRO_GP

#includes data from Joe 9/28/21:
#Qiime V1 and Qiime V2 outputs for the 442 samples that include controls
#Joe's code from the first paper in first100_code
#metadata with some errors repaired by Hai
#primer sequences from Christine - 27F: AGAGTTTGATCMTGGCTCAG, 336R: ACT GCT GCS YCC CGT AGG AGT CT
#updated for new alignments by Lijun June 2021

#open data, make MRexperiment file: https://www.bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
library(metagenomeSeq)
otu_tables <- read.csv("../data/filtered_oct2021/PIH_CSF_445Samples_Qiime2_dada2_442samples_OTU_Oct2021_MSKFiltered.txt", stringsAsFactors = F)
rownames(otu_tables) <- otu_tables$X
otu_tables <- otu_tables[,-1]
# #remove human contamination
# contam <- read.csv("../data/filtered_oct2021/pih_proteo_no_match_class_10122021_threeTaxRemoved.csv", stringsAsFactors = F)
# otu_tables <- otu_tables[!rownames(otu_tables) %in% contam$center,]

taxa = read.delim("../data/chiime2_20210615/PIH_CSF_445Samples_Qiime2_dada2_442samples_TAX_June2021.txt", stringsAsFactors=FALSE)
rownames(taxa) <- taxa$center
taxa <- taxa[,-1]
taxa <- taxa[rownames(taxa) %in% rownames(otu_tables),]

metadata_paeni_alll <- read.csv("../data/Hydrocephalus_CONSHA020321.csv", stringsAsFactors = F) #400, this was the one used for the original submission, just need the "paeni_all" status from this sheet
for(i in 1:nrow(metadata)){
  metadata$paeni_all[i] <- metadata_paeni_alll$paeni_all[metadata_paeni_alll$StudyID==metadata$SNUM[i]]
}
table(metadata$paeni_all)
metadata$StudyID 
metadata$SNUM 
#add in new pathogen and CT data
pathogen_status <- read.delim("../data/pathogenDB_10MAY2021.txt", stringsAsFactors=F)
pathogen_to_merge <- pathogen_status[pathogen_status$ID %in% metadata$StudyID,]
table(pathogen_to_merge$Pathogen) #CMV and Paeni
cmv <- pathogen_to_merge[pathogen_to_merge$Pathogen=="CMV",]
cmv_csf_pcr <- cmv %>% filter(SampleType=="CSF") %>% filter(Assay=="qpcr_CMV")
cmv_blood_pcr <- cmv %>% filter(SampleType=="Blood") %>% filter(Assay=="qpcr_CMV")
# paeni_test <- pathogen_to_merge[pathogen_to_merge$Pathogen=="Paenibacillus",] #replace this with new data new_PCR
# table(paeni_test$Assay)
# paeni_csf_pcr <- paeni_test %>% filter(Assay=="qpcr_genus")
# paeni_csf_pcr <- paeni_csf_pcr[!duplicated(paeni_csf_pcr$ID),]

for(i in 1:nrow(metadata)){
  metadata$cmv_csf_pcr[i] <- cmv_csf_pcr$Result[cmv_csf_pcr$ID==metadata$StudyID[i]]
}
for(i in 1:nrow(metadata)){
  metadata$cmv_csf_pcrQ[i] <- cmv_csf_pcr$Quant[cmv_csf_pcr$ID==metadata$StudyID[i]]
}
for(i in 1:nrow(metadata)){
  metadata$cmv_blood_pcr[i] <- cmv_blood_pcr$Result[cmv_blood_pcr$ID==metadata$StudyID[i]]
}
for(i in 1:nrow(metadata)){
  metadata$cmv_blood_pcrQ[i] <- cmv_blood_pcr$Quant[cmv_blood_pcr$ID==metadata$StudyID[i]]
}

ct_meta <- read.csv("../data/CTscores_300.csv", stringsAsFactors = F)
for(i in 1:nrow(metadata)){
  if (metadata$StudyID[i] %in% ct_meta$StudyID){
    metadata$CTscore_new[i] <- ct_meta$Total[ct_meta$StudyID==metadata$StudyID[i]]
  } else{
    metadata$CTscore_new[i] <- metadata$Score[i]    
  }
}
table(metadata$CTscore_new)
#use this full metadata for the Table 1
#do this filtering first to make Supplemental Table with included samples
#metadata <- metadata %>% filter(StudyID %in% included_samples)
table(metadata$HYDRO_GP)
days <- metadata[!is.na(metadata$Age_at_Infection_days),]
mean(days$Age_at_Infection_days)
mean(days$Age_at_Infection_days[days$HYDRO_GP=="PIH"])
mean(days$Age_at_Infection_days[days$HYDRO_GP=="NPIH"])
sd(days$Age_at_Infection_days)
sd(days$Age_at_Infection_days[days$HYDRO_GP=="PIH"])
sd(days$Age_at_Infection_days[days$HYDRO_GP=="NPIH"])
# wilcox.test(metadata$Age_at_Infection_days[metadata$HYDRO_GP=="PIH"], metadata$Age_at_Infection_days[metadata$HYDRO_GP=="NPIH"])
age <- metadata$Age_at_Infection_days[metadata$HYDRO_GP=="PIH"]
age <- age[!is.na(age)]
mean(age)
sd(age)
mean(metadata$Age_Sampling_Days)
mean(metadata$Age_Sampling_Days[metadata$HYDRO_GP=="PIH"])
mean(metadata$Age_Sampling_Days[metadata$HYDRO_GP=="NPIH"])
sd(metadata$Age_Sampling_Days)
sd(metadata$Age_Sampling_Days[metadata$HYDRO_GP=="PIH"])
sd(metadata$Age_Sampling_Days[metadata$HYDRO_GP=="NPIH"])
wilcox.test(metadata$Age_Sampling_Days[metadata$HYDRO_GP=="PIH"], metadata$Age_Sampling_Days[metadata$HYDRO_GP=="NPIH"])
interval <- metadata$DOI_to_Sampling_Days[metadata$HYDRO_GP=="PIH"]
interval <- interval[!is.na(interval)]
mean(interval)
sd(interval)
table(metadata$Gender)
table(metadata$Gender[metadata$HYDRO_GP=="PIH"])
table(metadata$Gender[metadata$HYDRO_GP=="NPIH"])
fisher.test(c(91,118), c(82,109))
mean(metadata$WBC_blood.1)
mean(metadata$WBC_blood.1[metadata$HYDRO_GP=="PIH"])
mean(metadata$WBC_blood.1[metadata$HYDRO_GP=="NPIH"])
sd(metadata$WBC_blood.1)
sd(metadata$WBC_blood.1[metadata$HYDRO_GP=="PIH"])
sd(metadata$WBC_blood.1[metadata$HYDRO_GP=="NPIH"])
wilcox.test(metadata$WBC_blood.1[metadata$HYDRO_GP=="PIH"], metadata$WBC_blood.1[metadata$HYDRO_GP=="NPIH"])
mean(metadata$CSF_cells)
mean(metadata$CSF_cells[metadata$HYDRO_GP=="PIH"])
mean(metadata$CSF_cells[metadata$HYDRO_GP=="NPIH"])
sd(metadata$CSF_cells)
sd(metadata$CSF_cells[metadata$HYDRO_GP=="PIH"])
sd(metadata$CSF_cells[metadata$HYDRO_GP=="NPIH"])
wilcox.test(metadata$CSF_cells[metadata$HYDRO_GP=="PIH"], metadata$CSF_cells[metadata$HYDRO_GP=="NPIH"])
mean(metadata$Hemaglobin.1)
mean(metadata$Hemaglobin.1[metadata$HYDRO_GP=="PIH"])
mean(metadata$Hemaglobin.1[metadata$HYDRO_GP=="NPIH"])
sd(metadata$Hemaglobin.1)
sd(metadata$Hemaglobin.1[metadata$HYDRO_GP=="PIH"])
sd(metadata$Hemaglobin.1[metadata$HYDRO_GP=="NPIH"])
wilcox.test(metadata$Hemaglobin.1[metadata$HYDRO_GP=="PIH"], metadata$Hemaglobin.1[metadata$HYDRO_GP=="NPIH"])
table(metadata$HIV.1)
table(metadata$HIV.1[metadata$HYDRO_GP=="PIH"])
table(metadata$HIV.1[metadata$HYDRO_GP=="NPIH"])
fisher.test(c(204,5), c(187,4))
table(metadata$cmv_blood_pcr)
table(metadata$cmv_blood_pcr[metadata$HYDRO_GP=="PIH"])
table(metadata$cmv_blood_pcr[metadata$HYDRO_GP=="NPIH"])
fisher.test(c(149,48), c(144,38))
library(fmsb)
oddsratio(52, 209-52, 37, 191-1-37)
table(metadata$cmv_csf_pcr)
table(metadata$cmv_csf_pcr[metadata$HYDRO_GP=="PIH"])
table(metadata$cmv_csf_pcr[metadata$HYDRO_GP=="NPIH"])
fisher.test(c(26,183), c(1,189))
oddsratio(26, 183, 1, 189)
table(metadata$CTscore_new)
x <- table(metadata$CTscore_new[metadata$HYDRO_GP=="PIH"])  %>% as.data.frame()
x$cohort <- "PIH"
y <- table(metadata$CTscore_new[metadata$HYDRO_GP=="NPIH"])  %>% as.data.frame()
y$cohort <- "NPIH"
library(MASS)
dat <- rbind(x,y)
dat$Var1 <- as.factor(dat$Var1)
dat$cohort <- as.factor(dat$cohort)
dat
dat[10,] <- c(3, 0, "NPIH")
pom <- polr(Var1 ~ cohort, data=dat)
summary(pom)
ctable <- coef(summary(pom))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ctable
#OR >=2 vs <2
oddsratio(129, 209-129, 10, 191-10)

# #and again for supplemental table of PIH based on paenibacillus status 
pmetadata <- metadata %>% filter(HYDRO_GP=="PIH")
pmetadata$infection_to_sample <- pmetadata$Age_at_Infection_days-pmetadata$Age_Sampling_Days
ppmetadata <- pmetadata %>% filter(paeni_all=="Y")
npmetadata <- pmetadata %>% filter(paeni_all=="N")
mean(ppmetadata$Age_at_Infection_days[!is.na(ppmetadata$Age_at_Infection_days)])
mean(npmetadata$Age_at_Infection_days[!is.na(npmetadata$Age_at_Infection_days)])
sd(ppmetadata$Age_at_Infection_days[!is.na(ppmetadata$Age_at_Infection_days)])
sd(npmetadata$Age_at_Infection_days[!is.na(npmetadata$Age_at_Infection_days)])
wilcox.test(ppmetadata$Age_at_Infection_days[!is.na(ppmetadata$Age_at_Infection_days)], npmetadata$Age_at_Infection_days[!is.na(npmetadata$Age_at_Infection_days)])
mean(ppmetadata$infection_to_sample[!is.na(ppmetadata$infection_to_sample)])
mean(npmetadata$infection_to_sample[!is.na(npmetadata$infection_to_sample)])
sd(ppmetadata$infection_to_sample[!is.na(ppmetadata$infection_to_sample)])
sd(npmetadata$infection_to_sample[!is.na(npmetadata$infection_to_sample)])
wilcox.test(ppmetadata$infection_to_sample[!is.na(ppmetadata$infection_to_sample)], npmetadata$infection_to_sample[!is.na(npmetadata$infection_to_sample)])
mean(ppmetadata$Age_Sampling_Days)
mean(npmetadata$Age_Sampling_Days)
sd(ppmetadata$Age_Sampling_Days)
sd(npmetadata$Age_Sampling_Days)
wilcox.test(ppmetadata$Age_Sampling_Days, npmetadata$Age_Sampling_Days)
table(ppmetadata$Gender)
table(npmetadata$Gender)
fisher.test(c(45,53), c(45,65))
mean(ppmetadata$WBC_blood.1)
mean(npmetadata$WBC_blood.1)
sd(ppmetadata$WBC_blood.1)
sd(npmetadata$WBC_blood.1)
wilcox.test(ppmetadata$WBC_blood.1, npmetadata$WBC_blood.1)
mean(ppmetadata$CSF_cells)
mean(npmetadata$CSF_cells)
sd(ppmetadata$CSF_cells)
sd(npmetadata$CSF_cells)
wilcox.test(ppmetadata$CSF_cells, npmetadata$CSF_cells)
mean(ppmetadata$Hemaglobin.1)
mean(npmetadata$Hemaglobin.1)
sd(ppmetadata$Hemaglobin.1)
sd(npmetadata$Hemaglobin.1)
wilcox.test(ppmetadata$Hemaglobin.1, npmetadata$Hemaglobin.1)
table(ppmetadata$HIV)
table(npmetadata$HIV)
fisher.test(c(95,1), c(99,4))
table(ppmetadata$cmv_blood_pcr)
table(npmetadata$cmv_blood_pcr)
fisher.test(c(77,21), c(79,31))
table(ppmetadata$cmv_csf_pcr)
table(npmetadata$cmv_csf_pcr)
fisher.test(c(13,85), c(13,97))
table(pmetadata$CTscore_new)
x <- table(ppmetadata$CTscore_new)  %>% as.data.frame()
x$cohort <- "Y"
y <- table(npmetadata$CTscore_new)  %>% as.data.frame()
y$cohort <- "N"
t.test(x$Freq, y$Freq)
library(MASS)
dat <- rbind(x,y)
dat$Var1 <- as.factor(dat$Var1)
dat$cohort <- as.factor(dat$cohort)
dat
pom <- polr(Var1 ~ cohort, data=dat)
summary(pom)
ctable <- coef(summary(pom))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ctable

#now subset to those with data for metagenomeseq 
ipm_id <- read.csv("../data/metadata_092432020.csv", stringsAsFactors = F)
ipm_id <-  ipm_id %>% filter(IPM_ID %in% colnames(otu_tables)) #442, same as OTU
ipm_id <-  ipm_id %>% filter(!Sample_Type %in% c("Standard", "NTC"))
ipm_id <-  ipm_id %>% filter(ID %in% metadata$SNUM)
ipm_id <-  ipm_id %>% filter(!duplicated(ipm_id$ID))
metadata <-  metadata %>% filter(SNUM %in% ipm_id$ID)
metadata <- merge(metadata, ipm_id, by.x = "SNUM", by.y = "ID")
for(i in 1:nrow(metadata)){
  metadata$IPM_ID[i] <- ipm_id$IPM_ID[ipm_id$ID==metadata$SNUM[i]]
}
metadata$days_since_infection <- as.numeric(metadata$Age_Sampling_Days) - as.numeric(metadata$Age_at_Infection_days)
summary(metadata$days_since_infection)

# for(i in 1:nrow(metadata)){
#   metadata$paeni_csf_pcr[i] <- paeni_csf_pcr$Result[paeni_csf_pcr$ID==metadata$StudyID[i]]
# }
# for(i in 1:nrow(metadata)){
#   metadata$paeni_csf_pcrQ[i] <- paeni_csf_pcr$Quant[cmv_csf_pcr$ID==metadata$StudyID[i]]
# }
for(i in 1:nrow(metadata)){
  if(metadata$StudyID[i] %in% new_PCR$StudyID){
    metadata$paeni_csf_pcr[i] <- new_PCR$qpcr_genus[new_PCR$StudyID==metadata$StudyID[i]]
  } else{
    metadata$paeni_csf_pcr[i] <- "NA"
  }
}
table(metadata$paeni_csf_pcr)
for(i in 1:nrow(metadata)){
  if(metadata$StudyID[i] %in% new_PCR$StudyID){
    metadata$paeni_csf_pcrQ[i] <- new_PCR$quant_genus[new_PCR$StudyID==metadata$StudyID[i]]
  } else{
    metadata$paeni_csf_pcrQ[i] <- "NA"
  }
}
table(metadata$paeni_csf_pcrQ)

#rest of Table 1
metadata_pih <- metadata %>% filter(Hydrocephalus=="PIH")
mean(metadata_pih$Age_at_Infection_days[!is.na(metadata_pih$Age_at_Infection_days)])
sd(metadata_pih$Age_at_Infection_days[!is.na(metadata_pih$Age_at_Infection_days)])
mean(metadata_pih$days_since_infection[!is.na(metadata_pih$days_since_infection)])
sd(metadata_pih$days_since_infection[!is.na(metadata_pih$days_since_infection)])

metadata_x <- metadata[metadata$IPM_ID %in% colnames(otu_tables),]
colnames(otu_tables)[!colnames(otu_tables) %in% metadata_x$IPM_ID] #46 missing from the new dataframe?

#first remove mito reads
table(taxa$Family)
mito <- taxa[grep("Mitochondria", taxa$Family),]
taxa_x <- taxa[!rownames(taxa) %in% rownames(mito),]
#now can filter out contaminants
contam <- taxa[grep("TRUE", taxa$contaminantP0.5),]
taxa_x <- taxa_x[!rownames(taxa_x) %in% rownames(contam),]
#and now remove unknown bacteria
unk <-  taxa[is.na(taxa$Phylum),]
taxa_x <- taxa_x[!rownames(taxa_x) %in% rownames(unk),]
#and remove chloroplasts
chloro <-  taxa[taxa$Order=="Chloroplast",]
taxa_x <- taxa_x[!rownames(taxa_x) %in% rownames(chloro),]
otu_tables_x <- otu_tables %>% filter(rownames(otu_tables) %in% rownames(taxa_x))
otu_tables_x <- otu_tables_x[,colnames(otu_tables_x) %in% metadata_x$IPM_ID]
metadata_x <- metadata_x[match(colnames(otu_tables_x), metadata_x$IPM_ID),] %>% as.data.frame()
rownames(metadata_x) <- metadata_x$IPM_ID
rownames(otu_tables_x)[!rownames(otu_tables_x) %in% rownames(taxa_x)] #all otu are in taxa
rownames(taxa_x)[!rownames(taxa_x) %in% rownames(otu_tables_x)] #all taxa are in otu
colnames(otu_tables_x)[!colnames(otu_tables_x) %in% rownames(metadata_x)] #all otu samples are in metadata
rownames(metadata_x)[!rownames(metadata_x) %in% colnames(otu_tables_x)] #all metadata samples are in otu

#all dataframes must be in the same order
metadata_x <- metadata_x[colnames(otu_tables_x),]
phenotypeData =AnnotatedDataFrame(metadata_x)
phenotypeData
taxa_x <- taxa_x[row.names(otu_tables_x),]
OTUdata =AnnotatedDataFrame(taxa_x)
OTUdata
head(rownames(pData(OTUdata)))
featureNames(phenotypeData)
featureNames(OTUdata)

obj =newMRexperiment(otu_tables_x,phenoData=phenotypeData,featureData=OTUdata)
obj
pData(obj)$Cohort = factor(pData(obj)$Hydrocephalus,level=c("NPIH","PIH"))
# save(obj, file="MRexp_16S_400_5dec22.rda")

otu_input <- otu_tables_x
otu_input$OTU <- rownames(otu_tables_x)
# write.table(taxa_x, file="taxa_input_5dec22.tsv", row.names = T, col.names = T, quote = F, sep = "\t")
# write.table(otu_tables_x, file="otu_input_5dec22.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
# write.table(metadata_x, file="pheno_input_5dec22.tsv", row.names = T, col.names = T, quote = F, sep = "\t")

obj_paeni<-which(fData(obj)$Genus=="Paenibacillus")
paeni_obj<-obj[obj_paeni,]
dim(fData(paeni_obj))

#normalize and filter
load("MRexp_16S_400_5dec22.rda")
head(libSize(obj))
normFactors(obj) <-rnorm(ncol(obj))
head(normFactors(obj))
head(normFactors(obj))
obj2 <- filterData(obj, present = 10, depth = 1)
p =cumNormStatFast(obj2) #some samples had zero features left, need to remove
obj2 =cumNorm(obj2, p = p)
table(obj2$Hydrocephalus)
# save(obj, file="normalized_MRexp_16S_400_5dec22.rda")

obj2_paeni<-which(fData(obj2)$Genus=="Paenibacillus")
paeni_obj2<-obj2[obj2_paeni,]
dim(fData(paeni_obj2)) #19 from 421

#table 1 of included
table(metadata_x$Hydrocephalus)
table(metadata_x$CTscore_new)
metadata_x_pih <- metadata_x %>% filter(Hydrocephalus=="PIH")
table(metadata_x_pih$CTscore_new)
metadata_x_npih <- metadata_x %>% filter(Hydrocephalus=="NPIH")
table(metadata_x_npih$CTscore_new)
t.test(as.numeric(metadata_x_pih$CTscore_new), as.numeric(metadata_x_npih$CTscore_new))

#plan to use microbiome explorer: https://github.com/zoecastillo/microbiomeExplorer
# BiocManager::install("zoecastillo/microbiomeExplorer")
#
# BiocManager::install(version='devel')
# BiocManager::install("microbiomeExplorer")

library(microbiomeExplorer)
load("normalized_MRexp_16S_400_5dec22.rda")
# runMicrobiomeExplorer() #can launch separate interface
#adapt example: https://github.com/zoecastillo/microbiomeExplorer/blob/master/vignettes/exploreMouseData.Rmd
meData <- filterMEData(obj2,minpresence = 1, minfeats = 10, minreads = 100)
# save(meData, file="filtered_normalized_MRexp_16S_400_5dec22.rda")

## Data QC: Before starting an analysis, it is recommend to review the results of the sequencing experiment and perform quality control. Some samples may have lower number of features than expected or an overall low number of reads. There are several ways to filter samples considered not useful for downstream analyses. Multiple QC plots can be generated, including those showing the number of unique features in each sample as a barplot or in a scatterplot against the number of reads. These can be colored by specific phenotypes stored in the pData slot of the MRexperiment. In addition, histograms show the overall distribution of feature and read frequencies.

makeQCPlot(meData, col_by = "Cohort",
           log = "both", #none as other option
           filter_feat = 101,
           filter_read = 511,
           allowWebGL = FALSE)
makeQCPlot(meData, col_by = "Seq_Batch",
           log = "both", #none as other option
           filter_feat = 101,
           filter_read = 511,
           allowWebGL = FALSE)
makeQCPlot(meData, col_by = "Extract_Batch",
           log = "both", #none as other option
           filter_feat = 101,
           filter_read = 511,
           allowWebGL = FALSE)
makeQCPlot(meData, col_by = "LP_Batch",
           log = "both", #none as other option
           filter_feat = 101,
           filter_read = 511,
           allowWebGL = FALSE)
plotlySampleBarplot(meData,
                    col_by = "Cohort")
plotlySampleBarplot(meData,
                    col_by = "Seq_Batch")
plotlySampleBarplot(meData,
                    col_by = "paeni_csf_pcr")

## Data Filtering and Subsetting: Within the application, three different sliders can be used to adjust quantitave restrictions on the data. The user can require a feature to be present in a minimum number of samples, and they can require a sample to have a minimum number of features or reads. Subsetting via a specific phenotype is also possible. This provides the option to exclude certain samples or limit the analysis to a subset of the data.

# meData <- filterMEData(obj,minpresence = 1, minfeats = 10, minreads = 100) #baseline 100, 500
table(meData$Cohort)
table(meData$Hydrocephalus)
table(meData$HYDRO_GP)
table(meData$CTscore_new)
included_samples <- meData$StudyID

## Normalization: Normalization allows for the user to account for library size differences before analysis. Certain app features are restricted if not done (e.g. percentage is unavailable if not normalized). Differential abundance testing also requires normalization, which we perform silently if the user does not choose to do so. The two available methods included in the package are based on either calculating proportions or by using cumulative sum scaling (CSS), Paulson, et al. Nat Meth 2013.

meData <- normalizeData(meData,norm_method = "Proportion")
table(meData$Cohort)

## Phenotable alteration: A new phenotable can easily be added to the MRexperiment. The application offers users similar options by allowing the extension or replacement of the pheno data via the load & filter section or by small in-app modifications of the table via the phenotype section. Here, the user can combine two columns to create a new phenotype based on a concatenation of the data. This new phenotype can then be used in the analysis. We provide the option to adjust data types, such as converting a declared numeric variable a factor, etc. This is typically not necessary, but may be useful or more appropriate, eg. in an example where cages are numeric. Finally, the user can simplify phenodata by selecting only necesssary columns.  Changes made need to be saved in order to update the underlying data structure for analysis.

new_pheno <- interaction(pData(meData)[,c("Hydrocephalus","cmv_csf_pcr")])
mutatedRows <- row.names(pData(meData))
mutatedData <- dplyr::mutate(pData(meData), "group_cmv" = new_pheno)
row.names(mutatedData) <- mutatedRows
meData <- addPhenoData(meData,mutatedData)

new_pheno <- interaction(pData(meData)[,c("paeni_csf_pcr", "Hydrocephalus")])
mutatedRows <- row.names(pData(meData))
mutatedData <- dplyr::mutate(pData(meData), "group_paeni" = new_pheno)
row.names(mutatedData) <- mutatedRows
meData <- addPhenoData(meData,mutatedData)

## Feature table alteration: The feature table provides an overview of the taxonomy annotation associated with each unique feature. Depending on the results, you might see many empty cells at specific taxonomy levels, eg. species or strain. These empty cells indicate that a feature could not be successfully assigned to a specific genus or strain. These can be left as are, marked as unknown or roll down the taxonomy. In this case, the most specific taxonomy annotation available for each feature is pushed down to more specific levels, e.g. "unknown_Firmicutes". The advantage of this method lies in the way analyses are done on the data. The user chooses a specific taxonomy level to analyze and all available feature counts are summarized down to unique features at this particular level. Two features with missing annotation cannot be distinguished, but an "unknown_Firmicutes" can be distinguished from an "unknown_Clostridiales". Any changes need to first be assigned and then saved in order to update the underlying data structure for analysis.

bufcolnames <- names(fData(meData))
df <- as.data.frame(t(apply(fData(meData),1, rollDownFeatures)))
names(df) <- bufcolnames
meData <- addFeatData(meData,df)

# Analysis: The analysis workflow within the application is split into five different sections: intra sample, inter sample, correlation, differential and longitudinal. Each section will be described below. All visualizations are implemented using the plotly R package which provides basic interactivity, including zooming or panning via its modebar. In addition, the user can export the plot in its current state (i.e. showing specific user interactions) as a svg file using the camera icon of the modebar.

## Aggregation: Before any analysis is possible, the user needs to aggregate the data down to a specific feature level. The available levels can be restricted via the code in global.R. Once this is completed, the analysis sections will be enable for use. Alternatively, the user can choose to add an analysis to a report by clicking the "Report" button which will recreate the visualizations as described in the reports section of this document.

# aggDat <- aggFeatures(meData, level = "Species")
aggDat <- aggFeatures(meData, level = "Genus")
#aggDat <- aggFeatures(meData, level = "Family")

## Intra-Sample Analysis: Intra-sample analysis contains functions focus on investigating the microbial composition within a sample or a group of samples. Different functions are available to visualize the relative abundance of top features, the abundance of a specific feature as well as the alpha diversity within the sample. Within the application, one common set of input elements is used to generate all visualization. 

### Relative Abundance: Relative abundance shows the most abundant feature in a barplot summarized by a user-defined variable across the x-asis. In addition, the user can choose to facet by phenotypes, adjust the number of features to show, switch between showing total numbers (Reads) and normalized value (if normalized), and modify the overall plot width. Clicking on a specific feature in the plot, automatically opens a feature abundance plot for this feature.

plotAbundance(aggDat,
              level = "genus",
              x_var = "Hydrocephalus",
              facet1 = NULL,
              facet2 = NULL,
              ind = 1:10,
              plotTitle = "Top 10 feature percentage at genus level",
              ylab = "Percentage")

plotAbundance(aggDat,
              level = "genus",
              x_var = "Hydrocephalus",
              facet1 = NULL,
              facet2 = NULL,
              ind = 1:20,
              plotTitle = "Top 20 feature percentage at genus level",
              ylab = "Percentage")

plotAbundance(aggDat,
              level = "genus",
              x_var = "group_paeni",
              facet1 = NULL,
              facet2 = NULL,
              ind = 1:10,
              plotTitle = "Top 10 feature percentage at genus level",
              ylab = "Percentage")

plotAbundance(aggDat,
              level = "genus",
              x_var = "StudyID",
              facet1 = "Hydrocephalus",
              facet2 = NULL,
              ind = 1:10,
              plotTitle = "Top 10 feature percentage at genus level",
              ylab = "Percentage")

aggDat2 <- aggDat[sort(aggDat["Paenibacillus"], decreasing = T)]
plotAbundance(aggDat,
              level = "genus",
              x_var = "StudyID",
              facet1 = "Hydrocephalus",
              facet2 = NULL,
              ind = 1:10,
              plotTitle = "Top 10 feature percentage at genus level",
              ylab = "Percentage")

# #now individual samples
cnts <- MRcounts(aggDat)
paeni_counts_new <- cnts["Paenibacillus",] %>% as.data.frame()
# paeni_counts_new$sum <- colSums(cnts)
# paeni_counts_new$proportion <- paeni_counts_new$./paeni_counts_new$sum
# paeni_counts_new$ID <- colnames(cnts)
# paeni_counts_new <- paeni_counts_new[order(paeni_counts_new$., decreasing=TRUE),]
# 
# saggDat = aggDat[,paeni_counts_new$ID]
# plotAbundance(saggDat,
#               level = "genus",
#               x_var = "StudyID",
#               facet1 = "Cohort",
#               facet2 = NULL,
#               ind = 1:10,
#               plotTitle = "Top 10 feature percentage at genus level",
#               ylab = "Percentage")

saggDat <- aggFeatures(meData, level = "Species")
scnts <- MRcounts(saggDat)
spaeni_counts_new <- scnts["popilliae",] %>% as.data.frame()

### Feature abundance: The feature abundance plot shows the individual abundance of a specific feature either as a boxplot or a categorical scatterplot depending on the x-axis variable chosen. The user can choose to employ a $\text{log}_2$ scale, define plot width and decide wether to show individual sample points or not. Feature plots can be opened by selecting a specific feature in the input section or by clicking on a feature in the relative abundance plot.

plotSingleFeature(aggDat,
                  x_var = "StudyID",
                  ind = 1:10,
                  plotTitle = "Percentage of Paenibacillus",
                  facet1 = "Cohort",
                  facet2 = NULL,
                  feature = "Paenibacillus",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Paenibacillus",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "Paenibacillus",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Alphaproteobacteria",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "unknown_Alphaproteobacteria",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Bacillus",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "Bacillus",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Staph",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "Staphylococcus",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Marinomonas",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "Marinomonas",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

plotSingleFeature(aggDat,
                  x_var = "Cohort",
                  ind = 1:10,
                  plotTitle = "Percentage of Psychromonas",
                  facet1 = NULL,
                  facet2 = NULL,
                  feature = "Psychromonas",
                  ylab = "Percentage",
                  log = TRUE,
                  showPoints = TRUE)

### Alpha Diversity: Alpha diversity is a measure of the complexity or diversity within a particular sample, eg. habitat or area. Alpha diversity is computed by functions in the vegan package and is visualized as a boxplot using the same input definitions by feature and relative abundance. The user can choose to color and thus split the boxes by a phenotype and set the overall plot width. Multiple diversity measures are offered with Shannon diversity provided as the default. We suggest users read up on the various measures and understand the differences in interpretation and nuances. Shannon diversity in particular measures how evenly the microbes are distributed in a sample and is defined by the following relationship where $p_i$, is the proportion of an individual feature:  $H = -\sum_{i=1}^T p_i \ln p_i$.  
plotAlpha(aggDat,
          level = "species",
          index = "shannon",
          x_var = "Hydrocephalus",
          facet1 = NULL,
          facet2 = NULL,
          col_by = "Hydrocephalus",
          plotTitle = "Shannon diversity index at species level")

plotAlpha(aggDat,
          level = "genus",
          index = "shannon",
          x_var = "Hydrocephalus",
          facet1 = NULL,
          facet2 = NULL,
          col_by = "Hydrocephalus",
          plotTitle = "Shannon diversity index at genus level")

#remove paeni positive
aggDat_paeniNeg <- filterByPheno(MRobj = aggDat,
              rm_phenovalues = list("paeni_csf_pcr" = c(1)))

plotAlpha(aggDat_paeniNeg,
          level = "genus",
          index = "shannon",
          x_var = "Hydrocephalus",
          facet1 = NULL,
          facet2 = NULL,
          col_by = "Hydrocephalus",
          plotTitle = "Shannon diversity index at genus level")

# plotAlpha(aggDat,
#           level = "genus",
#           index = "shannon",
#           x_var = "group_paeni",
#           facet1 = NULL,
#           facet2 = NULL,
#           col_by = "group_paeni",
#           plotTitle = "Shannon diversity index at genus level")

## Inter Sample Analysis: Inter-sample analyses focus on differences between samples or groups of samples via feature heatmaps and beta diversity calculations. Please note that these functions can be computationally intensive if there are many samples and a low aggregation level is chosen.

### Beta Diversity: Beta diversities are the measures of the complexity of communities between samples, as compared to within a sample (alpha diversity). Calculating beta diversity first requires the computation of a pairwise distance or similarity matrix. The user can select between different measures offered via the vegan package with Bray being the suggested default selection for microbiome analysis. We suggest users read up on the various measures and understand the differences in interpretation and nuances. Principal component analysis, a dimension reduction method, is subsequently performed on the chosen distance matrix and visualized in a scatter plot. The user has the option to choose the principal components to display, add coloring and confidence ellipses based on a phenotype, define the shape based on a phenotype and adjust both the point size as well as the overall plot width. PERMANOVA (permutational multivariate analysis of variance), from the vegan package is offered via the application (command line users will need to run this function independently and pass the results to the plotting function). Conceptually, a PERMANOVA analysis lets the user statistically determine if the centroids of a dissimilarity or distance matrix differ between groups of samples. Optionally, the user can select an phenotype as well as a strata variable with the results being shown, both within the visualization as well as in a table below it.

distMat <- computeDistMat(aggDat, "bray")
pcaVals <- calculatePCAs(distMat, 
                         c("PC1", "PC2"))
plotBeta(aggDat,
         dist_method = "bray",
         pcas = pcaVals,
         dim = c("PC1", "PC2"),
         col_by = "Hydrocephalus",
         shape_by = NULL,
         plotTitle = "Bray-Curtis diversity at genus level",
         pt_size = "6",
         # plotText = "R2: 0.478; Pr(>F): 0.002",
         confInterval = 0.95,
         allowWebGL = FALSE)

distMat <- computeDistMat(aggDat_paeniNeg, "bray")
pcaVals <- calculatePCAs(distMat, 
                         c("PC1", "PC2"))

plotBeta(aggDat_paeniNeg,
         dist_method = "bray",
         pcas = pcaVals,
         dim = c("PC1", "PC2"),
         col_by = "Hydrocephalus",
         shape_by = NULL,
         plotTitle = "Bray-Curtis diversity at genus level",
         pt_size = "6",
         # plotText = "R2: 0.478; Pr(>F): 0.002",
         confInterval = 0.95,
         allowWebGL = FALSE)

### Heatmap: The heatmap offers another view on differences and similarities between the samples in a dataset. The user can either choose specific features or show the top 50 features sorted either by variance, Fano factor or median absolute deviation (MAD). The visualization is done with heatmaply which in turns relies on plotly to render the heatmap. The same options to interact with the plot are thus available. Once rendered, the user can change the number of features to include, turn of log scale, and add annotation to both rows (phenotypes) and columns (higher taxonomy levels) of the heatmap. It is recommended to not use the heatmap functionality in datasets with many samples (5000+) as this can be quite slow to render.

plotHeatmap(aggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 25,
            col_by = c("Cohort"),
            row_by = "",
            plotTitle = "Top 25 features sorted by Variance at genus level")

plotHeatmap(aggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 50,
            col_by = c("CMV_Blood_or_csf"),
            row_by = "",
            plotTitle = "Top 50 features sorted by Variance at genus level")

plotHeatmap(aggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 50,
            col_by = c("Score"),
            row_by = "",
            plotTitle = "Top 50 features sorted by Variance at genus level")

#now make single heatmap with cohort, CMV and CT score across header
plotHeatmap(aggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 25,
            col_by = c("Cohort", "CMV_Blood_or_csf", "CTscore_new"),
            row_by = "",
            plotTitle = "Top 25 features sorted by Variance at genus level")

#again at species level
saggDat <- aggFeatures(meData, level = "Species")
plotHeatmap(saggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 25,
            col_by = c("Cohort", "CMV_Blood_or_csf", "CTscore_new"),
            row_by = "",
            plotTitle = "Top 25 features sorted by Variance at species level")

#improved colors
library("RColorBrewer")

heatpheno <- pData(saggDat)
new_row_names <- c("Paenibacillus popilliae", "Paenibacillus thiaminolyticus", "Staphlococcus epidermidis", "Rickettsiales spp", "Mycobacterium mucogenicum", "Streptococcus spp", "Streptococcus thermophilus", "Corynebacterium tuberculostearicum", "Staphylococcus capitis", "Mnemiopsis leidyi", "Tepidimonas spp", "Cutibacterium spp", "Acinetobacter bereziniae", "Corynebacterium kroppenstedtii", "Staphlococcus spp", "Staphlococcus aurea", "Lactobacillus iners", "Clavibacter michiganensis", "Acinetobacter spp", "Pelomonas spp", "Lactobacillus johnsonii", "Haemophilus parainfluenzae", "Staphylococcus hominis", "Massilia spp", "Alphaproteobacteria spp")
paeni_colors=gsub(0, "grey", heatpheno$paeni_csf_pcr)
paeni_colors=gsub(1, "darkred", paeni_colors)
cohort_colors=gsub("NPIH", "green", heatpheno$Hydrocephalus)
cohort_colors=gsub("PIH", "darkgreen", cohort_colors)
CT_colors=gsub("NA", "white", heatpheno$CTscore_new)
CT_colors=gsub(0, "darkblue", CT_colors)
CT_colors=gsub(1, "lightblue", CT_colors)
CT_colors=gsub(2, "lightgreen", CT_colors)
CT_colors=gsub(3, "orange", CT_colors)
CT_colors=gsub(4, "red", CT_colors)
clab=cbind(paeni_colors,cohort_colors,CT_colors) %>% as.data.frame()
colnames(clab)=c("Paeni_Status","Hydrocephalus_Status","CT_Score")
heatmapCols = colorRampPalette(brewer.pal(9, "Blues"))(5)
main_title="Top 25 species sorted by variance"

plotMRheatmap(saggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
          Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=paeni_colors,
          density.info="none", trace="none", main=main_title, labCol="", labRow=new_row_names, col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")
legend("bottomleft",legend=c("Paeni_Positive","Paeni_Negative","PIH","NPIH","N/A","0","1","2","3","4"), fill=c("darkred","grey","darkgreen","green","white","darkblue","lightblue","lightgreen","orange","red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

plotMRheatmap(saggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
              Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=cohort_colors,
              density.info="none", trace="none", main=main_title, labCol="", labRow=new_row_names, col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")

plotMRheatmap(saggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
              Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=CT_colors,
              density.info="none", trace="none", main=main_title, labCol="", labRow=new_row_names, col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")

#combine these in graphic software. now make genus version
heatpheno <- pData(aggDat)
new_row_names <- c("Paenibacillus popilliae", "Paenibacillus thiaminolyticus", "Staphlococcus epidermidis", "Rickettsiales spp", "Mycobacterium mucogenicum", "Streptococcus spp", "Streptococcus thermophilus", "Corynebacterium tuberculostearicum", "Staphylococcus capitis", "Mnemiopsis leidyi", "Tepidimonas spp", "Acinetobacter bereziniae", "Cutibacterium spp", "Corynebacterium kroppenstedtii", "Staphlococcus spp", "Staphlococcus aurea", "Lactobacillus iners", "Clavibacter michiganensis", "Acinetobacter spp", "Pelomonas spp", "Lactobacillus johnsonii", "Haemophilus parainfluenzae", "Staphylococcus hominis", "Massilia spp", "alphaproteobacteria spp")
paeni_colors=gsub(0, "grey", heatpheno$paeni_csf_pcr)
paeni_colors=gsub(1, "darkred", paeni_colors)
cohort_colors=gsub("NPIH", "green", heatpheno$Hydrocephalus)
cohort_colors=gsub("PIH", "darkgreen", cohort_colors)
CT_colors=gsub("NA", "white", heatpheno$CTscore_new)
CT_colors=gsub(0, "darkblue", CT_colors)
CT_colors=gsub(1, "lightblue", CT_colors)
CT_colors=gsub(2, "lightgreen", CT_colors)
CT_colors=gsub(3, "orange", CT_colors)
CT_colors=gsub(4, "red", CT_colors)
clab=cbind(paeni_colors,cohort_colors,CT_colors) %>% as.data.frame()
colnames(clab)=c("Paeni_Status","Hydrocephalus_Status","CT_Score")
heatmapCols = colorRampPalette(brewer.pal(9, "Blues"))(5)
main_title="Top 25 genera sorted by variance"

plotMRheatmap(aggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
              Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=paeni_colors,
              density.info="none", trace="none", main=main_title, labCol="", col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")
legend("bottomleft",legend=c("Paeni_Positive","Paeni_Negative","PIH","NPIH","N/A","0","1","2","3","4"), fill=c("darkred","grey","darkgreen","green","white","darkblue","lightblue","lightgreen","orange","red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

plotMRheatmap(aggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
              Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=cohort_colors,
              density.info="none", trace="none", main=main_title, labCol="", col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")

plotMRheatmap(aggDat, n=25, na.rm = TRUE, dendrogram="column", margins=c(6,12),
              Rowv=FALSE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE,ColSideColors=CT_colors,
              density.info="none", trace="none", main=main_title, labCol="", col=(heatmapCols), colCol = paeni_colors, key.title="Species Counts")

#and at OTU level
xaggDat <- aggFeatures(meData, level = "centerSeq")
plotHeatmap(xaggDat,
            features = NULL,
            log = TRUE,
            sort_by = "Variance",
            nfeat = 50,
            col_by = c("Cohort", "CMV_Blood_or_csf", "CTscore_new"),
            row_by = "",
            plotTitle = "Top 50 features sorted by Variance at species level")

## Correlation: Correlation allows the user to visualize the relationship between either two features or a feature and a numeric phenotype in a scatterplot enhanced with a linear regression statistic. Faceting and/or coloring by phenotypes is available in both correlation plots. The user is asked to choose between three different methods to aid in the evaluation of the association: Spearman (default), Pearson or Kendall.

cf <- corrFeature(aggDat,
                  feat1 = "Paenibacillus",
                  feat2 = "Staphylococcus",
                  log = TRUE,
                  facet1 = "Hydrocephalus",
                  facet2 = NULL,
                  method = "pearson",
                  plotTitle = "Pearson correlation of Paenibacillus vs Streptococcus split by Group",
                  col_by = "Hydrocephalus",
                  allowWebGL = FALSE) 

cf

## Differential abundance: Differential abundance (DA) analysis is focused on testing the null hypothesis that the mean or mean ranks between groups are the same for a specific feature. DA analysis can help detect changes in feature abundance across two or more different levels of a phenotype. Four different methods can be chosen via the application: DESeq2, Kruskal-Wallis, limma, or a zero-inflated log normal model. DESeq2 and limma are widely used methods for comparisons in microarray and RNA-sequencing data which can easily be adapted for microbiome data. Kruskal-Wallis is a non-parametric test for any differences in distribution between groups. The zero-inflated log normal model is implemented in the metagenomeSeq package to account for zero-inflation in microbiome data. Typically, DESeq2 would is used with small (<=25) sample sizes. The results will be displayed in an interactive table (DT) within the application and the user can open feature plots showing the specific levels by clicking on a row of interest.

diffResults <- runDiffTest(aggDat,
                           level = "Genus",
                           phenotype = "Hydrocephalus",
                           phenolevels = c("PIH", "NPIH"),
                           method = "limma")
diffResults$adj.P.Val
#estimate 0 p-values. Degrees if freedom is the smaller of the two group sizes -1
table(meData$Cohort)
diffResults$newPup <- pt(diffResults$t, 182, lower.tail = FALSE)
diffResults$new_adjPup <- p.adjust(diffResults$newPup, method="fdr")
diffResults$newPdown <- pt(diffResults$t, 182, lower.tail = TRUE)
diffResults$new_adjPdown <- p.adjust(diffResults$newPdown, method="fdr")
for(i in 1:nrow(diffResults)){
  diffResults$new_adjP[i] <- min(diffResults$new_adjPup[i], diffResults$new_adjPdown[i])
}
table(head(diffResults))
range(diffResults$log2FoldChange)
range(log2(diffResults$padj+0.001))

diffResults$genus <- gsub("unknown_", "", diffResults$Genus)
diffResults$BF_pval <- p.adjust(diffResults$new_adjP, method = "bonferroni")
diffResults$sig <- as.factor(diffResults$new_adjP<0.05)
diffResults$sig_color <- diffResults$sig
diffResults$sig_color <- gsub("TRUE", 'red', diffResults$sig_color)
diffResults$sig_color <- gsub("FALSE", 'black', diffResults$sig_color)
library(ggplot2)
library(ggrepel)
ggplot(diffResults) +
  geom_point(aes(x='PIH-NPIH', y=-log10(new_adjP)), color=diffResults$sig_color) +
  geom_text_repel(data=subset(diffResults, new_adjP<0.05), aes(label=genus, x='PIH-NPIH', y=-log10(new_adjP)), point.padding = NA, max.overlaps = 15) +
  # xlim(c(-11,31)) +
  # ylim(c(0,3)) +
  ggtitle('16S Genus') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 

# write.table(diffResults, file = "diff_abundance_5dec22.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#and now at the species level
diffsResults <- runDiffTest(saggDat,
                           level = "Species",
                           phenotype = "Hydrocephalus",
                           phenolevels = c("PIH", "NPIH"),
                           method = "limma")

table(head(diffsResults))
diffsResults$newPup <- pt(diffsResults$t, 182, lower.tail = FALSE)
diffsResults$new_adjPup <- p.adjust(diffsResults$newPup, method="fdr")
diffsResults$newPdown <- pt(diffsResults$t, 182, lower.tail = TRUE)
diffsResults$new_adjPdown <- p.adjust(diffsResults$newPdown, method="fdr")
for(i in 1:nrow(diffsResults)){
  diffsResults$new_adjP[i] <- min(diffsResults$new_adjPup[i], diffsResults$new_adjPdown[i])
}
range(diffsResults$log2FoldChange)
range(log2(diffsResults$padj+0.001))

diffsResults$sig <- as.factor(as.numeric(diffsResults$new_adjP<0.05))
diffsResults$sig_color <- diffsResults$sig
diffsResults$sig_color <- gsub("TRUE", 'red', diffsResults$sig_color)
diffsResults$sig_color <- gsub("FALSE", 'black', diffsResults$sig_color)

fdata <- fData(saggDat)
fdata <- fdata[!duplicated(fdata$Species),]
fdata <- fdata[!is.na(fdata$Species),]
fdata$Genus <- gsub(" 1", " ", fdata$Genus)
diffsResults[!diffsResults$Species %in% fdata$Species,]
for(i in 1:nrow(diffsResults)){
  if(diffsResults$Species[i] %in% fdata$Species){
    diffsResults$Genus[i] <- fdata$Genus[fdata$Species==diffsResults$Species[i]]
  } else(
    diffsResults$Genus[i] <- "NA"
  )
}
diffsResults$name <- paste(diffsResults$Genus, diffsResults$Species, sep = " ")

library(ggplot2)
library(ggrepel)
ggplot(diffsResults) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj+0.001)), color=diffsResults$sig_color) +
  geom_text_repel(data=subset(diffsResults, padj<0.05), aes(label=name, x=log2FoldChange, y=-log10(padj+0.001)), point.padding = NA, max.overlaps = 15) +
  # xlim(c(-11,31)) +
  # ylim(c(0,3)) +
  ggtitle('16S Species') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) 
# write.table(diffsResults, file = "diff_abundance_species_5dec22.txt", quote = F, sep = "\t", row.names = F, col.names = T)

## Longitudinal: Longitudinal analysis allows the user to generate feature plots with more control over the data shown within the plot. For a specific feature, the user can choose a phenotype and specific levels of that phenotype to show in the plot. The chosen order of the levels will be kept within the visualization which allows sorting by specific dates or tissues among other things. If desired and available, the user can choose a specific phenotype to summarize on which will then be connected by lines across the different levels. The resulting visualization is interactive and the user can then select and color specific IDs within the plot.

plotLongFeature(aggDat,
                x_var = "Age_at_Infection_days",
                id_var = "StudyID",
                plotTitle = "Abundance of Paenibacillus",
                feature = "Paenibacillus",
                ylab = "Percentage",
                log = FALSE,
                x_levels = c(0,1,2,3,4,5,6,7,10,14,21,30))


plotLongFeature(aggDat,
                x_var = "Season_Birth",
                id_var = "StudyID",
                plotTitle = "Abundance of Paenibacillus",
                feature = " Paenibacillus",
                ylab = "Reads",
                log = TRUE,
                x_levels = c("Rainy","Dry"))

cnts <- MRcounts(aggDat)
cnts["Paenibacillus",]
pni <- cnts["Paenibacillus",] %>% as.data.frame()
for(i in 1:nrow(pni)){
  pni$Age_at_Infection_days[i] <- metadata_x$Age_at_Infection_days[metadata_x$IPM_ID==rownames(pni)[i]]
}
for(i in 1:nrow(pni)){
  pni$Age_Sampling_Days[i] <- metadata_x$Age_Sampling_Days[metadata_x$IPM_ID==rownames(pni)[i]]
}
pni$log <- log(pni$.,2)
pni_age <- pni[!is.na(pni$Age_Sampling_Days),]
cor.test(pni_age$., pni_age$Age_Sampling_Days)
pni_age <- pni[!is.na(pni$Age_at_Infection_days),]
cor.test(pni_age$., pni_age$Age_at_Infection_days)
pni_age$days_since_infection <- pni_age$Age_Sampling_Days-pni_age$Age_at_Infection_days
cor.test(pni_age$., pni_age$days_since_infection)

# order_list <- as.numeric(metadata_x$days_since_infection)
# order_list <- order_list[!order_list ==-2]
# order_list <- order_list[!is.na(order_list)]
# order_list <- order_list[order(order_list)]
# order_list <- unique(order_list)
# plotLongFeature(aggDat,
#                 x_var = "days_since_infection",
#                 id_var = "StudyID",
#                 plotTitle = "Abundance of Paenibacillus",
#                 feature = "Paenibacillus",
#                 ylab = "Percentage",
#                 log = FALSE,
#                 x_levels = order_list)
library(ggplot2)
ggplot(pni_age, aes(days_since_infection, log)) +geom_point() +xlab("Days Since Symptom Onset") + ylab("log2 Paenibacillus Counts")

season <- pData(obj)["Season_Birth"]


#next:batch evaluation (remake old figures), sample characterization (stacked vars, heat map, alpha/beta, differential abundance, compare 1st 100 paeni with >50 counts here + ROC curve), mixomics RNA+16S+lab values (binary/continuous only; 75 training with 5-folds, 25 test, 300 validation balanced on paeni, csf cell count, age), unknown bacteria (blast sequences), clinical associaton (CT score, CSF cell count)

#now rolldown the unmapped genus: unknown_Bacteria, unknown_Enterobacteriaceae, unknown_OD1, unknown_OPB56, unknown_Rickettsiales
featrowuB <- fData(aggDat)[" unknown_k__Bacteria",]
rollDownFeatures(featrowuB)

featrowuE <- fData(obj)["unknown_p__OD1",]
rollDownFeatures(featrowuE)

#now make abundance percentage graph in order of Paeni counts for top 10
cnts <- MRcounts(aggDat)
top10 <- cnts[1:10,] %>% as.data.frame()
other <- cnts[11:nrow(cnts),]
other_sum <- colSums(other) %>% as.data.frame()
top10_full <- rbind(top10, t(other_sum))
rownames(top10_full)[11] <- "Other"
total_sum <- colSums(top10_full) %>% as.data.frame()

top10_proportion <- top10_full
for(i in 1:ncol(top10_proportion)){
  top10_proportion[,i] <- top10_proportion[,i]/total_sum[i,1]
}
top10_proportion <- top10_proportion[,!colSums(top10_proportion)<1]
paeni_proportion <- t(top10_proportion["Paenibacillus",]) %>% as.data.frame()
paeni_proportion$ID <- rownames(paeni_proportion)
paeni_proportion <- paeni_proportion[order(paeni_proportion$Paenibacillus, decreasing = T),]

top10_proportion <- top10_proportion[,rownames(paeni_proportion)]
top10_proportion <- top10_proportion[c("Paenibacillus", "Cutibacterium", rownames(top10_proportion)[3:11]),]

#check that is it 100% so that sum below = n samples
check_sum <- colSums(top10_proportion) 
sum(check_sum)#yes, works
  
bar_data <- unlist(top10_proportion) %>% as.data.frame()
colnames(bar_data) <- "Proportion"
library(stringr)
bar_data$Sample <- str_sub(rownames(bar_data), 1, str_length(rownames(bar_data))-1)
included_genera <- rownames(top10_proportion)
included_genera <- gsub("Corynebacterium 1", "Corynebacterium", included_genera)
included_genera <- gsub("Prevotella 7", "Prevotella", included_genera)
bar_data$Genus <- c(rep(included_genera, 376))


for(i in 1:nrow(bar_data)){
  if(bar_data$Sample[i] %in% metadata_x_npih$IPM_ID){
    bar_data$cohort[i] <- "NPIH"
  } else if(bar_data$Sample[i] %in% metadata_x_pih$IPM_ID){
    bar_data$cohort[i] <- "PIH"
  } else {
    bar_data$cohort[i] <- "NA"    
  }
}
 
for(i in 1:nrow(bar_data)){
  if(bar_data$cohort[i] == "NA"){
    bar_data$Sample[i] <- str_sub(bar_data$Sample[i], 1, str_length(bar_data$Sample[i])-1)
  }
}

for(i in 1:nrow(bar_data)){
  if(bar_data$Sample[i] %in% metadata_x_npih$IPM_ID){
    bar_data$cohort[i] <- "NPIH"
  } else if(bar_data$Sample[i] %in% metadata_x_pih$IPM_ID){
    bar_data$cohort[i] <- "PIH"
  } else {
    bar_data$cohort[i] <- "NA"    
  }
}

table(bar_data$cohort)
bar_data_pih <- bar_data[bar_data$cohort=="PIH",]
bar_data_npih <- bar_data[bar_data$cohort=="NPIH",]

library(ggplot2)
ggplot(bar_data_pih, aes(fill=factor(Genus, levels=included_genera), y=Proportion, x=factor(Sample, level=unique(c(Sample))))) +
  scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus, PIH Cohort", x = "Patient", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")
#linkage cases: 3323 (NS)—>2341 (PIH) = IPM3164 (column 18); 3334—>2346 (PIH) = IPM3169 (column 32)

ggplot(bar_data_npih, aes(fill=factor(Genus, levels=included_genera), y=Proportion, x=factor(Sample, level=unique(c(Sample))))) + 
  # scale_fill_manual(values=c("firebrick1", "khaki1", "lemonchiffon3", "lightblue1", "mediumaquamarine", "mistyrose", "darkseagreen1", "lightsteelblue3", "ghostwhite")) +
  scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus, NPIH Cohort", x = "Patient", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")

#now sum by Paeni/non-Paeni and remake bar graph with same colors
bar_data_pih_paeni <- bar_data_pih %>% filter(Sample %in% metadata_pih$IPM_ID[metadata_pih$paeni_csf_pcr==1])
length(unique(bar_data_pih_paeni$Sample))
bar_data_pih_NOpaeni <- bar_data_pih %>% filter(Sample %in% metadata_pih$IPM_ID[metadata_pih$paeni_csf_pcr==0])
length(unique(bar_data_pih_NOpaeni$Sample))

bar_data_npih_paeni <- bar_data_npih %>% filter(Sample %in% metadata_x_npih$IPM_ID[metadata_x_npih$paeni_csf_pcr==1])
length(unique(bar_data_npih_paeni$Sample))
bar_data_npih_NOpaeni <- bar_data_npih %>% filter(Sample %in% metadata_x_npih$IPM_ID[metadata_x_npih$paeni_csf_pcr==0])
length(unique(bar_data_npih_NOpaeni$Sample))

bar_data_pih_paeni_ave <- aggregate(bar_data_pih_paeni$Proportion, list(bar_data_pih_paeni$Genus), mean)
bar_data_pih_paeni_ave$cohort <- "PIH_Paeni"
bar_data_pih_paeni_ave$hydrocephalus <- "PIH"

bar_data_pih_NOpaeni_ave <- aggregate(bar_data_pih_NOpaeni$Proportion, list(bar_data_pih_NOpaeni$Genus), mean)
bar_data_pih_NOpaeni_ave$cohort <- "PIH_NOPaeni"
bar_data_pih_NOpaeni_ave$hydrocephalus <- "PIH"

bar_data_npih_paeni_ave <- aggregate(bar_data_npih_paeni$Proportion, list(bar_data_npih_paeni$Genus), mean)
bar_data_npih_paeni_ave$cohort <- "NPIH_Paeni"
bar_data_npih_paeni_ave$hydrocephalus <- "NPIH"

bar_data_npih_NOpaeni_ave <- aggregate(bar_data_npih_NOpaeni$Proportion, list(bar_data_npih_NOpaeni$Genus), mean)
bar_data_npih_NOpaeni_ave$cohort <- "NPIH_NOPaeni"
bar_data_npih_NOpaeni_ave$hydrocephalus <- "NPIH"

bar_data_ave <- rbind(bar_data_pih_paeni_ave, bar_data_pih_NOpaeni_ave, bar_data_npih_paeni_ave, bar_data_npih_NOpaeni_ave)

ggplot(bar_data_ave, aes(fill=factor(Group.1, levels=included_genera), y=x, x=factor(cohort, level=c("NPIH_NOPaeni", "PIH_NOPaeni", "NPIH_Paeni", "PIH_Paeni")))) + scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus", x = "Patient", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")
length(unique(bar_data_npih_NOpaeni$Sample))
length(unique(bar_data_pih_NOpaeni$Sample))
length(unique(bar_data_npih_paeni$Sample))
length(unique(bar_data_pih_paeni$Sample))

ggplot(bar_data_ave, aes(fill=factor(Group.1, levels=included_genera), y=x, x=factor(hydrocephalus, level=c("NPIH", "PIH")))) + scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus", x = "Cohort", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")

#and overall cohorts
ggplot(bar_data_ave, aes(fill=factor(Group.1, levels=included_genera), y=x, x=factor(cohort, level=c("NPIH_NOPaeni", "PIH_NOPaeni", "NPIH_Paeni", "PIH_Paeni")))) + scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus", x = "Patient", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")
#finaly by CT score <2 or >=2
bar_data_pih_highCT <- bar_data_pih %>% filter(Sample %in% metadata_pih$IPM_ID[metadata_pih$CTscore_new>1])
length(unique(bar_data_pih_highCT$Sample))
bar_data_pih_lowCT <- bar_data_pih %>% filter(Sample %in% metadata_pih$IPM_ID[metadata_pih$CTscore_new<2])
length(unique(bar_data_pih_lowCT$Sample))

bar_data_npih_highCT <- bar_data_npih %>% filter(Sample %in% metadata_x_npih$IPM_ID[metadata_x_npih$CTscore_new>1])
length(unique(bar_data_npih_highCT$Sample))
bar_data_npih_lowCT <- bar_data_npih %>% filter(Sample %in% metadata_x_npih$IPM_ID[metadata_x_npih$CTscore_new<2])
length(unique(bar_data_npih_lowCT$Sample))

bar_data_pih_highCT_ave <- aggregate(bar_data_pih_highCT$Proportion, list(bar_data_pih_highCT$Genus), mean)
bar_data_pih_highCT_ave$cohort <- "PIH_highCT"

bar_data_pih_lowCT_ave <- aggregate(bar_data_pih_lowCT$Proportion, list(bar_data_pih_lowCT$Genus), mean)
bar_data_pih_lowCT_ave$cohort <- "PIH_lowCT"

bar_data_npih_highCT_ave <- aggregate(bar_data_npih_highCT$Proportion, list(bar_data_npih_highCT$Genus), mean)
bar_data_npih_highCT_ave$cohort <- "NPIH_highCT"

bar_data_npih_lowCT_ave <- aggregate(bar_data_npih_lowCT$Proportion, list(bar_data_npih_lowCT$Genus), mean)
bar_data_npih_lowCT_ave$cohort <- "NPIH_lowCT"

bar_data_ave <- rbind(bar_data_pih_highCT_ave, bar_data_pih_lowCT_ave, bar_data_npih_highCT_ave, bar_data_npih_lowCT_ave)

ggplot(bar_data_ave, aes(fill=factor(Group.1, levels=included_genera), y=x, x=factor(cohort, level=c("NPIH_lowCT", "PIH_lowCT", "NPIH_highCT", "PIH_highCT")))) + scale_fill_manual(values=c("#FF3300", "#FFCC66", "#FFFF33", "#99FF33", "#339933", "#33CCCC", "#006666", "#3333CC", "#000033", "#CCCCCC", "snow")) +
  labs(title = "16S Proportion by Genus", x = "Patient", y = "16S Proportion", fill = "Genus") +
  geom_bar(position="stack", stat="identity")

####Paeni counts, use 307.5 as cutoff for positive. can't use this because is not normalized
cnts <- MRcounts(aggDat)
paeni_counts <- cnts[rownames(cnts) =="Paenibacillus",] %>% as.data.frame()
paeni <- taxa_x[grep("Paenibacillus", taxa_x$Genus),]
paeni_otu <- otu_tables_x %>% filter(rownames(otu_tables_x) %in% rownames(paeni))
paeni_counts <- colSums(paeni_otu)  %>% as.data.frame()
for(i in 1:nrow(paeni_counts)){
  paeni_counts$cohort[i] <- metadata_x$Hydrocephalus[metadata_x$IPM_ID==rownames(paeni_counts)[i]]
}
paeni_counts$id <- paste(paeni_counts$cohort, rownames(paeni_counts), sep = ":")
paeni_counts <- paeni_counts[order(paeni_counts$id),]

barplot(paeni_counts$., names.arg = paeni_counts$id, cex.axis = 0.8, cex.names = 0.5, xlab = "", ylab = "Paenibacillus Counts", las=2, col = as.factor(paeni_counts$cohort))
paeni_counts$positive <- paeni_counts$.>54.5
table(paeni_counts$positive)
paeni_counts$IPM_ID <- rownames(paeni_counts)
# save(paeni_counts, file = "paeni_status_new_otu_5dec22.rda")

og_cohort <- read.csv("../first100_code/16s_paeni_positive_04222019.csv", stringsAsFactors = F)
colnames(paeni_counts) <- c("new_paeni_counts", "cohort", "id", "new_paeni_status", "IPM_ID")
og_cohort <- merge(og_cohort, paeni_counts, by="IPM_ID", all.x = T, all.y = F)
og_cohort_pih <- og_cohort %>% filter(HYDROCEPHALUS=="PIH")
og_cohort_npih <- og_cohort %>% filter(HYDROCEPHALUS=="NPIH")
table(og_cohort_pih$PAENI_POSITIVE)
table(og_cohort_pih$PAENI_POSITIVE[og_cohort_pih$new_paeni_status==T])
table(og_cohort_pih$PAENI_POSITIVE[og_cohort_pih$new_paeni_status==F])
table(og_cohort_npih$PAENI_POSITIVE[og_cohort_npih$new_paeni_status==T])
table(og_cohort_npih$PAENI_POSITIVE[og_cohort_npih$new_paeni_status==F])

#now scatter plot to compare counts
paeni_100 <- readRDS("../first100_code/paeniData_filtered_paeni.rds")
aggDat_paeni_100 <- aggFeatures(paeni_100, level = "Genus")
cnts_100 <- MRcounts(aggDat_paeni_100)
paeni_counts_100 <- cnts_100[rownames(cnts_100) =="g__Paenibacillus",] %>% as.data.frame()
colnames(paeni_counts_100) <- "OG_Paeni_Counts"
paeni_counts_100$IPM_ID <- rownames(paeni_counts_100)
og_cohort <- merge(og_cohort, paeni_counts_100, by="IPM_ID", all.x = T, all.y = F)
ggplot(og_cohort, aes(OG_Paeni_Counts, new_paeni_counts, color=HYDROCEPHALUS)) + geom_point() +xlab("Original 16S Paenibacillus Counts") +ylab("New 16S Paenibacillus Counts") + xlim(0,5E+5) + ylim(0,5E+5) + geom_abline(intercept = -482.2243, slope = 0.5846)
lm(og_cohort$new_paeni_counts~og_cohort$OG_Paeni_Counts)

