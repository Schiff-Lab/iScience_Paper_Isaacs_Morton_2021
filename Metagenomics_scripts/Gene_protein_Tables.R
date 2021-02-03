library(dplyr)
library(ggplot2)
library(readxl)
library(phyloseq)
library(reshape2)

#Proteomics

MD_Col <- rep("text", 50)
MD<-read_excel("PIHNPHMetaData.xlsx", col_types=MD_Col )

Prot_Col <- rep("text", 143)
Prot<-read_excel("ProteomicsAlbertInitMeasurements.xlsx", col_types=Prot_Col )

Sig_Gene_col<- rep("text", 1)
Sig_Gen<-read_excel("ProteinsInBothRNAandProteinAlbertSara.xlsx", col_types=Sig_Gene_col )

##NormalizeToMedianAndNumericTransformation
rownames(Prot)<-Prot$Gene_ID
Protn<-Prot$Gene_ID
Prot<-Prot%>%dplyr::select(-Gene_ID)
rownames(Prot)<-Protn
ProtN <- mutate_all(Prot, function(x) as.numeric(as.character(x)))
rownames(ProtN)<-rownames(Prot)

#total1 = median(colSums(ProtN))
#standf1 = function(x, t=total1) (t * (x / sum(x)))
##Do apply on the row hence the 1 if 2 column but we want row then the function
#ProtNt = as.data.frame(t(apply(ProtN,1,standf1)))

##TakeOutProtSeenBothinMRNAandProt

#ProtNt$Symbol<-rownames(ProtNt)
ProtN$Symbol<-rownames(ProtN)
JProt<-left_join(Sig_Gen,ProtN,by="Symbol")
MJprot<-melt(JProt)
MJprot$IPM_ID<-MJprot$variable

MJprotD<-left_join(MJprot,MD,by="IPM_ID")
##DrawBoxPlot
MJprotDF <- mutate_all(MJprotD, function(x) as.character(x))

MJprotDF <- filter(MJprotD, Hydrocephalus == "NPIH" | Hydrocephalus == "PIH")

MJprotDF2<-filter(MJprotDF,Hydrocephalus=="NPIH")
NPIHAv<-unique(MJprotDF2%>%select(IPM_ID,Paeni_status))

ProtNcol<-ProtN$Symbol
ProtN<-ProtN%>%select(-Symbol)
tProtNt<-as.data.frame(t(ProtN))
colnames(tProtNt)<-ProtNcol
tProtNt$IPM_ID<-rownames(tProtNt)

LJNPIH<-left_join(NPIHAv,tProtNt,by="IPM_ID")
symbMJprotDF2<-unique(MJprotDF2$Symbol)

LJNPIH2<-LJNPIH[,symbMJprotDF2]
LJNPIH2_IPM<-LJNPIH$IPM_ID
rownames(LJNPIH2)<-LJNPIH2_IPM
LJNPIH2<-data.frame(t(LJNPIH2))
k<-rownames(LJNPIH2)
k<-data_frame(k)
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)
rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)


#PIH-Paeni_PIH
MJprotDF3<-filter(MJprotDF,Hydrocephalus=="PIH")
PIHAv<-MJprotDF3%>%select(Symbol,IPM_ID,Paeni_status)
MJprotDF4<-filter(MJprotDF,Paeni_status=="NonPaeni_PIH")
NPIHAv<-unique(MJprotDF4%>%select(IPM_ID,Paeni_status))

LJNPIH<-left_join(NPIHAv,tProtNt,by="IPM_ID")
symbMJprotDF2<-unique(MJprotDF2$Symbol)

LJNPIH2<-LJNPIH[,symbMJprotDF2]
LJNPIH2_IPM<-LJNPIH$IPM_ID
rownames(LJNPIH2)<-LJNPIH2_IPM
LJNPIH2<-data.frame(t(LJNPIH2))
k<-rownames(LJNPIH2)
k<-data_frame(k)
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)
rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)
LJNPIH3_f<-LJNPIH3%>%select(Gene_Mean,std)

#PIH-NonPaeni_PIH
MJprotDF3<-filter(MJprotDF,Hydrocephalus=="PIH")
PIHAv<-MJprotDF3%>%select(Symbol,IPM_ID,Paeni_status)
MJprotDF5<-filter(MJprotDF,Paeni_status=="Paeni_PIH")
NPIHAv<-unique(MJprotDF5%>%select(IPM_ID,Paeni_status))

LJNPIH<-left_join(NPIHAv,tProtNt,by="IPM_ID")
symbMJprotDF2<-unique(MJprotDF2$Symbol)

LJNPIH2<-LJNPIH[,symbMJprotDF2]
LJNPIH2_IPM<-LJNPIH$IPM_ID
rownames(LJNPIH2)<-LJNPIH2_IPM
LJNPIH2<-data.frame(t(LJNPIH2))
k<-rownames(LJNPIH2)
k<-data_frame(k)
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)
rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)
LJNPIH3_f<-LJNPIH3%>%select(Gene_Mean,std)

#RNASeq
load("Albert&Sara/Data for mercedeh/eSet.rda")
MD<-read.table("16s_paeni_positive_04222019.txt",sep="\t",header = TRUE)
MDRNA<-eSet@phenoData@data
MD$ID<-as.character(MD$ID)
MDF<-left_join(MD,MDRNA,by="ID")
MDFin<-MDF%>%select(IPM_ID.y,HYDROCEPHALUS,PAENI_POSITIVE)
names(MDFin)<-c("IPM_ID","Hydrocephalus","Paeni_status")

RNAname<-eSet@featureData@data
RNAnameS<-RNAname%>%select(gene_name)
RNAnameS$RNA_Ens<-rownames(RNAnameS)
names(RNAnameS)<-c("Symbol","RNA_Ens")
RNA<-eSet@assayData[["exprs"]]

Sig_Gene_col<- rep("text", 1)
Sig_Gen<-read_excel("ProteinsInBothRNAandProteinAlbertSara.xlsx", col_types=Sig_Gene_col )

RNAnameS2<-left_join(Sig_Gen,RNAnameS,by="Symbol")
RNA_Ens<-RNAnameS2$RNA_Ens
RNA<-RNA[RNA_Ens,]
RNA2<-as.data.frame(RNA)
RNA2$RNA_Ens<-RNA_Ens
RNAnameS3<-as.data.frame(RNAnameS2)
RNA3<-left_join(RNA2,RNAnameS3,by="RNA_Ens")
rownames(RNA3)<-RNA3$Symbol
colRNA3<-rownames(RNA3)
RNA4<-RNA3%>%select(-Symbol,-RNA_Ens)


#total1 = median(colSums(RNAN))
#standf1 = function(x, t=total1) (t * (x / sum(x)))
##Do apply on the row hence the 1 if 2 column but we want row then the function
#RNANt = as.data.frame(t(apply(RNAN,1,standf1)))

##TakeOutRNASeenBothinMRNAandRNA

#RNANt$Symbol<-rownames(RNANt)
RNA4t<-t(RNA4)
RNA4td<-as.data.frame(RNA4t)
RNA4td$IPM_ID<-rownames(RNA4td)
MJRNAD<-left_join(RNA4td,MDFin,by="IPM_ID")
##DrawBoxPlot
MJRNADF <- mutate_all(MJRNAD, function(x) as.character(x))

MJRNADF <- filter(MJRNAD, Hydrocephalus == "NPIH" | Hydrocephalus == "PIH")

MJRNADF2<-filter(MJRNADF,Hydrocephalus=="NPIH")
NPIHAv<-unique(MJRNADF2%>%select(IPM_ID))



LJNPIH<-left_join(NPIHAv,RNA4td,by="IPM_ID")
LJNPIH2<-LJNPIH

LJNPIH2<-data.frame(t(LJNPIH2))
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
LJNPIHA<-header.true(LJNPIH2)
LJNPIH2<-LJNPIHA
k<-rownames(LJNPIH2)
k<-data_frame(k)
#
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)

rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)
LJNPIH3_f<-LJNPIH3%>%select(Gene_Mean,std)


#NonPIH-Paeni_PIH
MJRNADF3<-filter(MJRNADF,Hydrocephalus=="PIH")
PIHAv<-MJRNADF3%>%select(IPM_ID,Paeni_status)
MJRNADF43<-filter(MJRNADF,Paeni_status=="NonPaeni_PIH")
NPIHAv<-unique(MJRNADF43%>%select(IPM_ID))



LJNPIH<-left_join(NPIHAv,RNA4td,by="IPM_ID")
LJNPIH2<-LJNPIH

LJNPIH2<-data.frame(t(LJNPIH2))
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
LJNPIHA<-header.true(LJNPIH2)
LJNPIH2<-LJNPIHA
k<-rownames(LJNPIH2)
k<-data_frame(k)
#
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)

rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)
LJNPIH3_f<-LJNPIH3%>%select(Gene_Mean,std)

#PIH-NonPaeni_PIH
MJRNADF3<-filter(MJRNADF,Hydrocephalus=="PIH")
PIHAv<-MJRNADF3%>%select(IPM_ID,Paeni_status)
MJRNADF42<-filter(MJRNADF,Paeni_status=="Paeni_PIH")
NPIHAv<-unique(MJRNADF42%>%select(IPM_ID))



LJNPIH<-left_join(NPIHAv,RNA4td,by="IPM_ID")
LJNPIH2<-LJNPIH

LJNPIH2<-data.frame(t(LJNPIH2))
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
LJNPIHA<-header.true(LJNPIH2)
LJNPIH2<-LJNPIHA
k<-rownames(LJNPIH2)
k<-data_frame(k)
#
LJNPIH3<-sapply(LJNPIH2,function(x){as.numeric(as.character(x))})
LJNPIH3<-as.data.frame(LJNPIH3)

rownames(LJNPIH3)<-k$k
LJNPIH3$Gene_Mean<-rowMeans(LJNPIH3)
##StandarErrorForTheMean
std <- function(x) sd(x)/sqrt(length(x))
#Apply for matrix bc I am indexing (matrix type obj)-1 is row
LJNPIH3$std<-apply(LJNPIH3,1,std)
LJNPIH3_f<-LJNPIH3%>%select(Gene_Mean,std)

