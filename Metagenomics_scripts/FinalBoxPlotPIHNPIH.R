library(dplyr)
library(ggplot2)
library(readxl)
library(phyloseq)
library(reshape2)

MD_Col <- rep("text", 50)
MD<-read_excel("PIHNPHMetaData.xlsx", col_types=MD_Col )

Prot_Col <- rep("text", 143)
Prot<-read_excel("ProteomicsAlbertInitMeasurements.xlsx", col_types=Prot_Col )


Sig_Gene_col<- rep("text", 1)
Sig_Gen<-read_excel("ProteinsInBothRNAandProteinAlbertSara.xlsx", col_types=Sig_Gene_col )

##NormalizeToMedianAndNumericTransformation
rownames(Prot)<-Prot$Gene_ID
Prot<-Prot%>%dplyr::select(-Gene_ID)

ProtN <- mutate_all(Prot, function(x) as.numeric(as.character(x)))
rownames(ProtN)<-rownames(Prot)

total1 = median(colSums(ProtN))
standf1 = function(x, t=total1) (t * (x / sum(x)))
##Do apply on the row hence the 1 if 2 column but we want row then the function
ProtNt = as.data.frame(t(apply(ProtN,1,standf1)))

##TakeOutProtSeenBothinMRNAandProt

ProtNt$Symbol<-rownames(ProtNt)

JProt<-left_join(Sig_Gen,ProtNt,by="Symbol")
MJprot<-melt(JProt)
MJprot$IPM_ID<-MJprot$variable

MJprotD<-left_join(MJprot,MD,by="IPM_ID")
##DrawBoxPlot
MJprotDF <- mutate_all(MJprotDF, function(x) as.character(x))

MJprotDF <- filter(MJprotD, Hydrocephalus == "NPIH" | Hydrocephalus == "PIH")

ggplot(MJprotDF,
       aes(Symbol, value,fill= Hydrocephalus)) +
  geom_boxplot(size = 1, alpha = 0.8, ) +
  scale_y_log10() +#ylim(0.0002, 1e-08)+
  #geom_jitter(size = 0.5)+
  ggtitle("ProtExpPIH_NPIH") +
  xlab("Protein")+ylab("Normalized Counts")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c("#9c954d","#b067a3"))


ggplot(MJprotDF,
       aes(Symbol, value,fill= Paeni_status)) +
  geom_boxplot(size = 1, alpha = 0.8, ) +
  scale_y_log10() +#ylim(0.0002, 1e-08)+
  #geom_jitter(size = 0.5)+
  ggtitle("ProtExpPIH_NPIH") +
  xlab("Protein")+ylab("Normalized Counts")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c("#e41a1c","#377eb8","#4daf4a"))
#PaeniStatusWithColorForAlbert
