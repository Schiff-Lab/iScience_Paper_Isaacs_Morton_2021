#SaraAlbertProject
load("eSetNF.rda")
#ExpressionData
ExpD<-(eSetNF@assayData$exprs)
#PhenotypeData
Phenotype<-eSetNF@phenoData@data
#Gene
Gene<-eSetNF@featureData@data
d<-density(ExpD)
plot(d)

library(dplyr)
library(plyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(readxl)
library(reshape)

DiffExp<-read.table("NS_TOTAL_SIGNIFICANT.txt",header=T,
                     sep="\t" ,as.is=F) # as.is=F because we need them factors!

#SideTracPlanToOnlyGetThoseWithabove2FCforNow
DiffExpFoldChange2<-DiffExp[(DiffExp$logFC > 2 | DiffExp$logFC < -2 ), ]

#SideTracPlanToOnlyGetThoseWithabove1FCforNow
DiffExpFoldChange1<-DiffExp[(DiffExp$logFC > 1 | DiffExp$logFC < -1 ), ]


write.table(DiffExpFoldChange2,"GeneExpabs_FC2.txt", sep="\t")
write.table(DiffExpFoldChange1,"GeneExpabs_FC1.txt", sep="\t")


ExpDC<-ExpD
ExpDC<-as.data.frame(ExpDC)
ExpDC <- tibble::rownames_to_column(ExpDC,"Symbol")
ExpDC$Symbol <- sub("\\..*","\\",ExpDC$Symbol)


FC_ExpD<-left_join(ExpDC,DiffExp, by ="Symbol")
rownames(FC_ExpD) <- FC_ExpD[,1]

FC_ExpD<-FC_ExpD %>% select (-Symbol,-genes)
FC_ExpD<-na.omit(FC_ExpD)
FC_ExpD <- FC_ExpD[ which(FC_ExpD$FDR <0.05),]
#FC_ExpD <- FC_ExpD[ which(abs(FC_ExpD$logFC) > 2),]

FC_ExpD$Between0_2<-0
FC_ExpD$Between3_6<-0
FC_ExpD$Between6_10<-0
FC_ExpD$Above10<-0

FC_ExpDM<-as.matrix(FC_ExpD)

#Mean

for (i in 1:nrow(FC_ExpDM))
{ ##This is a matrix so $ will not work for accessing the column use ""
  if (mean(FC_ExpDM[i,1:98]) < 3 & mean(FC_ExpDM[i,1:98]) > 0 ){
    FC_ExpDM[i,"Between0_2"] <-mean(FC_ExpDM[i,1:98])
  }
  else if (mean(FC_ExpDM[i,1:98]) < 6 & mean(FC_ExpDM[i,1:98]) > 2 ){
    FC_ExpDM[i, "Between3_6"] <- mean(FC_ExpDM[i,1:98])
  }
  else if (mean(FC_ExpDM[i,1:98]) < 10 & mean(FC_ExpDM[i,1:98]) > 5){
    FC_ExpDM[i,"Between6_10"] <- mean(FC_ExpDM[i,1:98])
  }
  else{
    FC_ExpDM[i,"Above10"] <- mean(FC_ExpDM[i,1:98])
  }
}

#Quantile 75%

for (i in 1:nrow(FC_ExpDM))
{ ##This is a matrix so $ will not work for accessing the column use ""
  if (quantile(FC_ExpDM[i,1:98],probs =0.75) < 3 & quantile(FC_ExpDM[i,1:98],probs =0.75) > 0 ){
    FC_ExpDM[i,"Between0_2"] <- quantile(FC_ExpDM[i,1:98],probs =0.75)
  }
  else if (quantile(FC_ExpDM[i,1:98],probs =0.75) < 6 & quantile(FC_ExpDM[i,1:98],probs =0.75) > 2 ){
    FC_ExpDM[i, "Between3_6"] <- quantile(FC_ExpDM[i,1:98],probs =0.75)
  }
  else if (quantile(FC_ExpDM[i,1:98],probs =0.75) < 10 & quantile(FC_ExpDM[i,1:98],probs =0.75) > 5){
    FC_ExpDM[i,"Between6_10"] <- quantile(FC_ExpDM[i,1:98],probs =0.75)
  }
  else{
    FC_ExpDM[i,"Above10"] <- quantile(FC_ExpDM[i,1:98],probs =0.75)
  }
}

#Qunatile 25%
for (i in 1:nrow(FC_ExpDM))
{ ##This is a matrix so $ will not work for accessing the column use ""
  if (quantile(FC_ExpDM[i,1:98],probs =0.25) < 3 & quantile(FC_ExpDM[i,1:98],probs =0.25) > 0 ){
    FC_ExpDM[i,"Between0_2"] <- quantile(FC_ExpDM[i,1:98],probs =0.25)
  }
  else if (quantile(FC_ExpDM[i,1:98],probs =0.25) < 6 & quantile(FC_ExpDM[i,1:98],probs =0.25) > 2 ){
    FC_ExpDM[i, "Between3_6"] <- quantile(FC_ExpDM[i,1:98],probs =0.25)
  }
  else if (quantile(FC_ExpDM[i,1:98],probs =0.25) < 10 & quantile(FC_ExpDM[i,1:98],probs =0.25) > 5){
    FC_ExpDM[i,"Between6_10"] <- quantile(FC_ExpDM[i,1:98],probs =0.25)
  }
  else{
    FC_ExpDM[i,"Above10"] <- quantile(FC_ExpDM[i,1:98],probs =0.25)
  }
}

FC_ExpDM<-as.data.frame(FC_ExpDM)

FC_ExpDM1 <- filter(FC_ExpDM, Between0_2 > 0)
FC_ExpDM2 <- filter(FC_ExpDM, Between3_6 > 0)
FC_ExpDM3 <- filter(FC_ExpDM, Between6_10 > 0)
FC_ExpDM4 <- filter(FC_ExpDM, Above10 > 0)


##Through Mean
par(mfrow=c(2,2))
plot(FC_ExpDM1$logFC,FC_ExpDM1$Between0_2, main="Mean Number of Reads
     Per Transcript vs LFC group 0-2", ylab="Mean Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM2$logFC,FC_ExpDM2$Between3_6,main="Mean Number of Reads
     Per Transcript vs LFC group 3-6", ylab="Mean Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM3$logFC,FC_ExpDM3$Between6_10,main="Mean Number of Reads
     Per Transcript vs LFC group 6-10", ylab="Mean Of Reads Per Transcript",xlab="FoldChange Per Transcript")


par(mfrow=c(2,2))
plot(log10(FC_ExpDM1$PValue),FC_ExpDM1$Between0_2, main="Mean Number of Reads
     Per Transcript vs P_Value group 0-2", ylab="Mean Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")
plot(log10(FC_ExpDM2$PValue),FC_ExpDM2$Between3_6,main="Mean Number of Reads
     Per Transcript vs P_Value group 3-6", ylab="Mean Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")
plot(log10(FC_ExpDM3$PValue),FC_ExpDM3$Between6_10,main="Mean Number of Reads
     Per Transcript vs P_Value group 6-10", ylab="Mean Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")

#Through Quantile 75%

FC_ExpDM<-as.data.frame(FC_ExpDM)
par(mfrow=c(2,2))
plot(FC_ExpDM$logFC,FC_ExpDM$Between0_2, main="75% Quantile Number of Reads
     Per Transcript vs LFC group 0-2", ylab="75% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Between3_6,main="75% Quantile Number of Reads
     Per Transcript vs LFC group 3-6", ylab="75% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Between6_10,main="75% Quantile Number of Reads
     Per Transcript vs LFC group 6-10", ylab="75% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Above10,main="75% Quantile Number of Reads Per Transcript
     vs LFC group above 10", ylab="75% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")



#Through Quantile 25%

FC_ExpDM<-as.data.frame(FC_ExpDM)
par(mfrow=c(2,2))
plot(FC_ExpDM$logFC,FC_ExpDM$Between0_2, main="25% Quantile Number of Reads
     Per Transcript vs LFC group 0-2", ylab="25% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Between3_6,main="25% Quantile Number of Reads
     Per Transcript vs LFC group 3-6", ylab="25% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Between6_10,main="25% Quantile Number of Reads
     Per Transcript vs LFC group 6-10", ylab="25% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")
plot(FC_ExpDM$logFC,FC_ExpDM$Above10,main="25% Quantile Number of Reads Per Transcript
     vs LFC group above 10", ylab="25% Quantile Of Reads Per Transcript",xlab="FoldChange Per Transcript")


par(mfrow=c(2,2))
plot(log10(FC_ExpDM$PValue),FC_ExpDM$Between0_2, main="25% Quantile Number of Reads
     Per Transcript vs P_Value group 0-2", ylab="25% Quantile Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")
plot(log10(FC_ExpDM$PValue),FC_ExpDM$Between3_6,main="25% Quantile Number of Reads
     Per Transcript vs P_Value group 3-6", ylab="25% Quantile Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")
plot(log10(FC_ExpDM$PValue),FC_ExpDM$Between6_10,main="25% Quantile Number of Reads
     Per Transcript vs P_Value group 6-10", ylab="25% Quantile Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")
plot(log10(FC_ExpDM$PValue),FC_ExpDM$Above10,main="25% Quantile Number of Reads
Per Transcript vs P_Value group above 10", ylab="25% Quantile Of Reads Per Transcript",xlab="Log 10 P_Value Per Transcript")


library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels


#q_value<0.05
Significance=ifelse(FC_ExpDM$FDR<0.01, "padj<0.01", "Not Sig")
volc = ggplot(FC_ExpDM, aes(logFC, -log10(PValue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=Significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "orange")) +
  ggtitle("Volcano Plot High Versus Medium Groups Gene Expression") #e.g. 'Volcanoplot DESeq2'
volc+geom_text_repel(data=filter(FC_ExpDM, FDR<0.05), aes(label=gene),
                     segment.size  = 0.05,size=2.3) #adding text for the top 20 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
volc
