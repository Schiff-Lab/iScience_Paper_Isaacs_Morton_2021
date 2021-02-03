library(plotrix)
library(reshape2)
library(dplyr)
library(ggplot2)
library(data.table)
library(readxl)
library(ggalluvial)
library(stringr)

colAllu <- rep("text", 13)##All  text
Alluvi<-read_excel("FinalBonferroniGOEnrichment.xlsx",col_types=colAllu)

Alluvi1<-strsplit(gsub(" ", "", as.character(Alluvi$Genes)), split = ",")
Alluvi2<-data.frame(Term = rep(Alluvi$Term, sapply(Alluvi1, length)), Genes = unlist(Alluvi1))

Alluvi2<-Alluvi2%>% dplyr::rename(Ensemble=Genes)
Alluvi2<-as.data.frame(Alluvi2)


Connector1<-read.table("Ensemble_List.txt", header=T,
                      sep="\t", as.is=F, quote = "")

colnames(Connector1)<-c("Symbol","Ensemble")

joinAlluvi2Conn1<-left_join(Alluvi2,Connector1,by="Ensemble")

joinAlluvi2Conn2<-left_join(joinAlluvi2Conn1,Alluvi,by="Term")

joinAlluvi2Conn3<-joinAlluvi2Conn2 %>%
  select(Category,Symbol,Term,Ensemble,Bonferroni)%>%
  mutate(Category=replace(Category,Category=="GOTERM_BP_DIRECT", "Biological_Process")) %>%
  mutate(Category=replace(Category,Category=="GOTERM_MF_DIRECT", "Molecular_Function")) %>%
  mutate(Category=replace(Category,Category=="GOTERM_CC_DIRECT", "Cellular_Component")) %>%
  as.data.frame()
ProtRNA<-read.table("OutPutFileForInitPathwayAnalysisAllGenesUpANDDown_MJM.txt", header=T,
                       sep="\t", as.is=F, quote = "")


ProtRNA<-ProtRNA%>% dplyr::rename(Symbol=Genes)
ProtRNA<-ProtRNA %>%
  mutate(Symbol=str_replace(Symbol,"_Protein", "")) %>%
  mutate(Symbol=str_replace(Symbol,"_RNA", "")) %>%
           as.data.frame()

ProtRNA$Status<- 0

for (i in 1:nrow(ProtRNA))
{
  if (ProtRNA$Protein[i]=="yes" & ProtRNA$RNA[i]=="yes" ){
    ProtRNA$Status[i] <- "Both"
  }
  else if (ProtRNA$Protein[i]=="yes" & ProtRNA$RNA[i]=="no" ){
    ProtRNA$Status[i] <- "Protein"
  }
  else if (ProtRNA$Protein[i]=="no" & ProtRNA$RNA[i]=="yes" ){
    ProtRNA$Status[i] <- "RNA"
  }
}

joinAlluvi2Conn4<-left_join(joinAlluvi2Conn3,ProtRNA,by="Symbol")
#joinAlluvi2Conn4[is.na(joinAlluvi2Conn4)] <- "Not_Enriched"

joinAlluvi2Conn4<-joinAlluvi2Conn4 %>% dplyr::select(Protein,RNA,Term,Symbol,Status)%>% dplyr::filter(!is.na(Protein))


alluFin<-unique(joinAlluvi2Conn4)
is_alluvia_form(as.data.frame(alluFin), axes = 1:3, silent = TRUE)

alluFin$Status <- as.factor(alluFin$Status)
ggplot(alluFin,
       aes(x = Term, stratum = Status, alluvium = Symbol,
           fill = Status, label = Status)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Assay Status Versus Go-Term")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Number of Genes")+xlab("Go Term")+
  scale_fill_manual(values=c("#83a444","#9170c7","#cb547e") )


write.table(alluFin,"alluFin.txt",sep="\t",row.names = FALSE,quote=FALSE)

