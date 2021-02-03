library(dplyr)
library(readxl)

options(stringsAsFactors=FALSE)
#df <- data.frame(pathway=c("BLA_1", "BLA_2"), genes=c("G1;G2;G3", "G10;G1"))
fKidCol <- rep("text", 2)##All  text
Kids<-read_excel("Reactome_Pathways.xlsx",col_types=fKidCol)
Kids<-as.data.frame(Kids)
awesome <- lapply(1:nrow(Kids), function(i) data.frame(pathway=Kids[i, "pathway"],
                                                       gene=unlist(strsplit(Kids[i, "gene"], ';')))) %>% bind_rows()


geneLCol <- rep("text", 2)##All  text
geneL<-read_excel("Gens_List.xlsx",col_types=geneLCol)


#LJ

Fin<-left_join(awesome,geneL,by="gene")


#Plot
library(plotrix)
library(reshape2)
library(dplyr)
library(ggplot2)
library(data.table)
library(readxl)
library(ggalluvial)
library(stringr)

alluFin<-unique(Fin)
is_alluvia_form(as.data.frame(Fin), axes = 1:3, silent = TRUE)
alluFin$Status <- as.factor(Fin$Symbol)

ggplot(Fin,
       aes(x = pathway, stratum = Symbol, alluvium = gene,
           fill = Symbol, label = Symbol)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Assay Status Versus Pathway")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Number of Genes")+xlab("Reactom Pathway Analysis")+
  scale_fill_manual(values=c("#83a444","#cb547e","blue") )




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

