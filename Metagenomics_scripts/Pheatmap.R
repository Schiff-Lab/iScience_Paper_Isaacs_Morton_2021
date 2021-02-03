#PheatmapSaraAndAlbert
library("pheatmap")
library("RColorBrewer")
library("reshape2")
library("dplyr")
library(tidyverse)


##StartFromTheOutputFileMadeInSara_Albert.R File
Reactome<-read.table("OutputProtLFC1RNALFC1MJ.txt", header=T,
               sep="\t", as.is=F, quote = "") # as.is=F because we need them factors!#****Also use quote="For words like alzymer's with quote

#ReactomeC <- rep("text", 16)##All 7 text
#Reactome<-read_excel("eexa.xlsx", col_types=ReactomeC )
###ForPresentInBoth
Reactome2<-Reactome %>% filter(RNA =="yes"& Protein=="yes")
View(unique(Reactome2$Genes))
View(unique(Reactome2$Pathwayname))

#SubReactomeP <- Reactome[grep("_Protein",Reactome$Genes), ]
#SubReactomeR <- Reactome[grep("_RNA",Reactome$Genes), ]
#SubReactomeC<-rbind(SubReactomeP,SubReactomeR)
#SubReactomeC<-SubReactomeC%>%select(Genes, Pathwayname,logFC_concat,P_Value_concat,FDR_concat,Protein,RNA)
#SubReactomeC<-as.data.frame(SubReactomeC)
#WideSubReactome<-dcast(SubReactomeC,Pathwayname~Genes,value.var = "logFC_concat")

WideSubReactome<-dcast(Reactome2, Pathwayname~Genes, value.var ='logFC_concat')

# encode NAs as a small negative value, so we can color it differently (i.e. as missing)
#WideSubReactome[is.na(WideSubReactome)] <- 0.00001
#WideSubReactome[is.na(WideSubReactome)] <- -1.5
#WideSubReactome[is.na(WideSubReactome)] <- 0.00001
##ForJoe'sCommentOnRedandGrennOnly
WideSubReactome[ WideSubReactome>0.45 ] <- 1
WideSubReactome[ WideSubReactome< -1] <- -1
WideSubReactome[is.na(WideSubReactome)] <- -1.5

##
row.names(WideSubReactome)<-WideSubReactome[,1]
WideSubReactome <-WideSubReactome %>% select (-Pathwayname)
names(WideSubReactome) <- gsub(x = names(WideSubReactome), pattern = "_Protein", replacement = "_P")
names(WideSubReactome) <- gsub(x = names(WideSubReactome), pattern = "_RNA", replacement = "_R")

##CreateAnnotationFile

AnnotFile<-Reactome2%>%select(Genes, Protein)
names(AnnotFile)<-str_replace_all(names(AnnotFile), c("Protein" = "Assay" ))

#AnnotFile<-AnnotFile %>%
 # mutate(Assay = str_replace(Assay, "yes" , "Protein"))
#AnnotFile<-AnnotFile %>%
 # mutate(Assay = str_replace(Assay, "no" , "RNA"))
AnnotFile<-AnnotFile %>%
   mutate(Assay = str_replace(Assay, "yes" , "RNA_Protein"))
AnnotFile<-as.data.frame(unique(AnnotFile$Genes))
###PatchWorkNeedToremakeThispartOfScript
##ConditionalForBothCategory(StillNeedsWorkButForNowPrintsFALSE)

rename_col_by_position <- function(df, position, new_name) {
  new_name <- enquo(new_name)
  new_name <- quo_name(new_name)
  select(df, !! new_name := !! quo(names(df)[[position]]), everything())
}
AnnotFile <- rename_col_by_position(AnnotFile,1,"Genes")
AnnotFile$Assay<-0


for (i in 1:nrow(AnnotFile))
{
  if (grepl("^.*Protein.*$", AnnotFile$Genes[i], ignore.case = T)) {
      AnnotFile$Assay[i] <- "Protein"
      }

  else if (grepl("^.*RNA.*$", AnnotFile$Genes[i], ignore.case = T)) {
    AnnotFile$Assay[i] <- "RNA"
  }
}
AnnotFile$Genes <- gsub("_Protein", "_P", AnnotFile$Genes)
AnnotFile$Genes <- gsub("_RNA", "_R", AnnotFile$Genes)

AnnotFile<-(distinct(AnnotFile,Genes, .keep_all = TRUE))
rownames(AnnotFile)<-AnnotFile$Genes
AnnotFile <-AnnotFile %>% select (-Genes)
AnnotFile$Assay <- factor(AnnotFile$Assay, levels = c("Protein", "RNA"))
#AnnotFile$Assay <- factor(AnnotFile$Assay, levels = c("RNA_Protein"))

library(stringr)

my_colour =  list(
  Assay = c(Protein = "#5977ff", RNA = "purple"))

#my_colour =  list(
 #  Assay = c(RNA_Protein="green"))


#my_colour =  list(
  #Assay = c(Protein = "#5977ff", RNA = "#f74747"))
# make up some colors (set lowest value to specific color, i.e. grey80)
#color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdYlGn"))(200))[50:200])
color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdGn"))(100)))
breaks <- c(0.00002, 0.00001, seq(0, 3.2, length.out=99))
#Default Hyrarchical Clustering
pheatmap(WideSubReactome,annotation_col = AnnotFile,
         annotation_colors = my_colour,
         breaks=breaks, col=color,cluster_cols = F,
         fontsize_row = 7,fontsize_col = 7,
         cellwidth=6,cellheight=6,main="Heatmap LogFold-Change Paeni VS Non-Paeni")


###NoteThisIsforDownregulatedGenes###UsedForFinalFigureJoeWanted
color <- c("white","darkgreen", rev(colorRampPalette(colors=brewer.pal(1, "Reds"))(3)))
breaks <- c(-1.4, -1.5, seq(-1, 4, length.out=2))
#Default Hyrarchical Clustering ####FinalOneUsedForThePaper

newnames <- lapply(
  rownames(WideSubReactome),
  function(x) bquote(bold(.(x))))

pheatmap(WideSubReactome,annotation_col = AnnotFile,
         annotation_colors = my_colour,
         breaks=breaks, col=color,cluster_cols = F,
         fontsize_row = 10,fontsize_col = 7,
         cellwidth=8,cellheight=12,labels_row = as.expression(newnames),
         main=" Heatmap Log2-Fold-Change Paeni VS Non-Paeni")
##UseFor
pheatmap(WideSubReactome,annotation_col = AnnotFile,
         annotation_colors = my_colour,
         col=color,
         fontsize_row = 7,fontsize_col = 7,labels_row = make_bold_names,
         cellwidth=6,cellheight=9)

###ForAllHeatmap
color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdBu"))(100)))
breaks <- c(-5.35,-5.4 , seq(-5.5, 3.2, length.out=101))
#Default Hyrarchical Clustering
pheatmap(WideSubReactome,annotation_col = AnnotFile,
         annotation_colors = my_colour,
         breaks=breaks, col=color,
         fontsize_row = 7,fontsize_col = 7,
         cellwidth=6,cellheight=9)

##ForAll
pheatmap(WideSubReactome,annotation_col = AnnotFile,
         annotation_colors = my_colour,
         breaks=breaks, col=color,
         fontsize_row = 7,fontsize_col = 7
         )

##Do this for Protein Only
Reactome<-read.table("OutPutFileForInitPathwayAnalysis_MJM.txt",header=T,
                     sep="\t" ,as.is=F) # as.is=F because we need them factors!
ReactomeP<-Reactome %>% filter(Protein=="yes")
ReactomeP<-ReactomeP %>% filter(Presence=="One")

SubReactomeP<-ReactomeP%>%select(Genes, Pathwayname,logFC_concat)

WideSubReactomeP<-dcast(SubReactomeP,Pathwayname~Genes,value.var = "logFC_concat")
# encode NAs as a small negative value, so we can color it differently (i.e. as missing)
WideSubReactomeP[is.na(WideSubReactomeP)] <- 0.00001
row.names(WideSubReactomeP)<-WideSubReactomeP[,1]
WideSubReactomeP <-WideSubReactomeP %>% select (-Pathwayname)

##CreateAnnotationFile

AnnotFileP<-Reactome%>%select(Genes, Protein)
names(AnnotFileP)<-str_replace_all(names(AnnotFileP), c("Protein" = "Assay" ))
AnnotFileP<-AnnotFileP %>%
  mutate(Assay = str_replace(Assay, "yes" , "Protein"))

AnnotFileP<-(distinct(AnnotFileP,Genes, .keep_all = TRUE))
rownames(AnnotFileP)<-AnnotFileP$Genes
AnnotFileP <-AnnotFileP %>% select (-Genes)
AnnotFileP$Assay <- factor(AnnotFileP$Assay, levels = c("Protein"))
my_colourP =  list(
  Assay = c(Protein = "#5977ff"))
# make up some colors (set lowest value to specific color, i.e. grey80)
#color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdYlGn"))(200))[50:200])
color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdYlGn"))(100)))
breaks <- c(0.00002, 0.00001, seq(0, 1.9, length.out=101))
#Default Hyrarchical Clustering
pheatmap(WideSubReactomeP,annotation_col = AnnotFileP,
         annotation_colors = my_colourP,
         breaks=breaks, col=color,
         fontsize_row = 7,fontsize_col = 7,
         cellwidth=6,cellheight=9)

###ForRNAOnly
Reactome<-read.table("OutPutFileForInitPathwayAnalysis_MJM.txt",header=T,
                     sep="\t" ,as.is=F) # as.is=F because we need them factors!
ReactomeR<-Reactome %>% filter(RNA=="yes")
ReactomeR<-ReactomeR %>% filter(Presence=="One")
SubReactomeR<-ReactomeR%>%select(Genes, Pathwayname,logFC_concat)

WideSubReactomeR<-dcast(SubReactomeR,Pathwayname~Genes,value.var = "logFC_concat")
# encode NAs as a small negative value, so we can color it differently (i.e. as missing)
WideSubReactomeR[is.na(WideSubReactomeR)] <- 0.00001
row.names(WideSubReactomeR)<-WideSubReactomeR[,1]
WideSubReactomeR <-WideSubReactomeR %>% select (-Pathwayname)

##CreateAnnotationFile

AnnotFileR<-ReactomeR%>%select(Genes, RNA)
names(AnnotFileR)<-str_replace_all(names(AnnotFileR), c("RNA" = "Assay" ))

AnnotFileR<-AnnotFileR %>%
  mutate(Assay = str_replace(Assay, "yes" , "RNA"))
AnnotFileR<-(distinct(AnnotFileR,Genes, .keep_all = TRUE))
rownames(AnnotFileR)<-AnnotFileR$Genes
AnnotFileR <-AnnotFileR %>% select (-Genes)
AnnotFileR$Assay <- factor(AnnotFileR$Assay, levels = c("RNA"))
my_colourR =  list(
  Assay = c(RNA = "#f74747"))
# make up some colors (set lowest value to specific color, i.e. grey80)
#color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdYlGn"))(200))[50:200])
color <- c("white", rev(colorRampPalette(colors=brewer.pal(11, "RdYlGn"))(100)))
breaks <- c(0.00002, 0.00001, seq(0, 2, length.out=101))
#Default Hyrarchical Clustering
pheatmap(WideSubReactomeR,annotation_col = AnnotFile,
         annotation_colors = my_colourR,
         breaks=breaks, col=color,
         fontsize_row = 7,fontsize_col = 7,
         cellwidth=6,cellheight=9)




