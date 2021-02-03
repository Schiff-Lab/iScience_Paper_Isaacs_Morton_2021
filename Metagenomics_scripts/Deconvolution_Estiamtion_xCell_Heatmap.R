library(yarn)
library(immunedeconv)
library(pheatmap)
library(dplyr)


  
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

es2 <- filterLowGenes(eSet,minSamples=18)
es2 <- es2[-which(is.na(fData(es2)$gene_name)),]
cnts <- exprs(es2)
rownames(cnts) <- fData(es2)$gene_name
width = (fData(es2)$end-fData(es2)$start)
tpm <- tpm3(counts = cnts,len = width)
xcell <- immunedeconv::deconvolute(tpm, 'xcell')
x2 <- xcell
xcell <- x2 %>% as.matrix
rownames(xcell) <- xcell[,1]
xcell = xcell[,-1]
for(i in seq(ncol(xcell))) xcell[,i] = as.numeric(as.character(xcell[,i]))

df <- array(NA,dim=c(nrow(xcell),ncol(xcell)))
for(i in seq(ncol(xcell))) df[,i] = as.numeric(as.character(xcell[,i]))
rownames(df) <- rownames(xcell);  colnames(df) <- colnames(xcell);
annodf <- data.frame(Paeni=factor(pData(es2)$new_paeni),
                     Hydrocephalus=pData(es2)$Hydrocephalus,row.names = colnames(df))
pheatmap(df,annotation_col = annodf)

###StartAfterKWTestHere

#Keep Only Those rows which are part of SigDeconvP
dfN <- df[rownames(df) %in% SigDeconvP, ]

pheatmap(dfN,annotation_col = annodf)

##MJ'sPrettyHeatmaps
annodf$Paeni<-as.character(annodf$Paeni)
annodfM<-annodf
annodfM<-annodfM %>%
  mutate(Paeni = sub("FALSE", "Negative", Paeni))

annodfM<-annodfM %>%
  mutate(Paeni = sub("TRUE", "Positive", Paeni))

rownames(annodfM)<-rownames(annodf)

my_colors <- list(
  Paeni = c(Negative = "#e41a1c", Positive = "#4daf4a"),
  Hydrocephalus = c(NPIH="#377eb8",PIH="orange")
)

library(viridis)
library(RColorBrewer)

pheatmap(dfN,annotation=annodfM,
         annotation_colors = my_colors,
         color= c("grey", rev(
           colorRampPalette(colors = rev(brewer.pal(11, "RdPu")))(100)
         )),
         fontsize_row = 8,
         fontsize_col=4,
         fontsize = 8,
         cellheight =14,
         cellwidth = 4,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")
##ThisIsForUnadjusted

#Keep Only Those rows which are part of SigDeconvP

##MJ'sPrettyHeatmaps
annodf$Paeni<-as.character(annodf$Paeni)
annodfM<-annodf
annodfM<-annodfM %>%
  mutate(Paeni = sub("FALSE", "Negative", Paeni))

annodfM<-annodfM %>%
  mutate(Paeni = sub("TRUE", "Positive", Paeni))

rownames(annodfM)<-rownames(annodf)

my_colors <- list(
  Paeni = c(Negative = "#e41a1c", Positive = "#4daf4a"),
  Hydrocephalus = c(NPIH="#377eb8",PIH="orange")
)

library(viridis)
library(RColorBrewer)

pheatmap(df,annotation=annodfM,
         annotation_colors = my_colors,
         color= c("grey", rev(
           colorRampPalette(colors = rev(brewer.pal(11, "RdPu")))(100)
         )),
         fontsize_row = 8,
         fontsize_col=4,
         fontsize = 8,
         cellheight =14,
         cellwidth = 4,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")

#PValueSigOnly

dfN <- df[rownames(df) %in% sigWOadjP, ]

pheatmap(dfN,annotation_col = annodf)

##MJ'sPrettyHeatmaps
annodf$Paeni<-as.character(annodf$Paeni)
annodfM<-annodf
annodfM<-annodfM %>%
  mutate(Paeni = sub("FALSE", "Negative", Paeni))

annodfM<-annodfM %>%
  mutate(Paeni = sub("TRUE", "Positive", Paeni))

rownames(annodfM)<-rownames(annodf)

my_colors <- list(
  Paeni = c(Negative = "#e41a1c", Positive = "#4daf4a"),
  Hydrocephalus = c(NPIH="#377eb8",PIH="orange")
)

library(viridis)
library(RColorBrewer)

pheatmap(dfN,annotation=annodfM,
         annotation_colors = my_colors,
         color= c("grey", rev(
           colorRampPalette(colors = rev(brewer.pal(11, "RdPu")))(100)
         )),
         fontsize_row = 8,
         fontsize_col=4,
         fontsize = 8,
         cellheight =14,
         cellwidth = 4,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean")
