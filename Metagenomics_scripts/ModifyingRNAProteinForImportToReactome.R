#ModifyiingProteinandRNAData
library(ggplot2)
library("readxl")
library("dplyr")
library(tidyverse)
library("forcats")
library("data.table")
library("stringr")
library("plyr")
library("reshape")



DiffExp<-read.table("NS_TOTAL_SIGNIFICANT.txt",header=T,
                    sep="\t" ,as.is=F) # as.is=F because we need them factors!

#SideTracPlanToOnlyGetThoseWithabove2FCforNow
DiffExpFoldChange2<-DiffExp[(DiffExp$logFC > 2 | DiffExp$logFC < -2 ), ]

write.table(DiffExpFoldChange2,"GeneExpabs_FC2.txt", sep="\t")

##Had Some error renaming with standard functions bc of ...1 name hence let's do this by position
rename_col_by_position <- function(df, position, new_name) {
  new_name <- enquo(new_name)
  new_name <- quo_name(new_name)
  select(df, !! new_name := !! quo(names(df)[[position]]), everything())
}

ProteomicsF <- rep("text", 7)##All 7 text
Proteomics<-read_excel("AllProteinFile.xlsx", col_types=ProteomicsF )
Proteomics <- rename_col_by_position(Proteomics,1,"Genes")
##WeirdFileitShould be Proteomics$logFC < -1 )somehow it is not getting the memo this works though
Proteomics2<-Proteomics[(Proteomics$logFC > 1 | Proteomics$logFC > -1 ), ]

write.table(Proteomics2,"Proteomics_AfterFilteringLFC2.txt", sep="\t")

