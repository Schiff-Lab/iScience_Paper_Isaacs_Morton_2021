
library(EnsDb.Hsapiens.v79)
library(PANTHER.db)


# 2. Convert from gene.symbol to ensembl.gene
geneSymbols <-  read.table("GoEnrichment/gene_list.txt",
                           sep="\t", as.is=F, quote = "")
geneSymbols<-as.vector(as.character(geneSymbols$V1))
geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

write.table(geneIDs2,"Ensemble_List.txt",quote=FALSE, sep="\t")
