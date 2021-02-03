#SaraAndAlbertDeconvolutionKWTest

#RunThisAfterInitDeconvScript
KWinit<-as.data.frame(t(df))
Annotation<-annodf
combined<-cbind(Annotation,KWinit)

testingKW <- function(x, cols=3:41) {
  return(sapply(cols, function(i) kruskal.test(x[, i] ~ Paeni , data = x)$p.value))
}

results <- as.data.frame(testingKW(combined))
colnames(results)<-c("KW_P")
colComb<-as.data.frame(t(colnames(combined)))
colComb<-colComb%>%dplyr::select(-V1,-V2)
colComb<-as.data.frame(t(colComb))
results$Pathway<-colComb$V1
results$p_adjP<-p.adjust(results$KW_P, method = "bonferroni")

combined<-combined%>%dplyr::select(-Paeni)

testingKWH <- function(x, cols=2:40) {
  return(sapply(cols, function(i) kruskal.test(x[, i] ~ Hydrocephalus , data = x)$p.value))
}

resultsH <- as.data.frame(testingKWH(combined))
colnames(resultsH)<-c("KW_H")
colCombH<-as.data.frame(t(colnames(combined)))
colCombH<-colCombH%>%dplyr::select(-V1)
colCombH<-as.data.frame(t(colCombH))
resultsH$Pathway<-colCombH$V1
resultsH$p_adjH<-p.adjust(resultsH$KW_H, method = "bonferroni")

JoinResults<-left_join(resultsH,results, by= "Pathway")
JoinResultsAll <- JoinResults[ , c("Pathway", "KW_H", "p_adjH", "KW_P" ,"p_adjP")]

write.table(JoinResultsAll, file = "DeconvKWTestForAllResults.txt",
            sep = "\t", col.names = TRUE,row.names = F)

Sig<-dplyr::filter(JoinResults,p_adjH<0.05 | p_adjP<0.05)
colnames(Sig)
SigDeconv <- Sig[ , c("Pathway", "KW_H", "p_adjH", "KW_P" ,"p_adjP")]

write.table(SigDeconv, file = "DeconvKWTest.txt",
            sep = "\t", col.names = TRUE,row.names = F)

SigDeconvP<-SigDeconv$Pathway

##PvalOnlySig
Sig2<-dplyr::filter(JoinResults,KW_H<0.05 | KW_P<0.05)

sigWOadjk<-Sig2[ , c("Pathway", "KW_H", "p_adjH", "KW_P" ,"p_adjP")]

sigWOadjP<-sigWOadjk$Pathway
#NowGoBackToHeatmapCode!!!
