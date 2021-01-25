filterLowGenes<-function(obj,minSamples,threshold=1){
  counts <- cpm(obj)
  keep <- rowSums(counts > threshold) >= minSamples
  obj <- obj[keep, ]
  obj
}