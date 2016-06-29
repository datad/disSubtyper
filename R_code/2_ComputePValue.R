kegg <- keggPathwayGraphs("hsa")
namsPath <- keggPathwayNames("hsa")
namsPath <- namsPath[names(kegg)]



############
#RANDOMIZATION
loadlib("bootstrap", bioc=FALSE)

#Number of random pathways 
Nrp <- 10000

#all genes of the data
totalNumberOfGenes <- dim(geneEData)[2]
allGenes <- colnames(geneEData)

#global for theta function
survivalDa <- survivalData

resultList <- rep(list(),length(kegg))

pValofCox <- rep(NA,length(kegg))

if(typeOfClustering == 1) #kmeans
{
  resultsFolder <- paste("output_kmeans_3_Colon_table2_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else if(typeOfClustering == 2){
  resultsFolder <- paste("output_SNF_3_Colon_table2_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else{
  resultsFolder <- paste("output_HC_3_Colon_table2_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}

dir.create(resultsFolder)


pdf( file=paste(resultsFolder, "plots.pdf", sep="")  ,  onefile=TRUE)  
for(i in seq_len(length(row.names(tablePaths)))){
  #get the size of the pathway of interest
  interst_Path <- kegg[[row.names(tablePaths)[i]]]

  interestingGenes <- nodes(interst_Path)
  interestingGenes <- geneSetFormat(interestingGenes)
  NgenesIp <- length(interestingGenes)
  #the meassured p-value
  indicesOfIP <- which(colnames(geneEData) %in% interestingGenes, arr.ind=TRUE)
  if(typeOfClustering == 1) #kmeans
  {
    crcPthpval <- thetaK(indicesOfIP)
    results <- resampling(x=1:totalNumberOfGenes,
                          sampleSize=NgenesIp,
                          nboot=Nrp, thetaK)
  }else if(typeOfClustering == 2){
    crcPthpval <- thetaSNF(indicesOfIP)
    results <- resampling(x=1:totalNumberOfGenes,
                          sampleSize=NgenesIp,
                          nboot=Nrp, thetaK)
  }else{
    crcPthpval <- thetaHC(indicesOfIP)
    results <- resampling(x=1:totalNumberOfGenes,
                          sampleSize=NgenesIp,
                          nboot=Nrp, thetaHC)
  }
  
  
  resultList[[i]] <- results
  hx <- results$thetastar
  FinalPCP = mean(hx <= crcPthpval)
  pValofCox[i] <- FinalPCP
  denshx <- density(hx, kernel="gaussian", from =0, to = 1, bw="SJ" )
  plot(denshx, main=paste("Kernel Density Estimation\n",
                          namsPath[row.names(tablePaths)[i]]))
  abline(v=FinalPCP,col=4)  
  legend("topright", paste("p-val=",FinalPCP))
}
dev.off()

tablePaths <- cbind(tablePaths,pValofCox)
View(tablePaths)

View(tablePaths[with(tablePaths, order(pValofCox,coxpValue)),] )


write.csv(tablePaths[with(tablePaths, order(pValofCox,coxpValue)),] , 
          file=paste(resultsFolder, "tablePaths.csv", sep="") )

save.image(file=paste(resultsFolder, "resultImage.RData", sep=""))

pathwaysTable(tablePaths)
