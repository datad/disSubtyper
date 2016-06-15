#k-means and SNF
#no normal
#rm(list=ls())
#cat("\014")

#typeOfClustering <- 1 #1=kmeans \ 2=SNF

#source("Testing/Colon_src/global_parameters.R")


# foldersData <- "output_kmeans_2_Colon_table1_Genes_Feb_15_2016/"
# tablePaths <- read.table(paste(foldersData, "dataTable3c.csv",sep =""), 
#                          sep=",", quote = '"', header = TRUE)
# row.names(tablePaths) <- as.character(tablePaths[,1])
# tablePaths <- tablePaths[,-1]
View(tablePaths)

kegg <- keggPathwayGraphs("hsa")
namsPath <- keggPathwayNames("hsa")
namsPath <- namsPath[names(kegg)]
############
#HERE THE RANDOMIZATION!
loadlib("bootstrap", bioc=FALSE)

#Number of random pathways 
Nrp <- 10000

#all genes of the data
totalNumberOfGenes <- dim(geneEData)[2]
allGenes <- colnames(geneEData)

#global for theta function
survivalDa = survivalData

resultList <- rep(list(),length(kegg))

pValofCox <- rep(NA,length(kegg))

if(typeOfClustering == 1) #kmeans
{
  resultsFolder <- paste("output_kmeans_3_Colon_table2_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else{
  resultsFolder <- paste("output_SNF_3_Colon_table2_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}

dir.create(resultsFolder)

pdf( file=paste(resultsFolder, "plots.pdf", sep="")  ,  onefile=TRUE)  
for(i in seq_len(length(row.names(tablePaths)))){
#get the size of the pathway of interest
  
  interst_Path <- kegg[[row.names(tablePaths)[i]]]
  #
  interestingGenes <- nodes(interst_Path)
  interestingGenes <- geneSetFormat(interestingGenes)
  NgenesIp <- length(interestingGenes)
  #the meassured p-value
  indicesOfIP <- which(colnames(geneEData) %in% interestingGenes, arr.ind=TRUE)
  if(typeOfClustering == 1) #kmeans
  {
    crcPthpval <- thetaK(indicesOfIP)
    #results <- bootstrap(1:NgenesIp,Nrp,theta)
    results <- resampling(x=1:totalNumberOfGenes,
                          sampleSize=NgenesIp,
                          nboot=Nrp, thetaK)
  }else{
    crcPthpval <- thetaSNF(indicesOfIP)
    results <- resampling(x=1:totalNumberOfGenes,
                          sampleSize=NgenesIp,
                          nboot=Nrp, thetaK)
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
