
tablePaths <- tablePaths[!is.na(tablePaths$pValofCox),]


ToptablePaths <-
  tablePaths[(tablePaths$pValofCox < 0.05 &
                tablePaths$coxpValue < 0.05),]


namPaths <- row.names(ToptablePaths)
unionGenes <- c()

for (i in seq_len(length(namPaths))) {
  iGenes <- nodes(kegg[[namPaths[i]]])
  unionGenes <- union(unionGenes, iGenes)
}

length(unionGenes)

unionGenes <- unionGenes
genesOfInterest <- unionGenes
genesOfInterest <- gsub("hsa:","",genesOfInterest)
genesOfInterest <- lapply(genesOfInterest, find_sym)
genesOfInterest <- unlist(genesOfInterest[!is.na(genesOfInterest)])
#Number of mRNA on pathway
genesOfInterest <-
  intersect(colnames(geneEData), as.character(genesOfInterest))
genesOfInterest


length(genesOfInterest)



kegg <- keggPathwayGraphs("hsa")
NgenesIp <- length(genesOfInterest)


#all genes of the data
totalNumberOfGenes <- dim(geneEData)[2]
allGenes <- colnames(geneEData)

############
#HERE THE RANDOMIZATION!

#Number of random pathways
Nrp <- 10000

#the meassured p-value
indicesOfIP <-
  which(colnames(geneEData) %in% genesOfInterest, arr.ind = TRUE)

genesOfInterest2 <- allGenes[indicesOfIP]


if (typeOfClustering == 1)
{     #kmeans
  groups <- getGroupsFromGenesK(genesOfInterest2)
  crcPthpval <- thetaK(indicesOfIP)
}else if (typeOfClustering == 2){
  groups <- getGroupsFromGenesSNF(genesOfInterest2)
  crcPthpval <- thetaSNF(indicesOfIP)
} else {
  groups <<- getGroupsFromGenesK(genesOfInterest2)
  crcPthpval <- thetaK(indicesOfIP)
}

disease = "colon"

############
#Writing the output
if (typeOfClustering == 1)
{    #kmeans
  resultsFolder <- paste("output_kmeans_4_",disease,"_union_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep = '')
  dir.create(resultsFolder)
  mainTitle <- "b) Survival curve, k-means, selected genes"
  pdf(
    file = paste(resultsFolder, "Kmeans_union_kaplan.pdf", sep = ""), onefile =
      TRUE, width = 9, height = 7
  )
  plot.SurvivalK3 (groups, mainTitle, survivalData,"bottomright",
                   colorsP = c("red","black", "blue"))
  dev.off()
}else if (typeOfClustering == 2){
  resultsFolder <- paste("output_SNF_4_",disease,"_union_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep = '')
  
  dir.create(resultsFolder)
  mainTitle <- "b) Survival curve, SNF, selected genes"
  pdf(
    file = paste(resultsFolder, "SNF_union_kaplan.pdf", sep = ""),  onefile =
      TRUE, width = 9, height = 7
  )
  plot.SurvivalK3 (groups, mainTitle, survivalData,"bottomright",
                   colorsP = c("red","blue","black"))
  dev.off()
} else {
  resultsFolder <- paste("output_HC_4_",disease,"_union_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep = '')
  
  dir.create(resultsFolder)
  mainTitle <- "b) Survival curve, hierarchical, selected genes"
  pdf(
    file = paste(resultsFolder, "HC_union_kaplan.pdf", sep = ""),  onefile =
      TRUE, width = 9, height = 7
  )
  plot.SurvivalK3 (groups, mainTitle, survivalData,"bottomright",
                   colorsP = c("red","blue","black"))
  dev.off()
}



save.image(file = paste(resultsFolder, "resultsImage.RData", sep = ""))


write.csv(
  genesOfInterest, row.names = FALSE,
  file = paste(resultsFolder, "genesOfInterest.csv", sep = "")
)


write(toLatex(sessionInfo(), locale = FALSE),
      file = paste(resultsFolder, "session.tex", sep = ""))

pathwaysTable(ToptablePaths)
finalTable(resultsFolder)

write.csv(ToptablePaths,
          file = paste(resultsFolder, "topPathws.csv", sep = ""))



save(genesOfInterest,groups, file=paste(resultsFolder, "resultGeneList.RData", sep=""))

