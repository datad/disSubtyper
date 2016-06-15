source("global_parameters.R")
typeOfClustering <- 3 #1=kmeans \ 2=SNF  \ 3=HC

source("survival_analysis.R")


#define a minimum threshold in number of node in the pathway
#threshold
minNnodes <- 5

resultList <- rep(list(),length(kegg))
kaplanq <- rep(NA,length(kegg))
coxScore <- rep(NA,length(kegg))
coxpValue <- rep(NA,length(kegg))
sizePath <- rep(0,length(kegg))
smallestC <- rep(NA,length(kegg))
allNormal <- rep(NA,length(kegg))



#Writing the output


if(typeOfClustering == 1) #kmeans
{
  resultsFolder <- paste("output_kmeans_2_Colon_table1_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else{
  resultsFolder <- paste("output_SNF_2_Colon_table1_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}


dir.create(resultsFolder)

pdf( file=paste(resultsFolder, "plots.pdf", sep="")  ,  onefile=TRUE)  
#compare the 149 pathways
for(i in seq_len(length(kegg))){
  if(length(nodes(kegg[[i]])) >= minNnodes){
    sizePath[i] <- length(nodes(kegg[[i]]))
    if(typeOfClustering == 1) #kmeans
    {
      resultList[[i]] <- as.list(clusterByPathwayK(kegg[[i]], survivalData))
    }else{
      resultList[[i]] <- as.list(clusterByPathwaySNF(kegg[[i]], survivalData))
    }
    k <-  resultList[[i]]$kaplan
    kaplanq[i] <-k$chisq
    coxScore[i] <- resultList[[i]]$coxScore
    coxpValue[i] <- resultList[[i]]$coxPval
    sizeC1 <- resultList[[i]]$sizGroup1
    sizeC2 <- resultList[[i]]$sizGroup2
    sizeC3 <- resultList[[i]]$sizGroup3
    smallestC <- min(sizeC1,sizeC2,sizeC3)  
    #writing plot
    mainTitle <- as.character(namsPath[i])
    plot.SurvivalK3 (resultList[[i]]$groups, mainTitle, survivalData,"bottomright")
  }
}

dataTable <- data.frame(pathway=namsPath, coxScore, coxpValue, sizePath, smallestC)
boxplot(dataTable$coxpValue, main = "Cox pvalue of the pathways")
dev.off()

tablePaths <- dataTable[order(coxpValue),]

View(tablePaths )
write.csv(tablePaths, 
          file=paste(resultsFolder, "dataTable3c.csv", sep="") )
save.image(file=paste(resultsFolder, "resultImage.RData", sep=""))

