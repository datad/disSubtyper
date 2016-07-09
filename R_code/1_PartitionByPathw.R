#author: Diana Diaz (dmd@wayne.edu)



###############
## SNF parameters
K = 20;##number of neighbors, usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
NIT = 10; ###Number of Iterations, usually (10~20)


#GENERAL PARAMETERS
C = 3  #number of clusters for colon cancer


libs<-c("survival_analysis.R","Packages.R", "permutation.R", "clustering.R", "Plots.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')

sapply(libs, function(u) {source(u)})

loadls("SNFtool survival ROntoTools org.Hs.eg.db")

data("geneEData")
data("survivalData")


typeOfClustering <- 3 #1=kmeans \ 2=SNF  \ 3=HC



resultList <- rep(list(),length(kegg))
kaplanq <- rep(NA,length(kegg))
coxScore <- rep(NA,length(kegg))
coxpValue <- rep(NA,length(kegg))
sizePath <- rep(0,length(kegg))
smallestC <- rep(NA,length(kegg))
allNormal <- rep(NA,length(kegg))



if(typeOfClustering == 1) #kmeans
{
  resultsFolder <- paste("output_kmeans_2_Colon_table1_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else if(typeOfClustering == 2){
  resultsFolder <- paste("output_SNF_2_Colon_table1_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}else{
  resultsFolder <- paste("output_HC_2_Colon_table1_Genes_",
                         format(Sys.time(), "%b_%d_%Y"),"/", sep='')
}


dir.create(resultsFolder)

pdf( file=paste(resultsFolder, "plots.pdf", sep="")  ,  onefile=TRUE)  
for(i in seq_len(length(kegg))){
  if(length(nodes(kegg[[i]])) >= minNnodes){
    sizePath[i] <- length(nodes(kegg[[i]]))
    if(typeOfClustering == 1) #kmeans
    {
      resultList[[i]] <- as.list(clusterByPathwayK(kegg[[i]], survivalData))
    }else if(typeOfClustering == 2){
      resultList[[i]] <- as.list(clusterByPathwaySNF(kegg[[i]], survivalData))
    }else{
      resultList[[i]] <- as.list(clusterByPathwayHC(kegg[[i]], survivalData))
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
