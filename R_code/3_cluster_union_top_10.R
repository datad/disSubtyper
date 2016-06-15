# jun 12 2016 ------------------------------
## TODO
# rm(list=ls())
# cat("\014")


source("R/Colon_src/global_parameters.R")
source("R/forPaper/tables_fdr.R")

source("D:/My_Docs/Programming_projects/R/Utils/utilitary/Packages.R")

loadlib("bootstrap", bioc=FALSE)

###############
# #from previous EXPERIMENT
# #save(genesOfInterest,FinalPCP,results, file=paste(resultsFolder, "resultGeneList.RData", sep=""))
# iFolder <-
#   "initialValidation/output_kmeans_4_Colon_union_Genes_Feb_16_2016/"
# load(file = paste(iFolder, "resultGeneList.RData" , sep = ''))
###############

#TODO manually change
# if(typeOfClustering == 1) #kmeans
# {
#   foldersData <- "initialValidation/output_kmeans_3_Colon_table2_Genes_Feb_15_2016/"
#
# }else{
#   foldersData <- "output_SNF_3_Colon_table2_Genes_Dec_15_2015/"
# }

# tablePaths <- read.table(paste(foldersData, "tablePaths.csv",sep =""),
#                          sep=",", quote = '"', header = TRUE)
# row.names(tablePaths) <- as.character(tablePaths[,1])
# tablePaths <- tablePaths[,-1]

# View(tablePaths)

#correct p-vlaues
#cannot all goes to 1
#p.adjust(tablePaths$pValofCox, method="fdr")


 
###########

###find pwathways with the same groups
# for(i in seq_len(length(resultList))){
#   groups <- resultList[[i]]$groups
#   if(all(groupGE == groups))
#     cat(i, " , ")
# }
#
# for(i in seq_len(length(resultList))){
#   groups <- resultList[[i]]$groups
#   csizesg <- c(length(groups[groups==1]), length(groups[groups==2]), length(groups[groups==3]),
#                +               length(groups[groups==4]) )
#   if(all(csizes %in% csizesg))
#     cat(i, " , ")
# }


#90 pathways get the same clustering that using all genes, so what is the improvement
#length(which(dataTable$sizeC2 == 4))


# resultListi <- clusterByPathway(kegg[[1]])
# for (i in seq_len(length(groupsKEGG))){
#   comparison[i] <- (groupsKEGG[i]==groupSNF[i])
# }
#
# which(!comparison)
# Only 6 out of 92 are different
# Is this normal for all the pathways?
# Repeat with all the pathways
###########

#findUnionGenes <-function()
#{
  tablePaths <- tablePaths[!is.na(tablePaths$pValofCox),]
  # head(tablePaths)
  
  ToptablePaths <-
    tablePaths[(tablePaths$pValofCox < 0.05 &
                  tablePaths$coxpValue < 0.05),]
  
  ToptablePaths <<- ToptablePaths[1:TOP,]
  
  #View(ToptablePaths)
  
  
  namPaths <- row.names(ToptablePaths)
  unionGenes <- c()
  
  for (i in seq_len(length(namPaths))) {
    iGenes <- nodes(kegg[[namPaths[i]]])
    unionGenes <- union(unionGenes, iGenes)
  }
  
  length(unionGenes)
  
  unionGenes <- unionGenes
  #get the size of the pathway of interest
  genesOfInterest <- unionGenes
  genesOfInterest <- gsub("hsa:","",genesOfInterest)
  genesOfInterest <- lapply(genesOfInterest, find_sym)
  genesOfInterest <- unlist(genesOfInterest[!is.na(genesOfInterest)])
  #Number of mRNA on pathway
  genesOfInterest <-
    intersect(colnames(geneEData), as.character(genesOfInterest))
  genesOfInterest
#}

  
  
  length(genesOfInterest)


#computeCoxPPofFinalSet <- function(disease, genesOfInterest)
#{
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
    
    #     results <- resampling(
    #       x = 1:totalNumberOfGenes,
    #       sampleSize = NgenesIp,
    #       nboot = Nrp, thetaK
    #     )
    
  }else{
    groups <- getGroupsFromGenesSNF(genesOfInterest2)
    crcPthpval <- thetaSNF(indicesOfIP)
    
    #     results <- resampling(
    #       x = 1:totalNumberOfGenes,
    #       sampleSize = NgenesIp,
    #       nboot = Nrp, thetaSNF
    #     )
    
  }
  
  
  #   hx <- results$thetastar
  #  
  #   FinalPCP <- (1 + sum(hx <= crcPthpval)) / (1 + totalNumberOfGenes)
  
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
  }else{
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
  }
  
  
  
  
  
  #recoloring
  #Try
  #by renaming levels of a factor
  # library(plyr)
  # groups <- revalue(groups2, c("2"="3", "3"="2"))
  # groups2 <- groups
  
  #OR manually
  # levels(groups)[levels(groups)=="2"] <- "4"
  # levels(groups)[levels(groups)=="3"] <- "2"
  # levels(groups)[levels(groups)=="4"] <- "3"
  
  #finally, if nothing works change the colors paprameters

  
  
  
  #   pdf(file = paste(resultsFolder, "density_plots.pdf", sep = "")  ,  onefile =
  #         TRUE)
  
  #   denshx <-
  #     density(
  #       hx, kernel = "gaussian", from = 0, to = 1, bw = "SJ"
  #     )
  #   plot(denshx, main = "Kernel Density Estimation\n Union of genes from top pathways")
  #   abline(v = FinalPCP,col = 4)
  #   abline(v = crcPthpval,col = 4)
  #   legend("topright", paste(
  #     "p-val=",round(FinalPCP,5),
  #     "\n Cox p-val=", round(crcPthpval,5)
  # ))
  
  # dev.off()
  
  #Save the whole workspace
  # set save defaults using option:
  #options(save.defaults = list(ascii = TRUE, safe = FALSE))
  save.image(file = paste(resultsFolder, "resultsImage.RData", sep = ""))
  
  
  write.csv(
    genesOfInterest, row.names = FALSE,
    file = paste(resultsFolder, "genesOfInterest.csv", sep = "")
  )
  
  #Write session
  write(toLatex(sessionInfo(), locale = FALSE),
        file = paste(resultsFolder, "session.tex", sep = ""))
  
  pathwaysTable(ToptablePaths)
  finalTable(resultsFolder)
  #TODO MANUALLY COPY THE VERSION OF KEEG
  #   cat("MANUALLY COPY THE VERSION OF KEEG")
  #   write("kegg version ", file = paste(resultsFolder, "keggVersion.tex", sep =
  #                                         ""))
  #   kpg <- keggPathwayGraphs("hsa")
  
  
  write.csv(ToptablePaths,
           file = paste(resultsFolder, "topPathws.csv", sep = ""))
  
  

  save(genesOfInterest,groups, file=paste(resultsFolder, "resultGeneList.RData", sep=""))
  # iFolder <-
  #   "initialValidation/output_kmeans_4_Colon_union_Genes_Feb_16_2016/"
  # load(file = paste(iFolder, "resultGeneList.RData" , sep = ''))
  
#}
