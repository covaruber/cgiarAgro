clustPheno <- function(
    phenoDTfile=NULL,
    fieldinst= NULL, # fieldinst= "-"
    traits =NULL, # traits <- unique(mydata$trait);traits
    kMax=30,
    kPicked=5, # to be updated by user while running
    method="wss", # "silhouette"  "gap_stat"
    verbose=FALSE
){

  id <- paste( paste("clu",idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "clu"
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  if (is.null(traits)) stop("No traits specified.")
  if (is.null(fieldinst)) stop("No fieldinst specified.")
  # library(caret) # should go into a different cgiar library, has too many dependencies
  # library('RANN')
  # library(factoextra) # Begin Partition Around Medoid Clustering
  # library(NbClust)
  # library(cluster)

  mydata=phenoDTfile$predictions
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),]
  mydataWide <- reshape(mydata[,c("geno","trait","predictedValue")], direction = "wide", idvar = "geno",
                        timevar = "trait", v.names = "predictedValue", sep= "_")
  colnames(mydataWide) <- gsub("predictedValue_","",colnames(mydataWide))
  # impute missing values
  preProcValues <- caret::preProcess(mydataWide[,traits], method = c("knnImpute", "center", "scale"))
  WsProc <- predict(preProcValues, mydataWide)
  # build dissimilarity matrix
  dsm <- cluster::daisy(WsProc[,traits], metric = "euclidean", stand = FALSE)
  # To create WSS plot to identify number of clusters
  res1 <- factoextra::fviz_nbclust(WsProc[,traits], pam, method = method, k.max = kMax) + theme_classic()
  res1$labels
  head(res1$data)
  # end of clustering
  pamRes <- cluster::pam(WsProc[,traits], kPicked)
  # Add cluster assignments to data
  WsClust <- cbind(WsProc, cluster = pamRes$cluster)
  # then bring back to long format again
  mydataRes <- merge(mydata, WsClust[,c("geno","cluster")], by="geno", all.x = TRUE)
  mydataRes$genoCode <- mydataRes$cluster
  mydataRes$cluster <- NULL

  ##############################

  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	phenoDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )

  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }

  phenoDTfile$predictions <- mydataRes
  phenoDTfile$metadata <-  db.params

  return(phenoDTfile)
}
