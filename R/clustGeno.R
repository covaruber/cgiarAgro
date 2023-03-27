clustGeno <- function(
    markerDTfile=NULL,
    kMax=30,
    kPicked=5, # to be updated by user while running
    method="wss", # "silhouette"  "gap_stat"
    verbose=FALSE
){

  id <- paste( paste("mcl",cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "mcl"
  if (is.null(markerDTfile)) stop("No input genotypic data file specified.")

  mydataWide=markerDTfile$cleaned$M
  # impute missing values
  preProcValues <- caret::preProcess(mydataWide, method = c("knnImpute", "center", "scale"))
  WsProc <- predict(preProcValues, mydataWide)
  # build dissimilarity matrix
  dsm <- cluster::daisy(WsProc, metric = "euclidean", stand = FALSE)
  # To create WSS plot to identify number of clusters
  res1 <- factoextra::fviz_nbclust(WsProc, cluster::pam, method = method, k.max = kMax)
  # end of clustering
  pamRes <- cluster::pam(WsProc, kPicked)
  # Add cluster assignments to data
  predictionsBind <- data.frame(analysisId=id, pipeline=NA,trait="cluster",genoCode= pamRes$cluster,
                                geno=rownames(mydataWide),mother="BASE_MOTHER",father="BASE_FATHER",
                                genoType=NA,genoYearOrigin=1,genoYearTesting=1, fieldinst=NA,
                                predictedValue=NA,stdError=NA,rel=NA,stage=NA)

  ##############################

  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	markerDTfile$id,
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

  markerDTfile$predictions <- predictionsBind
  markerDTfile$metadata <-  db.params
  markerDTfile$id <- id
  markerDTfile$clust <- res1$data
  return(markerDTfile)
}
