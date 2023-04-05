clustPheno <- function(
    phenoDTfile=NULL,
    fieldinst= NULL, # fieldinst= "-"
    traits =NULL, # traits <- unique(mydata$trait);traits
    wideBy="trait", # column to be used for reshaping
    kMax=30,
    kPicked=5, # to be updated by user while running
    method="wss", # "silhouette"  "gap_stat"
    verbose=FALSE
){

  id <- paste( paste("clu",cgiarPIPE::idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "clu"
  if (is.null(phenoDTfile)) stop("No input phenotypic data file specified.")
  if (is.null(traits)) stop("No traits specified.")
  if (is.null(fieldinst)) stop("No fieldinst specified.")
  if(length(fieldinst) > 1 & length(traits) >1){
    stop("We can only handle multiple traits for one field or one trait for multiple fields.")
  }
  if((wideBy == "trait") & length(traits) < 2){
    stop("To cluster by traits we need more than one trait.")
  }
  if((wideBy == "fieldinst") & length(fieldinst) < 2){
    stop("To cluster by fieldinst we need more than one fieldinst.")
  }

  mydata=phenoDTfile$predictions
  mydata <- mydata[which(mydata$fieldinst %in% fieldinst),]
  mydata <- mydata[which(mydata$trait %in% traits),]
  mydataWide <- reshape(mydata[,c("geno",wideBy,"predictedValue")], direction = "wide", idvar = "geno",
                        timevar = wideBy, v.names = "predictedValue", sep= "_")
  colnames(mydataWide) <- gsub("predictedValue_","",colnames(mydataWide))
  # impute missing values
  levelsOfWideBy <- unique(mydata[,wideBy])
  preProcValues <- caret::preProcess(mydataWide[,levelsOfWideBy], method = c("knnImpute", "center", "scale"))
  WsProc <- predict(preProcValues, mydataWide)
  # build dissimilarity matrix
  dsm <- cluster::daisy(WsProc[,levelsOfWideBy], metric = "euclidean", stand = FALSE)
  # To create WSS plot to identify number of clusters
  res1 <- factoextra::fviz_nbclust(WsProc[,levelsOfWideBy], cluster::pam, method = method, k.max = kMax)
  # end of clustering
  pamRes <- cluster::pam(WsProc[,levelsOfWideBy], kPicked)
  # Add cluster assignments to data
  WsClust <- cbind(WsProc, cluster = pamRes$cluster)

 # Summary statistics of each cluster per Trait
  des <- psych::describeBy(WsClust[-c(1,length(WsClust))], group=WsClust$cluster, skew=FALSE, mat=T, digits=2)
  des <- tibble::rownames_to_column(des, "Traits") %>% dplyr::select(-vars,-item) %>%
    dplyr::rename(Cluster=group1)

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

  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")

  phenoDTfile$predictions <- mydataRes[,predcols]
  phenoDTfile$metadata <-  db.params
  phenoDTfile$id <- id
  phenoDTfile$clust <- res1$data
  phenoDTfile$des <- des
  return(phenoDTfile)
}
