clustSTPGA <- function(
    markerDTfile=NULL,
    kMaxPop=30,
    kPicked=50, # to be updated by user while running
    nElite=5,
    nIterations=10,
    method="CDMEAN", # "silhouette"  "gap_stat"
    verbose=FALSE
){
  
  id <- paste( paste("ops",cgiarPIPE::idGenerator(5,5),sep=""), markerDTfile$idOriginal, sep = "_")
  type <- "ops"
  if (is.null(markerDTfile)) stop("No input genotypic data file specified.")
  
  mydataWide=markerDTfile$cleaned$M
  kPicked <- min(c(kPicked, nrow(mydataWide))) # take the smalles, selected size or number of indiviuals
  
  svdWide<-svd(mydataWide, nu=5, nv=5)
  PC50Wide<-mydataWide%*%svdWide$v
  colnames(PC50Wide) <- paste0("PC",1:ncol(PC50Wide))
  res1 <- STPGA::GenAlgForSubsetSelection(P=PC50Wide,
                              Candidates=rownames(PC50Wide),
                              Test=rownames(PC50Wide),  
                              InitPop=NULL, # a list of initial solutions
                              ntoselect=kPicked,
                              npop=kMaxPop, # number of solutions at each iteration
                              nelite=nElite, # number of solutions selected as elite parents which will generate the next set of solutions
                              mutprob=.5, # probability of mutation for each generated solution.
                              mutintensity = rpois(1,4), # mean of the poisson variable that is used to decide the number of mutations for each cross.
                              niterations=nIterations,
                              minitbefstop=5, # number of iterations before stopping if no change is observed in criterion value.
                              tabumemsize = 2, # Number of generations to hold in tabu memory.
                              plotiters=FALSE,
                              lambda=1e-9,errorstat=method, mc.cores=1)
 
  # Add cluster assignments to data
  pcaList <- list(); colnamespca <- colnames(PC50Wide); counter1 <- 1
  nPC <- min(c(ncol(PC50Wide),3))
  for(iCol in 1:nPC){
    pcaList[[counter1]] <- data.frame(analysisId=id, pipeline="unknown",
                                      trait=paste(colnamespca[iCol], sep="-"), geno=rownames(PC50Wide),
                                      mother="unknown",father="unknown",
                                      genoType="general",fieldinst="across", predictedValue=PC50Wide[,iCol],
                                      stdError=1, rel=1, genoYearOrigin=1, genoYearTesting=1, stage="rank1",
                                      genoCode= 1:nrow(PC50Wide)
    ); counter1 <- counter1+1
  }
  predictionsBind <- do.call(rbind,pcaList)
  predictionsBind$genoType[which(predictionsBind$geno %in% res1$`Solution with rank 1`)]="subset"
  
  predcols <- c("analysisId", "pipeline","trait","genoCode","geno","mother","father","genoType","genoYearOrigin",
                "genoYearTesting", "fieldinst","predictedValue","stdError","rel","stage")
  predictionsBind <- predictionsBind[,predcols]
  
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
  
  res <- list(cleaned=NA, predictions=predictionsBind, modeling=NA, metrics=NA,
              metadata=db.params, id=id, idOriginal=markerDTfile$idOriginal)
  
  return(res)
}
