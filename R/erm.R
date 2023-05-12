erm <- function(
    weatherDTfile= NULL,
    vars=NULL,
    verbose=FALSE
){
  
  id <- paste( paste("erm",cgiarPIPE::idGenerator(5,5),sep=""), weatherDTfile$idOriginal, sep = "_")
  type <- "erm"
  
  ############################
  # loading the dataset
  if (is.null(weatherDTfile)) stop("No input weather data file specified.")
  mydata <- weatherDTfile$cleaned #readRDS(file.path(wd,"files_cleaned",paste0(weatherDTfile)))
  mydata <- unique(mydata)
  
  if(is.null(vars)){
    vars <- c(
      "wind_speed_10m:ms", # "wind_dir_10m:d",
      "wind_gusts_10m_24h:ms", # "wind_gusts_10m_1h:ms",
      "t_2m:C", "t_max_2m_24h:C", "t_min_2m_24h:C", 
      "msl_pressure:hPa",
      "precip_24h:mm", # "precip_1h:mm",
      "relative_humidity_1000hPa:p",
      # "sunrise:sql","sunset:sql",
      "dew_point_2m:C",
      "uv:idx"
    )
  }
  vars <- gsub(":","",vars)
  
  wide <- reshape(df[,c("fieldinst","validdate",vars)], direction = "wide", idvar = "fieldinst",
                  timevar = "validdate", v.names = vars, sep= "_")
  X <- wide[,-1]
  rownames(X) <- wide[,1]
  X <- scale(X, center = TRUE, scale = TRUE)
  Acor <- cor(t(X))
  Acov <- cov(t(X))
  
  ## save as metrics
  Acor[lower.tri(Acor)] <- NA
  Cdf <- as.data.frame(as.table(Acor)) # converts a matrix in a data frame
  Cdf <- Cdf[which(!is.na(Cdf$Freq)),]
  ##
  Acov[lower.tri(Acov)] <- NA
  Gdf <- as.data.frame(as.table(Acov)) # converts a matrix in a data frame
  Gdf <- Gdf[which(!is.na(Gdf$Freq)),]
  # write pipeline metrics
  pm <- data.frame(value=Gdf$Freq,  stdError=NA,
                   fieldinst=paste(Gdf$Var1, Gdf$Var2, sep="###"),  trait=paste(vars, collapse = "::"),
                   analysisId=id, method="cov",
                   traitUnits=NA, parameter="covE",
                   pipeline=NA,
                   stage =NA
  )
  pm2 <- data.frame(value=Cdf$Freq,  stdError=NA,
                    fieldinst=paste(Cdf$Var1, Cdf$Var2, sep="###"),  trait=paste(vars, collapse = "::"),
                    analysisId=id, method="cor",
                    traitUnits=NA, parameter="corE",
                    pipeline=NA,
                    stage =NA
  )
  pms <- rbind(pm,pm2)
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile =	weatherDTfile$id,
    markerbooks	= NA,  markerDataFile =	NA,
    timePoint = NA,  season =	NA,  location =	NA,
    country	= NA,  trial	= NA,  design =	NA,
    geno = NA,  rep	= NA,  block =	NA,
    rowcoord =	NA,  colcoord = NA,
    stage = NA
  )
  # write the values used for cleaning to the modeling database
  # write predictions
  # write pipeline metrics
  
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
    cat(paste("Your results will be available in the files_cleaned folder under such id \n"))
  }
  result <- list(metrics=pms, predictions=NA, modeling=NA, metadata=db.params,
                 cleaned=Acor, outliers=NA, desire=NA, id=id, idOriginal=weatherDTfile$idOriginal,
                 metadataFieldinst=NA)
  return(result)#paste("grm done:",id))
}