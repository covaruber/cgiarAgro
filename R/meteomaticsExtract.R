meteomaticsExtract <- function(phenoDTfile= NULL, verbose=FALSE){
  
  id <- paste( paste("wea",cgiarPIPE::idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "wea"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("clp",phenoDTfile$id)) == 0){stop("Index can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  
  ############################
  # loading the dataset
  if(is.null(phenoDTfile$metadataFieldinst)){stop("There's no metadata for this file. Likely belongs to an older version of the cgiarPIPE package. Please match the columns of your data again.", call. = FALSE)}
  mydata <- phenoDTfile$metadataFieldinst # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  ############################
  ## index calculation
  
  badLats <- which(mydata$latitude == 1)
  badLons <- which(mydata$longitude == 1)
  if(length(badLats) > 0  | length(badLons) > 0){
    stop("You have provided latitudes and longitudes with value equal to 1. Please provide valid latitudes and longitudes to extract weather data.", call. = FALSE)
  }
  ###########################
  ## extract data
  username <- 'irri_covarrubias'
  password <- 'wG55H37mxF'
  
  l1 <- as.list(mydata$latitude)
  l2 <- as.list(mydata$longitude)
  coordinates <- mapply(c, l1, l2, SIMPLIFY=FALSE)
  
  # startdate0 <- as.Date(min(mydata$plantingDate))
  # enddate0 <- as.Date(max(mydata$harvestingDate))
  # startdate1 <- min(c(startdate0,enddate0))
  # enddate1 <- max(c(startdate0,enddate0))
  # ## Specify the start and end date
  # startdate <- ISOdatetime(year = as.integer(strftime(startdate1, '%Y')),
  #                          month = as.integer(strftime(startdate1, '%m')),
  #                          day = as.integer(strftime(startdate1, '%d')),
  #                          hour = 00, min = 00, sec = 00, tz = 'UTC')
  # enddate <- ISOdatetime(year = as.integer(strftime(enddate1, '%Y')),
  #                        month = as.integer(strftime(enddate1, '%m')),
  #                        day = as.integer(strftime(enddate1, '%d')) + 1,
  #                        hour = 00, min = 00, sec = 00, tz = 'UTC')
  
  startdate <- ISOdatetime(year = as.integer(strftime(lubridate::today(), '%Y')),
                           month = as.integer(strftime(lubridate::today(), '%m')),
                           day = as.integer(strftime(lubridate::today(), '%d')),
                           hour = 00, min = 00, sec = 00, tz = 'UTC')
  enddate <- ISOdatetime(year = as.integer(strftime(lubridate::today(), '%Y')),
                         month = as.integer(strftime(lubridate::today(), '%m')),
                         day = as.integer(strftime(lubridate::today(), '%d')) + 1,
                         hour = 00, min = 00, sec = 00, tz = 'UTC')
  
  ## Specify a time interval
  interval <- "PT1H"
  
  ## Specify the parameter(s) of interest
  parameters <- list("t_2m:C", "dew_point_2m:C", "relative_humidity_1000hPa:p", "precip_1h:mm")
  
  ## Call the MeteomaticsRConnector::query_time_series() function
  wdata0 <- MeteomaticsRConnector::query_time_series(coordinates, startdate, enddate, interval, parameters,
                          username, password)
  colnames(wdata) <- gsub(":","",colnames(wdata))
  
  myDates <- lapply(strsplit(as.character(wdata$validdate), split=" "), function(x){x[1]})
  myTimes <- lapply(strsplit(as.character(wdata$validdate), split=" "), function(x){x[2]})
  myDatesDf <- as.data.frame(do.call(rbind, lapply(myDates, function(x){strsplit(x,"-")[[1]]}) ) )
  myHoursDf <- as.data.frame(do.call(rbind, lapply(myTimes, function(x){strsplit(x,":")[[1]]}) ) )
  wdataDf <- cbind(myDatesDf,myHoursDf,wdata0)
  colnames(wdataDf)[1:6] <- c("year","month","day","hour","minute","second")
  ## come up with aggregated data for each fieldinst
  wdata2 <- merge(wdataDf,mydata, by.x = c("lat","lon"), by.y = c("latitude","longitude"), all.x = TRUE)
  myFormula <- paste("cbind(", paste(c( "t_2mC", "dew_point_2mC", "relative_humidity_1000hPap", "precip_1hmm"), collapse = ","),")~fieldinst")
  aggregatedWdata <- aggregate(as.formula(myFormula), data=wdata2, FUN=mean)
  ## merge environmental summaries to metadata
  metadataFieldinst <- merge(mydata,aggregatedWdata, by="fieldinst", all.x = TRUE )
  
  ## write the parameters to the parameter database
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
  ## write the values used for cleaning to the modeling database
  mod <- data.frame(
    trait = NA,
    traitLb = NA,
    traitUb = NA,
    outlierCoef = NA,
    analysisId = id,
    analysisType = type,
    fixedModel = NA,
    randomModel = NA,
    residualModel = NA,
    h2Threshold = NA
  )
  
  # write pipeline metrics
  pm <- data.frame(value=NA,
                   stdError=NA,
                   fieldinst="none",
                   trait="all",
                   analysisId=id, method= "terra",
                   traitUnits=NA,
                   stage = NA,
                   parameter= "many",
                   pipeline=NA
  )
  
  if(verbose){
    cat(paste("Your analysis id is:",id,"\n"))
  }
  result <- list(metrics=pm, predictions=NA, modeling=mod, metadata=db.params,
                 cleaned=wdataDf, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
  
}
