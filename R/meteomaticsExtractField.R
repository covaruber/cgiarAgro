meteomaticsExtractField <- function(phenoDTfile= NULL, verbose=FALSE, interval="PT12H"){
  
  id <- paste( paste("wea",cgiarPIPE::idGenerator(5,5),sep=""), phenoDTfile$idOriginal, sep = "_")
  type <- "wea"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("clp",phenoDTfile$id)) == 0){stop("This function can only be used in matched phenotypic files",call. = FALSE)}
  
  ############################
  # loading the dataset
  if(is.null(phenoDTfile$metadataFieldinst)){stop("There's no metadata for this file. Likely belongs to an older version of the cgiarPIPE package. Please match the columns of your data again.", call. = FALSE)}
  mydata <- phenoDTfile$metadataFieldinst # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  # mydata$latitude <- sample(seq(-90,90, 180/length(mydata$latitude)), length(mydata$latitude))
  # mydata$longitude <- sample(seq(-180,180, 360/length(mydata$longitude)), length(mydata$longitude))
  # mydata$plantingDate <- Sys.Date()
  # mydata$harvestingDate <- Sys.Date()
  badLats <- which(mydata$latitude == 1000)
  badLons <- which(mydata$longitude == 1000)
  if((length(badLons) > 0) | (length(badLats) > 0) ){
    stop("Some coordinates are not valid. Please correct and try again", call. = FALSE)
  }
  
  mydata$latitude <- round(mydata$latitude,2)
  mydata$longitude <- round(mydata$longitude,2)
  ############################
  ## index calculation
  
  ###########################
  ## extract data
  username <- 'irri_covarrubias'
  password <- 'wG55H37mxF'
  
  l1 <- as.list(mydata$latitude)
  l2 <- as.list(mydata$longitude)
  coordinates <- mapply(c, l1, l2, SIMPLIFY=FALSE)
  
  startdate0 <- as.Date(min(mydata$plantingDate))
  enddate0 <- as.Date(max(mydata$harvestingDate))
  startdate1 <- min(c(startdate0,enddate0))
  enddate1 <- max(c(startdate0,enddate0))
  ## Specify the start and end date
  startdate <- ISOdatetime(year = as.integer(strftime(startdate1, '%Y')),
                           month = as.integer(strftime(startdate1, '%m')),
                           day = as.integer(strftime(startdate1, '%d')),
                           hour = 00, min = 00, sec = 00, tz = 'UTC')
  enddate <- ISOdatetime(year = as.integer(strftime(enddate1, '%Y')),
                         month = as.integer(strftime(enddate1, '%m')),
                         day = as.integer(strftime(enddate1, '%d')) + 1,
                         hour = 00, min = 00, sec = 00, tz = 'UTC')
  
  ## Specify the parameter(s) of interest
  parameters <- list("wind_speed_10m:ms", # "wind_dir_10m:d",
                     "wind_gusts_10m_24h:ms", # "wind_gusts_10m_1h:ms",
                     "t_2m:C", "t_max_2m_24h:C", "t_min_2m_24h:C", 
                     "msl_pressure:hPa",
                     "precip_24h:mm", # "precip_1h:mm",
                     "relative_humidity_1000hPa:p",
                     # "sunrise:sql","sunset:sql",
                     "dew_point_2m:C",
                     "uv:idx"
  )
  
  ## Call the MeteomaticsRConnector::query_time_series() function
  
  start <- seq(1,length(coordinates),300)
  end <- start-1; end <- end[-1]; end <- c(end, length(coordinates)); end <- unique(end)
  
  wdataList <- list()
  for(i in 1:length(start)){
    iCoords <- coordinates[start[i]:end[i]]
    wdataList[[i]] <- MeteomaticsRConnector::query_time_series(iCoords, startdate, enddate, interval, parameters,
                                                               username, password)
  }
  wdata0 <- do.call(rbind, wdataList)
  # wdata0 <- MeteomaticsRConnector::query_time_series(coordinates, startdate, enddate, interval, parameters, username, password)
  colnames(wdata0) <- gsub(":","",colnames(wdata0))
  
  #######################
  
  myDates <- lapply(strsplit(as.character(wdata0$validdate), split=" "), function(x){x[1]})
  myTimes <- lapply(strsplit(as.character(wdata0$validdate), split=" "), function(x){x[2]})
  myDatesDf <- as.data.frame(do.call(rbind, lapply(myDates, function(x){strsplit(x,"-")[[1]]}) ) )
  myHoursDf <- as.data.frame(do.call(rbind, lapply(myTimes, function(x){strsplit(x,":")[[1]]}) ) )
  wdataDf <- cbind(myDatesDf,myHoursDf,wdata0)
  colnames(wdataDf)[1:6] <- c("year","month","day","hour","minute","second")
  ## come up with aggregated data for each fieldinst
  wdata2 <- merge(wdataDf,mydata, by.x = c("lat","lon"), by.y = c("latitude","longitude"), all.x = TRUE)
  myFormula <- paste("cbind(", paste(gsub(":","",unlist(parameters)), collapse = ","),")~fieldinst")
  aggregatedWdata <- aggregate(as.formula(myFormula), data=wdata2, FUN=function(x){mean(x, na.rm=TRUE)})
  ## merge environmental summaries to metadata
  metadataFieldinst <- merge(mydata,aggregatedWdata, by="fieldinst", all.x = TRUE )
  # growing degree days
  
  xmax <- ifelse(wdataDf$t_max_2m_24hC > 30, 30, wdataDf$t_max_2m_24hC)
  xmin <- ifelse(wdataDf$t_min_2m_24hC < 10, 10, wdataDf$t_min_2m_24hC)
  gdd <- ((xmax + xmin)/2) - 10
  wdataDf$gdd <- gdd * 30 # average degree days times 30 days
  wdataDf$cumgdd <- cumsum(wdataDf$gdd)
  
  
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
