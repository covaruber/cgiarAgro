nasapowerExtractSingleField <- function(phenoDTfile= NULL, verbose=FALSE, interval="daily"){
  
  id <- paste( paste("wea",cgiarPIPE::idGenerator(5,5),sep=""), "nasapower", sep = "_")
  type <- "wea"
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(length(grep("clp",phenoDTfile$id)) == 0){stop("Index can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  
  ############################
  # loading the dataset
  if(is.null(phenoDTfile$metadataFieldinst)){stop("There's no metadata for this file. Likely belongs to an older version of the cgiarPIPE package. Please match the columns of your data again.", call. = FALSE)}
  mydata <- phenoDTfile$metadataFieldinst # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydata$latitude <- round(mydata$latitude,2)
  mydata$longitude <- round(mydata$longitude,2)
  
  ###########################
  ## extract data
  mydata$plantingDate <- as.Date(mydata$plantingDate)
  mydata$harvestingDate <- as.Date(mydata$harvestingDate)
  
  ## Specify the parameter(s) of interest
  parameters <- c(
    "RH2M", # The ratio of actual partial pressure of water vapor to the partial pressure at saturation, expressed in percent
    "U2M", # The estimate of the eastward wind average speed for winds blowing 2 meters above the surface of the earth
    "T2M", # The average air (dry bulb) temperature at 2 meters above the surface of the earth
    "T2M_MAX","T2M_MIN", 
    "WS2M", # "The average of wind speed at 2 meters above the surface of the earth."
    "PRECTOTCORR", # The bias corrected average of total precipitation at the surface of the earth in water mass
    "PS", # The average of surface pressure at the surface of the earth.
    "GWETTOP", # The percent of soil moisture a value of 0 indicates a completely water-free soil and a value of 1 indicates a completely saturated soil
    "CLOUD_AMT", # "The average percent of cloud amount during the temporal period."
    "PW", # The total atmospheric water vapor contained in a vertical column of the atmosphere.
    "ALLSKY_SFC_LW_UP", # The upward thermal infrared irradiance under all sky conditions.
    "AIRMASS" # "Air Mass according to Perez et al. 1990."
    
    # "QV2M", # The ratio of the mass of water vapor to the total mass of air at 2 meters (g water/kg total air).
    # "EVLAND", # The evaporation over land at the surface of the earth.
    # "T2MDEW", # The dew/frost point temperature at 2 meters above the surface of the earth.
    # "PRECSNOLAND", # The snow precipitation only over land at the surface of the earth.
    # "ALLSKY_SFC_UVB", #  The ultraviolet B (UVB 280nm-315nm) irradiance under all sky conditions.
    # "SG_DAY_HOURS", # The number of hours where solar irradiance is present;
  )
  # extract weather data
  wdataList <- list()
  for(i in 1:nrow(mydata)){
    print(i)
    wdataList[[i]] <- nasapower::get_power(
      community = "ag",
      lonlat = c(mydata$longitude[i], mydata$latitude[i]),
      pars = parameters,
      dates = c(mydata$plantingDate[i], mydata$harvestingDate[i]),
      temporal_api = interval
    )
    wdataList[[i]]$fieldinst <- mydata$fieldinst[i]
  }
  wdataDf <- do.call(rbind, wdataList)
  colnames(wdataDf) <- gsub(":","",colnames(wdataDf))
  colnames(wdataDf)[1:7] <- c("lat","lon","year","month","day","dayOfYear", "date")
  missing <- apply(wdataDf,2,function(x){length(which(is.na(x)))/length(x)})
  keep <- which(missing < .8)
  parameters <- intersect(parameters,names(keep))
  wdataDf <- wdataDf[,keep]
  ## come up with aggregated data for each fieldinst
  mydata <- mydata[,setdiff(colnames(mydata), setdiff( colnames(wdataDf), "fieldinst"))]
  wdata2 <- merge(wdataDf,mydata, by.x = c("lat","lon","fieldinst"), by.y = c("latitude","longitude","fieldinst"), all.x = TRUE)
  myFormula <- paste("cbind(", paste(gsub(":","",unlist(parameters)), collapse = ","),")~fieldinst")
  aggregatedWdata <- aggregate(as.formula(myFormula), data=wdata2, FUN=function(x){mean(x, na.rm=TRUE)})
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
                 cleaned=wdata2, outliers=NA, desire=NA, id=id, idOriginal=phenoDTfile$idOriginal,
                 metadataFieldinst=metadataFieldinst
  )
  return(result)
  
}
