nasapowerExtractRegion <- function(west= NULL, east=NULL, south=NULL, north=NULL, space=4.5,
                                   initialDateVal=NULL, finalDateVal=NULL,
                                   verbose=FALSE, interval="daily"){

  id <- paste( paste("wea",cgiarPIPE::idGenerator(5,5),sep=""), "nasapowerRegion", sep = "_")
  type <- "wea"
  if(is.null(west) | is.null(east) | is.null(south) | is.null(north)){stop("Please provide the coordinates for all points", call. = FALSE)}
  if(space > 4.5){stop("grids cannot be larger than 4.5 x 4.5 degrees", call. = FALSE)}
  ############################
  # loading the dataset

  x1 <- sort(c(west,east), decreasing = FALSE)
  xs <- seq(x1[1],x1[2],space) # west > east breakpoints (0.5 degrees)
  y1 <- sort(c(south,north), decreasing = FALSE)
  ys <- seq(y1[1],y1[2],space) # south > north breakpoints (0.5 degrees)

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

  wdataList <- list(); counter <- 1
  for(i in 1:(length(xs)-1)){
    for(j in 1:(length(ys)-1)){
      print(counter)
      wdataList[[counter]] <- nasapower::get_power(
        community = "ag",
        lonlat = c(xs[i], ys[j],  xs[i+1], ys[j+1]), # west, south, east, north
        pars = parameters,
        dates = c(initialDateVal, finalDateVal),
        temporal_api = interval
      )
      counter <- counter + 1
    }
  }
  wdataDf <- do.call(rbind, wdataList)

  colnames(wdataDf) <- gsub(":","",colnames(wdataDf))
  colnames(wdataDf)[1:7] <- c("latitude","longitude","year","month","day","dayOfYear", "date")
  missing <- apply(wdataDf,2,function(x){length(which(is.na(x)))/length(x)})
  keep <- which(missing < .8)
  parameters <- intersect(parameters,names(keep))
  wdataDf <- wdataDf[,keep]

  ## come up with aggregated data for each fieldinst
  mydata <- unique(wdataDf[,c("latitude","longitude")]);# colnames(mydata) <- c("latitude","longitude")
  mydata$fieldinst <- paste0("L",1:nrow(mydata))

  wdata2 <- merge(wdataDf,mydata, by.x = c("latitude","longitude"), by.y = c("latitude","longitude"), all.x = TRUE)
  myFormula <- paste("cbind(", paste(gsub(":","",unlist(parameters)), collapse = ","),")~fieldinst")
  aggregatedWdata <- aggregate(as.formula(myFormula), data=wdata2, FUN=function(x){mean(x, na.rm=TRUE)})
  ## merge environmental summaries to metadata
  metadataFieldinst <- merge(mydata,aggregatedWdata, by="fieldinst", all.x = TRUE )
  ## write the parameters to the parameter database
  db.params <- data.frame(
    analysisId	= id,
    analysisType =	type,
    fieldbooks	= NA,
    phenoDataFile ="none",
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
                 cleaned=wdata2, outliers=NA, desire=NA, id=id, idOriginal="none",
                 metadataFieldinst=metadataFieldinst
  )
  return(result)

}
