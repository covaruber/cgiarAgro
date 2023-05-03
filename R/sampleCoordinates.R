sampleCoordinates <- function(  lat=NULL, lon=NULL, nSample=100,
                              startDate=NULL, endDate=NULL){
  
  if(is.null(lat)){stop("Please provide the latitude vector", call. = FALSE)}
  if(is.null(lon)){stop("Please provide the longitude vector", call. = FALSE)}
  if(length(lat) != length(lon)){stop("Latitude and longitude vector should have equal length", call. = FALSE)}
  if(is.null(startDate)){startDate = Sys.Date()-(365*4)}
  if(is.null(endDate)){endDate = startDate-10}
  
  lat <- sort(lat)
  lon <- sort(lon)
  latSeqs <- seq( min(lat), max(lat), abs(min(lat)-max(lat))/sqrt(nSample) )
  lonSeqs <- seq( min(lon), max(lon), abs(min(lon)-max(lon))/sqrt(nSample) )
  coordGrid <- expand.grid(latSeqs,lonSeqs)
  
  #######
  
  metadata <- data.frame(fieldinst=paste0("L",1:nrow(coordGrid)), timepoint=1, location="A", country="B", 
                         latitude= coordGrid$Var1 ,longitude= coordGrid$Var2, altitude=1,
                         plantingDate=startDate, harvestingDate=endDate
  )
  
  res <- list(metadataFieldinst=metadata, id="clpw")
  return(res)
}