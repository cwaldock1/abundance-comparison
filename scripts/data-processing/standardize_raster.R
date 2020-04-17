
# function to standaridse rasters to a common mask (values of NA and 1)

standardize_raster <- function(masterMask, targetRaster){
  
  # check if extent and resolution are the same
  if(extent(masterMask) != extent(targetRaster)){
    targetRaster <- crop(targetRaster, bbox(masterMask))
    }
  
  if(sum(res(masterMask) != res(targetRaster))){
    targetRaster <- resample(targetRaster, masterMask)
  }
  
  manipulatedRaster <- targetRaster
  
  # Turns all potential projected area that does not have a projection to -999
  
  manipulatedRaster[masterMask[] == 1 & is.na(manipulatedRaster[])] <- -9999
  
  # Create objects to fill through IDW
  
  dt = rasterToPoints(manipulatedRaster)
  if(sum(dt[,3]==-9999)!=0){
  knowndt <- data.frame(dt[-which(dt[,3]==-9999),])
  unknowndt <- data.frame(dt[which(dt[,3]==-9999),])
  coordinates(knowndt) <- ~x+y
  unknowndt[,3] <- NA
  coordinates(unknowndt) <- ~x+y
  names(unknowndt) <- 'value'
  names(knowndt) <- 'value'
  
  # Perform inverse distance weighting
  
  idwmodel <- idw(value~1, knowndt, unknowndt,
                  maxdist = 10, idp = 5)
  
  # Fill the target raster that are -999 (defined above)
  
  manipulatedRaster[manipulatedRaster[]==-9999] <- idwmodel$var1.pred
  
  }
  
  # Switch off all areas that are not target in the mask (target = 1)
  
  manipulatedRaster[masterMask[]!=1 | is.na(masterMask[])] <- NA
  
  return(manipulatedRaster)
  
}
