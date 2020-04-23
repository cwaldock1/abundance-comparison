### script to create a standardised global grid 

# libraries ----
lib_vect <- c('tidyverse', 'raster', 'rgdal', 'gdalUtils', 'rgeos')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# empty raster ----

# create empty extent object
global_grid <- raster(extent(c(-180, 180, -90, 90)), 
                      crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0', 
                      res = 0.025)

# set values to 1
global_grid[] <- 1

# test plot
plot(global_grid)

# read in RLS bounding box ----

# projection
Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

load(file = 'data/rls_abun_modelling_data_v2.RData')
rls_xy     <- data.frame(SiteLongitude = rls_abun$SiteLongitude, SiteLatitude = rls_abun$SiteLatitude) %>% unique()
rls_points <- SpatialPoints(cbind(rls_abun$SiteLongitude, rls_abun$SiteLatitude) %>% unique(), proj4string = crs(Proj))

# depth ----
depth <- raster('/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019.nc')
plot(depth)

# tuncate values
depth_clamp <- clamp(depth, lower = -200, upper = 50, useValues = FALSE)
# save the depth clamp
writeRaster(depth_clamp, filename = '/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019_depth_-200_50.tif')

# apply resampling from gdal as faster than raster
src_dataset <- system.file("/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019_depth_-200_50.tif", package="gdalUtils")

gdalwarp(srcfile = '/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019_depth_-200_50.tif',
         dstfile = '/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019_depth_-200_50_resampled_0.05.tif',
         t_srs ='+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0',
         tr    = c(0.05, 0.05), # resampling resolution
         r     = 'near', # resampling method
         output_Raster=TRUE,
         overwrite=F,
         verbose=TRUE)

plot(raster('/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019_depth_-200_50_resampled_0.05.tif'))

# coastlines ----

# read in the coastlines from world hires GSHHS_l_L1_disolve.shp
coastlines <- readOGR(dsn = '/Volumes/RF-env-data/reef-futures/env-data/coastal-shape-files/GSHHS_l_L1_disolve.shp')

coastlines_buffered <- gBuffer(coastlines, width = 0.25)

# coral reefs and buffer ----

# Read in reef spatial polygon dataset
ReefPoly <- readOGR(dsn = path.expand('/Volumes/RF-env-data/reef-futures/env-data/WCMC_ReefData/01_Data/WCMC008_CoralReef2010_Py_v3.shp'))

# Create a subset of smaller polygons and see if these get stuck in buffer too. 

Areas_Poly <- do.call(rbind, lapply(ReefPoly@polygons, function(x){x@area}))
SubsetOut <- round(0.9*length(Areas_Poly)):length(Areas_Poly)

# Function to disaggregate large polygons as causes issues with memory. 

ConvertPolygonToBuffer <- function(Polygon, ID, i){
  Polygon_disaggregated <- sp::disaggregate(Polygon)
  Polygon_disaggregated_simple <- gSimplify(Polygon_disaggregated, tol=0.1, topologyPreserve=T)
  Polygon_disaggregated_simple_Buffer <- gBuffer(Polygon_disaggregated_simple, byid = T, width = 0.25)
  Polygon_disaggregated_simple_Buffer_aggregate <- aggregate(Polygon_disaggregated_simple_Buffer)
  Polygon_disaggregated_simple_Buffer_aggregate <- SpatialPolygonsDataFrame(Polygon_disaggregated_simple_Buffer_aggregate, data = data.frame(ID=ID[i]))
  return(Polygon_disaggregated_simple_Buffer_aggregate)
}

# Run buffer function over all

result <- list()
SmallPolygonIDs <- order(Areas_Poly)[-SubsetOut]
for(i in 1:length(ReefPoly)){
  result[[i]] <- ConvertPolygonToBuffer(ReefPoly[i,], ID = SmallPolygonIDs, i=i); 
  print(i/length(ReefPoly))
  }
result <- do.call("rbind", result)

# Run buffer function over each polygon (and feature set). 

result_BIG <- list()
BigPolygonIDs <- order(Areas_Poly)[SubsetOut]
for(i in 1:length(BigPolygonIDs)){
  result_BIG[[i]] <- ConvertPolygonToBuffer(ReefPoly[BigPolygonIDs[i],], ID = BigPolygonIDs, i=i); print(i/length(BigPolygonIDs))
  }
result_BIG_bind <- do.call("rbind", result_BIG)

# comile outputs

reef_polygonsAll <- rbind(result_BIG_bind, result)
reef_polygonsAll_V1 <- aggregate(reef_polygonsAll)
reef_polygonsAll_V1 <- as(reef_polygonsAll_V1, 'SpatialPolygonsDataFrame')

# save to folder on hard-drive

writeOGR(reef_polygonsAll_V1, "/Volumes/RF-env-data/reef-futures/env-data/WCMC_ReefData/BufferedPolygons", "Reef_AsBufferedPolygon_0.25deg", driver="ESRI Shapefile")

# read in for faster fastness
coral <- readOGR(dsn = '/Volumes/RF-env-data/reef-futures/env-data/WCMC_ReefData/BufferedPolygons/Reef_AsBufferedPolygon_0.25deg.shp')

# created masked global raster ----

# crop coastlines

coastlines          <- crop(coastlines, bbox(rls_points) + c(-1, -1, 1, 1))
coastlines_buffered <- crop(coastlines_buffered, bbox(rls_points) + c(-1, -1, 1, 1))

# crop coral 

coral          <- crop(coral, bbox(rls_points) + c(-1, -1, 1, 1))

# create rls sites buffer
crs(rls_points) <- NULL
rls_polyons <- buffer(rls_points, 25000)

# create union of buffered polygons

union_reef <- gUnion(gUnion(coastlines_buffered, coral), rls_polyons)
plot(union_reef)

# mask land

union_reef_buffered_mask <- mask(aus_grid, union_reef)
final_reef_area_aus      <- mask(union_reef_buffered_mask, coastlines, inverse = T)

plot(final_reef_area_aus)
points(rls_points)

# edit mask errors manually ----


# save the environmental mask layers ---- 
writeRaster(final_reef_area_aus, 'data/rls-aus-grid.nc', overwrite = T)
final_reef_area_aus <- raster('')

