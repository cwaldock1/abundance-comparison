# script to match up rls sites with environmental rasters on a standardised grid 

# custom function to fill NAs with nearest neighbour ----
# This set of functions  replaces all NAs that are found in the raster with a nearest neighbour. Very slow. 
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[raster::which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}

FillSiteNAs <- function(StackLayer, # Environmental stack layer
                        RLS_Sites   # Coordiantes as a spatial points database.  
){
  
  # Do first extraction that finds the NA values inside the subset (extract(StackLayer, RLS_Sites))
  # Find the coordinates at the values
  xy <- RLS_Sites@coords[is.na(extract(StackLayer, RLS_Sites)),]
  
  # Match to the nearest raster cell using a pre-defined function 
  NA_coords <- sample_raster_NA(r = StackLayer, xy = xy)
  StackLayer[cellFromXY(StackLayer, as.matrix(xy))] <- NA_coords
  return(StackLayer)
  
}

# load packages ----
# pca methods requires special install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")

# load/install other packages for script
lib_vect <- c('tidyverse', 'rgdal', 'raster', 'maptools', 'summarytools')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load in rls sites for matching ----

# projection
Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# load data
load(file = 'data/rls_abun_modelling_data.RData')
rls_xy     <- data.frame(SiteCode = rls_abun$SiteCode, SiteLongitude = rls_abun$SiteLongitude, SiteLatitude = rls_abun$SiteLatitude)
rls_points <- SpatialPoints(cbind(rls_abun$SiteLongitude, rls_abun$SiteLatitude), proj4string = crs(Proj))

# load in and standardise environmental layers ----


# species distribution models will be build from the best resolution data available rather than comforming to a single grid
# this is because we are not evaluating the spatial distribution of the species distributions but instead just modelling the parameters

# depth ----
# extract depth
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/ 
# new data for 2019. 
depth <- raster("raw-data/GEBCO_2019/GEBCO_2019.nc", CRS = Proj)
plot(depth)
rls_xy$Depth_GEBCO <- raster::extract(depth, rls_points)

# bio-oracle-v2 ----
# stack and extract bio-oracle v2 
bio_orc <- list.files('raw-data/bio-oracle-v2', full.names = T)
bio_orc_stack <- stack(bio_orc)
crs(bio_orc_stack) <- Proj
bio_orc_stack <- crop(bio_orc_stack, bbox(rls_points) + c(-1, -1, 1, 1))

# rename values 
names(bio_orc_stack) <- c('calcite.mean', 'current.velocity.mean', 'o2.min', 
                          'iron.mean', 'nitrate.mean', 'pH', 'phosphate', 
                          'PP.mean', 'salinity.mean', 'silicate.mean', 
                          'sst.max', 'sst.min', 'sst.mean')

# plot distribution of variables
hist(log10(bio_orc_stack$calcite.mean))
hist(log10(bio_orc_stack$current.velocity.mean))
hist(bio_orc_stack$o2.min)
hist(log10(bio_orc_stack$iron.mean))
hist(log10(bio_orc_stack$nitrate.mean))
hist(log(bio_orc_stack$pH))
hist(bio_orc_stack$phosphate)
hist(log10(bio_orc_stack$PP.mean))
hist(bio_orc_stack$salinity.mean)
hist(log(bio_orc_stack$silicate.mean))
hist(bio_orc_stack$sst.max)
hist(bio_orc_stack$sst.min)
hist(bio_orc_stack$sst.mean)

# transform variables
bio_orc_stack$calcite.mean <- calc(bio_orc_stack$calcite.mean, log10)
bio_orc_stack$current.velocity.mean <- calc(bio_orc_stack$current.velocity.mean, log10)
bio_orc_stack$iron.mean <- calc(bio_orc_stack$iron.mean, log10)
bio_orc_stack$nitrate.mean <- calc(bio_orc_stack$nitrate.mean, log10)
bio_orc_stack$pH <- calc(bio_orc_stack$pH, log10)
bio_orc_stack$PP.mean <- calc(bio_orc_stack$PP.mean, log10)
bio_orc_stack$silicate.mean <- calc(bio_orc_stack$silicate.mean, log)

# plot all bio-orc distributions
dev.off(); plot(bio_orc_stack)

# find NAs
for(i in 1:length(names(bio_orc_stack))){
  
  bio_orc_stack[[i]] <- FillSiteNAs(StackLayer = bio_orc_stack[[i]], RLS_Sites = rls_points)
  
}

# extract
bio_orc_ext <- data.frame(extract(bio_orc_stack, rls_points))

# check for NAs
sum(is.na(bio_orc_ext))

# check distributions
view(dfSummary(bio_orc_ext))

# estimate principal components using robust pca
bio_orc_ext <- scale(bio_orc_ext, scale = T, center = T)
resRobSvd <- pcaMethods::pca(bio_orc_ext, 
                             method = "robustPca", 
                             center = F, 
                             scale = NULL, nPcs = dim(bio_orc_ext)[2])

resRobSvd # 6 loadings appear important
Loadings      <- resRobSvd@loadings[,1:6]
Loadings <- reshape2::melt(Loadings)

ggplot(Loadings) + 
  geom_bar(aes(x = Var1, y = value, fill = value), stat = 'identity') +
  theme_bw() + 
  facet_wrap(~Var2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 1') + xlab(NULL) 

slplot(resRobSvd, scoresLoadings = c(T,T))
biplot(resRobSvd)
plotPcs(resRobSvd)

# get values to model with
bio_orc_ext$robPCA_1 <- as.numeric(resRobSvd@scores[,1])
bio_orc_ext$robPCA_2 <- as.numeric(resRobSvd@scores[,2]) 
bio_orc_ext$robPCA_3 <- as.numeric(resRobSvd@scores[,3])
bio_orc_ext$robPCA_4 <- as.numeric(resRobSvd@scores[,4]) 
bio_orc_ext$robPCA_5 <- as.numeric(resRobSvd@scores[,5])
bio_orc_ext$robPCA_6 <- as.numeric(resRobSvd@scores[,6]) 

# bind to rls
rls_xy <- cbind(rls_xy, bio_orc_ext)


