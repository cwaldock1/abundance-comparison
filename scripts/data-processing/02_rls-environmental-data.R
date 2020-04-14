# script to match up rls sites with environmental rasters on a standardised grid 

# custom function to fill NAs with nearest neighbour ----
# This set of functions  replaces all NAs that are found in the raster with a nearest neighbour. Very slow. 
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, # edit the [1,]
        FUN = function(xy) r@data@values[raster::which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}

FillSiteNAs <- function(StackLayer, # Environmental stack layer
                        RLS_Sites   # Coordiantes as a spatial points database.  
){
  
  # Do first extraction that finds the NA values inside the subset (extract(StackLayer, RLS_Sites))
  # Find the coordinates at the values
  xy <- RLS_Sites@coords[is.na(extract(StackLayer, RLS_Sites)),]
  
  # Match to the nearest raster cell using a pre-defined function 
  NA_coords <- force(sample_raster_NA(r = StackLayer, xy = xy))
  StackLayer[cellFromXY(StackLayer, as.matrix(xy))] <- NA_coords
  return(StackLayer)
  
}

# load packages ----
# pca methods requires special install
if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")
BiocManager::install("pcaMethods")}

# load/install other packages for script
lib_vect <- c('tidyverse', 'rgdal', 'raster', 'maptools', 'summarytools')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load in rls sites for matching ----

# projection
Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# load data
load(file = 'data/rls_abun_modelling_data_v2.RData')
rls_xy     <- data.frame(SiteLongitude = rls_abun$SiteLongitude, SiteLatitude = rls_abun$SiteLatitude) %>% unique()
rls_points <- SpatialPoints(cbind(rls_abun$SiteLongitude, rls_abun$SiteLatitude) %>% unique(), proj4string = crs(Proj))

# load in and standardise environmental layers ----

# species distribution models will be build from the best resolution data available rather than comforming to a single grid
# this is because we are not evaluating the spatial distribution of the species distributions but instead just modelling the parameters

# depth ----
# extract depth
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/ 
# new data for 2019. 
depth <- raster("raw-data/GEBCO_2019/GEBCO_2019.nc", CRS = Proj)
# plot(depth)
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
  print(i)
  bio_orc_stack[[i]] <- force(FillSiteNAs(StackLayer = bio_orc_stack[[i]], RLS_Sites = rls_points))
}

# extract
bio_orc_ext <- data.frame(extract(bio_orc_stack, rls_points))

# check for NAs
sum(is.na(bio_orc_ext))

# check distributions
view(dfSummary(bio_orc_ext))

# estimate principal components using robust pca
bio_orc_ext_scaled <- scale(bio_orc_ext, scale = T, center = T)
resRobSvd <- pcaMethods::pca(bio_orc_ext_scaled, 
                             method = "robustPca", 
                             center = F, 
                             scale = NULL, 
                             nPcs = dim(bio_orc_ext)[2])

resRobSvd # 6 loadings appear important and explain 75% variation, and all more than 1%
Loadings <- resRobSvd@loadings[,1:5]
Loadings <- reshape2::melt(Loadings)

pdf(file = 'figures/robustPCA-env-covariates.pdf', width = 7.5, height = 7.5)
ggplot(Loadings) + 
  geom_bar(aes(x = Var1, y = value, fill = value), stat = 'identity') +
  theme_bw() + 
  facet_wrap(~Var2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Loading weight') + xlab(NULL) 
dev.off()

#slplot(resRobSvd, scoresLoadings = c(T,T))
biplot(resRobSvd)
plotPcs(resRobSvd)

# get values to model with
bio_orc_ext$robPCA_1 <- as.numeric(resRobSvd@scores[,1])
bio_orc_ext$robPCA_2 <- as.numeric(resRobSvd@scores[,2]) 
bio_orc_ext$robPCA_3 <- as.numeric(resRobSvd@scores[,3])
bio_orc_ext$robPCA_4 <- as.numeric(resRobSvd@scores[,4]) 
bio_orc_ext$robPCA_5 <- as.numeric(resRobSvd@scores[,5])
#bio_orc_ext$robPCA_6 <- as.numeric(resRobSvd@scores[,6]) 

# bind to rls
rls_xy <- cbind(rls_xy, bio_orc_ext)

# msec ----

# read in various layers
human_pop  <- list.files(path='raw-data/MSECData_Yeager2017',pattern='nc$',full.names = T)[c(3)]
msec  <- list.files(path='raw-data/MSECData_Yeager2017',pattern='nc$',full.names = T)[c(13,14)]

# stack and rotate layers
preds_MSEChuman <- stack(human_pop)
preds_MSEChuman <- preds_MSEChuman$X2015
preds_MSEC      <- stack(msec)
preds_MSEChuman <- rotate(preds_MSEChuman)
preds_MSEC      <- rotate(preds_MSEC)
crs(preds_MSEChuman)     <- Proj
crs(preds_MSEC)          <- Proj

# stack for simplicity
MSEC <- stack(preds_MSEChuman, preds_MSEC)

# crop
MSEC <- crop(MSEC, bbox(rls_points) + c(-1, -1, 1, 1))

# rename
names(MSEC) <- c('human_pop_2015_50km', 'reef_area_200km', 'wave_energy_mean')

# check distributions
#hist(log10(MSEC$human_pop_2015_50km+1))
#hist(log10(MSEC$reef_area_200km + 1))
#hist(log10(MSEC$wave_energy_mean + 1))

# change values
MSEC$human_pop_2015_50km <- calc(MSEC$human_pop_2015_50km, function(x) scale(log10(x+1)))
MSEC$reef_area_200km     <- calc(MSEC$reef_area_200km, function(x) scale(log10(x+1)))
MSEC$wave_energy_mean    <- calc(MSEC$wave_energy_mean, function(x) scale(log10(x+1)))

# perform extractions
for(i in 1:length(names(MSEC))){
  print(i)
  MSEC[[i]] <- force(FillSiteNAs(StackLayer = MSEC[[i]], RLS_Sites = rls_points))
}

# extract points
MSEC_ext <- data.frame(extract(MSEC, rls_points))

# check for NAs
table(is.na(MSEC_ext))

# bind to rls_xy object
rls_xy <- cbind(rls_xy, MSEC_ext)

# check for any NAs
table(is.na(rls_xy))

# transform, scale and centre all non-pca variables ---- 
rls_xy$Depth_GEBCO[which(rls_xy$Depth_GEBCO > 0)] <- 0 # lose the differentiation between land and sea.
hist(log10(round(abs(rls_xy$Depth_GEBCO+1))))
rls_xy$Depth_GEBCO_transformed <- log10(round(abs(rls_xy$Depth_GEBCO)+1))

# check for any NAs
table(is.na(rls_xy))

# save rls_xy ----
save(rls_xy, file = 'data/rls_covariates.RData')





