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
  xy <- RLS_Sites@coords[is.na(raster::extract(StackLayer, RLS_Sites)),]
  
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

# load data into xy coords

rls_xy <- do.call(rbind, lapply(list.files('data/rls_all_basic', full.names = T), function(x){
  # read in all the individual files
  ddata <- readRDS(x)
  ddata_all <- do.call(rbind, ddata)
  data.frame(SiteLongitude = ddata_all$SiteLongitude, SiteLatitude = ddata_all$SiteLatitude) %>% 
    unique()})) %>% 
  unique()

# convert to points

rls_points <- SpatialPoints(rls_xy, proj4string = crs(Proj))

# load in and standardise environmental layers ----

# species distribution models will be build from the best resolution data available rather than comforming to a single grid
# this is because we are not evaluating the spatial distribution of the species distributions but instead just modelling the parameters

# depth ----
# extract depth
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/ 
# new data for 2019. 
depth <- raster("/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019.nc", CRS = Proj)
# plot(depth)
rls_xy$Depth_GEBCO <- raster::extract(depth, rls_points)


# read in rasters developed for reef-futures ----


# read in an stack the rasters

rf_covs <- stack(list.files('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/raw-data/marine_rasters_reef_futures', 
           recursive = T, full.names = T))

# rename the rasters

names(rf_covs) <- gsub('.nc', '', lapply(str_split(list.files('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/raw-data/marine_rasters_reef_futures', 
                             recursive = T, full.names = F), '/'), function(x){x[[2]]}))

# plot to test

#plot(rf_covs)

# remove temperature variables and handle seperately for cross validation

rf_covs <- rf_covs[[-grep('sst', names(rf_covs))]]

# ensure projection is standard

crs(rf_covs) <- Proj

# crop to area of interest

rf_covs <- crop(rf_covs, bbox(rls_points) + c(-1, -1, 1, 1))

# plot distribution of all varibles
view(summarytools::dfSummary(rasterToPoints(rf_covs)))

hist(log10(rf_covs$chl_max))
hist(log10(rf_covs$chl_mean))
hist(log10(rf_covs$chl_min))

hist(log10(rf_covs$dhw_max))
hist(log10(rf_covs$dhw_mean))
hist(log10(rf_covs$dhw_min+1))

hist(log10(rf_covs$NPP_max))
hist(log10(rf_covs$NPP_mean))
hist(log10(rf_covs$NPP_min))

hist(rf_covs$pH_max)
hist(rf_covs$pH_mean)
hist(rf_covs$pH_min)

hist(rf_covs$SSS_max)
hist(rf_covs$SSS_mean)
hist(rf_covs$SSS_min)


# transform variables

log10_vars <- c('chl_max', 'chl_mean', 'chl_min', 
                'dhw_max', 'dhw_mean', 'dhw_min',
                'NPP_max', 'NPP_mean', 'NPP_min')

for(i in 1:length(log10_vars)){
  rf_covs[[log10_vars[i]]][] <- log10(rf_covs[[log10_vars[i]]]+1)[]
}


# find NAs (shouldn't be any because of mask standardisation)

for(i in 1:length(names(rf_covs))){
  print(i)
  rf_covs[[i]] <- force(FillSiteNAs(StackLayer = rf_covs[[i]], RLS_Sites = rls_points))
}

# extract
rf_covs_ext <- data.frame(raster::extract(rf_covs, rls_points))

# check for NAs

table(is.na(rf_covs_ext))

# check distributions
view(dfSummary(rf_covs_ext))

# estimate principal components using robust pca
resRobSvd <- pcaMethods::pca(rf_covs_ext, 
                             method = "robustPca", 
                             center = T,
                             scale = 'uv', 
                             nPcs = dim(rf_covs_ext)[2])

saveRDS(resRobSvd, file = 'data/rls_resRobSvd.rds')

resRobSvd # 3 loadings appear important and explain ~75% variation, and all more than 1%

# robustPca calculated PCA
# Importance of component(s):
#               PC1    PC2     PC3      PC4     PC5      PC6    PC7     PC8     PC9    PC10     PC11     PC12     PC13     PC14    PC15
# R2            0.4789 0.2163 0.08279 -0.08121 -0.1416 -0.03306 0.0508 -0.0568 -0.0405 -0.2314 -0.05585 -0.26783 -0.01642 -0.07162 -0.1746
# Cumulative R2 0.4789 0.6952 0.77796  0.69675  0.5552  0.52210 0.5729  0.5161  0.4756  0.2442  0.18837 -0.07946 -0.09588 -0.16750 -0.3421
# 15 	Variables
# 3790 	Samples
# 0 	NAs ( 0 %)
# 15 	Calculated component(s)
# Data was mean centered before running PCA 
# Data was scaled before running PCA 
# Scores structure:
#   [1] 3790   15
# Loadings structure:
#   [1] 15 15

Loadings <- resRobSvd@loadings[,1:3]
Loadings <- reshape2::melt(Loadings)

pdf(file = 'figures/robustPCA-rls_env-covariates.pdf', width = 7.5, height = 7.5)
ggplot(Loadings) + 
  geom_bar(aes(x = Var1, y = value, fill = value), stat = 'identity') +
  theme_bw() + 
  facet_wrap(~Var2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'none', aspect.ratio = 1, 
        panel.grid = element_blank())+
  scale_fill_viridis_c() + 
  ylab('Loading weight') + xlab(NULL) 
dev.off()

biplot(resRobSvd)

# get values to model with
rf_covs_ext$robPCA_1 <- as.numeric(resRobSvd@scores[,1])
rf_covs_ext$robPCA_2 <- as.numeric(resRobSvd@scores[,2]) 
rf_covs_ext$robPCA_3 <- as.numeric(resRobSvd@scores[,3])

# bind to rls
rls_xy <- cbind(rls_xy, rf_covs_ext)

# sst (doing this seperately because will be the axis of change) ----

sst <- stack(list.files('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/raw-data/marine_rasters_reef_futures/sst', 
                            recursive = T, full.names = T))

# rename the rasters

names(sst) <- gsub('.nc', '', lapply(str_split(list.files('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/raw-data/marine_rasters_reef_futures/sst', 
                                                              recursive = T, full.names = F), '/'), function(x){x[[1]]}))

# double check projection is correct

crs(sst) <- Proj

# crop to area of interest

sst <- crop(sst, bbox(rls_points) + c(-1, -1, 1, 1))

# double check for NAs

for(i in 1:length(names(sst))){
  print(i)
  sst[[i]] <- force(FillSiteNAs(StackLayer = sst[[i]], RLS_Sites = rls_points))
}

# extract

sst_ext <- data.frame(raster::extract(sst, rls_points))

# check for NAs

table(is.na(sst_ext))

# check distributions

view(dfSummary(sst_ext))

# bind to rls

rls_xy <- cbind(rls_xy, sst_ext)


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
MSEC$reef_area_200km     <- calc(MSEC$reef_area_200km,     function(x) scale(log10(x+1)))
MSEC$wave_energy_mean    <- calc(MSEC$wave_energy_mean,    function(x) scale(log10(x+1)))

# perform extractions
for(i in 1:length(names(MSEC))){
  print(i)
  MSEC[[i]] <- force(FillSiteNAs(StackLayer = MSEC[[i]], RLS_Sites = rls_points))
}

# extract points
MSEC_ext <- data.frame(raster::extract(MSEC, rls_points))

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


# create RLS raster stack of env. variables for projection ----

# read in robust PCA model

resRobSvd <- readRDS(file = 'data/rls_resRobSvd.rds')

# read in global mask create for reef-futures project in reef-futures-data-processing project

rls_aus_grid <- raster('data/rls-aus-grid.nc')

# load in function
source('scripts/data-processing/functions/standardize_raster.R')

# create raster list
rls_rasters <- c(as.list(rf_covs), 
                 list(depth, 
                      sst, 
                      MSEC))

# perform standardizations
rls_layers <- lapply(1:length(rls_rasters), 
                                 function(x){ 
                                   print(x) 
                                   standardize_raster(masterMask   = rls_aus_grid, 
                                                      targetRaster = rls_rasters[[x]], 
                                                      initial_constraint = T, 
                                                      maxdist = 2)}
)

# stack standarized rasters 

rls_standardized <- stack(rls_layers)

# rename 

raster_names <- do.call(c, lapply(rls_rasters, names))
names(rls_standardized) <- raster_names 

# save standardized raters 

dir.create('data/rls_rasters/')
for(i in 1:nlayers(rls_standardized)){
  print(i)
  writeRaster(rls_standardized[[i]],
              paste0('data/rls_rasters/', names(rls_standardized)[i], '.grd'),
              bylayer = T,
              overwrite=T)
  
}

# apply PCA parameters to climatic variables ----

# get raster

rls_climate  <- rls_standardized[[names(rf_covs)]]
rls_all_xy <- rls_climate[,1:2]

# check identical ordering of variables

identical(colnames(resRobSvd@completeObs), names(rls_climate))

# predict rasters from the robust PCA
rls_raster_predictions <- rls_aus_grid
resRobSVD_raster <- list()
for(i in 1:3){
  empty_grid <- rls_aus_grid
  predicted_values <- predict(resRobSvd, as.matrix(rls_climate), npcs = 3)$scores[,i]
  predicted_values[is.na(rls_aus_grid[])] <- NA
  resRobSVD_raster[[i]] <- setValues(empty_grid, predicted_values)
}
resRobSVD_raster <- stack(resRobSVD_raster)

# rename 

names(resRobSVD_raster) <- c('environmental-pc1', 'environmental-pc2', 'environmental-pc3')

# perform standardizations

rls_PCA_stack_standardized <- lapply(1:nlayers(resRobSVD_raster), 
                     function(x){ 
                       print(x) 
                       standardize_raster(masterMask   = rls_aus_grid, 
                                          targetRaster = resRobSVD_raster[[x]], 
                                          initial_constraint = T, 
                                          maxdist = 2)}
)

# write pcas to rasters

dir.create('data/rls_rasters/')
for(i in 1:nlayers(rls_PCA_stack)){
  
  writeRaster(rls_PCA_stack[[i]],
              paste0('data/rls_rasters/', names(rls_PCA_stack)[i], '.grd'),
              bylayer = T,
              overwrite=T)
  
}


# stack all of interest together, rescale and save ----

# select only covariates used in modelling
load('data/rls_covariates.RData')
covariates = rls_xy[c('SiteLongitude', 
                      'SiteLatitude',
                      "human_pop_2015_50km", 
                      "reef_area_200km", 
                      "wave_energy_mean", 
                      "Depth_GEBCO_transformed",
                      'robPCA_1', 
                      'robPCA_2', 
                      'robPCA_3',
                      'sst_mean')]

# stack all the rasters of interest 
rls_all_covariates_stacked <- stack(rls_standardized[['Elevation.relative.to.sea.level']], 
                                    resRobSVD_raster, 
                                    rls_standardized[[c('sst_mean')]])

# rename depth
names(rls_all_covariates_stacked)[which(names(rls_all_covariates_stacked) == 'Elevation.relative.to.sea.level')] <- 'Depth_GEBCO_transformed'

# transform depth
rls_all_covariates_stacked[['Depth_GEBCO_transformed']][which(rls_all_covariates_stacked[['Depth_GEBCO_transformed']][] > 0)] <- 0 # lose the differentiation between land and sea.
rls_all_covariates_stacked[['Depth_GEBCO_transformed']] <- log10(round(abs(rls_all_covariates_stacked[['Depth_GEBCO_transformed']])+1))

# rename 
names(rls_all_covariates_stacked) <- c("Depth_GEBCO_transformed",'robPCA_1', 'robPCA_2', 'robPCA_3','sst_mean')

# take the layers that are scaled and transformed
# apply rescaling to the covariates
vars <- c("Depth_GEBCO_transformed",
          'robPCA_1', 
          'robPCA_2', 
          'robPCA_3',
          'sst_mean')

# convert to points
aus_covariates <- data.frame(rasterToPoints(rls_all_covariates_stacked))
names(aus_covariates)[-c(1,2)] <- vars

# transform scale for predictions to the same as in models
aus_covariates_2 <- do.call(cbind, lapply(1:length(vars), function(x) {
  
  # get the values to rescale by
  cov_x <- covariates[vars[x]]
  scale_cov_x <- scale(cov_x)
  (aus_covariates[vars[x]] - attr(scale_cov_x,"scaled:center")) / attr(scale_cov_x,"scaled:scale")
  
}))

# combine all
all_rls_scaled <- cbind(aus_covariates[c('x', 'y')], aus_covariates_2, rasterToPoints(rls_standardized[[c('human_pop_2015_50km', 'reef_area_200km', 'wave_energy_mean')]])[,-c(1,2)])

names(all_rls_scaled)[1:2] <- c('SiteLongitude', 'SiteLatitude')

# save output
saveRDS(all_rls_scaled, file = 'data/rls_spatial_projection_data.rds')


