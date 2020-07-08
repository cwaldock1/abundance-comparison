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
lib_vect <- c('tidyverse', 'rgdal', 'raster', 'maptools', 'summarytools', 'pcaMethods', 'gstat', 'rnaturalearth')
#lib_vect <- c('tidyverse', 'rgdal', 'raster', 'maptools', 'summarytools', 'pcaMethods', 'gstat')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load in bbs sites for matching ----

# projection
Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
Proj_bbs <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0 +units=m"

# load data
bbs_xy <- do.call(rbind, lapply(list.files('data/bbs_all_basic', full.names = T), function(x){
  # read in all the individual files
  ddata <- readRDS(x)
  ddata_all <- do.call(rbind, ddata)
  data.frame(SiteLongitude = ddata_all$SiteLongitude, SiteLatitude = ddata_all$SiteLatitude) %>% 
    unique()})) %>% 
  unique()

# convert to points

bbs_points <- SpatialPoints(bbs_xy, proj4string = crs(Proj))

# load in and standardise environmental layers ----

# species distribution models will be build from the best resolution data available rather than comforming to a single grid
# this is because we are not evaluating the spatial distribution of the species distributions but instead just modelling the parameters

# elevation ----
# extract depth
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/ 
# new data for 2019. 
elevation <- raster("/Volumes/RF-env-data/reef-futures/env-data/GEBCO_2019/GEBCO_2019.nc", CRS = Proj)

# plot(depth)
bbs_xy$Elevation_GEBCO <- raster::extract(elevation, bbs_points)

elevation_scale_object <- scale(bbs_xy$Elevation_GEBCO, scale = T, center = T)

bbs_xy$Elevation_GEBCO <- as.numeric(elevation_scale_object)

# MERRAclim ----

# stack MERRAclim

stack_00s <- stack(list.files('/Volumes/RF-env-data/reef-futures/env-data/MERRAclim/5m_mean_00s/', full.names = T))
stack_90s <- stack(list.files('/Volumes/RF-env-data/reef-futures/env-data/MERRAclim/5m_mean_90s/', full.names = T))

# estimate the mean of the stacked objects

mean_stack <- (stack_00s + stack_90s)/2

names(mean_stack) <- gsub('_00s|X','', names(stack_00s))

# extract the sites from all locations

MERRAclim_values <- data.frame(extract(mean_stack, bbs_points))

# add names

names(MERRAclim_values) <- gsub('_00s|X','', names(stack_00s))

# back-transform those that have been multiple and divided by 100, 100,000

temp_vars <- paste0('bio', as.character(1:11), '_', collapse = '|')
temp_vars <- grep(temp_vars, paste0(names(MERRAclim_values), '_'), value = T)
temp_vars <- substr(temp_vars, 1, nchar(temp_vars)-1)

hum_vars <- paste0('bio', as.character(12:19), '_', collapse = '|')
hum_vars <- grep(hum_vars, paste0(names(MERRAclim_values), '_'), value = T)
hum_vars <- substr(hum_vars, 1, nchar(hum_vars)-1)

# divide by multiples

temp_stack <- mean_stack[[which(gsub('_00s|X','', names(mean_stack)) %in% temp_vars)]]/10
hum_stack  <- mean_stack[[which(gsub('_00s|X','', names(mean_stack)) %in% hum_vars)]]/100000
mean_stack <- stack(temp_stack, hum_stack)

MERRAclim_values[temp_vars] <- MERRAclim_values[temp_vars]/10
MERRAclim_values[hum_vars] <- MERRAclim_values[hum_vars]/100000

# reorder data 

MERRAclim_values <- MERRAclim_values[c(temp_vars[c(1,4:11,2,3)], 
                                       hum_vars)]

# check summary
view(dfSummary(MERRAclim_values))

# see table 1 in Vega et al. 2017
# bio1  = Annual mean temperature
# bio2  = Mean diurnal range temperature
# bio3  = Isothermality (bio2/bio7)
# bio4  = Temperature seasonality (sd*100)
# bio5  = Maximum temperature of warmest month
# bio6  = Minimum temperature of coldest month 
# bio7  = Temperature annual range (bio5-bio6)
# bio8  = Mean temperature of most humid quarter
# bio9  = Mean temperature of least humid quarter
# bio10 = Mean temperature of warmest quarter
# bio11 = Mean temperature of coldest quarter
# bio12 = Annual mean specific humidity
# bio13 = Specific humidity of the most humid month
# bio14 = Specific humidity of the least humid month
# bio15 = Specific humidity seasonality (CV)
# bio16 = Specific humidity of most humid quarter
# bio17 = Specific humidity of least humid quarter
# bio18 = Specific humidity of warmest quarter
# bio19 = Specific humidity of colder quarter

# estimate principal components using robust pca
resRobSvd <- pcaMethods::pca(MERRAclim_values, 
                             method = "robustPca", 
                             center = T, 
                             scale = 'uv', 
                             nPcs = dim(MERRAclim_values)[2])

resRobSvd # suggests 3 PCs lead to increasing variance explained.
#robustPca calculated PCA
#Importance of component(s):
#  PC1    PC2     PC3     PC4     PC5      PC6     PC7     PC8     PC9     PC10     PC11    PC12     PC13     PC14       PC15    PC16    PC17     PC18    PC19
#  R2            0.6486 0.1656 0.06458 -0.1724 -0.1084 0.008925 -0.1478 -0.1046 -0.1123 -0.02403 -0.24732 -0.2452 -0.04901 -0.06155 -0.0003978 -0.2343 -0.2503 -0.07001 -0.1246
#  Cumulative R2 0.6486 0.8142 0.87879  0.7064  0.5980 0.606945  0.4592  0.3546  0.2423  0.21826 -0.02907 -0.2742 -0.32323 -0.38478 -0.3851768 -0.6195 -0.8698 -0.93980 -1.0644
#  19 	Variables
#  4602 	Samples
#  0 	NAs ( 0 %)
#  19 	Calculated component(s)
#  Data was mean centered before running PCA 
#  Data was scaled before running PCA 
#  Scores structure:
#  [1] 4602   19
#  Loadings structure:
#   [1] 19 19


Loadings <- resRobSvd@loadings[,1:3]
Loadings <- reshape2::melt(Loadings)

pdf(file = 'figures/robustPCA-bbs-env-covariates.pdf', width = 7.5, height = 7.5)
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

#slplot(resRobSvd, scoresLoadings = c(T,T))
biplot(resRobSvd)
#plotPcs(resRobSvd)

# get values to model with
MERRAclim_values$robPCA_1 <- as.numeric(resRobSvd@scores[,1])
MERRAclim_values$robPCA_2 <- as.numeric(resRobSvd@scores[,2]) 
MERRAclim_values$robPCA_3 <- as.numeric(resRobSvd@scores[,3])

# bind to rls
bbs_xy <- cbind(bbs_xy, MERRAclim_values)

# Gridded Population of the World (GPW), v4 ----

# read in human population density

human_pop <- raster('/Volumes/RF-env-data/reef-futures/env-data/Gridded Population of the World (GPW), v4/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11_totpop_2pt5_min_nc/gpw_v4_population_density_adjusted_rev11_2pt5_min.nc')

# log transform

human_pop <- scale(log10(human_pop+1))

# check plot

plot(human_pop)
hist(human_pop[])

# perform extractions

human_pop <- force(FillSiteNAs(StackLayer = human_pop, RLS_Sites = bbs_points))

# rename 
names(human_pop) <- 'human_pop'

# extract points
human_pop_ext <- data.frame(human_pop = extract(human_pop, bbs_points))

# check for NAs
table(is.na(human_pop_ext))

# bind to rls_xy object
bbs_xy <- cbind(bbs_xy, human_pop_ext)

# primary forest data ----

# read in raster

primary_forest <- raster('/Volumes/RF-env-data/reef-futures/env-data/Global_30s_resolution_land_use_for_2005/PRI_2005/PRI_1km_2005_0ice.bil')

# change crs

crs(primary_forest) <- crs(human_pop)

# rename 

names(primary_forest) <- 'primary_forest'

# extract

primary_forest_ext <- data.frame(primary_forest = extract(primary_forest, bbs_points))

# find NAs

table(is.na(primary_forest_ext$primary_forest))
primary_forest_ext[is.na(primary_forest_ext$primary_forest),1] <- as.numeric(extract(primary_forest, 
          bbs_points[is.na(primary_forest_ext$primary_forest)], buffer = 5000, fun = mean))

# distribution of variable

hist(primary_forest_ext$primary_forest)

# bind to bbs

bbs_xy <- cbind(bbs_xy, primary_forest_ext)

# transform, scale and centre all non-pca variables ----

# check for any NAs
table(is.na(bbs_xy))

# save rls_xy ----
save(bbs_xy, file = 'data/bbs_covariates.RData')

# create stacked raster of key variables for projections ----

# load in mapped polygons with maptools

usa_polygon <- Rgshhs("/Volumes/RF-env-data/reef-futures/env-data/gshhg-bin-2/gshhs_i.b", 
                      xlim = extent(bbs_points)[1:2]+360,
                      ylim = extent(bbs_points)[3:4]+c(-10, 10), 
                      level = 2, 
                      minarea = 10000, 
                      shift = T, 
                      verbose = TRUE, 
                      no.clip = F, 
                      properly = F, 
                      avoidGEOS=F, 
                      checkPolygons=T)

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  rgeos::gIntersection(shp, b_poly, byid = TRUE)
}

usa_polygon <- gClip(usa_polygon$SP, bbox(bbs_points))

# plot as a test of resolution and minimum area

plot(usa_polygon)

# manually select dimensions of human pop as base raster to resample to at a 0.08° scale

base_raster <- human_pop
base_raster[!is.na(base_raster)] <- 0
base_raster <- crop(base_raster, extent(bbs_points))
plot(base_raster)

# rasterize the polygons of usa

base_raster <- rasterize(usa_polygon, base_raster, field = 1)

points(bbs_points)

masterMask   <- base_raster
targetRaster <- human_pop

# load in function
source('scripts/data-processing/functions/standardize_raster.R')

plot(standardize_raster(base_raster, human_pop))

# list all the layers to standardize

bbs_rasters <- c(as.list(mean_stack), # merraclim
                 list(
                   elevation,  
                   human_pop,
                   primary_forest
                   ))

# perform standardisation

bbs_layers <- parallel::mclapply(1:length(bbs_rasters), 
                                 mc.cores = 5, 
                                 function(x){ 
                                   print(x) 
                                   standardize_raster(base_raster, bbs_rasters[[x]])}
                                 )

# stack standarized rasters 

bbs_standardized <- stack(bbs_layers)

# rename 

raster_names <- do.call(c, lapply(bbs_rasters, names))
names(bbs_standardized) <- raster_names 

# transform as for extracted values (only necessary for elevation)
ele_centre <- attr(elevation_scale_object, 'scaled:center')
ele_scale  <- attr(elevation_scale_object, 'scaled:scale')

bbs_standardized$Elevation.relative.to.sea.level <- ((bbs_standardized$Elevation.relative.to.sea.level) - ele_centre) / ele_scale

# save standardized raters 

dir.create('data/bbs_rasters/')
for(i in 1:nlayers(bbs_standardized)){
  
  writeRaster(bbs_standardized[[i]],
              paste0('data/bbs_rasters/', names(bbs_standardized)[i], '.grd'),
              bylayer = T,
              overwrite=T)
  
}

# apply PCA parameters to merraclim variables ----

# get raster

bbs_merra  <- rasterToPoints(bbs_standardized[[grep('_bio', names(bbs_standardized))]])
bbs_all_xy <- rasterToPoints(bbs_standardized[[grep('_bio', names(bbs_standardized))]])[,c(1,2)]

# rename to match

colnames(bbs_merra) <- gsub('X', '', colnames(bbs_merra))
bbs_merra <- bbs_merra[,na.omit(match(colnames(MERRAclim_values), colnames(bbs_merra)))]

# perform scaling to exact same dimensions of 

raster_predictions <- predict(resRobSvd, bbs_merra, pcs = 3)

# convert PCAs to raster stacks

bbs_PCA_stack <- stack(lapply(1:3, function(x) {rasterFromXYZ(cbind(bbs_all_xy, raster_predictions$scores[,x]))}))

# rename 

names(bbs_PCA_stack) <- c('merra-pc1', 'merra-pc2', 'merra-pc3')

# save pca rasters too 

dir.create('data/bbs_rasters/')
for(i in 1:nlayers(bbs_PCA_stack)){
  
  writeRaster(bbs_PCA_stack[[i]],
              paste0('data/bbs_rasters/', names(bbs_PCA_stack)[i], '.grd'),
              bylayer = T,
              overwrite=T)
  
}


# read in spatial datasets as a large matrix ----

# list all spatial datasets 
spatial_datasets <- list.files('data/bbs_rasters', pattern = 'gri', full.names = T)

# read in example raster and convert to matrix
spatial_datasets_raster <- lapply(spatial_datasets, function(x){
  
  ras <- raster(x)
  crs(ras) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
  ras_points <- data.frame(rasterToPoints(ras))
  return(ras_points)
  
})

# ensure x and y are identical (previous code seems to have added 0.00001 to some of the xy values in conversions)
spatial_datasets_raster <- lapply(spatial_datasets_raster, function(x){
  x$x <- spatial_datasets_raster[[1]]$x 
  x$y <- spatial_datasets_raster[[1]]$y
  return(x)
})

# apply together
spatial_datasets_raster_all <- Reduce(left_join, spatial_datasets_raster)

# check for NAs 
lapply(spatial_datasets_raster_all, function(x) table(is.na(x)))

# select only covariates used in modelling
load('data/bbs_covariates.RData')
covariates = bbs_xy[c('SiteLongitude', 
                      'SiteLatitude',
                      "robPCA_1", 
                      "robPCA_2", 
                      "robPCA_3", 
                      'primary_forest')]

covariates_usa <- spatial_datasets_raster_all[c('x', 
                                                'y',
                                                "merra.pc1",
                                                "merra.pc2", 
                                                "merra.pc3", 
                                                'primary_forest')]

names(covariates_usa) <- names(covariates)

# apply rescaling to the covariates
vars <- c("robPCA_1", "robPCA_2", "robPCA_3", 'primary_forest')
covariates_usa_scaled <- do.call(cbind, lapply(1:length(vars), function(x) {
  
  # get the values to rescale by
  cov_x <- covariates[vars[x]]
  scale_cov_x <- scale(cov_x)
  (covariates_usa[vars[x]] - attr(scale_cov_x,"scaled:center")) / attr(scale_cov_x,"scaled:scale")
  
}))

# combine all
all_bbs_scaled <- cbind(covariates_usa[c('SiteLongitude', 'SiteLatitude')],
                        spatial_datasets_raster_all['human_pop'], 
                        spatial_datasets_raster_all['Elevation.relative.to.sea.level'], 
                        covariates_usa_scaled)[c('SiteLongitude', 
                                               'SiteLatitude',
                                               "robPCA_1", 
                                               "robPCA_2", 
                                               "robPCA_3", 
                                               "human_pop",
                                               'primary_forest', 
                                               'Elevation.relative.to.sea.level')]

names(all_bbs_scaled)[which(names(all_bbs_scaled) == 'Elevation.relative.to.sea.level')] <- 'Elevation_GEBCO'

# save RDS
saveRDS(all_bbs_scaled, file = 'data/bbs_spatial_projection_data.rds')


