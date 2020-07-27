# script to crop the spatial projections to the convex hulls

# load packages and functions ----
lib_vect <- c('tidyverse', 'sp', 'rgeos', 'raster', 'rnaturalearth', 'parallel')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# rls crops ----

# insert empty vector of abundance and occurrence
rls_xy <- na.omit(readRDS('data/rls_spatial_projection_data.rds'))
rls_xy <- rls_xy %>% dplyr::select(SiteLatitude, SiteLongitude)
rls_xy$richness <- 0
rls_xy$occupancy_rate <- 0
rls_xy$abundance_log <- 0

# define dataset 
dataset <- 'rls'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections/rls', full.names=T)

# create output directory
dir.create('results/spatial_projections_cropped/rls', recursive = T)

# loop to crop
mclapply(1:length(sp_proj_files), mc.cores = 10, function(i){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])[[1]]
  sp_proj <- cbind(rls_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  
  # read in convex hull for the species
  sp_ch <- readRDS(paste0('data/species_range_convex_hull/', dataset, '/', strsplit(sp_proj_files, '/')[[i]][4]))
  
  # convert to spatial points data frame and crop to convex hull
  pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(sp_proj$SiteLongitude, sp_proj$SiteLatitude))),
                                     sp_proj[,c('occupancy_rate', 'abundance_log')],
                                     proj4string = sp_ch@proj4string)
  crs(pointsDF) <- sp_ch@proj4string
  sp_proj[is.na(over(pointsDF, sp_ch)), c('occupancy_rate', 'abundance_log')] <- 0
  
  saveRDS(sp_proj[-c(1,2)], paste0('results/spatial_projections_cropped/rls/', strsplit(sp_proj_files, '/')[[i]][4]))

})


# bbs crops ----

# insert empty vector of abundance and occurrence
bbs_xy <- na.omit(readRDS('data/bbs_spatial_projection_data.rds'))
bbs_xy <- bbs_xy %>% dplyr::select(SiteLatitude, SiteLongitude)
bbs_xy$richness <- 0
bbs_xy$occupancy_rate <- 0
bbs_xy$abundance_log <- 0

# define dataset 
dataset <- 'bbs'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections/bbs', full.names=T)

# create output directory
dir.create('results/spatial_projections_cropped/bbs', recursive = T)

# create USA clipping
world <- ne_countries(scale = "large", returnclass = "sf") 
usa <- subset(world, admin == "United States of America")
usa_2 <- crop(as_Spatial(usa), extent(c(-130, -60, 25, 50)))

# crop to USA map (all usa is in buffer)
all_points <- SpatialPoints(coordinates(cbind(bbs_xy$SiteLongitude, bbs_xy$SiteLatitude)),
                                   proj4string = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
NA_vector <- is.na(sp::over(all_points, usa_2))
NA_vector <- NA_vector[,1]

# loop to crop
mclapply(1:length(sp_proj_files), mc.cores = 10, function(i){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])[[1]]
  sp_proj <- cbind(bbs_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  
  # read in convex hull for the species
  sp_ch <- readRDS(paste0('data/species_range_convex_hull/', dataset, '/', strsplit(sp_proj_files, '/')[[i]][4]))
  
  # convert to spatial points data frame and crop to convex hull
  pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(sp_proj$SiteLongitude, sp_proj$SiteLatitude))),
                                     sp_proj[,c('occupancy_rate', 'abundance_log')],
                                     proj4string = sp_ch@proj4string)
  crs(pointsDF) <- sp_ch@proj4string
  sp_NAs <- is.na(over(pointsDF, sp_ch))
  # combine with clipping for contiguous USA
  all_NAs <- ifelse((sp_NAs + NA_vector) != 0, TRUE, FALSE)
  sp_proj[all_NAs, c('occupancy_rate', 'abundance_log')] <- NA
  
  # also clip and refine analysis to contiguous USA
  saveRDS(sp_proj[-c(1,2)], paste0('results/spatial_projections_cropped/bbs/', strsplit(sp_proj_files, '/')[[i]][4]))
  
})

