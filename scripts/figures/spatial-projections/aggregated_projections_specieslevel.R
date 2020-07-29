# script to spatiall map aggregated abundances


# load packages and functions ----
lib_vect <- c('tidyverse', 'sp', 'rgeos', 'raster', 'rnaturalearth', 'gridExtra', 'ggplot2')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# source functions with partial dependence plots
source('scripts/figures/functions/spatial_projection_functions.R')


# create output directories 
dir.create('results/spatial_projections_specieslevel')
dir.create('figures/spatial_projections_specieslevel/aggregated_distributions/rls', recursive = T)
dir.create('figures/spatial_projections_specieslevel/aggregated_distributions/bbs', recursive = T)

# rls aggregation ----

# insert empty vector of abundance and occurrence
rls_xy <- na.omit(readRDS('data/rls_spatial_projection_data.rds'))
rls_xy <- rls_xy %>% dplyr::select(SiteLatitude, SiteLongitude)
rls_xy$richness <- 0
rls_xy$occupancy_rate <- 0
rls_xy$abundance <- 0

# define dataset 
dataset <- 'rls'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections_cropped/rls', full.names=T)

# filter to only species with good model performance based on TSS (occurrence) and spearmans rank (abundance)
high_performance_species <- readRDS('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/results/high_performance_species.RDS')

sp_proj_files <- sp_proj_files[grep(paste(gsub(' ','_',high_performance_species), collapse = '|'), sp_proj_files, gsub(' ','_',high_performance_species))]

for(i in 1:length(unique(sp_proj_files))){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])
  sp_proj <- cbind(rls_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  sp_proj$presence = ifelse(sp_proj$occupancy_rate == 0, 0, 1)
  
  # rescale between 0 and 1
  occupancy_rate <- as.numeric(sp_proj$occupancy_rate/1000)
  abundance      <- as.numeric(exp(sp_proj$abundance_log/1000)-1)
  
  occupancy_rate_rescaled <- (occupancy_rate-min(occupancy_rate,na.rm = T))/(max(occupancy_rate,na.rm = T)-min(occupancy_rate,na.rm=T))
  abundance_rescaled <- (abundance-min(abundance,na.rm = T))/(max(abundance,na.rm = T)-min(abundance,na.rm=T))
  
  # aggregate richness, occupancy and abundance throughout loop
  rls_xy[c('occupancy_rate')] <- rls_xy[c('occupancy_rate')] + occupancy_rate_rescaled
  rls_xy[c('abundance')] <- rls_xy[c('abundance')] + abundance_rescaled
  rls_xy[c('richness')] <- rls_xy[c('richness')] + (sp_proj[c('presence')])
  
}

# save output
saveRDS(rls_xy, 'results/spatial_projections_specieslevel/rls_aggregated.RDS')

rls_xy <- readRDS('results/spatial_projections_specieslevel/rls_aggregated.RDS')


# plot rls aggregation ----

# read in aggregated properties
rls_abun_occ <- readRDS('results/spatial_projections_specieslevel/rls_aggregated.RDS')

# average to a species level
rls_abun_occ$occupancy_rate <- rls_abun_occ$occupancy_rate / rls_abun_occ$richness
rls_abun_occ$abundance      <- rls_abun_occ$abundance / rls_abun_occ$richness

# create a spatial buffer of 2° around all sites within which to investigate biodiversity patterns
load("data/rls_covariates.RData")
rls_xy = rls_xy[c('SiteLongitude', 
                  'SiteLatitude')]
rls_xy_sp <- SpatialPoints(rls_xy)
crs(rls_xy_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
rls_xy_sp_buffer <- gBuffer(rls_xy_sp, width = 2)
save(rls_xy_sp_buffer, file = 'figures/spatial_projections/rls_polygons.RData')


# crop to buffer
pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(rls_abun_occ$SiteLongitude, rls_abun_occ$SiteLatitude))),
                                   rls_abun_occ[,c('richness', 'occupancy_rate', 'abundance')],
                                   proj4string = rls_xy_sp_buffer@proj4string)
crs(pointsDF) <- rls_xy_sp_buffer@proj4string
rls_abun_occ[is.na(over(pointsDF, rls_xy_sp_buffer)), names(rls_abun_occ)] <- NA
rls_abun_occ <- na.omit(rls_abun_occ)

# create map of countries
world <- ne_countries(scale = "large", returnclass = "sf") 

# ipcc colour scale
red2  <- rgb(103, 0, 31, maxColorValue = 255)
red1  <- rgb(253, 219, 199, maxColorValue = 255)
mid   <- rgb(247, 247, 247, maxColorValue = 255)
blue1 <- rgb(209,229,240, maxColorValue = 255)
blue2 <- rgb(5, 48, 97, maxColorValue = 255)

# create simple linear model between two and examine residuals
lm_rls <- lm(abundance~occupancy_rate, data = rls_abun_occ)
rls_abun_occ$abundance_residual <- resid(lm_rls)

# get the spearmans rank correlation between values
rls_cor <- cor.test(rls_abun_occ$occupancy_rate, rls_abun_occ$abundance, method = 'spearman')
# Spearman's rank correlation rho
# 
# data:  rls_abun_occ$occupancy_rate and rls_abun_occ$abundance
# S = 1.8126e+15, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2523325 

png(filename = 'figures/spatial_projections_specieslevel/aggregated_distributions/rls/spatial_maps.png', width = 1200, height = 5000, res = 300)
grid.arrange(
  ggplot(data = world) +
    geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = occupancy_rate)) + 
    geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
    xlim(range(rls_abun_occ$SiteLongitude)) + 
    ylim(range(rls_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(data = world) +
    geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance)) + 
    geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
    xlim(range(rls_abun_occ$SiteLongitude)) + 
    ylim(range(rls_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  # add plot of residuals between abundance and occurrence
  ggplot(data = world) +
    geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance_residual)) + 
    geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
    xlim(range(rls_abun_occ$SiteLongitude)) + 
    ylim(range(rls_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(-max(rls_abun_occ$abundance_residual), -0.01, -0.005, 
                                    0, 
                                    0.005, 0.01, max(rls_abun_occ$abundance_residual)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(data = rls_abun_occ) +
    geom_point(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate, fill = abundance_residual), alpha = 0.5, pch = 21, stroke=0 ) + 
    stat_smooth(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate), method = 'lm', se = F, colour = 'black') + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(-max(rls_abun_occ$abundance_residual), -0.01, -0.005, 
                                    0, 
                                    0.005, 0.01, max(rls_abun_occ$abundance_residual)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle(label = bquote('rho' == .(signif(rls_cor$estimate,2))~','~'\n'~'p <' ~ .(0.001))) + 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          legend.position = 'none', 
          plot.title = element_text(vjust = -8, hjust = 0.05, size = 10)) + 
    ylab('species mean abundance') + 
    xlab('species mean occupancy probability'),
  
  ncol = 1
  
)
dev.off()



## perform seperate random forests for positive residuals and negative residuals!
## this will be a cool way to tell when we expect high abundance than occurrence (competitive advantage - more niche space?)
# and high occurrence than abundance (forced niche partitioning with compromise for abundance)

# random forests for rls data ----

# estimate partial dependence plots for spatial projections

# covariate data
rls_covariates  <- readRDS('data/rls_spatial_projection_data.rds')
# human, reef and wave energy are scaled in the 02_rls-environmental-data script. 
rls_covariates[c("Depth_GEBCO_transformed",
                 'robPCA_1', 
                 'robPCA_2', 
                 'robPCA_3',
                 'sst_mean')] <- as.numeric(scale(rls_covariates[c("Depth_GEBCO_transformed",
                                                                   'robPCA_1', 
                                                                   'robPCA_2', 
                                                                   'robPCA_3',
                                                                   'sst_mean')]))
rls_covariates <- rename(rls_covariates,
                         'depth/elevation' = 'Depth_GEBCO_transformed', 
                         'human' = 'human_pop_2015_50km', 
                         'sst'   =  'sst_mean',
                         'wave energy' = 'wave_energy_mean',
                         'reef area'   = 'reef_area_200km', 
                         'climate PC1' = 'robPCA_1', 
                         'climate PC2' = 'robPCA_2', 
                         'climate PC3' = 'robPCA_3')
rls_cov_names <- c('climate PC1', 'climate PC2', 'climate PC3', 'human', 'sst', 'wave energy', 'reef area', 'depth/elevation')

# projection data
rls_abun_occ_positive_residual <- rls_abun_occ %>% filter(abundance_residual > 0)
rls_abun_occ_positive_residual <- rls_abun_occ_positive_residual[sample(1:nrow(rls_abun_occ_positive_residual), 5000),]
rls_abun_occ_negative_residual <- rls_abun_occ %>% filter(abundance_residual < 0)
rls_abun_occ_negative_residual <- rls_abun_occ_negative_residual[sample(1:nrow(rls_abun_occ_negative_residual), 5000),]

# directories and names
base_dir <- 'figures/spatial_projections_specieslevel/pdp_varimp'
dataset  <- 'rls' 

# run pdp functions
partial_dependence_plots(covariates = rls_covariates,
                         cov_names  = rls_cov_names, 
                         projections = rls_abun_occ_positive_residual,
                         base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
                         dataset = 'rls', 
                         name = 'positive_abundance',
                         option = ifelse(dataset == 'rls', 'viridis', 3))

partial_dependence_plots(covariates = rls_covariates,
                         cov_names  = rls_cov_names, 
                         projections = rls_abun_occ_negative_residual,
                         base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
                         dataset = 'rls', 
                         name = 'negative_abundance',
                         option = ifelse(dataset == 'rls', 'viridis', 3))

partial_dependence_plots(covariates = rls_covariates,
                         cov_names  = rls_cov_names, 
                         projections = rbind(rls_abun_occ_positive_residual, rls_abun_occ_negative_residual),
                         base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
                         dataset = 'rls', 
                         name = 'all_residuals',
                         option = ifelse(dataset == 'rls', 'viridis', 3), 
                         N = 6)

# Mean of squared residuals: 0.2372335
# % Var explained: 93.02



# bbs aggregation ----

# insert empty vector of abundance and occurrence
bbs_xy <- na.omit(readRDS('data/bbs_spatial_projection_data.rds'))
bbs_xy <- bbs_xy %>% dplyr::select(SiteLatitude, SiteLongitude)
bbs_xy$richness <- 0
bbs_xy$occupancy_rate <- 0
bbs_xy$abundance <- 0

# define dataset 
dataset <- 'bbs'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections_cropped/bbs', full.names=T)

# filter to only species with good model performance based on TSS (occurrence) and spearmans rank (abundance)
high_performance_species <- readRDS('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/results/high_performance_species.RDS')

sp_proj_files <- sp_proj_files[grep(paste(gsub(' ','_',high_performance_species), collapse = '|'), sp_proj_files, gsub(' ','_',high_performance_species))]

for(i in 1:length(unique(sp_proj_files))){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])
  sp_proj <- cbind(bbs_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  sp_proj$presence = ifelse(sp_proj$occupancy_rate == 0, 0, 1)
  
  # deal with NAs
  sp_proj[is.na(sp_proj$occupancy_rate), c('occupancy_rate', 'abundance_log', 'presence')] <- 0
  
  # rescale between 0 and 1
  occupancy_rate <- as.numeric(sp_proj$occupancy_rate/1000)
  abundance      <- as.numeric(exp(sp_proj$abundance_log/1000)-1)
  
  occupancy_rate_rescaled <- (occupancy_rate-min(occupancy_rate,na.rm = T))/(max(occupancy_rate,na.rm = T)-min(occupancy_rate,na.rm=T))
  abundance_rescaled <- (abundance-min(abundance,na.rm = T))/(max(abundance,na.rm = T)-min(abundance,na.rm=T))
  
  # aggregate richness, occupancy and abundance throughout loop
  bbs_xy[c('occupancy_rate')] <- bbs_xy[c('occupancy_rate')] + occupancy_rate_rescaled
  bbs_xy[c('abundance')] <- bbs_xy[c('abundance')] + abundance_rescaled
  bbs_xy[c('richness')] <- bbs_xy[c('richness')] + (sp_proj[c('presence')])
  
  print(nrow(na.omit(bbs_xy)))
  
}

# save output
saveRDS(bbs_xy, 'results/spatial_projections_specieslevel/bbs_aggregated.RDS')


# plot bbs aggregation ----

# read in aggregated properties
bbs_abun_occ <- readRDS('results/spatial_projections_specieslevel/bbs_aggregated.RDS')

# remove empty grid cells
bbs_abun_occ <- bbs_abun_occ[-which(bbs_abun_occ$richness==0),]

# average to a species level
bbs_abun_occ$occupancy_rate <- bbs_abun_occ$occupancy_rate / bbs_abun_occ$richness
bbs_abun_occ$abundance      <- bbs_abun_occ$abundance / bbs_abun_occ$richness


# create a spatial buffer of 2° around all sites within which to investigate biodiversity patterns
load("data/bbs_covariates.RData")
bbs_xy = bbs_xy[c('SiteLongitude', 
                  'SiteLatitude')]
bbs_xy_sp <- SpatialPoints(bbs_xy)
crs(bbs_xy_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
bbs_xy_sp_buffer <- gBuffer(bbs_xy_sp, width = 2)
save(bbs_xy_sp_buffer, file = 'figures/spatial_projections/bbs_polygons.RData')
load('figures/spatial_projections/bbs_polygons.RData')

# create map of countries
require(sf)
world <- ne_countries(scale = "large", returnclass = "sf") 
usa <- subset(world, admin == "United States of America")
usa_2 <- crop(as_Spatial(usa), extent(c(-130, -60, 25, 50)))

# crop to USA map (all usa is in buffer)
pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(bbs_abun_occ$SiteLongitude, bbs_abun_occ$SiteLatitude))),
                                   bbs_abun_occ[,c('richness', 'occupancy_rate', 'abundance')],
                                   proj4string = bbs_xy_sp_buffer@proj4string)
crs(pointsDF) <- usa_2@proj4string
NA_vector <- is.na(sp::over(pointsDF, usa_2))
bbs_abun_occ[which(NA_vector[,1]), names(bbs_abun_occ)] <- NA
bbs_abun_occ <- na.omit(bbs_abun_occ)


# ipcc colour scale
red2  <- rgb(103, 0, 31, maxColorValue = 255)
red1  <- rgb(253, 219, 199, maxColorValue = 255)
mid   <- rgb(247, 247, 247, maxColorValue = 255)
blue1 <- rgb(209,229,240, maxColorValue = 255)
blue2 <- rgb(5, 48, 97, maxColorValue = 255)

# create simple linear model between two and examine residuals
lm_bbs <- lm(abundance~occupancy_rate, data = bbs_abun_occ)
bbs_abun_occ$abundance_residual <- resid(lm_bbs)

# get the spearmans rank correlation between values
bbs_cor <- cor.test(bbs_abun_occ$occupancy_rate, bbs_abun_occ$abundance, method = 'spearman')
# Spearman's rank correlation rho
# 
# data:  bbs_abun_occ$occupancy_rate and bbs_abun_occ$abundance
# S = 1.2537e+16, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2751153 

png(filename = 'figures/spatial_projections_specieslevel/aggregated_distributions/bbs/spatial_maps.png', width = 2200, height = 4500, res = 300)
grid.arrange(
  ggplot(data = st_as_sf(usa_2)) +
    geom_tile(data = bbs_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = occupancy_rate)) + 
    geom_sf(col = 'black', fill= 'transparent', lwd = 0.2) + 
    theme(legend.position = c(0.15, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(data = st_as_sf(usa_2)) +
    geom_tile(data = bbs_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance)) + 
    geom_sf(col = 'black', fill= 'transparent', lwd = 0.2) + 
    theme(legend.position = c(0.15, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  # add plot of residuals between abundance and occurrence
  ggplot(data = st_as_sf(usa_2)) +
    geom_tile(data = bbs_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance_residual)) + 
    geom_sf(col = 'black', fill= 'transparent', lwd = 0.2) + 
    theme(legend.position = c(0.15, 0.1)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(-max(bbs_abun_occ$abundance_residual), -0.01, -0.005, 
                                    0, 
                                    0.005, 0.01, max(bbs_abun_occ$abundance_residual)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal'), 
                         breaks = c(-0.04, 0, 0.04)) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(data = bbs_abun_occ) +
    geom_point(data = bbs_abun_occ, aes(y = abundance, x = occupancy_rate, fill = abundance_residual), alpha = 0.5, pch = 21, stroke=0 ) + 
    stat_smooth(data = bbs_abun_occ, aes(y = abundance, x = occupancy_rate), method = 'lm', se = F, colour = 'black') + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(-max(bbs_abun_occ$abundance_residual), -0.01, -0.005, 
                                    0, 
                                    0.005, 0.01, max(bbs_abun_occ$abundance_residual)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle(label = bquote('rho' == .(signif(bbs_cor$estimate,2))~','~'\n'~'p <' ~ .(0.001))) + 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          legend.position = 'none', 
          plot.title = element_text(vjust = -8, hjust = 0.03, size = 10)) + 
    xlab('summed occupancy probability') + 
    ylab('summed abundance'),
  
  nrow = 4, 
  ncol = 1
  
)
dev.off()



# random forests for bbs data ----

# estimate partial dependence plots for spatial projections

# covariate data
bbs_covariates  <- readRDS('data/bbs_spatial_projection_data.rds')
# scale variables before modelling
# elevation and human pop got scaled earlier in 02_bbs-environmental-data
bbs_covariates[c("robPCA_1", 
                 "robPCA_2", 
                 "robPCA_3", 
                 'primary_forest')] <- as.numeric(scale(bbs_covariates[c("robPCA_1", 
                                                                         "robPCA_2", 
                                                                         "robPCA_3", 
                                                                         'primary_forest')]))
bbs_covariates <- rename(bbs_covariates,
                         'depth/elevation' = 'Elevation_GEBCO', 
                         'human' = 'human_pop', 
                         'forest' = 'primary_forest', 
                         'climate PC1' = 'robPCA_1', 
                         'climate PC2' = 'robPCA_2', 
                         'climate PC3' = 'robPCA_3')
bbs_cov_names   <- c('climate PC1', 'climate PC2', 'climate PC3', 'human', 'forest', 'depth/elevation')

# projection data
bbs_abun_occ_positive_residual <- bbs_abun_occ %>% filter(abundance_residual > 0)
bbs_abun_occ_positive_residual <- bbs_abun_occ_positive_residual[sample(1:nrow(bbs_abun_occ_positive_residual), 5000),]
bbs_abun_occ_negative_residual <- bbs_abun_occ %>% filter(abundance_residual < 0)
bbs_abun_occ_negative_residual <- bbs_abun_occ_negative_residual[sample(1:nrow(bbs_abun_occ_negative_residual), 5000),]

# directories and names
base_dir <- 'figures/spatial_projections/pdp_varimp'
dataset  <- 'bbs' 

# run pdp functions
partial_dependence_plots(covariates = bbs_covariates,
                         cov_names  = bbs_cov_names, 
                         projections = bbs_abun_occ_positive_residual,
                         base_dir = 'figures/spatial_projections/pdp_varimp', 
                         dataset = 'bbs', 
                         name = 'positive_abundance',
                         option = ifelse(dataset == 'rls', 'viridis', 3))

partial_dependence_plots(covariates = bbs_covariates,
                         cov_names  = bbs_cov_names, 
                         projections = bbs_abun_occ_negative_residual,
                         base_dir = 'figures/spatial_projections/pdp_varimp', 
                         dataset = 'bbs', 
                         name = 'negative_abundance',
                         option = ifelse(dataset == 'rls', 'viridis', 3))

partial_dependence_plots(covariates = bbs_covariates,
                         cov_names  = bbs_cov_names, 
                         projections = rbind(bbs_abun_occ_positive_residual, bbs_abun_occ_negative_residual),
                         base_dir = 'figures/spatial_projections/pdp_varimp', 
                         dataset = 'bbs', 
                         name = 'all_residuals',
                         option = ifelse(dataset == 'rls', 'viridis', 3), 
                         N = 6)
# Mean of squared residuals: 6986.836
# % Var explained: 98.73





