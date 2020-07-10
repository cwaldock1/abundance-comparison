# script to spatiall map aggregated abundances


# load packages and functions ----
lib_vect <- c('tidyverse', 'sp', 'rgeos', 'raster', 'rnaturalearth', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)


# create output directories 
dir.create('figures/spatial_projections/aggregated_distributions/rls', recursive = T)
dir.create('figures/spatial_projections/aggregated_distributions/bbs', recursive = T)

#  ----

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

for(i in 1:length(unique(sp_proj_files))){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])
  sp_proj <- cbind(rls_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  sp_proj$presence = ifelse(sp_proj$occupancy_rate == 0, 0, 1)
  
  # aggregate richness, occupancy and abundance throughout loop
  rls_xy[c('occupancy_rate')] <- rls_xy[c('occupancy_rate')] + (sp_proj[c('occupancy_rate')]/1000)
  rls_xy[c('abundance')] <- rls_xy[c('abundance')] + (exp(sp_proj[c('abundance_log')]/1000)-1)
  rls_xy[c('richness')] <- rls_xy[c('richness')] + (sp_proj[c('presence')])

}

# save output
saveRDS(rls_xy, 'results/spatial_projections_cropped/rls_aggregated.RDS')



# plot rls ----

# read in aggregated properties
rls_abun_occ <- readRDS('results/spatial_projections_cropped/rls_aggregated.RDS')

# create a spatial buffer of 2° around all sites within which to investigate biodiversity patterns
load("data/rls_covariates.RData")
rls_xy = rls_xy[c('SiteLongitude', 
                      'SiteLatitude')]
rls_xy_sp <- SpatialPoints(rls_xy)
crs(rls_xy_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
rls_xy_sp_buffer <- gBuffer(rls_xy_sp, width = 2)

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
cor.test(rls_abun_occ$occupancy_rate, rls_abun_occ$abundance, method = 'spearman')
# data:  rls_abun_occ$occupancy_rate and rls_abun_occ$abundance
# S = 4.4861e+14, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8149706

png(filename = 'figures/spatial_projections/aggregated_distributions/rls/spatial_maps.png', width = 3000, height = 3000, res = 300)
grid.arrange(
ggplot(data = world) +
  geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = occupancy_rate)) + 
  geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
  xlim(range(rls_abun_occ$SiteLongitude)) + 
  ylim(range(rls_abun_occ$SiteLatitude)) + 
  theme(legend.position = c(0.2, 0.1)) + 
  viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
  ggtitle('(a) summed occupancy probability')+ 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        axis.title   = element_blank()),

ggplot(data = world) +
  geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance)) + 
  geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
  xlim(range(rls_abun_occ$SiteLongitude)) + 
  ylim(range(rls_abun_occ$SiteLatitude)) + 
  theme(legend.position = c(0.2, 0.1)) + 
  viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
  ggtitle('(b) summed abundance') + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        axis.title   = element_blank()),

# add plot of residuals between abundance and occurrence
ggplot(data = world) +
  geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance_residual)) + 
  geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
  xlim(range(rls_abun_occ$SiteLongitude)) + 
  ylim(range(rls_abun_occ$SiteLatitude)) + 
  theme(legend.position = c(0.2, 0.1)) + 
  scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                       values = c(-max(rls_abun_occ$abundance_residual), -500, -10, 
                                  0, 
                                  10, 500, max(rls_abun_occ$abundance_residual)),
                       rescaler = function(x,...) x,
                       name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
  ggtitle('(c) abundance residual') + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        axis.title   = element_blank()),

ggplot(data = rls_abun_occ) +
  geom_point(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate, fill = abundance_residual), alpha = 0.5, pch = 21, stroke=0 ) + 
  stat_smooth(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate), method = 'lm', se = F, colour = 'black') + 
  scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                       values = c(-max(rls_abun_occ$abundance_residual), -500, -10, 
                                  0, 
                                  10, 500, max(rls_abun_occ$abundance_residual)),
                       rescaler = function(x,...) x,
                       name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        legend.position = 'none') + 
  xlab('summed occupancy probability') + 
  ylab('abundance') + 
  ggtitle('(d) comparison of abundance and occurrence'),

nrow = 2

)
dev.off()






