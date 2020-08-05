# script to spatiall map aggregated abundances


# load packages and functions ----
lib_vect <- c('tidyverse', 'sp', 'rgeos', 'raster', 'rnaturalearth', 'gridExtra', 'ggplot2', 'Hmisc')
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
rls_xy$disconnect <- 0

# define dataset 
dataset <- 'rls'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections_cropped/rls', full.names=T)

# filter to only species with good model performance based on TSS (occurrence) and spearmans rank (abundance)
high_performance_species <- readRDS('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/results/high_performance_species.RDS')

sp_proj_files <- sp_proj_files[grep(paste(gsub(' ', '_', high_performance_species), collapse = '|'), sp_proj_files, gsub(' ','_',high_performance_species))]

# check which sites have species 
rich = 0
for(i in 1:length(sp_proj_files)){
  print(i)
  rich = rich + (readRDS(sp_proj_files[i])[,1]/1000)
}
non_0 <- which(rich>0)
set.seed(123)
sites <- sort(non_0[sample(1:length(non_0), 1000)])
assemblage_disconnect <- matrix(nrow = length(sites), ncol = length(sp_proj_files))

for(i in 1:length(unique(sp_proj_files))){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])
  sp_proj <- cbind(rls_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  sp_proj$presence = ifelse(sp_proj$occupancy_rate == 0, 0, 1)
  
  # convert back to actual numbers
  occupancy_rate <- as.numeric(sp_proj$occupancy_rate/1000)
  abundance      <- as.numeric(exp(sp_proj$abundance_log/1000)-1)
  
  # turn absences to 0
  occupancy_rate[sp_proj$presence == 0] <- 0
  abundance[sp_proj$presence == 0] <- 0
  
  # rescale to percentiles. This asks, what proportion of the data fall above or below a certain value.
  occupancy_rate[sp_proj$presence!=0] <- ecdf(occupancy_rate[sp_proj$presence!=0])(occupancy_rate[sp_proj$presence!=0])
  abundance[sp_proj$presence!=0]      <- ecdf(abundance[sp_proj$presence!=0])(abundance[sp_proj$presence!=0])
  
  # aggregate richness, occupancy and abundance throughout loop
  rls_xy[c('occupancy_rate')] <- rls_xy[c('occupancy_rate')] + occupancy_rate
  rls_xy[c('abundance')] <- rls_xy[c('abundance')] + abundance
  rls_xy[c('richness')] <- rls_xy[c('richness')] + (sp_proj[c('presence')])
  rls_xy[sp_proj$presence != 0, c('disconnect')] <- rls_xy[sp_proj$presence != 0, c('disconnect')] + (abundance[sp_proj$presence!=0] - occupancy_rate[sp_proj$presence!=0])

  # find assemblage disconnect
  assemblage_disconnect[,i] <- abundance[sites] - occupancy_rate[sites]
  # remove disconnect when species arent present
  assemblage_disconnect[sp_proj$presence[sites]==0,i] <- NA
  
}

# distribution of disconnect across all assemblages
hist(apply(assemblage_disconnect, 1, function(x) mean(x, na.rm=T)))

# save output
saveRDS(rls_xy, 'results/spatial_projections_specieslevel/rls_aggregated_V2.RDS')
saveRDS(assemblage_disconnect, 'results/spatial_projections_specieslevel/rls_assemblage_disconnect.RDS')


# plot rls aggregation ----

# read in aggregated properties
rls_abun_occ <- readRDS('results/spatial_projections_specieslevel/rls_aggregated_V2.RDS')

# average to a species level
rls_abun_occ$occupancy_rate <- rls_abun_occ$occupancy_rate / rls_abun_occ$richness
rls_abun_occ$abundance      <- rls_abun_occ$abundance / rls_abun_occ$richness
rls_abun_occ$disconnect      <- rls_abun_occ$disconnect / rls_abun_occ$richness

# create a spatial buffer of 2° around all sites within which to investigate biodiversity patterns
load("data/rls_covariates.RData")
rls_xy = rls_xy[c('SiteLongitude', 
                  'SiteLatitude')]
rls_xy_sp <- SpatialPoints(rls_xy)
crs(rls_xy_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
# rls_xy_sp_buffer <- gBuffer(rls_xy_sp, width = 2)
# save(rls_xy_sp_buffer, file = 'figures/spatial_projections/rls_polygons.RData')
load(file = 'figures/spatial_projections/rls_polygons.RData')

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
# S = 1.844e+15, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2394126 

png(filename = 'figures/spatial_projections_specieslevel/aggregated_distributions/rls/spatial_maps_V2.png', width = 1200, height = 5000, res = 300)
grid.arrange(
  ggplot(data = world) +
    geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = occupancy_rate)) + 
    geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
    xlim(range(rls_abun_occ$SiteLongitude)) + 
    ylim(range(rls_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal'), 
                                breaks = c(0.3, 0.5, 0.7)) +
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
    viridis::scale_fill_viridis(option = "magma", name = NULL, guide = guide_colourbar(direction = 'horizontal'), 
                                breaks = c(0.3, 0.5, 0.7)) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  # add plot of residuals between abundance and occurrence
  ggplot(data = world) +
    geom_tile(data = rls_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = disconnect)) + 
    geom_sf(col = 'black', fill= 'gray90', lwd = 0.2) + 
    xlim(range(rls_abun_occ$SiteLongitude)) + 
    ylim(range(rls_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(min(rls_abun_occ$disconnect), -0.05, -0.0001, 
                                    0, 
                                    0.0001, 0.05, max(rls_abun_occ$disconnect)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(rls_abun_occ) + 
    geom_histogram(data =rls_abun_occ, aes(x = disconnect, fill = ..x..)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(min(rls_abun_occ$disconnect), -0.05, -0.0001, 
                                    0, 
                                    0.0001, 0.05, max(rls_abun_occ$disconnect)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          legend.position = 'none', 
          plot.title = element_text(vjust = -8, hjust = 0.05, size = 10)), 
    
  # gplot(data = rls_abun_occ) +
  #  geom_point(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate, fill = disconnect), alpha = 0.5, pch = 21, stroke=0 ) + 
  #  stat_smooth(data = rls_abun_occ, aes(y = abundance, x = occupancy_rate), method = 'lm', se = F, colour = 'black') + 
  #  scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
  #                       values = c(-max(rls_abun_occ$abundance_residual), -0.01, -0.005, 
  #                                  0, 
  #                                  0.005, 0.01, max(rls_abun_occ$abundance_residual)),
  #                       rescaler = function(x,...) x,
  #                       name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
  #  ggtitle(label = bquote('rho' == .(signif(rls_cor$estimate,2))~','~'\n'~'p <' ~ .(0.001))) + 
  #  theme(panel.background = element_blank(), 
  #        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
  #        legend.position = 'none', 
  #        plot.title = element_text(vjust = -8, hjust = 0.05, size = 10)) + 
  #  ylab('species mean abundance') + 
  #  xlab('species mean occupancy probability'),
  
  ncol = 1
  
)
dev.off()



## perform seperate random forests for positive residuals and negative residuals!
## this will be a cool way to tell when we expect high abundance than occurrence (competitive advantage - more niche space?)
# and high occurrence than abundance (forced niche partitioning with compromise for abundance)

# plot rls species level within 1000 sample assemblages disconnection ----

rls_assemblage_disconnect <- readRDS('results/spatial_projections_specieslevel/rls_assemblage_disconnect.RDS')

# re-order all values and remove 0s to NAs
rls_assemblage_disconnect[is.na(rls_assemblage_disconnect)] <- 0
for(i in 1:ncol(rls_assemblage_disconnect)){
  rls_assemblage_disconnect[,i] <- sort(rls_assemblage_disconnect[,i])
  #rls_assemblage_disconnect[rls_assemblage_disconnect[,i] == 0,i] <- NA
}

# convert to long dataframe
rls_assemblage_disconnect <- data.frame(rls_assemblage_disconnect)
rls_assemblage_disconnect$community <- 1:nrow(rls_assemblage_disconnect)
rls_assemblage_disconnect <- pivot_longer(rls_assemblage_disconnect, cols = X1:X456)
rls_assemblage_disconnect_sorted <- rls_assemblage_disconnect %>% 
  # estimate ordered values and species per community
  group_by(community) %>% 
  nest() %>% 
  mutate(value_sort = purrr::map(data, ~sort(.$value, na.last = T)), 
         row_sort   = purrr::map(data, ~1:length(.$value))) %>% 
  unnest(value_sort, row_sort) %>% 
  ungroup() %>% 
  group_by(community) %>% 
  nest() %>% 
  # estimate mean value overall for a community
  mutate(mean_value = purrr::map(data, ~mean(.$value_sort))) %>% 
  unnest(mean_value) %>% 
  ungroup() %>% 
  unnest(data)

col_low_rls <- viridis::viridis(10, option =7)[2]
col_high_rls <- viridis::viridis(10, option = 7)[9]

png(filename = 'figures/spatial_projections_specieslevel/assemblage_disconnect/rls_assemblage_disconnect.png', width = 1200, height = 1200, res = 300)
ggplot(rls_assemblage_disconnect_sorted) + 
  geom_line(aes(x = row_sort, y = value_sort, group = community, col = mean_value), lwd = 0.2, alpha = 0.5) + 
  geom_line(data = rls_assemblage_disconnect_sorted %>% 
              filter(community %in% round(seq(1, 1000, length.out=50))), 
            aes(x = row_sort, y = value_sort, group = community, col = mean_value), lwd = 0.7) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        legend.position = 'none') + 
  scale_colour_gradientn(colours = c( viridis::viridis(10, option =7)[2],  viridis::viridis(10, option = 7)[4], 'gray25', 
                                      'black', 'gray25',  viridis::viridis(10, option = 7)[6],  viridis::viridis(10, option = 7)[9]),
                         values = c(min(rls_assemblage_disconnect_sorted$mean_value), -0.25, -0.01, 
                                    0, 
                                    0.01, 0.25, max(rls_assemblage_disconnect_sorted$mean_value)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) + 
  ylab('species disconnect') + 
  xlab('rank order')
dev.off()

# summary statistics of assemblage disconnect

# what proportion of assemblages have mean disconnects < - or + certain percentage points
rls_assemblage_disconnect %>% 
  group_by(community) %>% 
  filter(value != 0) %>% 
  do(mean_value = mean(.$value)) %>% 
  dplyr::select(mean_value, community) %>% 
  unnest(mean_value) %>% 
  unique() %>% 
  mutate(mean_value = abs(mean_value)) %>%
  do(percentage_0.05 = nrow(.[.$mean_value > 0.05,])/1000,
     percentage_0.1 = nrow(.[.$mean_value > 0.1,])/1000, 
     percentage_0.25 = nrow(.[.$mean_value > 0.25,])/1000, 
     percentage_0.5 = nrow(.[.$mean_value > 0.5,])/1000) %>% 
  unnest()
# percentage_0.05 percentage_0.1 percentage_0.25 percentage_0.5
# <dbl>          <dbl>           <dbl>          <dbl>
# 0.718          0.547           0.292          0.066


# is every community sensitive?
# what percentage of species in a community are disconnected?
rls_wide_disconnect <- readRDS('results/spatial_projections_specieslevel/rls_assemblage_disconnect.RDS')

# convert species disconnects to absolute values
rls_wide_disconnect <- abs(rls_wide_disconnect)
percentage_assemblage_disconnected <- c()
for(i in 1:nrow(rls_wide_disconnect)){
  
  ass_dis <- rls_wide_disconnect[i,]
  richness <- sum(!is.na(ass_dis))
  more_than_25 <- sum(ass_dis > c(0.5), na.rm = T)
  percentage_assemblage_disconnected[i] <- more_than_25 / richness
  
}
hist(percentage_assemblage_disconnected)

disconnect_vectors <- sapply(seq(from = 0, to = 1, by = 0.01), function(x){ 
  percentage_assemblage_disconnected <- c()
  for(i in 1:nrow(rls_wide_disconnect)){
  ass_dis <- rls_wide_disconnect[i,]
  richness <- sum(!is.na(ass_dis))
  more_than_25 <- sum(ass_dis > c(x), na.rm = T)
  percentage_assemblage_disconnected[i] <- more_than_25 / richness
  }
  return(quantile(percentage_assemblage_disconnected, c(0.05, 0.5, 0.95)))
  })

disconnect_vectors_df  <- gather(data.frame(disconnect_vectors))
disconnect_vectors_df$quantiles   <- rownames(disconnect_vectors)
disconnect_vectors_df$break_point <- sort(rep(seq(from = 0, to = 1, by = 0.01), 3))
disconnect_vectors_df <- pivot_wider(disconnect_vectors_df, names_from = quantiles, values_from = value)


rls_disconnect_breakpoints <- ggplot(data = disconnect_vectors_df) + 
  geom_ribbon(aes(x = break_point, ymin=`5%`, ymax=`95%`), alpha=0.5, 
              fill = viridis::viridis(10,option = 7)[5]) + 
  geom_line(aes(x = break_point,
                y = `50%`)) + 
  geom_line(aes(x = break_point, y=`5%`), alpha=1, lty = 2) + 
  geom_line(aes(x = break_point, y=`95%`), alpha=1, lty= 2) + 
  #geom_segment(data = disconnect_vectors_df[disconnect_vectors_df$break_point %in% c(0.25, 0.5, 0.75),], 
  #             aes(x=break_point, y=`50%`, xend=c(0.25, 0.5, 0.75), yend=0), color=c("gray25", 'gray50', 'gray75'))  + 
  geom_segment(data = disconnect_vectors_df[disconnect_vectors_df$break_point %in% c(0.25, 0.5, 0.75),], 
               aes(x=break_point, y=`50%`, xend=c(0,0,0), yend=`50%`), color=c("gray25", 'gray50', 'gray75')) + 
  xlab('disconnection threshold') + 
  ylab('% species in assemblage disconnected') + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)
  
png(file = 'figures/spatial_projections_specieslevel/assemblage_disconnect/rls_breakpoint.png', width = 1200, height = 1200, res = 300)
rls_disconnect_breakpoints
dev.off()


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

# create stratified sample
rls_abun_occ_sampled <- rls_abun_occ %>% 
  mutate(sample_cut = cut2(.$disconnect, g = 5, onlycuts = F)) %>% 
  group_by(sample_cut) %>% 
  nest() %>% 
  mutate(samples = purrr::map(data, ~.[sample(1:nrow(.), 10000, replace = T),])) %>% 
  dplyr::select(-data) %>% 
  unnest(samples) %>% ungroup()

# projection data
rls_abun_occ_positive_residual <- rls_abun_occ_sampled %>% filter(abundance_residual > 0)
#rls_abun_occ_positive_residual <- rls_abun_occ_positive_residual[sample(1:nrow(rls_abun_occ_positive_residual), 5000),]
rls_abun_occ_negative_residual <- rls_abun_occ_sampled %>% filter(abundance_residual < 0)
#rls_abun_occ_negative_residual <- rls_abun_occ_negative_residual[sample(1:nrow(rls_abun_occ_negative_residual), 5000),]

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
                         response_column = 'disconnect',
                         projections = rls_abun_occ_sampled,
                         base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
                         dataset = 'rls', 
                         name = 'all_residuals_V2',
                         option = ifelse(dataset == 'rls', 'viridis', 3), 
                         N = 6)

# Mean of squared residuals: 0.0001413254
# % Var explained: 99.13



# bbs aggregation ----

# insert empty vector of abundance and occurrence
bbs_xy <- na.omit(readRDS('data/bbs_spatial_projection_data.rds'))
bbs_xy <- bbs_xy %>% dplyr::select(SiteLatitude, SiteLongitude)
bbs_xy$richness <- 0
bbs_xy$occupancy_rate <- 0
bbs_xy$abundance <- 0
bbs_xy$disconnect <- 0

# define dataset 
dataset <- 'bbs'

# species projections folder 
sp_proj_files <- list.files('results/spatial_projections_cropped/bbs', full.names=T)

# filter to only species with good model performance based on TSS (occurrence) and spearmans rank (abundance)
high_performance_species <- readRDS('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/results/high_performance_species.RDS')
sp_proj_files <- sp_proj_files[grep(paste(gsub(' ','_',high_performance_species), collapse = '|'), sp_proj_files, gsub(' ','_',high_performance_species))]

# check which sites have species 
rich = 0
for(i in 1:length(sp_proj_files)){
  print(i)
  read_file <- (readRDS(sp_proj_files[i])[,1]/1000)
  read_file[is.na(read_file)] <- 0
  rich = rich + read_file
}
non_0 <- which(rich>0)
set.seed(123)
sites <- sort(non_0[sample(1:length(non_0), 1000)])
assemblage_disconnect <- matrix(nrow = length(sites), ncol = length(sp_proj_files))

for(i in 1:length(unique(sp_proj_files))){
  
  print(i)
  # read in spatial projections and align to xy
  sp_proj <- readRDS(sp_proj_files[i])
  sp_proj <- cbind(bbs_xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
  sp_proj$presence = ifelse(sp_proj$occupancy_rate == 0, 0, 1)
  
  # deal with NAs
  sp_proj[is.na(sp_proj$occupancy_rate), c('occupancy_rate', 'abundance_log', 'presence')] <- 0
  
  # convert back to actual numbers
  occupancy_rate <- as.numeric(sp_proj$occupancy_rate/1000)
  abundance      <- as.numeric(exp(sp_proj$abundance_log/1000)-1)
  
  # turn absences to 0
  occupancy_rate[sp_proj$presence == 0] <- 0
  abundance[sp_proj$presence == 0] <- 0
  
  # rescale to percentiles. This asks, what proportion of the data fall above or below a certain value.
  occupancy_rate[sp_proj$presence!=0] <- ecdf(occupancy_rate[sp_proj$presence!=0])(occupancy_rate[sp_proj$presence!=0])
  abundance[sp_proj$presence!=0]      <- ecdf(abundance[sp_proj$presence!=0])(abundance[sp_proj$presence!=0])
  
  # aggregate richness, occupancy and abundance throughout loop
  bbs_xy[c('occupancy_rate')] <- bbs_xy[c('occupancy_rate')] + occupancy_rate
  bbs_xy[c('abundance')] <- bbs_xy[c('abundance')] + abundance
  bbs_xy[c('richness')] <- bbs_xy[c('richness')] + (sp_proj[c('presence')])
  bbs_xy[sp_proj$presence != 0, c('disconnect')] <- bbs_xy[sp_proj$presence != 0, c('disconnect')] + (abundance[sp_proj$presence!=0] - occupancy_rate[sp_proj$presence!=0])
  
  # find assemblage disconnect
  assemblage_disconnect[,i] <- abundance[sites] - occupancy_rate[sites]
  # remove disconnect when species arent present
  assemblage_disconnect[sp_proj$presence[sites]==0,i] <- NA

  print(nrow(na.omit(bbs_xy)))
  
}

# save output
# saveRDS(bbs_xy, 'results/spatial_projections_specieslevel/bbs_aggregated.RDS')
saveRDS(bbs_xy, 'results/spatial_projections_specieslevel/bbs_aggregated_V2.RDS') # v2 is with the disconnect estimated
saveRDS(assemblage_disconnect, 'results/spatial_projections_specieslevel/bbs_assemblage_disconnect.RDS')


# plot bbs aggregation ----

# read in aggregated properties
bbs_abun_occ <- readRDS('results/spatial_projections_specieslevel/bbs_aggregated_V2.RDS')

# remove empty grid cells
bbs_abun_occ <- bbs_abun_occ[-which(bbs_abun_occ$richness==0),]

# average to a species level
bbs_abun_occ$occupancy_rate <- bbs_abun_occ$occupancy_rate / bbs_abun_occ$richness
bbs_abun_occ$abundance      <- bbs_abun_occ$abundance / bbs_abun_occ$richness
bbs_abun_occ$disconnect     <- bbs_abun_occ$disconnect / bbs_abun_occ$richness


# create a spatial buffer of 2° around all sites within which to investigate biodiversity patterns
load("data/bbs_covariates.RData")
bbs_xy = bbs_xy[c('SiteLongitude', 
                  'SiteLatitude')]
bbs_xy_sp <- SpatialPoints(bbs_xy)
crs(bbs_xy_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
#bbs_xy_sp_buffer <- gBuffer(bbs_xy_sp, width = 2)
#save(bbs_xy_sp_buffer, file = 'figures/spatial_projections/bbs_polygons.RData')
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
# S = 1.8219e+16, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.05344556 

png(filename = 'figures/spatial_projections_specieslevel/aggregated_distributions/bbs/spatial_maps_V2.png', width = 2200, height = 4500, res = 300)
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
  
  ggplot(data = st_as_sf(usa_2)) +
    geom_tile(data = bbs_abun_occ, aes(y = SiteLatitude, x = SiteLongitude, fill = disconnect)) + 
    geom_sf(col = 'black', fill= 'transparent', lwd = 0.2) + 
    xlim(range(bbs_abun_occ$SiteLongitude)) + 
    ylim(range(bbs_abun_occ$SiteLatitude)) + 
    theme(legend.position = c(0.15, 0.1)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(min(bbs_abun_occ$disconnect), -0.05, -0.001, 
                                    0, 
                                    0.001, 0.05, max(bbs_abun_occ$disconnect)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()),
  
  ggplot(bbs_abun_occ) + 
    geom_histogram(data =bbs_abun_occ, aes(x = disconnect, fill = ..x..)) + 
    scale_fill_gradientn(colours = c(blue2, blue1, mid,  mid, mid, red1, red2),
                         values = c(min(bbs_abun_occ$disconnect), -0.05, -0.001, 
                                    0, 
                                    0.001, 0.05, max(bbs_abun_occ$disconnect)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) +
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          legend.position = 'none', 
          plot.title = element_text(vjust = -8, hjust = 0.05, size = 10)), 
  
  nrow = 4, 
  ncol = 1
  
)
dev.off()

# plot bbs species level within 1000 sample assemblages disconnection ----

bbs_assemblage_disconnect <- readRDS('results/spatial_projections_specieslevel/bbs_assemblage_disconnect.RDS')

# re-order all values and remove 0s to NAs
bbs_assemblage_disconnect[is.na(bbs_assemblage_disconnect)] <- 0
for(i in 1:ncol(bbs_assemblage_disconnect)){
  bbs_assemblage_disconnect[,i] <- sort(bbs_assemblage_disconnect[,i])
  #bbs_assemblage_disconnect[bbs_assemblage_disconnect[,i] == 0,i] <- 0
}

# convert to long dataframe
bbs_assemblage_disconnect <- data.frame(bbs_assemblage_disconnect)
bbs_assemblage_disconnect$community <- 1:nrow(bbs_assemblage_disconnect)
bbs_assemblage_disconnect <- pivot_longer(bbs_assemblage_disconnect, cols = X1:X1009)
bbs_assemblage_disconnect_sorted <- bbs_assemblage_disconnect %>% 
  # estimate ordered values and species per community
  group_by(community) %>% 
  nest() %>% 
  mutate(value_sort = purrr::map(data, ~sort(.$value, na.last = T)), 
         row_sort   = purrr::map(data, ~1:length(.$value))) %>% 
  unnest(value_sort, row_sort) %>% 
  ungroup() %>% 
  group_by(community) %>% 
  nest() %>% 
  # estimate mean value overall for a community
  mutate(mean_value = purrr::map(data, ~mean(.$value_sort))) %>% 
  unnest(mean_value) %>% 
  ungroup() %>% 
  unnest(data)

png(filename = 'figures/spatial_projections_specieslevel/assemblage_disconnect/bbs_assemblage_disconnect.png', width = 1200, height = 1200, res = 300)
ggplot(bbs_assemblage_disconnect_sorted) + 
  geom_line(aes(x = row_sort, y = value_sort, group = community, col = mean_value), lwd = 0.2, alpha = 0.5) + 
  geom_line(data = bbs_assemblage_disconnect_sorted %>% 
              filter(community %in% round(seq(1, 1000, length.out=50))), 
            aes(x = row_sort, y = value_sort, group = community, col = mean_value), lwd = 0.5) + 
  theme_bw() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1, 
        legend.position = 'none') + 
  scale_colour_gradientn(colours = c(viridis(10, option = 3)[2], viridis(10, option = 3)[4], 'gray25', 'black', 'gray25', viridis(10, option = 3)[6], viridis(10, option = 3)[9]),
                         values = c(min(bbs_assemblage_disconnect_sorted$mean_value), -0.25, -0.01, 
                                    0, 
                                    0.01, 0.25, max(bbs_assemblage_disconnect_sorted$mean_value)),
                         rescaler = function(x,...) x,
                         name = NULL, guide = guide_colourbar(direction = 'horizontal')) + 
  ylab('species disconnect') + 
  xlab('rank order')
dev.off()

# summary statistics of assemblage disconnect

# what proportion of assemblages have mean disconnects < - or + certain percentage points
bbs_assemblage_disconnect %>% 
  group_by(community) %>% 
  filter(value != 0) %>% 
  do(mean_value = mean(.$value)) %>% 
  dplyr::select(mean_value, community) %>% 
  unnest(mean_value) %>% 
  unique() %>% 
  mutate(mean_value = abs(mean_value)) %>%
  do(percentage_0.05 = nrow(.[.$mean_value > 0.05,])/1000,
     percentage_0.1 = nrow(.[.$mean_value > 0.1,])/1000, 
     percentage_0.25 = nrow(.[.$mean_value > 0.25,])/1000, 
     percentage_0.5 = nrow(.[.$mean_value > 0.5,])/1000) %>% 
  unnest()
#percentage_0.05 percentage_0.1 percentage_0.25 percentage_0.5
#<dbl>          <dbl>           <dbl>          <dbl>
# 0.864          0.725           0.395          0.123


# is every community sensitive?
# what percentage of species in a community are disconnected?
bbs_wide_disconnect <- readRDS('results/spatial_projections_specieslevel/bbs_assemblage_disconnect.RDS')

# convert species disconnects to absolute values
bbs_wide_disconnect <- abs(bbs_wide_disconnect)

disconnect_vectors <- sapply(seq(from = 0, to = 1, by = 0.01), function(x){ 
  percentage_assemblage_disconnected <- c()
  for(i in 1:nrow(bbs_wide_disconnect)){
    ass_dis <- bbs_wide_disconnect[i,]
    richness <- sum(!is.na(ass_dis))
    more_than_25 <- sum(ass_dis > c(x), na.rm = T)
    percentage_assemblage_disconnected[i] <- more_than_25 / richness
  }
  return(quantile(percentage_assemblage_disconnected, c(0.05, 0.5, 0.95)))
})

disconnect_vectors_df  <- gather(data.frame(disconnect_vectors))
disconnect_vectors_df$quantiles   <- rownames(disconnect_vectors)
disconnect_vectors_df$break_point <- sort(rep(seq(from = 0, to = 1, by = 0.01), 3))
disconnect_vectors_df <- pivot_wider(disconnect_vectors_df, names_from = quantiles, values_from = value)


bbs_disconnect_breakpoints <- ggplot(data = disconnect_vectors_df) + 
  geom_ribbon(aes(x = break_point, ymin=`5%`, ymax=`95%`), alpha=0.5, 
              fill = viridis::viridis(10,option = 3)[5]) + 
  geom_line(aes(x = break_point,
                y = `50%`)) + 
  geom_line(aes(x = break_point, y=`5%`), alpha=1, lty = 2) + 
  geom_line(aes(x = break_point, y=`95%`), alpha=1, lty= 2) + 
  #geom_segment(data = disconnect_vectors_df[disconnect_vectors_df$break_point %in% c(0.25, 0.5, 0.75),], 
  #             aes(x=break_point, y=`50%`, xend=c(0.25, 0.5, 0.75), yend=0), color=c("gray25", 'gray50', 'gray75'))  + 
  geom_segment(data = disconnect_vectors_df[disconnect_vectors_df$break_point %in% c(0.25, 0.5, 0.75),], 
               aes(x=break_point, y=`50%`, xend=c(0,0,0), yend=`50%`), color=c("gray25", 'gray50', 'gray75')) + 
  xlab('disconnection threshold') + 
  ylab('% species in assemblage disconnected') + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)

png(file = 'figures/spatial_projections_specieslevel/assemblage_disconnect/bbs_breakpoint.png', width = 1200, height = 1200, res = 300)
bbs_disconnect_breakpoints
dev.off()

# whats the mean disconnect for bbs ---- 

bbs_abun_occ %>% 
  filter(abundance_residual >= 0) %>% 
  do(mean = mean(.$abundance_residual), 
     median = median(.$abundance_residual),
     sd   = sd(.$abundance_residual)) %>% unnest()
# mean  median      sd
# <dbl>   <dbl>   <dbl>
# 0.0110 0.00883 0.00931

bbs_abun_occ %>% 
  filter(abundance_residual <= 0) %>% 
  do(mean = mean(.$abundance_residual), 
     median = median(.$abundance_residual),
     sd   = sd(.$abundance_residual)) %>% unnest()
# A tibble: 1 x 3
# mean   median      sd
# <dbl>    <dbl>   <dbl>
# -0.0105 -0.00856 0.00811

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

# create stratified sample
bbs_abun_occ_sampled <- bbs_abun_occ %>% 
  mutate(sample_cut = cut2(.$abundance_residual, g = 5, onlycuts = F)) %>% 
  group_by(sample_cut) %>% 
  nest() %>% 
  mutate(samples = purrr::map(data, ~.[sample(1:nrow(.), 10000, replace = T),])) %>% 
  dplyr::select(-data) %>% 
  unnest(samples) %>% ungroup()

# projection data
bbs_abun_occ_positive_residual <- bbs_abun_occ_sampled %>% filter(abundance_residual >= 0)
#bbs_abun_occ_positive_residual <- bbs_abun_occ_positive_residual[sample(1:nrow(bbs_abun_occ_positive_residual), 5000),]
bbs_abun_occ_negative_residual <- bbs_abun_occ_sampled %>% filter(abundance_residual <= 0)
#bbs_abun_occ_negative_residual <- bbs_abun_occ_negative_residual[sample(1:nrow(bbs_abun_occ_negative_residual), 5000),]

# directories and names
base_dir <- 'figures/spatial_projections/pdp_varimp'
dataset  <- 'bbs' 

# run pdp functions
# partial_dependence_plots(covariates = bbs_covariates,
#                          cov_names  = bbs_cov_names, 
#                          projections = bbs_abun_occ_positive_residual,
#                          base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
#                          dataset = 'bbs', 
#                          name = 'positive_abundance',
#                          option = ifelse(dataset == 'rls', 'viridis', 3))
# 
# partial_dependence_plots(covariates = bbs_covariates,
#                          cov_names  = bbs_cov_names, 
#                          projections = bbs_abun_occ_negative_residual,
#                          base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
#                          dataset = 'bbs', 
#                          name = 'negative_abundance',
#                          option = ifelse(dataset == 'rls', 'viridis', 3))

partial_dependence_plots(covariates = bbs_covariates,
                         cov_names  = bbs_cov_names, 
                         response_column = 'disconnect',
                         projections = bbs_abun_occ_sampled,
                         base_dir = 'figures/spatial_projections_specieslevel/pdp_varimp', 
                         dataset = 'bbs', 
                         name = 'all_residuals',
                         option = ifelse(dataset == 'rls', 'viridis', 3), 
                         N = 6)
# Mean of squared residuals: 3.92872e-05
# % Var explained: 99.39





