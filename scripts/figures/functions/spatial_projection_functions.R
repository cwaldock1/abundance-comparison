
# function to map spatial projections ----

plot_distributions <- function(xy, 
                               species_name, 
                               dataset, 
                               save_dir){
  
  require(raster)
  require(sp)
  require(sf)
  require(rgeos)
  require(rgdal)
  require(rnaturalearth)
  require(viridis)
  
  # get raw species abundances
  species_data <- readRDS(list.files(paste0('data/', dataset, '_all_basic/'), pattern = paste0(species_name, '.RDS'), full.names = T))
  species_data <- species_data$fitting
  
  # get range convex hull
  species_hull <- readRDS(list.files(paste0('data/species_range_convex_hull/', dataset), pattern = paste0(species_name, '.RDS'), full.names = T))
  
  # get species projections
  species_projections <- readRDS(list.files('results/spatial_projections_cropped/', pattern = species_name, full.names = T, recursive = T))
  species_proj_2 <- species_projections/1000#species_projections[[1]]/1000
  species_proj_2$presence <- 1 #ifelse(species_proj_2$occupancy_rate > species_projections[[2]]$threshold.TSS.suitability, 1, 0) # perform threshold
  # add in xy
  species_proj_2 <- cbind(species_proj_2, na.omit(xy))
  # remove NAs 
  species_proj_2 <- na.omit(species_proj_2)
  
  # # convert to spatial points data frame
  # pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(species_proj_2$SiteLongitude, species_proj_2$SiteLatitude))),
  #                                    species_proj_2[,1:3],
  #                                    proj4string = species_hull@proj4string)
  # crs(pointsDF) <- species_hull@proj4string
  # species_proj_3 <- species_proj_2[!is.na(over(pointsDF, species_hull)), ]
  
  # get maximum and minimum of hull or ranges for xy lims
  ch_range <- rbind(matrix(extent(species_hull))[1:2,], matrix(extent(species_hull))[3:4,]) 
  sp_range <- unlist(rbind(range(species_proj_2$SiteLongitude), range(species_proj_2$SiteLatitude)))
  
  # minimum long
  min_long <- min(c(extent(species_hull)[1,], min(range(species_proj_2$SiteLongitude))))
  # maximum long
  max_long <- max(c(extent(species_hull)[2,], max(range(species_proj_2$SiteLongitude))))
  # minimum lat
  min_lat <- min(c(extent(species_hull)[3,], min(range(species_proj_2$SiteLatitude))))
  # maximum lat
  max_lat <- max(c(extent(species_hull)[4,], max(range(species_proj_2$SiteLatitude))))

  # create map of spatial countries
  world <- ne_countries(scale = "large", returnclass = "sf")  
  
  # create set of species level maps
  dir.create(paste0(save_dir, '/', dataset), recursive = T)
  png(paste0(paste0(save_dir, '/', dataset), '/', species_name, '.png'), height = 1000, width = 2000)
  grid.arrange(
    
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_polygon(data = species_hull, aes(x = long, y = lat), col = 'black', fill = 'light blue', alpha = 0.5) + 
    geom_tile(data = species_data %>% filter(Num > 0), aes(y = SiteLatitude, x = SiteLongitude, fill = NULL)) + 
    geom_point(data = species_data %>% filter(Num > 0), aes(y = SiteLatitude, x = SiteLongitude)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    ggtitle('survey presence') + 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_data %>% filter(Num == 0), aes(y = SiteLatitude, x = SiteLongitude, fill = NULL)) + 
    geom_point(data = species_data %>% filter(Num == 0), aes(y = SiteLatitude, x = SiteLongitude), col = 'red') + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    ggtitle('survey absence')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = as.factor(presence))) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    scale_fill_manual(values = c('white', 'black'), name = NULL) +
    ggtitle('presence - threshold')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = occupancy_rate)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, 
                                guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle('occupancy rate')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance_log)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, 
                                guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle('log abundance')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = exp(abundance_log)-1)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, 
                                guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle('abundance')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = abundance_log*occupancy_rate)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, 
                                guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle('log abundance * occupancy')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()),
  
  ggplot(data = world) +
    geom_sf(col = 'grey70', fill= 'gray70', lwd = 0.001) + 
    geom_tile(data = species_proj_2, aes(y = SiteLatitude, x = SiteLongitude, fill = (exp(abundance_log)-1)*occupancy_rate)) + 
    xlim(c(min_long, max_long)) + 
    ylim(c(min_lat, max_lat)) + 
    theme(legend.position = c(0.2, 0.1)) + 
    viridis::scale_fill_viridis(option = "magma", name = NULL, 
                                guide = guide_colourbar(direction = 'horizontal')) +
    ggtitle('abundance * occupancy')+ 
    theme(panel.background = element_blank(), 
          panel.border = element_rect(colour = 'black', fill = 'transparent'), 
          axis.title   = element_blank()), 
  
  ncol = 4
  
  )
  
  dev.off()
  
}

# function for species level spearmans rank correlations ----

#spatial_dir = 'results/spatial_projections/rls'
#xy = na.omit(readRDS('data/rls_spatial_projection_data.rds'))
#dataset = 'rls'
 
species_spearmans_spatial <- function(spatial_dir, 
                                      xy, 
                                      dataset){
  
  require(sp)
  require(rgeos)
  
  sp_correlations <- parallel::mclapply(list.files(spatial_dir, full.names = T), mc.cores = 10, function(x){
    
    print(x)
    
    # read in spatial projections and align to xy
    sp_proj <- readRDS(x)[[1]]
    sp_proj$species_name <- gsub('.RDS', '', str_split(x, '/')[[1]][4])
    sp_proj <- cbind(xy[c('SiteLongitude', 'SiteLatitude')], sp_proj)
    
    # read in convex hull for the species
    sp_ch <- readRDS(paste0('data/species_range_convex_hull/', dataset, '/', unique(sp_proj$species_name), '.RDS'))
    
    # convert to spatial points data frame and crop to convex hull
    pointsDF <- SpatialPointsDataFrame(SpatialPoints(coordinates(cbind(sp_proj$SiteLongitude, sp_proj$SiteLatitude))),
                                       sp_proj[,c('occupancy_rate', 'abundance_log')],
                                       proj4string = sp_ch@proj4string)
    crs(pointsDF) <- sp_ch@proj4string
    sp_proj_2 <- sp_proj[!is.na(over(pointsDF, sp_ch)), ]
    
    # estimate rank correlation coefficient for species
    cor_rank <- cor(sp_proj_2$occupancy_rate, sp_proj_2$abundance_log, method = 'spearman')
    
    return(data.frame(species_name = unique(sp_proj$species_name), 
               cor          = cor_rank))
    
    })
  
   sp_correlations_2 <- do.call(rbind, sp_correlations)
  
   return(sp_correlations_2)
   
}


# function to create partial dependence plots explaining species projections ----


partial_dependence_plots <- function(covariates,
         cov_names, 
         response_column,
         projections,
         base_dir, 
         dataset, 
         name,
         N = 3,
         option # viridis option
         ){
  
  require(randomForest)
  require(DALEX)
  # require(edarf) # had to be built from source
  
  # join together spatial projections and covariates
  comb_proj <- left_join(projections, covariates)
  
  # response data
  response_y <- data.frame(comb_proj[,which(colnames(comb_proj) == response_column)])[,1]
  
  # create formula for random forest
  rf_spatial <- randomForest(x = comb_proj[,which(colnames(comb_proj) %in% cov_names)], 
                             y = response_y, 
                             ntree = 1000, 
                             importance=T, 
                             keep.inbag=T)
  
  # combine predictions to dataset
  comb_proj <- comb_proj %>% 
    mutate(predicted = predict(rf_spatial))
  
  # exploring using DALEX for visualisation
  exp_rf_spatial <- explain(model = rf_spatial, 
                            data  = comb_proj[,which(colnames(comb_proj) %in% cov_names)], 
                            y = response_y)
  
  # calculate variable importance
  vip <- model_parts(exp_rf_spatial, 
                     loss_function = loss_root_mean_square, 
                     n_sample = 1000)
  vip_plot <- ggplot_build(plot(vip) +
                             ggtitle(NULL, ' ')+
                             ylab('root mean square loss') + 
                             scale_x_discrete(position = "top") +
                             theme_bw() + 
                             theme(panel.grid = element_blank(), 
                                   strip.background = element_blank(), 
                                   strip.text = element_blank(), 
                                   legend.position = 'none'))
  vip_plot$data[[1]]$linetype <- 'solid'
  vip_plot$data[[1]]$colour <- 'gray50'
  vip_plot$data[[2]]$colour <- viridis::viridis(10,option = option)[5]
  vip_plot$data[[3]]$colour <- 'gray75'
  vip_plot$data[[3]]$fill <- viridis::viridis(10,option = option)[5]
  vip_plot2 <- ggplot_gtable(vip_plot)
  vip_plot2 <- ggplotify::as.ggplot(vip_plot2)
  
  # select most important 3
  vars <- vip %>% 
    filter(!variable %in% c('_full_model_', '_baseline_')) %>% 
    group_by(variable) %>% 
    do(dropout_loss = mean(.$dropout_loss)) %>% 
    unnest() %>% 
    .[order(.$dropout_loss, decreasing = T),] %>% 
    .[1:N,] %>% 
    .$variable %>% 
    as.character()
  
  # estimate partial dependence
  pdp_rf_spatial <- model_profile(explainer = exp_rf_spatial,
                                  variables = vars, 
                                  N = 100)
  
  # create partial dependence plots
  pdp_plot <- plot(pdp_rf_spatial, geom = 'profiles') + 
    theme_bw() + 
    scale_colour_viridis_d(1, option = option, begin = 0.1, end = 0.9) + 
    theme(panel.grid = element_blank(), 
          strip.text = element_text(hjust = 0), 
          strip.background = element_blank(), 
          legend.position = 'none', 
          panel.background = element_rect(fill = 'transparent'), 
          panel.ontop = element_blank(), 
          plot.background = element_blank()) + 
    ylab('abundance residual') + 
    ggtitle(NULL, NULL)
  
  pdp_plot_build <- ggplot_build(pdp_plot)
  pdp_plot_build$data[[2]]$colour <- viridis::viridis(10,option = option)[5]
  pdp_plot_build <- ggplot_gtable(pdp_plot_build)
  pdp_plot2 <- ggplotify::as.ggplot(pdp_plot_build)
  
  
  # save outputs
  model_final_path <- paste0(base_dir, '/', dataset)
  dir.create(model_final_path, recursive = T)
  png(filename = paste0(model_final_path, '/', name, '_partial_dependence_plot.png'), width = 3000, height = ifelse(N !=3, 2500/2, 2500/4), res=300)
  grid.arrange(pdp_plot2 + xlab('standardized values'), vip_plot2, 
               nrow = 1, 
               widths = c(3,1))
  dev.off()
  
  return(rf_spatial)
  
}



