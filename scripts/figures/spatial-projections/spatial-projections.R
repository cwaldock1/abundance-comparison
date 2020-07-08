# produce folder of figures for species level spatial projections 

# load packages ----

lib_vect <- c('tidyverse', 'ggplot2', 'rnaturalearth', 'gridExtra', 'raster')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load in functions ----

source('scripts/figures/functions/spatial_projection_functions.R')

# load in rls spatial data ----

# rls all spatial information
rls_xy <- readRDS('data/rls_spatial_projection_data.rds')

# rls test species
#species_name <- 'Chromis_cinerascens'

# rls list all species
all_sp <- gsub('.RDS', '', list.files('results/spatial_projections/rls'))

# test function
lapply(1:length(all_sp), function(x)
  plot_distributions(xy           = rls_xy, 
                   species_name = all_sp[x], 
                   dataset      = 'rls', 
                   save_dir     = 'figures/spatial_projections/species_distributions'))

# bbs figures ----

# rls all spatial information
bbs_xy <- readRDS('data/bbs_spatial_projection_data.rds')

# bbs list all species
all_sp <- gsub('.RDS', '', list.files('results/spatial_projections/bbs'))

# test function
lapply(1:length(all_sp), function(x)
  plot_distributions(xy           = bbs_xy, 
                     species_name = all_sp[x], 
                     dataset      = 'bbs', 
                     save_dir     = 'figures/spatial_projections/species_distributions'))


