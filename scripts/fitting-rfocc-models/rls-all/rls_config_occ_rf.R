# script to fit each of the functions to species matrix

# load in packages ---- 
lib_vect <- c('tidyverse', 'raster')
sapply(lib_vect,require,character=TRUE)

# get trailing arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) 
print(i)

# load in data ----

# load in abundance data
rls_abun <- readRDS(list.files('data/rls_all_basic', full.names = T)[i])

# load in abundance data
abundance_input  <- rls_abun[[1]]
print(unique(abundance_input$TAXONOMIC_NAME))

# for testing scripts 
# abundance  <- abundance_input
# validation <- validation_input

# load in covariates
load("data/rls_covariates.RData")

# select covariates
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

# human, reef and wave energy are scaled in the 02_rls-environmental-data script. 
covariates[c("Depth_GEBCO_transformed",
             'robPCA_1', 
             'robPCA_2', 
             'robPCA_3',
             'sst_mean')] <- as.numeric(scale(covariates[c("Depth_GEBCO_transformed",
                                                           'robPCA_1', 
                                                           'robPCA_2', 
                                                           'robPCA_3',
                                                           'sst_mean')]))

# read in gridded covariate data for rls 
# not implemented in function yet

# setup number of bootstraps 
n_boots = 10

# set up dataset save and base-directory for file saves
dataset = 'rls'
base_dir = 'results/variable_importance'
spatial_dir = 'results/spatial_projections'

# spatial projection covariates
spatial_projections <- na.omit(readRDS('data/rls_spatial_projection_data.rds'))

# call script to call function to fit model
source('scripts/fitting-rfocc-models/fit-models/fit_occ_rf.R')
