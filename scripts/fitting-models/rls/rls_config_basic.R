# script to fit each of the functions to species matrix

# load in packages ---- 
lib_vect <- c('tidyverse')
sapply(lib_vect,require,character=TRUE)

# get trailing arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) 
print(i)

# load in data ----

# load in abundance data
load("data/rls_abun_modelling_data_V2.RData")
abundance_input  = rls_abun_fitting[[i]]    # change this in all the functions
validation_input = rls_abun_validation[[i]] # change this in all the functions
print(unique(abundance_input$TAXONOMIC_NAME))

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

# setup number of bootstraps 
n_boots = 10

# set up dataset save and base-directory for file saves
dataset = 'rls'
base_dir        = 'results/rls_basic'

source('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/scripts/fitting-models/fit_all_models.R')








