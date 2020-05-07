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
rls_abun <- readRDS(list.files('data/rls_all_CV', full.names = T)[i])

# load in abundance data
#load("data/rls_abun_modelling_data.RData")
#abundance_input  = rls_abun_fitting[[i]]    # change this in all the functions
#validation_input = rls_abun_validation[[i]] # change this in all the functions
#print(unique(abundance_input$TAXONOMIC_NAME))

# load in abundance data
abundance_input  <- rls_abun[[1]]
validation_input <- rls_abun[[2]]
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
base_dir        = 'results/rls_cv_all'

source('scripts/fitting-models/fit-models/fit_all_models.R')








