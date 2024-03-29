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
load("data/bbs_abun_modelling_data_oob_crossValidation.RData")
abundance_input  = bbs_abun_cv_fitting[[i]]    # change this in all the functions
validation_input = bbs_abun_cv_validation[[i]] # change this in all the functions
print(unique(abundance_input$TAXONOMIC_NAME))

# load in covariates
load("data/bbs_covariates.RData")

# select covariates
covariates = bbs_xy[c('SiteLongitude', 
                      'SiteLatitude',
                      "robPCA_1", 
                      "robPCA_2", 
                      "robPCA_3", 
                      "human_pop",
                      'primary_forest', 
                      'Elevation_GEBCO')]

# setup number of bootstraps 
n_boots = 10

# set up dataset save and base-directory for file saves
dataset = 'bbs'
base_dir        = 'results/bbs_oob_cv'

source('scripts/fitting-models/fit-models/fit_all_models.R')

