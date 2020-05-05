# script to fit each of the functions to species matrix

# load in packages ---- 

lib_vect <- c('tidyverse')
sapply(lib_vect,require,character=TRUE)


# find the species that are missing from suitability

suitability <- list.files('/Volumes/Simulation/conor/abundance-comparison/results/bbs_basic/suitability')

missing <- which(!grepl(gsub('.RData','',paste(suitability, collapse = '|')), abundance_key$TAXONOMIC_NAME))

i = missing[1]

# load in data ----

# load in abundance data
load("data/bbs_abun_modelling_data.RData")
abundance_input  = bbs_abun_fitting[[i]]    # change this in all the functions
validation_input = bbs_abun_validation[[i]] # change this in all the functions
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
dataset  = 'bbs'
base_dir = 'results/bbs_basic'

source('scripts/fitting-models/fit-models/fit_all_models.R')

