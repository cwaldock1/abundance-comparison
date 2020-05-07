# script to fit each of the functions to species matrix

# load in packages ---- 

lib_vect <- c('tidyverse')
sapply(lib_vect,require,character=TRUE)

# get trailing arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) 
print(i)

# load in data ----

bbs_abun <- readRDS(list.files('data/bbs_all_basic', full.names = T)[i])

# load in abundance data
#load("data/bbs_abun_modelling_data.RData")
#abundance_input  = bbs_abun_fitting[[i]]    # change this in all the functions
#validation_input = bbs_abun_validation[[i]] # change this in all the functions
#print(unique(abundance_input$TAXONOMIC_NAME))

# load in abundance data
abundance_input  <- bbs_abun[[1]]
validation_input <- bbs_abun[[2]]
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
base_dir = 'results/bbs_basic_all'

source('scripts/fitting-models/fit-models/fit_all_models.R')

