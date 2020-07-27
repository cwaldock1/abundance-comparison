# script to fit each of the functions to species matrix

# load in packages ---- 

lib_vect <- c('tidyverse', 'raster')
sapply(lib_vect,require,character=TRUE)

# get trailing arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) 
print(i)

# load in data ----

bbs_abun <- readRDS(list.files('data/bbs_all_basic', full.names = T)[i])

# load in abundance data
abundance_input  <- bbs_abun[[1]]
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

# scale variables before modelling
# elevation and human pop got scaled earlier in 02_bbs-environmental-data
covariates[c("robPCA_1", 
             "robPCA_2", 
             "robPCA_3", 
             'primary_forest')] <- as.numeric(scale(covariates[c("robPCA_1", 
                                                                  "robPCA_2", 
                                                                  "robPCA_3", 
                                                                  'primary_forest')]))

# read in gridded covariate data for bbs 
# not implemented in function yet

# setup number of bootstraps 
n_boots = 10

# set up dataset save and base-directory for file saves
dataset = 'bbs'
base_dir = 'results/variable_importance'
spatial_dir = 'results/spatial_projections'

# spatial projection covariates
spatial_projections <- na.omit(readRDS('data/bbs_spatial_projection_data.rds'))

source('scripts/fitting-rfocc-models/fit-models/fit_occ_rf.R')

# # load script for making species specific figures
# source('scripts/figures/functions/spatial_projection_functions.R')
# 
# plot_distributions(xy = readRDS('data/bbs_spatial_projection_data.rds'), 
#                    species_name = gsub(' ','_',unique(abundance_input$TAXONOMIC_NAME)),
#                    dataset = dataset,
#                    save_dir = 'figures/spatial_projections/species_distributions')
# 

