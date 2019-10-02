# script to fit each of the functions to species matrix

# load in packages ---- 
lib_vect <- c('tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)


# load in data ----

# load in abundance data
load("data/rls_abun_modelling_data_v2.RData")
abundance = rls_abun_fitting

# load in covariates
load("data/rls_covariates.RData")
covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
                      'Depth_GEBCO_transformed', 
                      "human_pop_2015_50km", "reef_area_200km", "wave_energy_mean", "Depth_GEBCO_transformed",
                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# ensure that all the focal variables are transformed to be scaled and 
# transformed if necessary

# run glm options ----
source('scripts/model-functions/glm_function.R')

# lm, log, no zero inflation
lapply(abundance, 
       FUN = function(x){
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = 'log', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2')})


# lm, log10, no zero inflation
lapply(abundance, 
       FUN = function(x){
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = 'log10', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2')})


# glm, poisson, no zero inflation
lapply(abundance, 
       FUN = function(x){
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2')})


# glm, nbinom, no zero inflation
lapply(abundance, 
       FUN = function(x){
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2')})


# glm, nbinom, no zero inflation
lapply(abundance, 
       FUN = function(x){
               glm_function(abundance = x,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'glm', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                            zi = F,     # option is F or T
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})


# glm, poisson, with zero inflation
lapply(abundance, 
       FUN = function(x){
         print(x)
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2')})


# glm, nbinom, with zero inflation
lapply(abundance, 
       FUN = function(x){
         tryCatch(
         glm_function(abundance = x,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = unique(x$TAXONOMIC_NAME), 
                      base_dir = 'results/test2'), 
         error = function(e) NA)})


# glm, nbinom, with zero inflation
lapply(abundance, 
       FUN = function(x){
               tryCatch(
                       glm_function(abundance = x,
                                    covariates = covariates,
                                    transformation = NA, # option is NA, log, log10
                                    model = 'glm', # option is lm or glm
                                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                                    zi = T,     # option is F or T
                                    species_name = unique(x$TAXONOMIC_NAME), 
                                    base_dir = 'results/test2'), 
                       error = function(e) NA)})


# run gam options ----
source('scripts/model-functions/gam_function.R')

# gam, log
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = 'log', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})

# gam, log10
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = 'log10', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})

# gam, raw, poisson
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})

# gam raw nbinom
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})

# gam raw tweedie
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})

# gam raw zip
lapply(abundance, 
       FUN = function(x){
               gam_function(abundance = x,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = unique(x$TAXONOMIC_NAME), 
                            base_dir = 'results/test2')})


# run random forest options ----
source('scripts/model-functions/rf_function.R')

# rf, raw, continuous
lapply(abundance, 
       FUN = function(x){
         rf_function(abundance = x, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     base_dir = 'results/test2')})

# rf, log, continuous
lapply(abundance, 
       FUN = function(x){
         rf_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     base_dir = 'results/test2')})

# rf, log, discrete
lapply(abundance, 
       FUN = function(x){
         rf_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     base_dir = 'results/test2')})

# rf, log10, continuous
lapply(abundance, 
       FUN = function(x){
         rf_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     base_dir = 'results/test2')})

# rf, log10, discrete
lapply(abundance, 
       FUN = function(x){
         rf_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     base_dir = 'results/test2')})



# run boosted regression trees options ----

# source function
source('scripts/model-functions/brt_function.R')

# brt, raw, continuous
lapply(abundance, 
       FUN = function(x){
         brt_function(abundance = x, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     n.cores = 10,
                     base_dir = 'results/test2')})

# brt, log, continuous
lapply(abundance, 
       FUN = function(x){
         brt_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     n.cores = 10,
                     base_dir = 'results/test2')})

# brt, log, discrete
lapply(abundance, 
       FUN = function(x){
         brt_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     n.cores = 10,
                     base_dir = 'results/test2')})

# need to run
# brt, log10, continuous
lapply(abundance, 
       FUN = function(x){
         brt_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     n.cores = 10,
                     base_dir = 'results/test2')})

# need to run
# brt, log10, discrete
lapply(abundance, 
       FUN = function(x){
         brt_function(abundance = x, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = unique(x$TAXONOMIC_NAME), 
                     n.cores = 10,
                     base_dir = 'results/test2')})











