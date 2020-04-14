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
load("data/rls_abun_modelling_data.RData")
abundance = rls_abun_fitting[,c(1:3, i+3)]
print(names(abundance))

# load in covariates
load("data/rls_covariates.RData")
covariates = rls_xy[c('SiteCode', 'SiteLongitude', 'SiteLatitude',
                      'human_pop_2015_50km', 'reef_area_200km', 'wave_energy_mean', 'Depth_GEBCO_transformed',
                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# ensure that all the focal variables are transformed to be scaled and 
# transformed if necessary

# run glm options ----
source('scripts/model-functions/glm_function.R')

# lm, log, no zero inflation
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = 'log', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits')


# lm, log10, no zero inflation
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = 'log10', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits')


# glm, poisson, no zero inflation
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits')


# glm, nbinom, no zero inflation
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits')


# glm, nbinom, no zero inflation
               glm_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'glm', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                            zi = F,     # option is F or T
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')


# glm, poisson, with zero inflation
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits')


# glm, nbinom, with zero inflation
         tryCatch(
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = as.character(names(abundance)[4]), 
                      base_dir = 'results/rls_modelfits'), 
         error = function(e) NA)


# glm, nbinom, with zero inflation
               tryCatch(
                       glm_function(abundance = abundance,
                                    covariates = covariates,
                                    transformation = NA, # option is NA, log, log10
                                    model = 'glm', # option is lm or glm
                                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                                    zi = T,     # option is F or T
                                    species_name = as.character(names(abundance)[4]), 
                                    base_dir = 'results/rls_modelfits'), 
                       error = function(e) NA)


# run gam options ----
source('scripts/model-functions/gam_function.R')

# gam, log
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = 'log', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')

# gam, log10
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = 'log10', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')

# gam, raw, poisson
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')

# gam raw nbinom
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')

# gam raw tweedie
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')

# gam raw zip
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[4]), 
                            base_dir = 'results/rls_modelfits')


# run random forest options ----
source('scripts/model-functions/rf_function.R')

# rf, raw, continuous
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     base_dir = 'results/rls_modelfits')

# rf, log, continuous
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     base_dir = 'results/rls_modelfits')

# rf, log, discrete
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     base_dir = 'results/rls_modelfits')

# rf, log10, continuous
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     base_dir = 'results/rls_modelfits')

# rf, log10, discrete
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     base_dir = 'results/rls_modelfits')



# run boosted regression trees options ----

# source function
source('scripts/model-functions/brt_function.R')

# brt, raw, continuous
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# brt, log, continuous
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# brt, log, discrete
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# need to run
# brt, log10, continuous
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# need to run
# brt, log10, discrete
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')











