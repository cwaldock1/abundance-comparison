# script to fit each of the functions to species matrix

# load in packages ---- 
lib_vect <- c('tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)


# load in data ----

# load in abundance data
load("data/rls_abun_modelling_data.RData")
abundance = rls_abun_fitting

# load in covariates
load("data/rls_covariates.RData")
covariates = rls_xy[c('SiteCode', 'SiteLongitude', 'SiteLatitude',
                      'Depth_GEBCO_transformed', 
                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# ensure that all the focal variables are transformed to be scaled and 
# transformed if necessary

# run glm options ----
source('scripts/model-functions/glm_function.R')

# lm, log, no zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = 'log', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits')})


# lm, log10, no zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = 'log10', # option is NA, log, log10
                      model = 'lm', # option is lm or glm
                      family = NA, # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits')})


# glm, poisson, no zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits')})


# glm, nbinom, no zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = F,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits')})


# glm, nbinom, no zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               glm_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'glm', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                            zi = F,     # option is F or T
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})


# glm, poisson, with zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         print(x)
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits')})


# glm, nbinom, with zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         tryCatch(
         glm_function(abundance = abundance,
                      covariates = covariates,
                      transformation = NA, # option is NA, log, log10
                      model = 'glm', # option is lm or glm
                      family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                      zi = T,     # option is F or T
                      species_name = as.character(names(abundance)[x+3]), 
                      base_dir = 'results/modelfits'), 
         error = function(e) NA)})


# glm, nbinom, with zero inflation
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               tryCatch(
                       glm_function(abundance = abundance,
                                    covariates = covariates,
                                    transformation = NA, # option is NA, log, log10
                                    model = 'glm', # option is lm or glm
                                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                                    zi = T,     # option is F or T
                                    species_name = as.character(names(abundance)[x+3]), 
                                    base_dir = 'results/modelfits'), 
                       error = function(e) NA)})


# run gam options ----
source('scripts/model-functions/gam_function.R')

# gam, log
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = 'log', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})

# gam, log10
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = 'log10', # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})

# gam, raw, poisson
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})

# gam raw nbinom
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})

# gam raw tweedie
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})

# gam raw zip
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
               gam_function(abundance = abundance,
                            covariates = covariates,
                            transformation = NA, # option is NA, log, log10
                            model = 'gam', # option is lm or glm
                            family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                            species_name = as.character(names(abundance)[x+3]), 
                            base_dir = 'results/modelfits')})


# run random forest options ----
source('scripts/model-functions/rf_function.R')

# rf, raw, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     base_dir = 'results/modelfits')})

# rf, log, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     base_dir = 'results/modelfits')})

# rf, log, discrete
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     base_dir = 'results/modelfits')})

# rf, log10, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     base_dir = 'results/modelfits')})

# rf, log10, discrete
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         rf_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     base_dir = 'results/modelfits')})



# run boosted regression trees options ----

# source function
source('scripts/model-functions/brt_function.R')

# brt, raw, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     n.cores = 10,
                     base_dir = 'results/modelfits')})

# brt, log, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     n.cores = 10,
                     base_dir = 'results/modelfits')})

# brt, log, discrete
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     n.cores = 10,
                     base_dir = 'results/modelfits')})

# need to run
# brt, log10, continuous
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     n.cores = 10,
                     base_dir = 'results/modelfits')})

# need to run
# brt, log10, discrete
lapply(1:(length(colnames(abundance))-3), 
       FUN = function(x){
         brt_function(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[x+3]), 
                     n.cores = 10,
                     base_dir = 'results/modelfits')})











