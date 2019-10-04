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
load("data/rls_abun_modelling_data_v2.RData")
abundance_input  = rls_abun_fitting[[i]]    # change this in all the functions
validation_input = rls_abun_validation[[i]] # change this in all the functions
print(unique(abundance_input$TAXONOMIC_NAME))

# load in covariates
load("data/rls_covariates.RData")
covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
                      'Depth_GEBCO_transformed', 
                      "human_pop_2015_50km", "reef_area_200km", "wave_energy_mean", "Depth_GEBCO_transformed",
                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# run glm options ----
source('scripts/model-functions/glm_function_boot.R')

# list of potential function arguments
#abundance = abundance, 
#validation = validation,
#covariates = covariates, 
#transformation = NA, # option is NA, log, log10
#model = NA, # option is lm or glm
#family = NA, # option is NA, poisson, or nbinom, or tweedie
#zi = F,     # option is F or T
#species_name = species_name, 
#n_bootstrap = 10,
#dataset = 'rls',
#base_dir        = 'results/rls',
#model_path      = 'model', 
#prediction_path = 'predictions'

# lm, log, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# lm, log10, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# glm, poisson, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'NA', # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# glm, nbinom, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'NA', # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# glm, tweedie, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'NA', # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# glm, poisson, with zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'NA', # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = T,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# glm, nbinom, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = covariates,
                    transformation = 'NA', # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = 10,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model', 
                    prediction_path = 'predictions'),
  error = function(e) NA)


# glm, tweedie, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = covariates,
                    transformation = 'NA', # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = 10,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model', 
                    prediction_path = 'predictions'), 
  error = function(e) NA)


# run gam options ----
source('scripts/model-functions/gam_function_boot.R')

# gam, log
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')

# gam, log10
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')

# gam, raw, poisson
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')

# gam raw nbinom
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')

# gam raw tweedie
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')

# gam raw zip
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = 10,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model', 
                  prediction_path = 'predictions')


# run random forest options ----
source('scripts/model-functions/rf_function_boot.R')

# abundance = abundance, 
# validation = validation,
# covariates = covariates, 
# transformation = NA, # option is NA, log, log10
# discrete = NA,       # option is T or F 
# species_name = species_name, 
# n_bootstrap = 10,
# dataset = 'rls',
# base_dir        = 'results/rls',
# model_path      = 'model', 
# prediction_path = 'predictions'

# rf, raw, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = NA, # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = 10,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model', 
                 prediction_path = 'predictions')

# rf, log, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = 10,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model', 
                 prediction_path = 'predictions')

# rf, log, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = 10,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model', 
                 prediction_path = 'predictions')

# rf, log10, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = 10,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model', 
                 prediction_path = 'predictions')

# rf, log10, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = 10,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model', 
                 prediction_path = 'predictions')



# run boosted regression trees options ----
stop('boosted regression trees not implemented on server')
# source function
source('scripts/model-functions/brt_function_boot.R')

# brt, raw, continuous
         brt_function_boot(abundance = abundance, 
                     covariates = covariates, 
                     transformation = NA, # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# brt, log, continuous
         brt_function_boot(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# brt, log, discrete
         brt_function_boot(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# need to run
# brt, log10, continuous
         brt_function_boot(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = F,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')

# need to run
# brt, log10, discrete
         brt_function_boot(abundance = abundance, 
                     covariates = covariates, 
                     transformation = 'log10', # option is NA, log, log10
                     discrete = T,       # option is T or F 
                     species_name = as.character(names(abundance)[4]), 
                     n.cores = 10,
                     base_dir = 'results/rls_modelfits')











