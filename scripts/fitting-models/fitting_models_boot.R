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
                      "human_pop_2015_50km", "reef_area_200km", "wave_energy_mean", "Depth_GEBCO_transformed",
                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5')]

# setup number of bootstraps 
n_boots = 10



# RUN MODELS FOR ABUNDANCE AND OCCUPANCY DATA ----
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
#n_bootstrap = n_boots,
#dataset = 'rls',
#base_dir        = 'results/rls',
#model_path      = 'model', 
#prediction_path = 'predictions_abunocc'

# lm, log, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# lm, log10, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# glm, poisson, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# glm, nbinom, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# glm, tweedie, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# glm, poisson, with zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = T,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# glm, nbinom, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = covariates,
                    transformation = NA, # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = n_boots,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model_abunocc', 
                    prediction_path = 'predictions_abunocc'),
  error = function(e) NA)


# glm, tweedie, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = covariates,
                    transformation = NA, # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = n_boots,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model_abunocc', 
                    prediction_path = 'predictions_abunocc'), 
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
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# gam, log10
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# gam, raw, poisson
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# gam raw nbinom
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# gam raw tweedie
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# gam raw zip
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# run random forest options ----
source('scripts/model-functions/rf_function_boot.R')

# abundance = abundance, 
# validation = validation,
# covariates = covariates, 
# transformation = NA, # option is NA, log, log10
# discrete = NA,       # option is T or F 
# species_name = species_name, 
# n_bootstrap = n_boots,
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
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc', 
                 prediction_path = 'predictions_abunocc')

# rf, log, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc', 
                 prediction_path = 'predictions_abunocc')

# rf, log, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc', 
                 prediction_path = 'predictions_abunocc')

# rf, log10, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc', 
                 prediction_path = 'predictions_abunocc')

# rf, log10, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc', 
                 prediction_path = 'predictions_abunocc')



# RUN MODELS FOR ABUNDANCE ONLY ----
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
#n_bootstrap = n_boots,
#dataset = 'rls',
#base_dir        = 'results/rls',
#model_path      = 'model', 
#prediction_path = 'predictions'

# lm, log, no zero inflation
glm_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')


# lm, log10, no zero inflation
glm_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')


# glm, poisson, no zero inflation
glm_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')


# glm, nbinom, no zero inflation
glm_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')


# glm, tweedie, no zero inflation
glm_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# run gam options ----
source('scripts/model-functions/gam_function_boot.R')

# gam, log
gam_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# gam, log10
gam_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# gam, raw, poisson
gam_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# gam raw nbinom
gam_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# gam raw tweedie
gam_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')


# run random forest options ----
source('scripts/model-functions/rf_function_boot.R')

# abundance = abundance, 
# validation = validation,
# covariates = covariates, 
# transformation = NA, # option is NA, log, log10
# discrete = NA,       # option is T or F 
# species_name = species_name, 
# n_bootstrap = n_boots,
# dataset = 'rls',
# base_dir        = 'results/rls',
# model_path      = 'model', 
# prediction_path = 'predictions'

# rf, raw, continuous
rf_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                 validation = validation_input,
                 covariates = covariates,
                 transformation = NA, # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abun', 
                 prediction_path = 'predictions_abun')

# rf, log, continuous
rf_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abun', 
                 prediction_path = 'predictions_abun')

# rf, log, discrete
rf_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abun', 
                 prediction_path = 'predictions_abun')

# rf, log10, continuous
rf_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abun', 
                 prediction_path = 'predictions_abun')

# rf, log10, discrete
rf_function_boot(abundance = abundance_input[which(abundance_input$Num != 0),],
                 validation = validation_input,
                 covariates = covariates,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abun', 
                 prediction_path = 'predictions_abun')



# RUN TWO STAGE MODELS ----
# run occurrence models - two stage ---- 
source('scripts/model-functions/occupancy_ensemble.R')

suitability <- occupancy_ensemble(abundance = abundance_input,
                                  validation = validation_input,
                                  covariates = covariates,
                                  species_name = unique(abundance_input$TAXONOMIC_NAME),
                                  n.cores=1) # remove the number of cores function when on server...

dir.create('results/rls/suitability', recursive = T)
save(suitability, file = paste0('results/rls/suitability', '/', unique(abundance_input$TAXONOMIC_NAME), '.RData'))

## run glm options - two stage ----
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
#n_bootstrap = n_boots,
#dataset = 'rls',
#base_dir        = 'results/rls',
#model_path      = 'model', 
#prediction_path = 'predictions_abunocc'

# lm, log, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# lm, log10, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'lm', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# glm, poisson, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# glm, nbinom, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# glm, tweedie, no zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                  zi = F,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# glm, poisson, with zero inflation
glm_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'glm', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie
                  zi = T,     # option is F or T
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls', 
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# glm, nbinom, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = suitability,
                    transformation = NA, # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'nbinom', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = n_boots,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model_abunocc_2stage', 
                    prediction_path = 'predictions_abunocc_2stage'),
  error = function(e) NA)


# glm, tweedie, with zero inflation
tryCatch(
  glm_function_boot(abundance = abundance_input,
                    validation = validation_input,
                    covariates = suitability,
                    transformation = NA, # option is NA, log, log10
                    model = 'glm', # option is lm or glm
                    family = 'tweedie', # option is NA, poisson, or nbinom, tweedie
                    zi = T,     # option is F or T
                    species_name = unique(abundance_input$TAXONOMIC_NAME), 
                    n_bootstrap = n_boots,
                    dataset = 'rls',
                    base_dir        = 'results/rls', 
                    model_path      = 'model_abunocc_2stage', 
                    prediction_path = 'predictions_abunocc_2stage'), 
  error = function(e) NA)


# run gam options - two stage ----
source('scripts/model-functions/gam_function_boot.R')

# gam, log
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# gam, log10
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log10', # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = NA, # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# gam, raw, poisson
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'poisson', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# gam raw nbinom
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'nbinom', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# gam raw tweedie
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'tweedie', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# gam raw zip
gam_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  model = 'gam', # option is lm or glm
                  family = 'zip', # option is NA, poisson, or nbinom, tweedie, zip
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')


# run random forest options - two stage ----
source('scripts/model-functions/rf_function_boot.R')

# abundance = abundance, 
# validation = validation,
# covariates = covariates, 
# transformation = NA, # option is NA, log, log10
# discrete = NA,       # option is T or F 
# species_name = species_name, 
# n_bootstrap = n_boots,
# dataset = 'rls',
# base_dir        = 'results/rls',
# model_path      = 'model', 
# prediction_path = 'predictions'

# rf, raw, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = suitability,
                 transformation = NA, # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc_2stage', 
                 prediction_path = 'predictions_abunocc_2stage')

# rf, log, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = suitability,
                 transformation = 'log', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc_2stage', 
                 prediction_path = 'predictions_abunocc_2stage')

# rf, log, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = suitability,
                 transformation = 'log', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc_2stage', 
                 prediction_path = 'predictions_abunocc_2stage')

# rf, log10, continuous
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = suitability,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = F,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc_2stage', 
                 prediction_path = 'predictions_abunocc_2stage')

# rf, log10, discrete
rf_function_boot(abundance = abundance_input,
                 validation = validation_input,
                 covariates = suitability,
                 transformation = 'log10', # option is NA, log, log10
                 discrete = T,       # option is T or F 
                 species_name = unique(abundance_input$TAXONOMIC_NAME), 
                 n_bootstrap = n_boots,
                 dataset = 'rls',
                 base_dir        = 'results/rls',
                 model_path      = 'model_abunocc_2stage', 
                 prediction_path = 'predictions_abunocc_2stage')



# RUNNING BOOSTED REGRESSION TREES AT END  ---- 
# run boosted regression trees options - ABUNDANCE_OCCURRENCE ----

# source function
source('scripts/model-functions/brt_function_boot_noJava.R')

# brt, raw, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'poisson',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# brt, log, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# brt, log, discrete
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# brt, log10, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')

# brt, log10, discrete
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc', 
                  prediction_path = 'predictions_abunocc')


# run boosted regression trees options - ABUNDANCE_ONLY ----

# source function
source('scripts/model-functions/brt_function_boot_noJava.R')

# brt, raw, continuous
brt_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = NA, # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'poisson',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# brt, log, continuous
brt_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# brt, log, discrete
brt_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# brt, log10, continuous
brt_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# brt, log10, discrete
brt_function_boot(abundance = abundance_input[which(abundance_input$Num!=0),],
                  validation = validation_input,
                  covariates = covariates,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abun', 
                  prediction_path = 'predictions_abun')

# run boosted regression trees options - two stage ----

# source function
source('scripts/model-functions/brt_function_boot_noJava.R')

# brt, raw, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = NA, # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'poisson',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# brt, log, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# brt, log, discrete
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial', # cannot have mulitinomial for the single variate models fitted with gbm
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# brt, log10, continuous
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = F,       # option is T or F 
                  family = 'gaussian',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

# brt, log10, discrete
brt_function_boot(abundance = abundance_input,
                  validation = validation_input,
                  covariates = suitability,
                  transformation = 'log10', # option is NA, log, log10
                  discrete = T,       # option is T or F 
                  family = 'multinomial',
                  species_name = unique(abundance_input$TAXONOMIC_NAME), 
                  n_bootstrap = n_boots,
                  dataset = 'rls',
                  base_dir        = 'results/rls',
                  model_path      = 'model_abunocc_2stage', 
                  prediction_path = 'predictions_abunocc_2stage')

stop()


