# function to fit glms 

# buildmer for model fitting and stepwise model selection of glmmTMB
# remotes::install_github("cvoeten/buildmer"); https://github.com/cvoeten/buildmer

#load("data/rls_abun_modelling_data_v2.RData")
#abundance = rls_abun_fitting[[1]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
# species_name <- as.character(names(abundance)[11])

# write model
# model = 'gam'

# transformation
# transformation = NA

# family
# family = 'poisson'
# family = 'tweedie'

# function to fit glms
gam_function <- function(abundance = abundance, 
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         model = NA, # option is gam
                         family = NA, # option is NA, poisson, or nbinom, tweedie or zip
                         species_name = species_name, 
                         base_dir = 'results/modelfits'){
  
  require(tidyverse)
  require(mgcv)
  
  if(!is.na(transformation) & !is.na(family)){stop('dont be naughty, you shouldnt fit transformations and error structures')}
  
  # filter to focal species
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
  }
    if(transformation == 'log10'){
      abundance$abundance <- log10(abundance$abundance+1)
    }}
  
  # response variable name
  response <- 'abundance~'
  
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  
  # create formula with new covariate names
  covNames_new <- names(covariates)
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  covNames_splines <- paste0(rep('s(', length(covNames_new)),  covNames_new, rep(', k = 3)', length(covNames_new)))
  covNames_combined <- paste0(covNames_splines, collapse = '+')
  
  # find the full model formula
  model_formula <- as.formula(paste0(response, covNames_combined))
  
  # MODEL SELECTION NOTE:
  #	Model selection is performed through 'Null space penalization' as there is no approriate way to perform 
  # stepwise model selection for gams https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/step.gam.html
  # select = TRUE: If this is TRUE then gam can add an extra penalty to each term so that it can be penalized 
  # to zero. This means that the smoothing parameter estimation that is part of fitting can completely remove 
  # terms from the model. If the corresponding smoothing parameter is estimated as zero then the extra penalty 
  # has no effect. Use gamma to increase level of penalization.
  
  # gaussian with log and log10 transformations
  if(is.na(family)){
    model_fit <- gam(model_formula,data=abundance,family=gaussian,select=TRUE, method = 'ML')
  }else{
  
  # poisson
  if(family == 'poisson'){
    model_fit <- gam(model_formula,data=abundance,family=poisson,select=TRUE, method = 'ML')
  }
  
  # nbinom
  if(family == 'nbinom'){
    model_fit <- gam(model_formula,data=abundance,family=nb,select=TRUE, method = 'ML')
  }
  
  # tweedie
  if(family == 'tweedie'){
    model_fit <- gam(model_formula,data=abundance,family=tw,select=TRUE, method = 'ML')
  }
  
  # zip
  if(family == 'zip'){
    model_fit <- gam(model_formula,data=abundance,family=ziP,select=TRUE, method = 'ML')
  }
  # ziplss doesn't appear to work...
  }
  if(!family %in% c(NA, 'poisson', 'nbinom', 'tweedie', 'zip')){stop('family must be one of: NA, poisson, nbinom, tweedie, zip')}

  # save output into appropriate folder system
  # transformation/model/family/zi/
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/', 
                      model , '/', 
                      if(is.na(family)){'gaussian'}else{family}
  )
  model_path <- paste0(base_dir, '/', model_dir)
  dir.create(model_path, recursive = T)
  save(model_fit, file = paste0(model_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
}



