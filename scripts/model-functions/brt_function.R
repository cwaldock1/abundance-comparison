# Function for fitting boosted regression tree abundance models 

# load in abundance data
#load("data/rls_abun_modelling_data.RData")
#abundance = rls_abun_fitting

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteCode', 'SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]
#covariates[,4] <- as.numeric(scale(log(abs(covariates[,4])), center = T))

# get species names
#species_name <- as.character(names(abundance)[8])

# discrete
#discrete <- T

# transformation
#transformation = 'log'

brt_function <- function(abundance = abundance, 
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         discrete = NA,       # option is T or F 
                         species_name = species_name, 
                         n.cores = 10,
                         base_dir = 'results/modelfits'){
  
  require(tidyverse)
  #require(gbm)
  #install.packages('h2o')
  require('h2o')
  #h2o.init(port = sample(c(1:10000), 1))
  # this makes sure that a port is opened
  while(is.logical(tryCatch(h2o.clusterInfo(), error = function(e)NA))){
    tryCatch(h2o.init(port = sample(c(1:10000), 1)), error = function(e) NA)
             }
  
  # filter to focal species
  abundance <- abundance[c('SiteCode', 'SiteLongitude', 'SiteLatitude', species_name)]
  names(abundance)[4] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
  }
    if(transformation == 'log10'){
      abundance$abundance <- log10(abundance$abundance+1)
    }}
  
  # convert to discrete values based on groupings as in Howard, C., Stephens, P. A., Pearce-Higgins, J. W., Gregory, R. D. & Willis, S. G. Improving species distribution models: the value of data on abundance. Methods Ecol. Evol. 5, 506–513 (2014).
  if(discrete == T){
    abundance$abundance <- round(abundance$abundance)     # round logged abundances
    abundance$abundance[abundance$abundance > 6] <- 6     # truncate abundances
    abundance$abundance <- as.factor(abundance$abundance) # turn into factors for random forests
  }
  
  # rename and get general names for covariates for generalism
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[4+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[3+i] <- paste0('cov', i)}
  covNames_new <- names(covariates) # brt forumla made with this object 
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  # create formula
  #response <- 'abundance ~'
  #covs     <- paste(covNames_new, collapse = '+')
  #brt_formula <- as.formula(paste(response, covs))
  
  # fit models with gbm packages
  # fit boosted regression tree with stochastic gradient boosting (i.e., bag.fraction != 1)
  #model_fit <- gbm(formula = brt_formula,
  #                 data = abundance, 
  #                 distribution = 'gaussian', 
  #                 n.trees = 10000,
  #                 interaction.depth = 3, 
  #                 shrinkage = 0.001,
  #                 bag.fraction = 0.75, 
  #                 cv.folds = 10, 
  #                 n.cores = n.cores)
  
  #gbm.mod.perf <- gbm.perf(model_fit, method = "cv", plot.it = F) 
  
  #model_fit <- gbm(formula = brt_formula,
  #                 data = abundance, 
  #                 distribution = 'gaussian', 
  #                 n.trees = gbm.mod.perf,
  #                 interaction.depth = 3, 
  #                 shrinkage = 0.001)
  
  # summary(model_fit, method = relative.influence, plotit = T)
  
  
  # fit models with h2o package (benefits include automated stopping, poisson errors, faster fitting and parameter selection)
  # create feature names
  y <- "abundance"
  x <- covNames_new
  train.h2o <- as.h2o(abundance[,-c(1:3)])
 
  # set up a hypergrid
  hyper_grid <- list(
    max_depth = c(5),
    min_rows = c(5, 10),
    learn_rate = c(0.001, 0.005),
    sample_rate = c(.5, .75))
  
  # fit model to parameter grid
  if(discrete == T & !is.na(transformation)){distribution <- 'multinomial'}
  if(discrete == F &  is.na(transformation)){distribution <- 'poisson'}
  if(discrete == F & !is.na(transformation)){distribution <- 'AUTO'}
  
  grid <- h2o.grid(
    algorithm = "gbm",
   # grid_id = "gbm_grid_new",
    x = x, 
    y = y, 
    training_frame = train.h2o,
    distribution = distribution,
    hyper_params = hyper_grid,
    learn_rate_annealing = 0.99,
    ntrees = 10000,
    nfolds = 10, 
    stopping_rounds = 10,
    stopping_tolerance = 0.01,
    seed = 123
  )
  
  # get the model grid
  grid_perf <- h2o.getGrid(
    grid_id = grid@grid_id, 
    sort_by = "r2", 
    decreasing = T
  )

  # obtain the best model
  best_model_id <- grid_perf@model_ids[[1]]
  best_model <- h2o.getModel(best_model_id)
  
  # create save directories model
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/',
                      'brt' , '/', 
                      if(discrete == T){'discrete/'}else{'continuous/'})
  model_path <- paste0(base_dir, '/', model_dir)
  dir.create(model_path, recursive = T)
  
  # old save
  #save(model_fit, file = paste0(model_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
  # new save
  h2o.saveModel(object=best_model, path=paste0(model_path, gsub(' ', '_', species_name)), force=TRUE)
  
  # example load
  #h2o.loadModel('results/modelfits/log/brth20/continuous/Acanthurus_lineatus/Grid_GBM_data.frame_sid_9c50_33_model_R_1569423010894_1_model_6')
  h2o.shutdown(prompt = F)
  
  }
  