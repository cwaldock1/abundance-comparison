# Function for fitting boosted regression tree abundance models 

# load in abundance data
#load("data/rls_abun_modelling_data_v2.RData")
#abundance = rls_abun_fitting[[1]]

# load validation set
#validation = rls_abun_validation[[1]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

# discrete
#discrete <- T

# transformation
#transformation = 'log'

brt_function_boot <- function(abundance = abundance, 
                         validation = validation,
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         discrete = NA,       # option is T or F 
                         species_name = species_name, 
                         n.cores=10,
                         n_bootstrap = 10,
                         dataset = 'rls',
                         base_dir        = 'results/rls',
                         model_path      = 'model', 
                         prediction_path = 'predictions'){
  
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
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  validation <- validation[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(validation)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  validation <- left_join(validation, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
    validation$abundance <- log(validation$abundance+1)
  }
    if(transformation == 'log10'){
      abundance$abundance <- log10(abundance$abundance+1)
      validation$abundance <- log10(validation$abundance+1)
    }}
  
  # convert to discrete values based on groupings as in Howard, C., Stephens, P. A., Pearce-Higgins, J. W., Gregory, R. D. & Willis, S. G. Improving species distribution models: the value of data on abundance. Methods Ecol. Evol. 5, 506–513 (2014).
  if(discrete == T){
    
    # cannot put into a ordered factor as in brt has to just be classes with a multinomial distribution
    # abundance
    abundance$abundance <- round(abundance$abundance)     # round logged abundances
    abundance$abundance[abundance$abundance > 6] <- 6     # truncate abundances
    abundance$abundance <- as.factor(abundance$abundance) # turn into factors for random forests
    
    # validation
    validation$abundance <- round(validation$abundance)     # round logged abundances
    validation$abundance[validation$abundance > 6] <- 6     # truncate abundances
    validation$abundance <- as.factor(validation$abundance) # turn into factors for random forests
    
  }
  
  if(discrete == T & is.na(transformation)){stop('bad behaviour - dont use discrete with no transformtion - too many classes')}
  
  # rename and get general names for covariates for generalism
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(validation)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
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
  
  
  ### Bootstrap all model input data to make the presences only 2* the absences where possible
  
  # outputs for bootstrap 
  abundance_boot        <- list() # bootstrapped abundance object
  model_fit             <- list() # model fit objects
  verification_observed <- list() # bootstrapped verification observtions for prediction objects
  validation_observed   <- list() # bootstrapped validation for predictions objects
  verification_predict  <- list() # bootstrapped predictions from verificaiton data
  validation_predict    <- list() # bootstrapped predictions from validation data 
  
  # for loop to create the bootstraps and fit the models 
  for(boot in 1:n_bootstrap){
    
    print(boot)
    
    ### CREATING BOOTSTRATS
    set.seed(123)
    seed_123 <- sample(1:1000, n_bootstrap, replace = F)
    
    # get abundance data
    abundance_only <- abundance[which(abundance$abundance != 0),]
    
    # get sample size for absence bootstrap
    n_subsample <- nrow(abundance[which(abundance$abundance != 0),])*2
    
    # get absence 
    set.seed(seed_123[[boot]])
    boot_absence <- abundance[sample(which(abundance$abundance == 0), n_subsample, replace = F),]
    
    # combine absence and presence
    abundance_boot[[boot]] <- rbind(abundance_only, boot_absence) # this also acts as the verification of the model
    
    
    ### FITTING BOOSTED REGRESSION TREES
    # fit models with h2o package (benefits include automated stopping, poisson errors, faster fitting and parameter selection)
    # create feature names
    y <- "abundance"
    x <- covNames_new
    train.h2o <- as.h2o(abundance_boot[[boot]][,-c(1:2)])
 
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
  
    # fitting boosted regression trees
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
      seed = 123)
  
    # get the model grid
    grid_perf <- h2o.getGrid(
      grid_id = grid@grid_id, 
      sort_by = "mse", 
      decreasing = T
      )

    # obtain the best model
    best_model_id <- grid_perf@model_ids[[1]]
    model_fit[[boot]] <- h2o.getModel(best_model_id)
    rm(grid, grid_perf, best_model_id, hyper_grid, train.h2o)
    
    
    ### PREDICTIONS
    # create objects to predict into
    verification.h20 <- as.h2o(abundance_boot[[boot]][,-c(1:3)])
    validation.h20 <- as.h2o(validation[,-c(1:3)])

    # predict data using verification
    verification_predict[[boot]] <- as.numeric(as.character(as.data.frame(h2o.predict(model_fit[[boot]], verification.h20))[,1]))
    verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
    
    # predict data using validation
    validation_predict[[boot]]   <- as.numeric(as.character(as.data.frame(h2o.predict(model_fit[[boot]], validation.h20))[,1]))
    validation_observed[[boot]] <- as.numeric(as.character(validation$abundance))
    
    
    # convert all values back to raw abundances
    if(is.na(transformation)){NULL}else{if(transformation == 'log10'){
      verification_observed[[boot]] <- 10^(verification_observed[[boot]])-1
      validation_observed[[boot]]   <- 10^(validation_observed[[boot]])-1
      verification_predict[[boot]]  <- 10^(verification_predict[[boot]])-1
      validation_predict[[boot]]    <- 10^(validation_predict[[boot]])-1
    }
      if(transformation == 'log'){
        verification_observed[[boot]] <- exp(verification_observed[[boot]])-1
        validation_observed[[boot]]   <- exp(validation_observed[[boot]])-1
        verification_predict[[boot]]  <- exp(verification_predict[[boot]])-1
        validation_predict[[boot]]    <- exp(validation_predict[[boot]])-1
      }
    } # end of else statement
  } # end of bootstrapping loop

  # MAKE PREDICTIONS OBJECT
  # extract boostrapped predictions
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  fitted_model = 'gbm', 
                                  family  = distribution, 
                                  transformation = transformation, 
                                  zi = NA,
                                  n_abundance = nrow(abundance_only), 
                                  n_absence   = nrow(abundance[which(abundance$abundance == 0),]),
                                  n_boot_absence = nrow(boot_absence), 
                                  
                                  # estimate mean predictions
                                  verification_observed_mean = list(apply(simplify2array(verification_observed), 1, mean)),# list(verification_observed), 
                                  verification_predict_mean = list(apply(simplify2array(verification_predict), 1, mean)),#list(verification_predict), 
                                  validation_observed_mean = list(apply(simplify2array(validation_observed), 1, mean)),#list(validation_observed),
                                  validation_predict_mean = list(apply(simplify2array(validation_predict), 1, mean)),
                                  
                                  # estimate median predictions
                                  verification_observed_median = list(apply(simplify2array(verification_observed), 1, median)),# list(verification_observed), 
                                  verification_predict_median = list(apply(simplify2array(verification_predict), 1, median)),#list(verification_predict), 
                                  validation_observed_median = list(apply(simplify2array(validation_observed), 1, median)),#list(validation_observed),
                                  validation_predict_median = list(apply(simplify2array(validation_predict), 1, median)),
                                  
                                  # the amount of variation caused by bootstrapping to random 0s
                                  sd_verification = mean(apply(simplify2array(verification_predict), 1, sd)), 
                                  sd_validation   = mean(apply(simplify2array(validation_predict), 1, sd)))

  
  ### SAVE
  # create save directories model
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/',
                      'brt' , '/', 
                      if(discrete == T){'discrete/'}else{'continuous/'})
  model_final_path <- paste0(base_dir, '/', model_path, '/', model_dir, gsub(' ', '_', species_name), '/')
  dir.create(model_final_path, recursive = T)
  for(boot in 1:n_bootstrap){
  h2o.saveModel(object=model_fit[[boot]], path=paste0(model_final_path, gsub(' ', '_', species_name), '_', boot), force=TRUE)
  }
  
  # example load
  #h2o.loadModel('results/modelfits/log/brth20/continuous/Acanthurus_lineatus/Grid_GBM_data.frame_sid_9c50_33_model_R_1569423010894_1_model_6')
  h2o.shutdown(prompt = F)
  
  
  # save prediction object
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_',
                      'brt' , '_', 
                      if(discrete == T){'discrete'}else{'continuous'})
  prediction_final_path <- paste0(base_dir, '/', prediction_path)
  dir.create(prediction_final_path, recursive = T)
  save(extracted_predictions, file = paste0(prediction_final_path, '/', model_dir, '_', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
} # end of function
  

