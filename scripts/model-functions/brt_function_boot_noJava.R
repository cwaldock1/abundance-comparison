# Function for fitting boosted regression tree abundance models 

brt_function_boot <- function(abundance = abundance, 
                         validation = validation,
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         discrete = NA,       # option is T or F 
                         family = NA, # must be multinomial for discrete transformations
                         species_name = NA, 
                         n.cores=1,
                         n_bootstrap = 10,
                         dataset = 'rls',
                         base_dir        = 'results/rls',
                         model_path      = 'model', 
                         prediction_path = 'predictions'){
  
  require(tidyverse)
  require(gbm)
  #install.packages('h2o')
  #require('h2o')
  #h2o.init(port = sample(c(1:10000), 1))
  # this makes sure that a port is opened
  #while(is.logical(tryCatch(h2o.clusterInfo(), error = function(e)NA))){
  #  tryCatch(h2o.init(port = sample(c(1:10000), 1)), error = function(e) NA)
  #           }
  
  tryCatch({
  
  # filter to focal species
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  validation <- validation[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(validation)[3] <- 'abundance'
  
  # estimate number of absences for if statements laters
  n_absences <- sum(abundance$abundance %in% 0)
  
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
    # abundance
    abundance$abundance <- if(n_absences == 0){ceiling(abundance$abundance)}else{round(abundance$abundance)}     # round logged abundances
    abundance$abundance[abundance$abundance > 6] <- 6     # truncate abundances
    abundance$abundance <- ordered(as.factor(abundance$abundance)) # turn into factors for random forests
    
    # validation
    validation$abundance <- if(n_absences == 0){ceiling(validation$abundance)}else{round(validation$abundance)}     # round logged abundances
    validation$abundance[validation$abundance > 6] <- 6     # truncate abundances
    validation$abundance <- ordered(as.factor(validation$abundance)) # turn into factors for random forests
  }
  
  if(discrete == T & is.na(transformation)){return(print('bad behaviour - dont use discrete with no transformtion - too many classes'))}
  if(length(unique(abundance$abundance)) == 1 | length(unique(validation$abundance)) == 1){return(print('bad behaviour - too few classes for discrete transformtion'))}
  
  # rename and get general names for covariates for generalism
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(validation)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  covNames_new <- names(covariates) # brt forumla made with this object 
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  # create formula
  response <- 'abundance ~'
  covs     <- paste(covNames_new, collapse = '+')
  brt_formula <- as.formula(paste(response, covs))
  
  ### Bootstrap all model input data to make the presences only 2* the absences where possible
  # outputs for bootstrap 
  abundance_boot        <- list() # bootstrapped abundance object
  model_fit             <- list() # model fit objects
  verification_observed <- list() # bootstrapped verification observtions for prediction objects
  validation_observed   <- list() # bootstrapped validation for predictions objects
  verification_predict  <- list() # bootstrapped predictions from verificaiton data
  validation_predict    <- list() # bootstrapped predictions from validation data 
  
  # for loop to create the bootstraps and fit the models
  n_bootstrap <- ifelse(n_absences > 0, n_bootstrap, 1)
  for(boot in 1:n_bootstrap){
    
    print(boot)
    
    ### CREATING BOOTSTRATS
    set.seed(123)
    seed_123 <- sample(1:1000, n_bootstrap, replace = F)
    
    # get abundance data
    abundance_only <- abundance[which(as.numeric(as.character(abundance$abundance)) > 0),]
    
    # get sample size for absence bootstrap
    n_subsample <- nrow(abundance[which(as.numeric(as.character(abundance$abundance)) > 0),])*2
    
    # get absence 
    set.seed(seed_123[[boot]])
    if(sum(abundance$abundance %in% 0) > 0){boot_absence <- abundance[sample(which(as.numeric(as.character(abundance$abundance)) == 0), n_subsample, replace = F),]}else{boot_absence <- NULL}
    
    # combine absence and presence
    abundance_boot[[boot]] <- rbind(abundance_only, boot_absence) # this also acts as the verification of the model
    
    ### FITTING BOOSTED REGRESSION TREES

    # boosted regression tree modelling in 'gbm'
    if(length(covNames_new) != 1){
    # fit models with gbm packages
    # fit boosted regression tree with stochastic gradient boosting (i.e., bag.fraction != 1)
    model_fit[[boot]] <- tryCatch(gbm(formula = brt_formula,
                                      data = abundance_boot[[boot]], 
                                      distribution = family, 
                                      n.trees = 10000,
                                      interaction.depth = 3, 
                                      shrinkage = 0.001,
                                      bag.fraction = 0.8, 
                                      cv.folds = 10, 
                                      n.cores = n.cores), error = function(e) NA)
     
    if(!is.na(model_fit[[boot]])){
      
    # selecting the best number of trees from cross validations
    gbm.mod.perf <- gbm.perf(model_fit[[boot]], method = "cv", plot.it = F) 
    
    # fit model to all data
    model_fit[[boot]] <- gbm(formula = brt_formula,
                             data = abundance_boot[[boot]], 
                             distribution = family, 
                             n.trees = gbm.mod.perf,
                             bag.fraction = 0.8,
                             interaction.depth = 3, 
                             shrinkage = 0.001)
    
    }else{ # if the optimal number of trees cannot be identified in gbm
      
      abundance_boot[[boot]]$abundance <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
      
      # find the best model using gbm
      times = 0
      gbm.mod.perf <- NULL
      while(is.null(gbm.mod.perf)){
      gbm.mod.perf <- tryCatch(dismo::gbm.step(data = data.frame(abundance_boot[[boot]]),
                                      gbm.x = 4:ncol(data.frame(abundance_boot[[boot]])),
                                      gbm.y = 3,
                                      family = 'poisson', # cannot have multinomial in this package, but also cannot fit models with single covariate and cross validation in gbm
                                      tree.complexity = 3,
                                      learning.rate = 0.001-(0.0001*times),
                                      n.folds = 10, 
                                      bag.fraction = 0.8, 
                                      plot.main = F)$n.trees, error = function(e) NULL)
      times <- times + 1
      if(times == 500){stop()}
      }
      
      abundance_boot[[boot]]$abundance <- as.ordered(abundance_boot[[boot]]$abundance)
      
      model_fit[[boot]] <- tryCatch(gbm(formula = brt_formula,
                                        data = abundance_boot[[boot]], 
                                        distribution = family, 
                                        n.trees = gbm.mod.perf,
                                        bag.fraction = 0.8,
                                        interaction.depth = 3, 
                                        shrinkage = 0.001), error = identity)
      
      if(class(model_fit[[boot]])[1] == 'simpleError'){
        verification_observed[[boot]] <- NA
        validation_observed[[boot]]   <- NA
        verification_predict[[boot]]  <- NA
        validation_predict[[boot]]    <- NA
        next
        }
    
    }
    
    }else{
      
      # the conditional statements occur here because
      # if the dataset is a 2nd stage model with one covariate then the multinomial throwns an error in the gbm function
      # note that we only reach this stage in the function if the above is true, so we can code conditional on all multinomial errors being converted to poisson for practicality
      if(family == 'multinomial'){family_conditional <- 'poisson'}else{family_conditional<-family}
      abundance_boot[[boot]]$abundance <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
      
      gbm.mod.perf <- dismo::gbm.step(data = data.frame(abundance_boot[[boot]]),
                      gbm.x = 4, 
                      gbm.y = 3,
                      family = family_conditional, # cannot have multinomial in this package, but also cannot fit models with single covariate and cross validation in gbm
                      tree.complexity = 3,
                      learning.rate = 0.001,
                      n.folds = 10, 
                      bag.fraction = 0.8, 
                      plot.main = F)$n.trees
      
      #abundance_boot[[boot]]$abundance <- if(family_conditional == 'multinomial'){as.ordered(abundance_boot[[boot]]$abundance)}else{as.numeric(abundance_boot[[boot]]$abundance)}
      
      model_fit[[boot]] <- gbm(formula = brt_formula,
                               data = tibble(abundance = abundance_boot[[boot]]$abundance,
                                             cov1      = abundance_boot[[boot]]$cov1), 
                               distribution = family_conditional, 
                               n.trees = gbm.mod.perf,
                               bag.fraction = 0.8,
                               interaction.depth = 3, 
                               shrinkage = 0.001)
      
    }
    
    ### PREDICTIONS
    
    # make discrete predictions
    if(discrete == T){
      
      predictions_boot_verification <- predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), n.trees = gbm.mod.perf, type = 'response')
      predictions_boot_validation   <- predict(model_fit[[boot]], validation, n.trees = gbm.mod.perf, type = 'response')
      if(length(covNames_new) != 1){ # this conditional is in place because the 2-stage abundance models do not use multinomial distributions
        predictions_boot_verification <- apply(predictions_boot_verification, 1,
                                             function(x){
                                               weights = x
                                               values = colnames(predictions_boot_verification)
                                               weighted.mean(as.numeric(values), weights)
                                             })
       predictions_boot_validation   <- apply(predictions_boot_validation, 1,
                                             function(x){
                                               weights = x
                                               values = colnames(predictions_boot_validation)
                                               weighted.mean(as.numeric(values), weights)
                                             })
      }
      
      verification_predict[[boot]]  <- as.numeric(as.character(predictions_boot_verification))
      verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
      validation_predict[[boot]]    <- as.numeric(as.character(predictions_boot_validation))
      validation_observed[[boot]]   <- as.numeric(as.character(validation$abundance))
      
    }else{

    # predict data using verification
    predictions_boot_verification <- predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), n.trees = gbm.mod.perf, type = 'response')
    predictions_boot_validation   <- predict(model_fit[[boot]], validation, n.trees = gbm.mod.perf, type = 'response')
    
    }
    
    # apply to bootstrapped estimates
    verification_predict[[boot]]  <- as.numeric(as.character(predictions_boot_verification))
    verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
    
    # predict data using validation
    validation_predict[[boot]]  <- as.numeric(as.character(predictions_boot_validation))
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
    }
  } # end of bootstrapping loop
  
  if(sum(is.na(verification_observed)) == 10){
    error_name <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_',
                        'brt' , '_', 
                        if(discrete == T){'discrete_'}else{'continuous_'}, 
                        family, 
                        '.txt')
    dir.create(base_dir, recursive = T)
    # save error message so can quickly diagnose
    fileConn<-file(paste0(base_dir, '/', error_name))
    writeLines('lack of data error', fileConn)
    close(fileConn)
    return(NA)}
  
  verification_observed <- verification_observed[!is.na(verification_observed)]
  validation_observed <- validation_observed[!is.na(validation_observed)]
  verification_predict <- verification_predict[!is.na(verification_predict)]
  validation_predict <- validation_predict[!is.na(validation_predict)]
  
  # MAKE PREDICTIONS OBJECT
  # extract boostrapped predictions
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  
                                  # model properties
                                  fitted_model = 'gbm', 
                                  only_abundance = if(sum(abundance$abundance %in% 0) > 0){F}else{T},
                                  family  = paste(ifelse(discrete == T, 'discrete', 'continuous'), family, sep = '_'), 
                                  transformation = transformation, 
                                  zi = NA,
                                  n_abundance = nrow(abundance_only), 
                                  n_absence   = if(sum(abundance$abundance %in% 0) > 0){nrow(abundance[which(abundance$abundance == 0),])}else{0},
                                  n_boot_absence = if(sum(abundance$abundance %in% 0) > 0){nrow(boot_absence)}else{0}, 
                                  
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
                                  sd_validation   = mean(apply(simplify2array(validation_predict), 1, sd)), 
                                  
                                  verification_locations = list(data.frame(SiteLongitude = abundance$SiteLongitude, SiteLatitude = abundance$SiteLatitude)), 
                                  validation_locations = list(data.frame(SiteLongitude = validation$SiteLongitude, SiteLatitude = validation$SiteLatitude)))

  
  ### SAVE
  # create save directories model
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/',
                      'brt' , '/', 
                      if(discrete == T){'discrete/'}else{'continuous/'}, 
                      family, '/')
  model_final_path <- paste0(base_dir, '/', model_path, '/', model_dir)
  dir.create(model_final_path, recursive = T)
  save(model_fit, file = paste0(model_final_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
  # create prediction object to save
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_',
                      'brt' , '_', 
                      if(discrete == T){'discrete'}else{'continuous'}, 
                      '_', family)
  prediction_final_path <- paste0(base_dir, '/', prediction_path)
  dir.create(prediction_final_path, recursive = T)
  save(extracted_predictions, file = paste0(prediction_final_path, '/', model_dir, '_', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
  }, 
  error = function(e) NA)
  
} # end of function
  

