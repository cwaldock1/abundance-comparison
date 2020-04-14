
# Function for fitting random forest abundance models 

rf_function_boot <- function(abundance = abundance, 
                             validation = validation,
                             covariates = covariates, 
                             transformation = NA, # option is NA, log, log10
                             discrete = NA,       # option is T or F 
                             species_name = species_name, 
                             n_bootstrap = 10,
                             dataset = 'rls',
                             base_dir        = 'results/rls',
                             model_path      = 'model', 
                             prediction_path = 'predictions'){
  
  require(tidyverse)
  require(randomForest)
  
  # filter to focal species
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  validation <- validation[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(validation)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  validation <- left_join(validation, covariates)
  
  # estimate number of absences for if statements laters
  n_absences <- sum(abundance$abundance %in% 0)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
    validation$abundance <- log(validation$abundance+1)
  }
    if(transformation == 'log10'){
      abundance$abundance <- log10(abundance$abundance+1)
      validation$abundance <- log10(validation$abundance+1)
    }}
  
  if(discrete == T & is.na(transformation)){stop('bad behaviour - dont use discrete with no transformtion - too many classes')}
  
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
  
  # rename and get general names for covariates for generalism
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(validation)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  covNames_new <- names(covariates) # randomForests take matrix which can be subset with this object 
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  
  
  
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
    ### CREATING BOOTSTRAPS
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
  

    ### FITTING MODELS 
    # fit the random forests
     model_fit[[boot]] <- randomForest(x = abundance_boot[[boot]][covNames_new], 
                                       y = abundance_boot[[boot]]$abundance, 
                                       ntree = 1000, 
                                       importance=FALSE)
     
     
     ### PREDICTIONS
     # make validations
     # predict data using verification
     if(discrete == T){ 
       
       # predictions for discrete data
       # predictions for abundance and occurrences 
         predictions_boot_verification <- predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'prob')
         predictions_boot_verification <- apply(predictions_boot_verification, 1,
                                                function(x){
                                                  weights = x
                                                  values = colnames(predictions_boot_verification)
                                                  weighted.mean(as.numeric(values), weights)
                                                  }
                                                )
         
         predictions_boot_validation   <- predict(model_fit[[boot]], validation, type = 'prob')
         predictions_boot_validation   <- apply(predictions_boot_validation, 1,
                                                function(x){
                                                  weights = x
                                                  values = colnames(predictions_boot_validation)
                                                  weighted.mean(as.numeric(values), weights)
                                                }
                                                )
         
         # apply to loop
         verification_predict[[boot]]  <- as.numeric(as.character(predictions_boot_verification))
         verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
         validation_predict[[boot]]    <- as.numeric(as.character(predictions_boot_validation))
         validation_observed[[boot]]   <- as.numeric(as.character(validation$abundance))
         
       # apply hurdle approach to probabilities
       # make verification predictions
       #predictions_boot_verification <- predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'prob')
       #if(length(predictions_boot_verification[1,]) <= 2){
       #   predictions_boot_verification <- apply(data.frame(predictions_boot_verification), 1, which.max)
       #}else{
       #predictions_boot_verification <- ifelse(predictions_boot_verification[,1] < 0.5, apply(data.frame(predictions_boot_verification[,-1]), 1, which.max), 1)
       #}
       
       # make validation predictions
       #predictions_boot_validation   <- predict(model_fit[[boot]], validation, type = 'prob')
       #if(length(predictions_boot_validation[1,]) <= 2){
       #   predictions_boot_validation <- apply(data.frame(predictions_boot_validation), 1, which.max)
       #}else{
       #   predictions_boot_validation <- ifelse(predictions_boot_validation[,1] < 0.5, apply(data.frame(predictions_boot_validation[,-1]), 1, which.max), 1)
       #}
       
       # Apply -1 because it is based on column number 
       # predict data using verification
       #if(n_absences != 0){verification_predict[[boot]] <- as.numeric(as.character(predictions_boot_verification))-1}else{
       #   verification_predict[[boot]] <- as.numeric(as.character(predictions_boot_verification))}
       #verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
       
       # predict data using validation
       #if(n_absences != 0){validation_predict[[boot]] <- as.numeric(as.character(predictions_boot_validation))-1}else{
       #   validation_predict[[boot]] <- as.numeric(as.character(predictions_boot_validation))}
       #validation_observed[[boot]] <- as.numeric(as.character(validation$abundance))
       
       
       }else{
     
         
     # predictions for continous data
     verification_predict[[boot]]  <- as.numeric(as.character(predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'response')))
     verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
     
     # predict data using validation
     validation_predict[[boot]]  <- as.numeric(as.character(predict(model_fit[[boot]], data.frame(validation), type = 'response')))
     validation_observed[[boot]] <- as.numeric(as.character(validation$abundance))
       
     }
     
     # convert all values back to raw abundances
     if(is.na(transformation)){NULL}else{if(transformation == 'log10'){
       verification_observed[[boot]] <- 10^(verification_observed[[boot]])-1
       verification_predict[[boot]]  <- 10^(verification_predict[[boot]])-1
       validation_observed[[boot]]   <- 10^(validation_observed[[boot]])-1
       validation_predict[[boot]]    <- 10^(validation_predict[[boot]])-1
     }
       if(transformation == 'log'){
         verification_observed[[boot]] <- exp(verification_observed[[boot]])-1
         verification_predict[[boot]]  <- exp(verification_predict[[boot]])-1
         validation_observed[[boot]]   <- exp(validation_observed[[boot]])-1
         validation_predict[[boot]]    <- exp(validation_predict[[boot]])-1
       }
     }
     
  } # end of bootstrap forloop
  
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  
                                  fitted_model = 'rf', 
                                  only_abundance = if(sum(abundance$abundance %in% 0) > 0){F}else{T},
                                  family  = ifelse(discrete == T, 'discrete', 'continuous'), 
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
  
  
  # save output into appropriate folder system
  # transformation/model/family/zi/
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/',
                      'rf' , '/', 
                      if(discrete == T){'discrete/'}else{'continuous/'})
  model_final_path <- paste0(base_dir, '/', model_path, '/', model_dir)
  dir.create(model_final_path, recursive = T)
  save(model_fit, file = paste0(model_final_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
  # create prediction object to save
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_',
                      'rf' , '_', 
                      if(discrete == T){'discrete'}else{'continuous'})
  prediction_final_path <- paste0(base_dir, '/', prediction_path)
  dir.create(prediction_final_path, recursive = T)
  save(extracted_predictions, file = paste0(prediction_final_path, '/', model_dir, '_', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
}
  
