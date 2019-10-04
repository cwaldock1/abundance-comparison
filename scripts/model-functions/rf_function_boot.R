
# Function for fitting random forest abundance models 

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

# discrete
#discrete <- T

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

# transformation
#transformation = 'log'

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
    abundance$abundance <- round(abundance$abundance)     # round logged abundances
    abundance$abundance[abundance$abundance > 6] <- 6     # truncate abundances
    abundance$abundance <- ordered(as.factor(abundance$abundance)) # turn into factors for random forests
    
    # validation
    validation$abundance <- round(validation$abundance)     # round logged abundances
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
  for(boot in 1:n_bootstrap){
    print(boot)
    ### CREATING BOOTSTRATS
    set.seed(123)
    seed_123 <- sample(1:1000, n_bootstrap, replace = F)
    
    # get abundance data
    abundance_only <- abundance[which(abundance$abundance > 0),]
    
    # get sample size for absence bootstrap
    n_subsample <- nrow(abundance[which(abundance$abundance > 0),])*2
    
    # get absence 
    set.seed(seed_123[[boot]])
    boot_absence <- abundance[sample(which(abundance$abundance == 0), n_subsample, replace = F),]
    
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
     verification_predict[[boot]]  <- as.numeric(as.character(predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'response')))
     verification_observed[[boot]] <- as.numeric(as.character(abundance_boot[[boot]]$abundance))
     
     # predict data using validation
     validation_predict[[boot]]  <- as.numeric(as.character(predict(model_fit[[boot]], data.frame(validation), type = 'response')))
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
     
  } # end of bootstrap forloop
  
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  fitted_model = 'rf', 
                                  family  = ifelse(discrete == T, 'discrete', 'continuous'), 
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
  
