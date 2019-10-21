# function to fit glms 

# buildmer for model fitting and stepwise model selection of glmmTMB
# remotes::install_github("cvoeten/buildmer"); https://github.com/cvoeten/buildmer

#load("data/rls_abun_modelling_data_v2.RData")
#abundance = rls_abun_fitting[[42]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

# write model
#model = 'gam'

# transformation
#transformation = NA

# family
#family = 'poisson'
# family = 'tweedie'

# load verification and validation data
# verification     = rls_abun_fitting[[1]][,c(2,3,5)] # this isn't needed as the verification is simply the data put into the model
#validation         = rls_abun_validation[[1]]


# function to fit glms
gam_function_boot <- function(abundance = abundance, 
                              validation = validation,
                              covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         model = NA, # option is gam
                         family = NA, # option is NA, poisson, or nbinom, tweedie or zip
                         species_name = species_name, 
                         n_bootstrap = 10,
                         dataset = 'rls',
                         base_dir        = 'results/rls',
                         model_path      = 'model', 
                         prediction_path = 'predictions'){
  
  require(tidyverse)
  require(mgcv)
  
  if(!is.na(transformation) & !is.na(family)){stop('dont be naughty, you shouldnt fit transformations and error structures')}
  
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
  
  # response variable name
  response <- 'abundance~'
  
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(validation)[3+i] <- paste0('cov', i)}
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
  
  
  ### Bootstrap all model input data to make the presences only 2* the absences where possible
  
  # outputs for bootstrap 
  abundance_boot        <- list() # bootstrapped abundance object
  model_fit             <- list() # model fit objects
  verification_observed <- list() # bootstrapped verification observtions for prediction objects
  validation_observed   <- list() # bootstrapped validation for predictions objects
  verification_predict  <- list() # bootstrapped predictions from verificaiton data
  validation_predict    <- list() # bootstrapped predictions from validation data 
  
  # for loop to create the bootstraps and fit the models 
  n_bootstrap <- ifelse(sum(abundance$abundance %in% 0) > 0, n_bootstrap, 1)
  for(boot in 1:n_bootstrap){
    print(boot)
    # CREATE BOOTSTRAPED DATA
    set.seed(123)
    seed_123 <- sample(1:1000, n_bootstrap, replace = F)
    
    # get abundance data
    abundance_only <- abundance[which(abundance$abundance > 0),]
    
    # get sample size for absence bootstrap
    n_subsample <- nrow(abundance[which(abundance$abundance > 0),])*2
    
    # get absence 
    set.seed(seed_123[[boot]])
    if(sum(abundance$abundance %in% 0) > 0){boot_absence <- abundance[sample(which(abundance$abundance == 0), n_subsample, replace = F),]}else{boot_absence <- NULL}
    
    # combine absence and presence
    abundance_boot[[boot]] <- rbind(abundance_only, boot_absence) # this also acts as the verification of the model
    
    # FIT THE MODELS
    #plot(gam(abundance ~ s(cov3, k = 3) + s(cov1, k = 3) + s(cov10, k = 3), data = abundance_boot[[boot]], family = gaussian, method = 'ML'))
    
    # test indiviual terms fail
    fail_cov <- !is.na(lapply(1:(length(abundance_boot[[boot]])-3), 
                              function(x) 
                                tryCatch(
                                  
                                  if(is.na(family)){
                                    gam(abundance ~ s(cov, k = 3),
                                             data = data.frame(abundance = abundance_boot[[boot]]$abundance,
                                                               cov = data.frame(abundance_boot[[boot]])[,x+3]), family = gaussian, method = 'ML')}else{
                                  if(family == 'poisson'){
                                    gam(abundance ~ s(cov, k = 3),
                                      data = data.frame(abundance = abundance_boot[[boot]]$abundance,
                                                        cov = data.frame(abundance_boot[[boot]])[,x+3]), family = poisson, method = 'ML')
                                    }
                                  if(family == 'nbinom'){
                                  gam(abundance ~ s(cov, k = 3),
                                      data = data.frame(abundance = abundance_boot[[boot]]$abundance,
                                                        cov = data.frame(abundance_boot[[boot]])[,x+3]), family = nb, method = 'ML')
                                    }
                                  if(family == 'tweedie'){
                                  gam(abundance ~ s(cov, k = 3),
                                      data = data.frame(abundance = abundance_boot[[boot]]$abundance,
                                                        cov = data.frame(abundance_boot[[boot]])[,x+3]), family = tw, method = 'ML')
                                    }
                                  if(family == 'tweedie'){
                                  gam(abundance ~ s(cov, k = 3),
                                      data = data.frame(abundance = abundance_boot[[boot]]$abundance,
                                                        cov = data.frame(abundance_boot[[boot]])[,x+3]), family = ziP, method = 'ML')
                                    }},error = function(e) NA)))
    covNames_combined_modified <- paste0(covNames_splines[fail_cov], collapse = '+')
    model_formula_modified <- as.formula(paste0(response, covNames_combined_modified))
    
    # gaussian with log and log10 transformations
    if(is.na(family)){
      
        model_fit[[boot]] <- gam(model_formula_modified, data=abundance_boot[[boot]], family=gaussian, select=TRUE, method = 'ML')
     
        }else{

        # poisson
        if(family == 'poisson'){
            model_fit[[boot]] <- gam(model_formula_modified,data=abundance_boot[[boot]],family=poisson,select=TRUE, method = 'ML')
          }
  
        # nbinom
        if(family == 'nbinom'){
            model_fit[[boot]] <- gam(model_formula_modified,data=abundance_boot[[boot]],family=nb,select=TRUE, method = 'ML')
          }
  
        # tweedie
        if(family == 'tweedie'){
            model_fit[[boot]] <- gam(model_formula_modified,data=abundance_boot[[boot]],family=tw,select=TRUE, method = 'ML')
          }
  
        # zip
        if(family == 'zip'){
             model_fit[[boot]] <- gam(model_formula_modified,data=abundance_boot[[boot]],family=ziP,select=TRUE, method = 'ML')
        }
        }# end of else statement

    if(!family %in% c(NA, 'poisson', 'nbinom', 'tweedie', 'zip')){stop('family must be one of: NA, poisson, nbinom, tweedie, zip')}

    
    # MAKE PREDICTIONS
    # predict data using verification
    verification_predict[[boot]]  <- as.numeric(predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'response'))
    verification_observed[[boot]] <- abundance_boot[[boot]]$abundance
    
    # predict data using validation
    validation_predict[[boot]]  <- as.numeric(predict(model_fit[[boot]], data.frame(validation), type = 'response'))
    validation_observed[[boot]] <- validation$abundance
    
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
    }# end of else statement
    
    
  } # end of bootstrapped for loop
    
  
  # create prediction object
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  
                                  fitted_model = 'gam', 
                                  only_abundance = if(sum(abundance$abundance %in% 0) > 0){F}else{T},
                                  family  = family, 
                                  transformation = transformation, 
                                  zi = ifelse(family == 'zip', T, F),
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
                      model , '/', 
                      if(is.na(family)){'gaussian'}else{family}
  )
  model_final_path <- paste0(base_dir, '/', model_path, '/', model_dir)
  dir.create(model_final_path, recursive = T)
  save(model_fit, file = paste0(model_final_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
  # save prediction output in the same file structure
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_', 
                      model , '_', 
                      if(is.na(family)){'gaussian'}else{family}
  )
  prediction_final_path <- paste0(base_dir, '/', prediction_path)
  dir.create(prediction_final_path, recursive = T)
  save(extracted_predictions, 
       file = paste0(prediction_final_path, '/', model_dir, '_', gsub(' ', '_', species_name), '.RData'), 
       recursive = T)
  
  
}



