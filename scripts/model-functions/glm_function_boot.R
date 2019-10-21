# function to fit glms 

# buildmer for model fitting and stepwise model selection of glmmTMB
# remotes::install_github("cvoeten/buildmer"); https://github.com/cvoeten/buildmer

# load in abundance data
#load("data/rls_abun_modelling_data_v2.RData")
#abundance = rls_abun_fitting[[1]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

# model
#model <- 'glm'

# transformation
#transformation = NA

# family
#family = 'poisson'
#family = 'tweedie'

# zi
#zi = F

# load verification and validation data
# verification     = rls_abun_fitting[[1]][,c(2,3,5)] # this isn't needed as the verification is simply the data put into the model
#validation         = rls_abun_validation[[1]]


# function to fit glms
glm_function_boot <- function(abundance = abundance, 
                              validation = validation,
                              covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         model = NA, # option is lm or glm
                         family = NA, # option is NA, poisson, or nbinom, or tweedie
                         zi = F,     # option is F or T
                         species_name = species_name, 
                         n_bootstrap = 10,
                         dataset = 'rls',
                         base_dir        = 'results/rls',
                         model_path      = 'model', 
                         prediction_path = 'predictions'){
  require(tidyverse)
  require(glmmTMB)
  require(buildmer)
  
  if(model == 'lm' & is.na(transformation)){stop('attempting to fit a lm to raw abundance data')}
  
  # filter to columns of interest
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  validation <- validation[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(validation)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance  <- left_join(abundance, covariates)
  validation <- left_join(validation, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance  <- log(abundance$abundance+1)
    validation$abundance <- log(validation$abundance+1)
    
  }
  if(transformation == 'log10'){
    abundance$abundance  <- log10(abundance$abundance+1)
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
  covNames_new_2 <- paste0('I(', covNames_new, '^2',')')
  covNames_combined <- paste0(c(covNames_new, covNames_new_2), collapse = '+')
  model_formula <- as.formula(paste0(response, covNames_combined))
  
  # build formula structure using buildmer
  model_tab <- tabulate.formula(model_formula)
  model_tab[-1, 'block'] <- rep(covNames_new, 2)
  
  
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
    
    
  # fit the models in the bootstrap
    
  if(model == 'lm'){
  
  # fit model with backwards stepwise regression
  model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                            data = abundance_boot[[boot]], 
                            dep = 'abundance') 
  
  # remove quadratic terms that are non-significant and refit model
  coef_table       <- coef(summary(model_fit[[boot]]))$cond               # get coefficient table
  coef_names       <- rownames(coef(summary(model_fit[[boot]]))$cond)[-1]
  coef_names_ns    <- coef_names[which(coef_table[-1,4] > 0.05)]  # get coefficients to remove
  
  # if there are non-significant I(cov^2) terms
  if(sum(grepl('I(', coef_names_ns, fixed = T)) > 0){
    model_terms <- coef_names[-which(coef_names %in% coef_names_ns[grep('I(', coef_names_ns, fixed = T)])] # get non-significant quadraties from table
    model_tab_2 <- model_tab[which(model_tab$term %in% model_terms),]    # subset model_tab by all significant terms
    model_tab <- rbind(model_tab[1,], model_tab_2) # put back intercept
    
  # refit model
  model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                            data = abundance_boot[[boot]], 
                            dep = 'abundance')

  }else{print('all terms are significant')}
  
  }
  
  if(model == 'glm'){
    
    
    # fit model with backwards stepwise regression
    if(!family %in% c('poisson', 'nbinom', 'tweedie')){stop('family for glm must be one of poisson or nbinom')}
    
    if(family == 'poisson'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                              data = abundance_boot[[boot]], 
                              dep = 'abundance', 
                              family = 'poisson')}
    
    if(family == 'nbinom'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                                      data = abundance_boot[[boot]], 
                                                      dep = 'abundance', 
                                                      family = nbinom1(link = 'log'))}
    
    if(family == 'tweedie'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                                     data = abundance_boot[[boot]], 
                                                     dep = 'abundance', 
                                                     family = tweedie(link = 'log'))}
    
    # remove quadratic terms that are non-significant and refit model
    coef_table       <- coef(summary(model_fit[[boot]]))$cond               # get coefficient table
    coef_names       <- rownames(coef(summary(model_fit[[boot]]))$cond)[-1]
    coef_names_ns    <- coef_names[which(coef_table[-1,4] > 0.05)]  # get coefficients to remove
    
    # if there are non-significant I(cov^2) terms
    if(sum(grepl('I(', coef_names_ns, fixed = T)) > 0){
      
    model_terms <- coef_names[-which(coef_names %in% coef_names_ns[grep('I(', coef_names_ns, fixed = T)])] # get non-significant quadraties from table
    model_tab_2 <- model_tab[which(model_tab$term %in% model_terms),]    # subset model_tab by all significant terms
    model_tab <- rbind(model_tab[1,], model_tab_2) # put back intercept
      
      # refit model
      if(family == 'poisson'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                                        data = abundance_boot[[boot]], 
                                                        dep = 'abundance', 
                                                        family = 'poisson')}
      
      if(family == 'nbinom'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                                       data = abundance_boot[[boot]], 
                                                       dep = 'abundance', 
                                                       family = nbinom1(link = 'log'))}
      
      if(family == 'tweedie'){model_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                                        data = abundance_boot[[boot]], 
                                                        dep = 'abundance', 
                                                        family = tweedie(link = 'log'))}
      
    }else{print('All terms are significant')}
    
  }
  
  
  
  # for zero-inflation take the selected model with no zero-inflation and remodel covariate set
  if(zi == T){
    
    # get the zi model structure
    model_structure <- formula(model_fit[[boot]]@model)
    zi_formula <- formula(paste('~', paste(attr(terms(model_structure), 'term.labels'), collapse = '+')))
    
    if(family == 'poisson'){model_fit[[boot]] <- glmmTMB(formula = model_structure,
                                                 zi = zi_formula, 
                                                 data = abundance_boot[[boot]], 
                                                 family = 'poisson')}
    
    if(family == 'nbinom'){model_fit[[boot]] <- glmmTMB(formula = model_structure,
                                                zi = zi_formula,
                                                data = abundance_boot[[boot]],                   
                                                family = nbinom1(link = 'log'))}
    
    if(family == 'tweedie'){model_fit[[boot]] <- glmmTMB(formula = model_structure,
                                                zi = zi_formula,
                                                data = abundance_boot[[boot]],                   
                                                family = tweedie(link = 'log'))}
    
    
  }
    
    
    
    # make validations
    # predict data using verification
    if(zi != T){verification_predict[[boot]]  <- predict(model_fit[[boot]]@model, data.frame(abundance_boot[[boot]]), type = 'response')}
    if(zi == T){verification_predict[[boot]]  <- predict(model_fit[[boot]], data.frame(abundance_boot[[boot]]), type = 'response')}
    verification_observed[[boot]] <- abundance_boot[[boot]]$abundance
      
    # predict data using validation
    if(zi != T){validation_predict[[boot]]  <- predict(model_fit[[boot]]@model, data.frame(validation), type = 'response')}
    if(zi == T){validation_predict[[boot]]  <- predict(model_fit[[boot]], data.frame(validation), type = 'response')}
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
      }

    
  } # end of bootstrapped for loop
  
  extracted_predictions <- tibble(dataset = dataset, 
                                  species_name = species_name, 
                                  
                                  fitted_model = 'glm', 
                                  only_abundance = if(sum(abundance$abundance %in% 0) > 0){F}else{T},
                                  family  = family, 
                                  transformation = transformation, 
                                  zi = zi,
                                  n_abundance    = nrow(abundance_only), 
                                  n_absence      = if(sum(abundance$abundance %in% 0) > 0){nrow(abundance[which(abundance$abundance == 0),])}else{0},
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
  
  # save model output into appropriate folder system
  # transformation/model/family/zi/
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/', 
                      model , '/', 
                      if(is.na(family)){'gaussian'}else{family}, '/', 
                      if(zi==F){'no_ZI'}else{'ZI'}
                      )
  model_final_path <- paste0(base_dir, '/', model_path, '/', model_dir)
  dir.create(model_final_path, recursive = T)
  save(model_fit, file = paste0(model_final_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)

  # save prediciton output in same file structure
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '_', 
                      model , '_', 
                      if(is.na(family)){'gaussian'}else{family}, '_', 
                      if(zi==F){'no_ZI'}else{'ZI'}
  )
  prediction_final_path <- paste0(base_dir, '/', prediction_path)
  dir.create(prediction_final_path, recursive = T)
  save(extracted_predictions, file = paste0(prediction_final_path, '/', model_dir, '_', gsub(' ', '_', species_name), '.RData'), recursive = T)

  }



