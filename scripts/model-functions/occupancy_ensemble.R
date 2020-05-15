
# here we want to use both the validation and the verification data because we are 
# simply producing a new dataset that has the same length as the entier set of surveys -> 
# otherwise we can't have covariates for the validation data
#load("data/rls_abun_modelling_data_v2.RData")
#abundance = do.call(rbind, rls_abun_list[[2]][[1]])
#abundance  = rls_abun_fitting[[1]]
#validation = rls_abun_validation[[1]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

occupancy_ensemble <- function(abundance = abundance, 
                               validation = validation,
                               covariates = covariates, 
                               species_name = species_name, 
                               n.cores = 10){
  require(tidyverse)
  require(glmmTMB)
  require(randomForest)
  require(mgcv)
  require(gbm)
  require(buildmer)
  
  # filter to columns of interest
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  validation <- validation[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(validation)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance  <- left_join(abundance, covariates)
  validation <- left_join(validation, covariates)
  
  # create occurrence frames
  occurrence <- abundance
  names(occurrence)[3] <- 'occurrence'
  occurrence$occurrence[which(occurrence$occurrence >=1)] <- 1
  val_occurrence <- validation
  names(val_occurrence)[3] <- 'occurrence'
  val_occurrence$occurrence[which(val_occurrence$occurrence >=1)] <- 1
  
  # check for constant values
  if(length(unique(occurrence$occurrence)) == 1 | length(unique(val_occurrence$occurrence)) == 1){return(print('bad behaviour - too few classes for discrete transformtion'))}
  
  # response variable name
  response <- 'occurrence~'
  
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(occurrence)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(val_occurrence)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  
  # get covariate names
  covNames_new <- names(covariates)
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  # create formula with new covariate names: GLM
  covNames_new_2 <- paste0('I(', covNames_new, '^2',')')
  covNames_glm <- paste0(c(covNames_new, covNames_new_2), collapse = '+')
  model_formula_glm <- as.formula(paste0(response, covNames_glm))
  model_tab <- tabulate.formula(model_formula_glm)
  model_tab[-1, 'block'] <- rep(covNames_new, 2)
  
  # create formula with new covariate names: GAM
  covNames_splines <- paste0(rep('s(', length(covNames_new)),  covNames_new, rep(', k = 3)', length(covNames_new)))
  covNames_gam <- paste0(covNames_splines, collapse = '+')
  model_formula_gam <- as.formula(paste0(response, covNames_gam))
  
  # create formula with new covariate names: BRT
  brt_covs     <- paste(covNames_new, collapse = '+')
  brt_formula <- as.formula(paste(response, brt_covs))
  
  
  ### Bootstrap all model input data to make the presences only 2* the absences where possible
  # outputs for bootstrap 
  occurrence_boot        <- list() # bootstrapped abundance object
  validation_ensemble    <- list()
  #abundance_boot        <- list() # bootstrapped abundance object
  #model_fit             <- list() # model fit objects
  verification_observed <- list() # bootstrapped verification observations for prediction objects
  verification_predict  <- list() # bootstrapped predictions from verificaiton data
  validation_observed <- list() # bootstrapped validation observations for prediction objects
  validation_predict  <- list() # bootstrapped predictions from verificaiton data
  glm_fit <- list() # empty model object
  gam_fit <- list() # empty model object
  rf_fit  <- list() # empty model object
  brt_fit <- list() # empty model object
  glm_predict_verification <- list() # empty predict object
  gam_predict_verification <- list() # empty predict object
  rf_predict_verification  <- list() # empty predict object
  brt_predict_verification <- list() # empty predict object
  glm_predict_validation <- list() # empty predict object
  gam_predict_validation <- list() # empty predict object
  rf_predict_validation  <- list() # empty predict object
  brt_predict_validation <- list() # empty predict object
  errors <- c()
  
  # for loop to create the bootstraps and fit the models 
  #n_bootstrap <- ifelse(sum(abundance$abundance %in% 0) > 0, n_bootstrap, 1)
  
  # get occurrence data to use in each iteration
  occurrence_only <- occurrence[which(occurrence$occurrence > 0),]
  abundance_only  <- abundance[which(abundance$abundance > 0),]
  
  # number of times more absences than presences, then split data into bootstrapped sets
  absence_occ  <- occurrence[which(occurrence$occurrence == 0),]
  #absence_abun <- abundance[which(abundance$abundance == 0),]
  nrow_times       <- ceiling(nrow(occurrence[which(occurrence$occurrence == 0),]) / nrow(abundance_only))
  absence_splits   <- split(1:nrow(occurrence[which(occurrence$occurrence == 0),]), 1:nrow_times)
  
  for(boot in 1:length(absence_splits)){
    
  ### SET UP BOOTSTRAP DATA

    # get absence 
    boot_absence_occ  <- absence_occ[absence_splits[[boot]],]
    #boot_absence_abun <- absence_abun[absence_splits[[boot]],]
    
    # combine absence and presence
    occurrence_boot[[boot]] <- rbind(occurrence_only, boot_absence_occ) # this also acts as the verification of the model
    #abundance_boot[[boot]]  <- rbind(abundance_only, boot_absence_abun) # this also acts as the verification of the model
    
  
  ### FIT GLM
    # fit model with backwards stepwise regression
    glm_fit[[boot]] <- buildglmmTMB(formula = model_tab,
                                      data = occurrence_boot[[boot]], 
                                      dep = 'occurrence', 
                                      family = 'binomial')
    
    # remove quadratic terms that are non-significant and refit model
    coef_table       <- coef(summary(glm_fit[[boot]]))$cond               # get coefficient table
    coef_names       <- rownames(coef(summary(glm_fit[[boot]]))$cond)[-1]
    coef_names_ns    <- coef_names[which(coef_table[-1,4] > 0.05)]  # get coefficients to remove
    
    # if there are non-significant I(cov^2) terms
    if(sum(grepl('I(', coef_names_ns, fixed = T)) > 0){
      model_terms <- coef_names[-which(coef_names %in% coef_names_ns[grep('I(', coef_names_ns, fixed = T)])] # get non-significant quadraties from table
      model_tab_2 <- model_tab[which(model_tab$term %in% model_terms),]    # subset model_tab by all significant terms
      model_tab_glm <- rbind(model_tab[1,], model_tab_2) # put back intercept
      
      # refit model
      glm_fit[[boot]] <- buildglmmTMB(formula = model_tab_glm,
                                        data = occurrence_boot[[boot]], 
                                        dep = 'occurrence',
                                        family = 'binomial')
      
      }else{print('all terms are significant')}

    # glm predictions
    glm_predict_verification[[boot]] <- predict(glm_fit[[boot]]@model, data.frame(occurrence_boot[[boot]]), type = 'response')
    glm_predict_validation[[boot]]   <- predict(glm_fit[[boot]]@model, val_occurrence, type = 'response')
    
    
  ### FIT GAM
    
    # test indiviual terms fail
    fail_cov <- !is.na(lapply(1:(length(occurrence_boot[[boot]])-3), 
                              function(x) 
                                tryCatch(
                      gam(occurrence ~ s(cov, k = 3),
                          data = data.frame(occurrence = occurrence_boot[[boot]]$occurrence,
                                           cov = data.frame(occurrence_boot[[boot]])[,x+3]),
                          family = binomial, 
                          method = 'ML'),
                                  error = function(e) NA)))
    
    covNames_combined_modified <- paste0(covNames_splines[fail_cov], collapse = '+')
    model_formula_modified <- as.formula(paste0(response, covNames_combined_modified))
    
    gam_fit[[boot]] <- gam(model_formula_modified, data=occurrence_boot[[boot]], family='binomial', select=TRUE, method = 'ML')
   
    gam_predict_verification[[boot]] <- as.numeric(predict(gam_fit[[boot]], data.frame(occurrence_boot[[boot]]), type = 'response'))
    gam_predict_validation[[boot]]   <- as.numeric(predict(gam_fit[[boot]], val_occurrence, type = 'response'))
    
    
  ### FIT RF
    rf_fit[[boot]] <- randomForest(x = occurrence_boot[[boot]][covNames_new], 
                                      y = occurrence_boot[[boot]]$occurrence, 
                                      ntree = 1000, 
                                      importance=FALSE)
    rf_predict_verification[[boot]] <- as.numeric(predict(rf_fit[[boot]], data.frame(occurrence_boot[[boot]]), type = 'response'))
    rf_predict_validation[[boot]]   <- as.numeric(predict(rf_fit[[boot]], val_occurrence, type = 'response'))
    
  ### FIT BRT
    # fit boosted regression tree with stochastic gradient boosting (i.e., bag.fraction != 1)
    brt_fit[[boot]] <- gbm(formula = brt_formula,
                             data = occurrence_boot[[boot]], 
                             distribution = 'bernoulli', 
                             n.trees = 10000,
                             interaction.depth = 3, 
                             shrinkage = 0.001,
                             bag.fraction = 0.75, 
                             cv.folds = 10, 
                             n.cores = n.cores)
    
    # selecting the best number of trees from cross validations
    gbm.mod.perf <- gbm.perf(brt_fit[[boot]], method = "cv", plot.it = F) 
    
    # fit model to all data
    brt_fit[[boot]] <- gbm(formula = brt_formula,
                             data = occurrence_boot[[boot]], 
                             distribution = 'bernoulli', 
                             n.trees = gbm.mod.perf,
                             interaction.depth = 3, 
                             shrinkage = 0.001)
    
    brt_predict_verification[[boot]] <- as.numeric(as.character(as.data.frame(predict(brt_fit[[boot]], 
                                                                                      data.frame(occurrence_boot[[boot]]), 
                                                                                      n.trees = gbm.mod.perf, 
                                                                                      type = 'response'))[,1]))
    
    brt_predict_validation[[boot]] <- as.numeric(as.character(as.data.frame(predict(brt_fit[[boot]], 
                                                                                      val_occurrence, 
                                                                                      n.trees = gbm.mod.perf, 
                                                                                      type = 'response'))[,1]))
    
    }
  
  # aggregate suitability estimates for each model type
  suitability_ensemble <- list()
    for(i in 1:length(absence_splits)){
      occurrence_boot[[i]]$suitability_ensemble <- rowMeans(simplify2array(list(glm_predict_verification[[i]],
                                                                                gam_predict_verification[[i]],
                                                                                rf_predict_verification[[i]],
                                                                                brt_predict_verification[[i]])))
      validation_ensemble[[i]]                  <- rowMeans(simplify2array(list(glm_predict_validation[[i]],
                                                                                gam_predict_validation[[i]],
                                                                                rf_predict_validation[[i]],
                                                                                brt_predict_validation[[i]])))
      
    }
  
  
  
  # combine into single suitability ensemble object
  suitability_ensemble <- do.call(rbind, occurrence_boot)[,c('SiteLongitude', 'SiteLatitude', 'suitability_ensemble')]
  
  # take average over locations when multiple locations are used in a loop
  ver_suitability_ensemble <- suitability_ensemble %>% 
    group_by(SiteLongitude, SiteLatitude) %>% 
    do(suitability_ensemble = mean(.$suitability_ensemble, na.rm=T)) %>% 
    unnest(suitability_ensemble) %>% 
    left_join(occurrence %>% dplyr::select(SiteLatitude, SiteLongitude), .)
  
  # create average over boots for validation dataset
  val_suitability_ensemble <- as_tibble(data.frame(validation[c('SiteLongitude', 'SiteLatitude')], 
                                                   suitability_ensemble = rowMeans(simplify2array(validation_ensemble))))
  
  # output the two suitability estimates
  # this is not circular because any sites that are contined in the validation dataset were not used to build the
  # initial models 
  suitability_ensemble <- rbind(ver_suitability_ensemble, val_suitability_ensemble) %>% 
    group_by(SiteLatitude, SiteLongitude) %>% 
    do(suitability_ensemble = mean(.$suitability_ensemble, na.rm=T)) %>% 
    unnest(c('suitability_ensemble'))
  
  # end of bootstrap
  return(suitability_ensemble)
      
  } #end of function
