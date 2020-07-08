
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

occ_rf <- function(abundance = abundance, 
                   covariates = covariates, 
                   species_name = species_name, 
                   spatial_projections, 
                   base_dir,
                   spatial_dir,
                   n_bootstrap,
                   n.cores = 10){
  
  require(tidyverse)
  require(randomForest)
  require(ecospat)
  
  ecoVersion <- packageVersion("ecospat")
  
  # filter to columns of interest for abundance and occurrence datasets
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  occurrence <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  names(occurrence)[3] <- 'occurrence'
  # modify occurrences to range between 0 and 1
  occurrence$occurrence[which(occurrence$occurrence >=1)] <- 1

  # remove absences and log transform abundance data
  abundance <- abundance[which(abundance$abundance != 0),]
  abundance$abundance <- log(abundance$abundance+1)
  
  # join together abundance, occurrence and covariates dataframes 
  abundance   <- left_join(abundance, covariates)
  occurrence  <- left_join(occurrence, covariates)
  
  # CHECK THAT I DON'T NEED THIS SECTION ANY LONGER
  # # response variable name
  # response <- 'occurrence~'
  # 
  # # rename covariates
  # covNames_org <- names(covariates)
  # covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  # for(i in 1:length(covNames_org)){names(occurrence)[3+i] <- paste0('cov', i)}
  # for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  # 
  # # get covariate names
  # covNames_new <- names(covariates)
  # covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  ### Bootstrap all model input data to make the presences only 2* the absences where possible
  # outputs for bootstrap 
  occurrence_boot <- list() # bootstrapped occurrence object
  rf_occ  <- list()
  imp_occ <- list()
  eval_occ <- list()
  occ_predictions <- list()
  occ_spatial_predictions <- list()
  errors <- c()
  
  # get occurrence data to use in each iteration
  presences_only <- occurrence[which(occurrence$occurrence > 0),]

  # number of times more absences than presences, then split data into bootstrapped sets
  absence_occ      <- occurrence[which(occurrence$occurrence == 0),]
  #nrow_times       <- ceiling(nrow(occurrence[which(occurrence$occurrence == 0),]) / nrow(presences_only))
  #absence_splits   <- split(1:nrow(occurrence[which(occurrence$occurrence == 0),]), 1:nrow_times)
  
  
  
  # RUN OCCURRENCE MODELS USING BOOTSTRAPS
  for(boot in 1:n_bootstrap){
    
    # set up seed
    set.seed(123)
    seed_123 <- sample(1:1000, n_bootstrap, replace = F)
    
    # get sample size for absence bootstrap
    n_subsample <- nrow(presences_only)*2
    
    # get absence 
    set.seed(seed_123[[boot]])
    if(sum(occurrence$occurrence %in% 0) > 0){boot_absence <- occurrence[sample(which(as.numeric(as.character(occurrence$occurrence)) == 0), n_subsample, replace = F),]}else{boot_absence <- NULL}
    
    # combine absence and presence
    occurrence_boot[[boot]] <- rbind(presences_only, boot_absence) # this also acts as the verification of the model
  
    # fit models
    rf_occ[[boot]] <- randomForest(x = occurrence_boot[[boot]][-c(1:3)], 
                                   y = as.factor(occurrence_boot[[boot]]$occurrence), 
                                   ntree = 1000, 
                                   importance=T)

    # get variable importance score for the occurrences
    imp_occ[[boot]] <- data.frame(importance(rf_occ[[boot]]))
    
    # estimate cut-off based on TSS
    eval_occ[[boot]] <- data.frame(predictions  = predict(rf_occ[[boot]], occurrence_boot[[boot]][-c(1:3)], type = 'prob')[,2], 
                                   observations = occurrence_boot[[boot]]['occurrence'])
    
    # predict on covariate values at sites with abundances
    occ_predictions[[boot]] <- predict(rf_occ[[boot]], 
                                       occurrence[which(occurrence$occurrence==1),names(covariates)], 
                                       type = 'prob')[,2]
    
    # make spatial predictions and average 
    occ_spatial_predictions[[boot]] <- predict(rf_occ[[boot]], 
                                               spatial_projections[,names(covariates)],
                                               type = 'prob')[,2]
  
  }
  
  # average variable importance across occurrence bootstraps
  occ_importance <- apply(do.call(cbind, lapply(imp_occ, function(x) x[,3])), 1, mean)
  
  # estimate best threshold for converting probability of presences to - presence-absence to constrain spatial projections DO THIS LATER
  obs_preds <- do.call(rbind, eval_occ)
  occ_eval <- ecospat.max.tss(obs_preds[,1], obs_preds[,2])
  
  # average over occupancy spatial predictions
  occ_spatial_predictions <- apply(do.call(cbind, occ_spatial_predictions), 1, mean)
  
  
  ### fit abundance models
  # model fitting
  rf_abun <- randomForest(x = abundance[-c(1:3)], 
                          y = abundance$abundance, 
                          ntree = 1000, 
                          importance=T)
  
  # get abundance variable importance
  imp_rf <- as.numeric(importance(rf_abun)[,1])
  
  # make predictions
  abun_predictions <- predict(rf_abun, occurrence[which(occurrence$occurrence==1),names(covariates)])
  
  # make spatial predictions 
  abun_spatial_predictions <- predict(rf_abun, spatial_projections[,names(covariates)])
  
  
  # extract occurrence predictions from occurrences
  occ_predictions <- apply(do.call(cbind, occ_predictions), 1, mean)
  
  ### CREATE OUTPUTS FOR SAVING 
  # create variable importance dataframe for a given species
  var_imp_comparison <- data.frame(species_name = species_name, 
             covariate = names(covariates)[-c(1,2)], 
             occurrence_imp = occ_importance, 
             abundance_imp = imp_rf, 
             occurrence_imp_rescaled = (occ_importance-min(occ_importance,na.rm = T))/(max(occ_importance,na.rm = T)-min(occ_importance,na.rm=T)),
             abundance_imp_rescaled  = (imp_rf-min(imp_rf,na.rm = T))/(max(imp_rf,na.rm = T)-min(imp_rf,na.rm=T)))
  
  var_imp_comparison$spearmans_rank = cor(var_imp_comparison$occurrence_imp_rescaled, 
                                          var_imp_comparison$abundance_imp_rescaled, 
                                          method = 'spearman')
  
  # create predictions output frame for correlation between occurrence and abundance
  predictions <- data.frame(species_name = species_name, 
                            occ_predictions  = occ_predictions, 
                            abun_predictions = abun_predictions, 
                            occ_predictions_scaled  = (occ_predictions-min(occ_predictions,na.rm = T))/(max(occ_predictions,na.rm = T)-min(occ_predictions,na.rm=T)), 
                            abun_predictions_scaled = (abun_predictions-min(abun_predictions,na.rm = T))/(max(abun_predictions,na.rm = T)-min(abun_predictions,na.rm=T)))

  predictions$spearman = cor(predictions$occ_predictions, 
                             predictions$abun_predictions, 
                             method = 'spearman')
  
  # save model outputs
  model_final_path <- paste0(base_dir, '/', dataset)
  dir.create(model_final_path, recursive = T)
  saveRDS(list(var_imp_comparison, predictions), 
          file = paste0(model_final_path, '/', gsub(' ', '_', species_name), '.RDS'))
  
  
  
  ### create object for spatial projections
  spatial_abundance_prediction <- data.frame(
    occupancy_rate = abs(round(occ_spatial_predictions*1000)), 
    abundance_log  = abs(round(abun_spatial_predictions*1000)))

  # save spatial projections
  spatial_path <- paste0(spatial_dir, '/', dataset)
  dir.create(spatial_path, recursive = T)
  saveRDS(list(spatial_abundance_prediction, TSS_threshold = ifelse(ecoVersion == 3.0, occ_eval[[2]][2,2], occ_eval[[3]])), 
            file = paste0(spatial_path, '/', gsub(' ', '_', species_name), '.RDS')
            )
  
  # save tss scores seperately
  #tss_path <- paste0(spatial_dir, '/TSS/', dataset)
  #dir.create(tss_path, recursive = T)
  #saveRDS(occ_eval, 
  #        file = paste0(tss_path, '/', gsub(' ', '_', species_name), '.RDS')
  #)
  
  
  } #end of function
