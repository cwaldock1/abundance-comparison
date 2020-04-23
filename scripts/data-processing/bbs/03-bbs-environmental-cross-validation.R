# setting up the cross validations based on temperature covariates

# Load packages ----
lib_vect <- c('tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# Load data to set up new cross validations ----

# load all objects
load('data/bbs_abun_modelling_data.RData')

# load covariates
load('data/bbs_covariates.RData')

# Join temperature data and set up new cross-validations ----

# join in cross validation covariate
bbs_abun_cv <- left_join(bbs_abun %>% select(-set), bbs_xy %>% dplyr::select(SiteLongitude, SiteLatitude, robPCA_1))

# lapply subsets 
bbs_abun_cv_2 <- lapply(1:length(unique(bbs_abun_cv$TAXONOMIC_NAME)), function(x){
  
  species_data <- bbs_abun_cv %>% 
    filter(TAXONOMIC_NAME == unique(bbs_abun_cv$TAXONOMIC_NAME)[x])
  
  all_absences  <- species_data %>% filter(Num == 0)  # absences
  all_presences <- species_data %>% filter(Num != 0)  # presences
  
  sample_size <- round(round(nrow(species_data)*0.8) / 2) # get sample size for cross validations
  
  # estimate 80th quantile in absences and presences
  absences_cv_quantile <- as.numeric(quantile(all_absences$robPCA_1, 0.8, na.rm = T))
  presences_cv_quantile <- as.numeric(quantile(all_presences$robPCA_1, 0.8, na.rm = T))
  
  absences_sample_fitting    <- all_absences[ -which(all_absences$robPCA_1  > absences_cv_quantile),]
  presence_sample_fitting    <- all_presences[-which(all_presences$robPCA_1 > presences_cv_quantile),]
  absences_sample_validation <- all_absences[ which(all_absences$robPCA_1  > absences_cv_quantile),]
  presence_sample_validation <- all_presences[which(all_presences$robPCA_1 > presences_cv_quantile),]
  
  fitting    <- rbind(absences_sample_fitting, presence_sample_fitting)
  validation <- rbind(absences_sample_validation, presence_sample_validation)
  
  fitting$set    <- 'fitting'
  validation$set <- 'validation'
  
  list_outputs <- list(list(fitting = fitting, validation = validation))
  names(list_outputs) <- unique(species_data$TAXONOMIC_NAME)
  
  return(list_outputs)
  
  }) # end of lapply

# create outputs to save - this is the random validation
bbs_abun_cv_fitting    <- lapply(bbs_abun_cv_2, function(x) x[[1]][[1]])
bbs_abun_cv_validation <- lapply(bbs_abun_cv_2, function(x) x[[1]][[2]])
bbs_abun_cv_list       <- bbs_abun_cv_2
bbs_abun_cv            <- rbind(do.call(rbind, bbs_abun_cv_fitting), do.call(rbind, bbs_abun_cv_validation))

# save objects
save(bbs_abun_cv, # whole data object as datafrane
     bbs_abun_cv_list, # data object as list
     bbs_abun_cv_fitting, # listed fitting data
     bbs_abun_cv_validation, # list validation data
     file = 'data/bbs_abun_modelling_data_oob_crossValidation.RData')





