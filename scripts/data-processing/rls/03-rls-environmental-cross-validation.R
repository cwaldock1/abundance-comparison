# setting up the cross validations based on temperature covariates

# Load packages ----
lib_vect <- c('tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# Load data to set up new cross validations ----

# load all objects
load('data/rls_abun_modelling_data_v2.RData')

# load covariates
load('data/rls_covariates.RData')

# 50 species: Join temperature data and set up new cross-validations ----

# join in sst_mean
rls_abun_sst <- left_join(rls_abun %>% select(-set), rls_xy %>% select(SiteLongitude, SiteLatitude, sst_mean))

# lapply subsets 
rls_abun_sst_2 <- lapply(1:length(unique(rls_abun_sst$TAXONOMIC_NAME)), function(x){
  
  species_data <- rls_abun_sst %>% 
    filter(TAXONOMIC_NAME == unique(rls_abun_sst$TAXONOMIC_NAME)[x])
  
  all_absences  <- species_data %>% filter(Num == 0)  # absences
  all_presences <- species_data %>% filter(Num != 0)  # presences
  
  sample_size <- round(round(nrow(species_data)*0.8) / 2) # get sample size for cross validations
  
  # estimate 80th quantile in absences and presences
  absences_sst_quantile <- as.numeric(quantile(all_absences$sst_mean, 0.8, na.rm = T))
  presences_sst_quantile <- as.numeric(quantile(all_presences$sst_mean, 0.8, na.rm = T))
  
  absences_sample_fitting    <- all_absences[ -which(all_absences$sst_mean  > absences_sst_quantile),]
  presence_sample_fitting    <- all_presences[-which(all_presences$sst_mean > presences_sst_quantile),]
  absences_sample_validation <- all_absences[ which(all_absences$sst_mean  > absences_sst_quantile),]
  presence_sample_validation <- all_presences[which(all_presences$sst_mean > presences_sst_quantile),]
  
  fitting    <- rbind(absences_sample_fitting, presence_sample_fitting)
  validation <- rbind(absences_sample_validation, presence_sample_validation)
  
  fitting$set    <- 'fitting'
  validation$set <- 'validation'
  
  list_outputs <- list(list(fitting = fitting, validation = validation))
  names(list_outputs) <- unique(species_data$TAXONOMIC_NAME)
  
  return(list_outputs)
  
  }) # end of lapply

# create outputs to save - this is the random validation
rls_abun_sst_fitting    <- lapply(rls_abun_sst_2, function(x) x[[1]][[1]])
rls_abun_sst_validation <- lapply(rls_abun_sst_2, function(x) x[[1]][[2]])
rls_abun_sst_list       <- rls_abun_sst_2
rls_abun_sst            <- rbind(do.call(rbind, rls_abun_sst_fitting), do.call(rbind, rls_abun_sst_validation))

# create directory
dir.create('data/rls_50_CV')

# save objects
save(rls_abun_sst, # whole data object as datafrane
     rls_abun_sst_list, # data object as list
     rls_abun_sst_fitting, # listed fitting data
     rls_abun_sst_validation, # list validation data
     file = 'data/rls_50_CV/rls_abun_modelling_data_sst_oob_crossValidation.RData')

# all species cross validation saves ----

all_sp <- list.files('data/rls_all_basic', full.names = T)

lapply(1:length(all_sp), function(x){
  print(x)
  if(file.exists(gsub('rls_all_basic', 'rls_all_CV', all_sp[x]))){return(NULL)}
  
  rls_abun <- do.call('rbind', readRDS(all_sp[x]))
  
  species_data <- left_join(rls_abun %>% select(-set), rls_xy %>% dplyr::select(SiteLongitude, SiteLatitude, sst_mean))
  
  all_absences  <- species_data %>% filter(Num == 0)  # absences
  all_presences <- species_data %>% filter(Num != 0)  # presences
  
  sample_size <- round(round(nrow(species_data)*0.8) / 2) # get sample size for cross validations
  
  # estimate 80th quantile in absences and presences
  absences_cv_quantile <- as.numeric(quantile(all_absences$sst_mean, 0.8, na.rm = T))
  presences_cv_quantile <- as.numeric(quantile(all_presences$sst_mean, 0.8, na.rm = T))
  
  absences_sample_fitting    <- all_absences[ -which(all_absences$sst_mean  > absences_cv_quantile),]
  presence_sample_fitting    <- all_presences[-which(all_presences$sst_mean > presences_cv_quantile),]
  absences_sample_validation <- all_absences[ which(all_absences$sst_mean  > absences_cv_quantile),]
  presence_sample_validation <- all_presences[which(all_presences$sst_mean > presences_cv_quantile),]
  
  fitting    <- rbind(absences_sample_fitting, presence_sample_fitting)
  validation <- rbind(absences_sample_validation, presence_sample_validation)
  
  fitting$set    <- 'fitting'
  validation$set <- 'validation'
  
  list_outputs <- list(fitting = fitting, validation = validation)
  
  # create output directory
  dir.create('data/rls_all_CV/')
  
  # save RDS
  saveRDS(list_outputs, paste0('data/rls_all_CV/', 
                               as.character(gsub(' ','_',unique(species_data$TAXONOMIC_NAME))),
                               '.RDS'))
}
)




