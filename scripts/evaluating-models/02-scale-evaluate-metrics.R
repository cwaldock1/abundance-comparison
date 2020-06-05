# extract evaluation at different spatial scales

# set up inputs
# get trailing arguments
args = commandArgs(trailingOnly=TRUE)
j <- as.numeric(args[1]) 
i <- as.numeric(args[2]) 
print(j)
print(i)


library(tidyverse)


# all functions for evaluating outputs
source('scripts/evaluating-models/functions/evaluation_functions.R')

### perform evaluation on bindings ----

bind_files <- list.files('results/predictions_all/bind', full.names = T)

spatial_scale <- c(0.1, 1, 5, 10, 20, 35, 50)

#for(j in 1:length(spatial_scale)){
  
  save_directory <- paste0('results/model_assessment_scale/', 'spatial_scale_', spatial_scale[j], '/')

#for(i in 1:length(bind_files)){
  
  cat(paste0('performing assessments for model set ', i))
  
  clean_files <- readRDS(bind_files[i])
  
  clean_files <- clean_levels(clean_files)
  
  model_type <- gsub('results/predictions_all/bind/|_all_predictions_|.rds','_', bind_files[i])
  model_type <- substring(model_type, 2, nchar(model_type)-1)
  
  cross_validation_2 <- gsub('results/predictions_all/bind/|_all_predictions_abun.rds|_all_predictions_abunocc_2stage.rds|_all_predictions_abunocc.rds','', bind_files[i])
  
  # perform validation (out the bag)
  model_assessment <- clean_files %>% 
        rowwise() %>% 
        do(metrics = abundance_assessment_metrics(predictions   = .$validation_predict_mean, 
                                                  observations  = .$validation_observed_mean, 
                                                  locations     = .$validation_locations, 
                                                  scale = spatial_scale[j])) %>% 
        unnest(metrics) %>% 
        bind_cols(clean_files[,1:13], .) %>% 
        mutate(cross_validation = cross_validation_2)
  
  dir.create(save_directory, recursive = T)
  
  saveRDS(model_assessment, file = paste0(save_directory, model_type, '.rds'))
  
#}

#}

