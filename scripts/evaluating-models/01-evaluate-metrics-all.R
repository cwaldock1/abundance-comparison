# extract evaluation 

library(tidyverse)

# all functions for evaluating outputs
source('scripts/evaluating-models/functions/evaluation_functions.R')

# # get folders of interest
# result_folders <- c(list.files('results', recursive = F, full.names = T, pattern = 'cv_all'), 
#                     list.files('results', recursive = F, full.names = T, pattern = 'basic_all'))
# 
# # get level 2 folders
# results_folders_2 <- sapply(result_folders, function(x){list.files(x, pattern = 'predictions_', full.names = T)})
# 
# for(folder in 1:length(results_folders_2)){
#   
#   # here in the future run iteratively for each folder of interest
#   predictions_folder_files <- lapply(results_folders_2[,folder], function(x){print(x); list.files(x, full.names = T)})
#   
#   # in the future, run this over each modelling subset for scalability
#   
#   bind_files <- lapply(predictions_folder_files, function(x){
#     bind_results_save(x,
#                       directory = 'results/predictions_all/bind/', 
#                       name = paste0(strsplit(x, '/')[[1]][2], '_',strsplit(x, '/')[[1]][3]))})
#   
#   # Clean data levels ----
#   
#   # some of the encodings output from the models aren't well match so match values here across modelling frameworks using the clean_levels functions
#   
#   clean_files <- lapply(bind_files, clean_levels)
#   
# }

# perform evaluation on bindings ----

bind_files <- list.files('results/predictions_all/bind', full.names = T)

for(i in 1:length(bind_files)){
  
  cat(paste0('performing assessments for model set ', i))
  
  clean_files <- readRDS(bind_files[i])
  
  model_type <- gsub('results/predictions_all/bind/|_all_predictions_|.rds','_', bind_files[i])
  model_type <- substring(model_type, 2, nchar(model_type)-1)
  
  cross_validation_2 <- gsub('results/predictions_all/bind/|_all_predictions_abun.rds|_all_predictions_abunocc_2stage.rds|_all_predictions_abunocc.rds','', bind_files[i])
  
  # perform validation (out the bag)
  model_assessment <- clean_files %>% 
        rowwise() %>% 
        do(metrics = abundance_assessment_metrics(predictions   = .$validation_predict_mean, 
                                                  observations  = .$validation_observed_mean)) %>% 
        unnest(metrics) %>% 
        bind_cols(clean_files[,1:13], .) %>% 
        mutate(cross_validation = cross_validation_2)
  
  dir.create('results/model_assessment_all/validation', recursive = T)
  
  saveRDS(model_assessment, file = paste0('results/model_assessment_all/validation/', model_type, '.rds'))
  
  # perform verification (in the bag)
  model_assessment <- clean_files %>% 
    rowwise() %>% 
    do(metrics = abundance_assessment_metrics(predictions   = .$verification_predict_mean, 
                                              observations  = .$verification_observed_mean)) %>% 
    unnest(metrics) %>% 
    bind_cols(clean_files[,1:13], .) %>% 
    mutate(cross_validation = cross_validation_2)
  
  dir.create('results/model_assessment_all/verification', recursive = T)
  
  saveRDS(model_assessment, file = paste0('results/model_assessment_all/verification/', model_type, '.rds'))
  
  
  }

#  clean_files <- rbind_all(clean_files)
#  
#  # Calculate assessment metrics ----
#  
#  model_assessment <- clean_files %>% 
#    rowwise() %>% 
#    do(metrics = abundance_assessment_metrics(predictions   = .$validation_predict_mean, 
#                                              observations  = .$validation_observed_mean)) %>% 
#    unnest(metrics) %>% 
#    bind_cols(clean_files[,1:13], .) %>% 
#    mutate(cross_validation = gsub('results/predictions/','', result_folders[folder]))
#  
#  # Attached abundance information into metrics ----
#  
#  model_assessment <- left_join(model_assessment, 
#                                readRDS(paste0('data/', if(model_assessment$dataset == 'bbs'){'bbs_species_properties.RDS'}else{'rls_species_properties.RDS'})) %>% 
#                                  rename(., species_name = TAXONOMIC_NAME))
#  
#  dir.create('results/model_assessment_all/', recursive = T)
#  
#  saveRDS(model_assessment, file = paste0('results/model_assessment_all/', strsplit(model_assessment$cross_validation, '/')[[1]][3], '.rds'))
#  
#}

