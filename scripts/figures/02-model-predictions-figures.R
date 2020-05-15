# script to produce plots of predicted vs. modelled abundance across all model frameworks 

# load in functions in external scripts
source('scripts/evaluating-models/functions/evaluation_functions.R')



# load in data files that contain predicted and known abundances ----

# results are stored on the server for now

predictions_all <- list.files('/Volumes/Simulation/conor/abundance-comparison/results/predictions_all/bind', full.names = T)

#for(model in 1:length(model_levels)){
  
# read in one set of predictions
model_predictions <- readRDS(predictions_all[model])
  
# clean the data
mean_predictions_all <- clean_levels(model_predictions[sample(1:nrow(model_predictions), 100),]) %>% 
  select(-family_grouped_simple, -family_grouped, -family, 
         -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
  
  # create new cross validation level
  mutate(cross_validation = ifelse(grepl('_cv_',predictions_all[model]), 'cv', 'basic')) %>% 
  
  # select final columns for this script
  select(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name,
         verification_observed_mean, verification_predict_mean, 
         validation_observed_mean,   validation_predict_mean)


# function needs to be by a species
input_data <- mean_predictions_all






names(model_predictions)
  
  # why are there more records in the validation datasets than in the verification datasets 
  test <- model_predictions %>% 
    filter(only_abundance==F, abundance_response == 'abunocc') %>% 
    filter(fitted_model == 'rf') %>% 
    filter(species_name == 'Accipiter herodias') %>% 
    rowwise() %>% 
    mutate(new_col = list(keep_non_0s(.$verification_observed_mean, .$verification_predict_mean)))
  
  test %>% 
    rowwise() %>% 
    mutate(metrics = abundance_assessment_metrics(observations = .$new_col[,1], 
                                                  predictions = .$new_col[,2]))
  
  unnest(metrics) %>% 
  bind_cols(clean_files[,1:13], .) %>% 
  mutate(cross_validation = cross_validation_2)    
  
  keep_non_0s <- function(observations, predictions){
    to_keep <- (observations[[1]]!=0)
    observations <- observations[[1]][to_keep]
    predictions  <- predictions[[1]][to_keep]
    return(cbind(observations, predictions))
  }
  
  keep_non_0s(test$verification_observed_mean, test$verification_predict_mean)
  
  test <- test %>% do(new_col =  keep_non_0s(.$verification_observed_mean, .$verification_predict_mean))
  
  test$new_col
  
  input_species <- readRDS('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/data/bbs_all_basic/Accipiter_herodias.RDS')
  table(input_species$fitting$Num>0)
  table(input_species$validation$Num>0)
  table(test$verification_observed_mean[[1]]>0)
  table(test$verification_predict_mean[[1]]>0)
  table(test$validation_observed_mean[[1]]>0)
  table(test$validation_predict_mean[[1]]>0)
  
}

