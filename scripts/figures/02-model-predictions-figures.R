# script to produce plots of predicted vs. modelled abundance across all model frameworks 

# load libraries ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)


# load in functions in external scripts ----

source('scripts/evaluating-models/functions/evaluation_functions.R')
source('scripts/figures/functions/model-performance-functions.R')
source('scripts/figures/functions/model-prediction-functions.R')

# select best fitted model for each model type based on a concensus metrics ----

# read and clean assessment data as in the script 01-model-performance-figures.R

all_assessments <- lapply(list.files('results/model_assessment', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# select only the columns to be used later 
all_assessments <- all_assessments %>% select(-family_grouped_simple, -family_grouped, -family, 
                                              -transformation, -n_absence, -n_boot_absence, -abundance_response_simple, 
                                              -mean_abundance, -frequency, -mean_abundance_perc, -frequency_perc, -n_abundance) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  select(dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name,
         Armse, Amae, Dintercept, Dslope, Dpearson, Dspearman, Psd, Pdispersion, Pr2)

# overall this code finds the best fitted model for each species within a type of model and cross-validation and dataset combination. Use this object to left join to the true predictions and create plots
# comparing each modelling frameworks overall predictions vs. observations. 
best_models <- all_assessments %>% 
  select(-Armse, -Psd) %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name, fitted_model, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model)


# plot observed vs. predicted abundances for all model types ----

# results are stored on the server for now

predictions_all <- list.files('/Volumes/Simulation/conor/abundance-comparison/results/predictions_all/bind', full.names = T)


for(model in 1:length(predictions_all)){

      # read in one set of predictions
      model_predictions <- readRDS(predictions_all[model])
        
      # clean the data
      mean_predictions_all <- clean_levels(model_predictions) %>% #clean_levels(model_predictions[sample(1:nrow(model_predictions), 500),]) %>% 
        select(-family_grouped_simple, -family_grouped, -family, 
               -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
        
        # create new cross validation level
        mutate(cross_validation = ifelse(grepl('_cv_',predictions_all[model]), 'cv', 'basic')) %>% 
        
        # select final columns for this script
        select(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name,
               verification_observed_mean, verification_predict_mean, 
               validation_observed_mean,   validation_predict_mean)
      
      
      # get seperate data for validations and verifications
      verification_data <- mean_predictions_all %>% 
        select(cross_validation, names(.)[!grepl('validation', names(.))]) %>%   
        group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
        do(verification_observed_mean = .$verification_observed_mean[[1]][.$verification_observed_mean[[1]] > 0], 
           verification_predict_mean = .$verification_predict_mean[[1]][.$verification_observed_mean[[1]] > 0]) %>% 
        ungroup()
      
      validation_data   <- mean_predictions_all %>% 
        select(names(.)[!grepl('verification', names(.))]) %>% 
        group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
        do(validation_observed_mean = .$validation_observed_mean[[1]][.$validation_observed_mean[[1]] > 0], 
           validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0]) %>% 
        ungroup()
      
      # edit plot level text to be shorter
      validation_data$plot_level <- ifelse(validation_data$abundance_response == 'abunocc_2stage', gsub('abunocc_2stage', 'a2stage', validation_data$plot_level), validation_data$plot_level)
      verification_data$plot_level <- ifelse(verification_data$abundance_response == 'abunocc_2stage', gsub('abunocc_2stage', 'a2stage', verification_data$plot_level), verification_data$plot_level)
      
      # create plots 
      observed_predicted_plot(input_data = validation_data, 
                              rescale = T, 
                              model_level = 'plot_level', # options are fitted_model or plot_level
                              directory   = 'figures/model-prediction-figures/validation/',
                              name        = unique(paste(validation_data$dataset, validation_data$cross_validation, validation_data$abundance_response, sep = '-')), 
                              width = 2000, 
                              height = 2000)
      
      
      observed_predicted_plot(input_data = verification_data, 
                              rescale = T, 
                              model_level = 'plot_level', # options are fitted_model or plot_level
                              directory   = 'figures/model-prediction-figures/verification/',
                              name        = unique(paste(verification_data$dataset, verification_data$cross_validation, verification_data$abundance_response, sep = '-')), 
                              width = 2000, 
                              height = 2000)
      
      
}


# plot observed and predicted for best fitted models across multiple dataframe types ----

predictions_all <- list.files('/Volumes/Simulation/conor/abundance-comparison/results/predictions_all/bind', full.names = T)

model_type <- c('bbs_basic_all', 'bbs_cv_all', 'rls_basic_all', 'rls_cv_all')

for(i in 1:length(model_type)){
  
  # get subsets
  predictions_subsets <- predictions_all[grep(model_type[i], predictions_all)]
  
  # read in subsets (easier on the memory..)
  cat('reading subset 1')
  subset1 <- readRDS(predictions_subsets[1])
  cat('reading subset 2')
  subset2 <- readRDS(predictions_subsets[2])
  cat('reading subset 3')
  subset3 <- readRDS(predictions_subsets[3])
  
  # bind together subsets
  model_predictions <- rbind(subset1, subset2, subset3)
  
  # clean the data
  mean_predictions_all <- clean_levels(model_predictions) %>% #clean_levels(model_predictions[sample(1:nrow(model_predictions), 500),]) %>% 
    select(-family_grouped_simple, -family_grouped, -family, 
           -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
    
    # create new cross validation level
    mutate(cross_validation = ifelse(grepl('_cv_',predictions_all[i*3]), 'oob_cv', 'basic')) %>% 
    
    # select final columns for this script
    select(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name,
           verification_observed_mean, verification_predict_mean, 
           validation_observed_mean,   validation_predict_mean)
  
  # join in and select models
  model_predictions_join <- left_join(best_models, mean_predictions_all)
  model_predictions_join <- model_predictions_join[!unlist(lapply(model_predictions_join$verification_observed_mean, is.null)),]
  
  # get seperate data for validations and verifications
  verification_data <- model_predictions_join %>% 
    select(cross_validation, names(.)[!grepl('validation', names(.))]) %>%   
    group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
    do(verification_observed_mean = .$verification_observed_mean[[1]][.$verification_observed_mean[[1]] > 0], 
       verification_predict_mean = .$verification_predict_mean[[1]][.$verification_observed_mean[[1]] > 0]) %>% 
    ungroup()
  
  validation_data   <- model_predictions_join %>% 
    select(names(.)[!grepl('verification', names(.))]) %>% 
    group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
    do(validation_observed_mean = .$validation_observed_mean[[1]][.$validation_observed_mean[[1]] > 0], 
       validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0]) %>% 
    ungroup()
  
  
  # create plots 
  observed_predicted_plot(input_data = verification_data, 
                          rescale = F, 
                          model_level = 'fitted_model', # options are fitted_model or plot_level
                          directory   = 'figures/model-prediction-figures/best_fitted_model/verification/',
                          name        = unique(paste(verification_data$dataset, verification_data$cross_validation, sep = '-')), 
                          width = 1000, 
                          height = 1000, 
                          upper_limit = 1.5,
                          nbins = 10)
  
  # create plots 
  observed_predicted_plot(input_data = validation_data, 
                          rescale = F, 
                          model_level = 'fitted_model', # options are fitted_model or plot_level
                          directory   = 'figures/model-prediction-figures/best_fitted_model/validation/',
                          name        = unique(paste(validation_data$dataset, validation_data$cross_validation, sep = '-')), 
                          width = 1000, 
                          height = 1000, 
                          upper_limit = 1.5, 
                          nbins = 10)
  

}
