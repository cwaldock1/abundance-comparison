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

all_assessments <- lapply(list.files('results/model_assessment_all/validation', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# select only the columns to be used later 
all_assessments <- all_assessments %>% select(-family_grouped_simple, -family_grouped, -family, 
                                              -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  select(dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name,
         Armse, Amae, Dintercept, Dslope, Dpearson, Dspearman, Psd, Pdispersion, Pr2, Evaluation_number, Evaluation_message) %>% 
  
  # change abundance_response
  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))


# overall this code finds the best fitted model for each species within a type of model and cross-validation and dataset combination. Use this object to left join to the true predictions and create plots
# comparing each modelling frameworks overall predictions vs. observations. 
best_models <- all_assessments %>% 
  select(-Armse, -Psd) %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 
                                                                              'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
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

# results are downloaded from the server

predictions_all <- list.files('results/predictions_all/bind', full.names = T)

model_type <- c('bbs_basic_all', 'bbs_cv_all', 'rls_basic_all', 'rls_cv_all')

for(model in 1:length(model_type)){

  
      # get subsets
      predictions_subsets <- predictions_all[grep(model_type[model], predictions_all)]
      
      # read in subsets (easier on the memory..)
      cat('reading subset 1')
      subset1 <- readRDS(predictions_subsets[1])
      cat('reading subset 2')
      subset2 <- readRDS(predictions_subsets[2])
      cat('reading subset 3')
      subset3 <- readRDS(predictions_subsets[3])
    
      # bind together subsets
      model_predictions <- rbind(subset1, subset2, subset3)
      
      rm(subset1, subset2, subset3)
      
      # clean the data
      mean_predictions_all <- clean_levels(model_predictions) %>% #clean_levels(model_predictions[sample(1:nrow(model_predictions), 500),]) %>% 
        select(-family_grouped_simple, -family_grouped, -family, 
               -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
        
        # create new cross validation level
        mutate(cross_validation = ifelse(grepl('_cv_',predictions_subsets)[1], 'cv', 'basic')) %>% 
        
        # select final columns for this script
        select(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name,
               verification_observed_mean, verification_predict_mean, 
               validation_observed_mean,   validation_predict_mean)
      
      
      rm(model_predictions)
      
      # get seperate data for validations and verifications
      validation_data   <- mean_predictions_all %>% 
        select(names(.)[!grepl('verification', names(.))]) %>% 
        group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
        do(validation_observed_mean = .$validation_observed_mean[[1]][.$validation_observed_mean[[1]] > 0], 
           validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0]) %>% 
        ungroup()
      
      # remove species with identical observations (cannot fit a model here)
      validation_data <- validation_data[-which(sapply(validation_data$validation_observed_mean, function(x) length(unique(x))==1)),]
      
      # edit plot level text to be shorter
      validation_data$plot_level <- ifelse(validation_data$abundance_response == 'abunocc_2stage', gsub('abunocc_2stage', 'ao-2', validation_data$plot_level), validation_data$plot_level)

      # change abundance_response
      validation_data <- validation_data %>%  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))
      
      # create plots 
      observed_predicted_plot(input_data = validation_data, 
                              rescale = T, 
                              model_level = 'plot_level', # options are fitted_model or plot_level
                              directory   = 'figures/model-prediction-figures/validation_all_models_plot_level/',
                              name        = unique(paste(validation_data$dataset, validation_data$cross_validation, sep = '-')), 
                              width = 8000, 
                              height = 8000, 
                              nbins = 20, 
                              upper_limit = 1.5)
      
      
      
}


# plot observed and predicted for best fitted models across multiple dataframe types ----

# overall this code finds the best fitted model for each species within a type of model and cross-validation and dataset combination. Use this object to left join to the true predictions and create plots
# comparing each modelling frameworks overall predictions vs. observations. 
best_models_overall <- all_assessments %>% 
  select(-Armse, -Psd) %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 
                                                                              'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model)

predictions_all <- list.files('results/predictions_all/bind', full.names = T)

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
    mutate(cross_validation = ifelse(grepl('_cv_',predictions_all[i*3]), 'cv', 'basic')) %>% 
    
    # select final columns for this script
    select(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name,
           verification_observed_mean, verification_predict_mean, 
           validation_observed_mean,   validation_predict_mean)

  # join in and select models
  model_predictions_join <- left_join(best_models_overall %>% 
                                        mutate(dataset = plyr::revalue(.$dataset, c(bbs_cv = 'bbs', rls_cv = 'rls'))), 
                                      mean_predictions_all)
  model_predictions_join <- model_predictions_join[!unlist(lapply(model_predictions_join$verification_observed_mean, is.null)),]
  
  # get validations
  validation_data   <- model_predictions_join %>% 
    select(names(.)[!grepl('verification', names(.))]) %>% 
    group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
    do(validation_observed_mean = .$validation_observed_mean[[1]][.$validation_observed_mean[[1]] > 0], 
       validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0]) %>% 
    ungroup()
  
  # remove species with identical observations (cannot fit a model here)
  if(length(which(sapply(validation_data$validation_observed_mean, function(x) length(unique(x))==1)))!=0){
    validation_data <- validation_data[-which(sapply(validation_data$validation_observed_mean, function(x) length(unique(x))==1)),]
    }
  
  # edit plot level text to be shorter
  validation_data$plot_level <- ifelse(validation_data$abundance_response == 'abunocc_2stage', gsub('abunocc_2stage', 'ao-2', validation_data$plot_level), validation_data$plot_level)
  
  # change abundance_response
  validation_data <- validation_data %>%  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))
  
  # testing why there are nulls
  table(sapply(validation_data$validation_observed_mean, is.null))
  table(sapply(validation_data$validation_predict_mean, is.null))
  
  
  # create plots 
   observed_predicted_plot(input_data = validation_data, 
                           rescale = T, 
                           model_level = 'plot_level', # options are fitted_model or plot_level
                           directory   = 'figures/model-prediction-figures/validation_best_models_plot_level',
                           name        = unique(paste(validation_data$dataset, validation_data$cross_validation, sep = '-')), 
                           width = 8000, 
                           height = 8000, 
                           nbins = 20, 
                           upper_limit = 1.5)
  
  observed_predicted_plot(input_data = validation_data, 
                          rescale = T, 
                          model_level = 'aggregated', # options are fitted_model or plot_level
                          directory   = 'figures/model-prediction-figures/validation_best_models_aggregated',
                          name        = unique(paste(validation_data$dataset, validation_data$cross_validation, sep = '-')), 
                          width = 2000, 
                          height = 2000, 
                          nbins = ifelse(unique(validation_data$dataset) == 'rls', 20, 30), 
                          upper_limit = ifelse(unique(validation_data$dataset) == 'rls', 1.5, 2), 
                          option = ifelse(unique(validation_data$dataset) == 'rls', 'viridis', 3))
  

}
