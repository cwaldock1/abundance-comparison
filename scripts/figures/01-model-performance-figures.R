# script to produce figures that evaluate model performance across different modelling proceedures

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# source functions ----

source('scripts/figures/functions/model-performance-functions.R')

colours = colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(4)

levels = c('glm', 'gam', 'gbm', 'rf')

metrics = c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2')

targets = list(Armse    = c(0, -20, 20), 
               Amae     = c(0, -20, 20), 
               Dintercept  = c(0, -5, 20), 
               Dslope      = c(1,  0, 2), 
               Dpearson    = c(1,  0, 1), 
               Dspearman   = c(1,  0, 1), 
               Psd         = c(0,  0, 100), 
               Pdispersion = c(1,  0, 20), 
               Pr2         = c(1,  0, 1))

# load in evaluation data ----

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




# full plots across all species and model scenarios for supporting materials ----

### REMOVE THIS PLOTTING PROCEEDURE AND FUNCTION AND ADAPT FOR EVALUATING SPECIES TRAITS AND ABUNDANCE EFFECTS ON MODEL PERFORMANCE
nested_assessments <- all_assessments %>% 
  group_by(dataset, 
           abundance_response, 
           cross_validation) %>% 
  nest()

lapply(1:nrow(nested_assessments), 
       function(x){
         all_model_plots(plot_data = nested_assessments$data[[x]], 
                         directory = paste0('figures/all_model_plots/', 
                                            nested_assessments$cross_validation[x], '/', 
                                            nested_assessments$abundance_response[x]))})

# plot of all metrics for each individual model type ----

nested_assessments <- all_assessments %>% 
  group_by(dataset, cross_validation_2) %>% 
  nest()

# plot_data <- nested_assessments$data[[1]] # for testing functions

# apply across both validation types
lapply(1:nrow(nested_assessments), 
       function(x){
         all_model_plots_v2(plot_data = nested_assessments$data[[x]], 
                            directory = paste0('figures/model-performance-figures/all_model_all_metric_combined/'), 
                            name =  paste0(nested_assessments$dataset[x], '_', nested_assessments$cross_validation_2[x]), 
                            height = 14, 
                            width = 14)})


# plot aggregated models within metrics types ----

nested_assessments <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'all',
                                     levels = levels, 
                                     metrics = metrics,
                                     targets = targets,
                                     colours = colours,
                                     directory  = 'figures/model-performance-figures/all_model_boxplots/', 
                                     name = nested_assessments$cross_validation_2[x], 
                                     width = 12, 
                                     height = 10)})

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     levels = levels, 
                                     metrics = metrics,
                                     targets = targets,
                                     colours = colours,
                                     response = 'abundance_response',
                                     directory  = 'figures/model-performance-figures/all_model_boxplots/', 
                                     name = paste0(nested_assessments$cross_validation_2[x], '--abundance_response'), 
                                     width = 12, 
                                     height = 10)})

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'fitted_model',
                                     levels = levels, 
                                     metrics = metrics,
                                     targets = targets,
                                     colours = colours,
                                     directory  = 'figures/model-performance-figures/all_model_boxplots/', 
                                     name = paste0(nested_assessments$cross_validation_2[x], '--fitted_model'), 
                                     width = 12, 
                                     height = 10)})


# plots of rescaled values comparing between models ----

all_assessments_relative <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

# aggreagte basic models
plot_all_aggregated(all_assessments_relative$data[[1]] %>% 
  group_by(dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
    mutate(dataset = all_assessments_relative$data[[1]]$dataset), 
  directory = 'figures/model-performance-figures/all_model_rescaled', 
  colours = colours,
  name = 'basic')

# aggregate oob_cv models
plot_all_aggregated(all_assessments_relative$data[[2]] %>% 
                      group_by(dataset) %>% 
                      nest() %>% 
                      mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.))) %>% 
                      .$metric_aggregation %>% 
                      do.call(rbind, .) %>% 
                      mutate(dataset = all_assessments_relative$data[[2]]$dataset), 
                    directory = 'figures/model-performance-figures/all_model_rescaled', 
                    colours = colours,
                    name = 'cv')

# plots of ranked models ----

nested_assessments <- split(all_assessments, f = all_assessments$cross_validation)

lapply(1:length(nested_assessments), 
       function(x){
         rank_plots(plot_data = nested_assessments[[x]], 
                    levels = levels, 
                    metrics = metrics,
                    targets = targets,
                    colours = colours,
                    directory  = 'figures/model-performance-figures/rank_plots/', 
                    name = unique(paste0(nested_assessments[[x]]$dataset, '_', nested_assessments[[x]]$cross_validation_2)), 
                    width = 12, 
                    height = 10)})


plot_data <- nested_assessments$data[[1]]






# plots of rescaled values comparing between cross-validation scenarios ----

# for cross validations we present the response ratio
cv_data <- all_assessments %>% 
  group_by(dataset, 
           abundance_response) %>% 
  nest() 

plot_data <- cv_data$data[[1]]
  
lapply(1:nrow(cv_data),
       function(x){
         plot_cv_comparison(plot_data = cv_data$data[[x]], 
                            directory  = paste0('figures/cv_comparisons/', cv_data$dataset[[x]],'/', cv_data$abundance_response[[x]]), 
                            width = 10, 
                            height = 10)})

plot_cv_comparison











