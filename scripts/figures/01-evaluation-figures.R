# script to produce figures

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# source functions ----

source('scripts/figures/functions/figure-functions.R')

# load in evaluation data ----

all_assessments <- lapply(list.files('results/model_assessment', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# full plots across all species and model scenarios for supporting materials ----

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

# full plots but combine across species within metrics types ----

nested_assessments <- all_assessments %>% 
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  group_by(cross_validation_2) %>% 
  nest()

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'all',
                                     directory  = 'figures/all_model_plots_median/', 
                                     name = nested_assessments$cross_validation_2[x], 
                                     width = 12, 
                                     height = 10)})

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'abundance_response',
                                     directory  = 'figures/all_model_plots_median/', 
                                     name = paste0(nested_assessments$cross_validation_2[x], '--abundance_response'), 
                                     width = 12, 
                                     height = 10)})

lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'fitted_model',
                                     directory  = 'figures/all_model_plots_median/', 
                                     name = paste0(nested_assessments$cross_validation_2[x], '--fitted_model'), 
                                     width = 12, 
                                     height = 10)})



# full plots on the robust PCA of metrics ----

plot_data = nested_assessments$data[[1]]


