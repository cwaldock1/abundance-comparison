# script to produce figures that evaluate model performance across different modelling proceedures
# and compare outputs across species

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# source functions and targets ----

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

# load in species' attributes ----

bbs_species <- readRDS('data/bbs_species_properties.RDS') %>% dplyr::rename(.,  species_name = TAXONOMIC_NAME)
rls_species <- readRDS('data/rls_species_properties.RDS') %>% dplyr::rename(.,  species_name = TAXONOMIC_NAME)

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


# standardise measures of model performance ----

all_assessments_relative <- all_assessments %>% 
  group_by(dataset, cross_validation_2, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data)

# fit model across all species ----

bbs_basic <- all_assessments_relative %>% 
  filter(cross_validation == 'bbs_basic') %>% 
  ungroup() %>% 
  left_join(., bbs_species) %>% 
  plot_trait_models(input_data = ., 
                    directory = paste0('figures/species-performance-figures/', unique(.$cross_validation), '/'), 
                    name = 'bbs_basic_species_performance')


all_assessments_relative %>% 
  group_by(cross_validation) %>% 
  left_join(., rbind(bbs_species, rls_species)) %>% 
  do(plot_trait_models(input_data = ., 
                       directory = paste0('figures/species-performance-figures/', unique(.$cross_validation), '/'), 
                       name = paste0(unique(.$dataset), '_', unique(.$cross_validation_2), '_species_performance')))




