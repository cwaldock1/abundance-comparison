# script to produce figures that evaluate model performance across different modelling proceedures
# and compare outputs across species

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'viridis')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# source functions and targets ----

source('scripts/figures/functions/model-performance-functions.R')
source('scripts/figures/functions/species-performance-functions.R')

colours = colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(4)

levels = c('glm', 'gam', 'gbm', 'rf')

metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')

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

all_assessments <- lapply(list.files('results/model_assessment_all/validation', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# select only the columns to be used later 
all_assessments <- all_assessments %>% select(-family_grouped_simple, -family_grouped, -family, 
                                              -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  select(dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name, n_abundance,
         Armse, Amae, Dintercept, Dslope, Dpearson, Dspearman, Psd, Pdispersion, Pr2, Evaluation_number, Evaluation_message) %>% 
  
  # change abundance_response
  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))

# find the best model for a species
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
  na.omit(.) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model) %>% 
  mutate(dataset = plyr::revalue(.$dataset, c(bbs_cv = 'bbs', rls_cv = 'rls')))


# standardise measures of model performance ----

all_assessments_relative <- all_assessments %>% 
  group_by(dataset, cross_validation_2, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, 
                                         ~aggregate_metrics(., 
                                                            metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data) %>% 
  ungroup()


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


# compare relative performance of best models across species and sampling properties ----

# get the best models and their assessment values
best_model_assessments <- left_join(best_models , 
                                    all_assessments_relative %>% mutate(cross_validation = .$cross_validation_2))

# check for NAs any species assessments if evaluations could not be assessed 
lapply(best_model_assessments, function(x) table(is.na(x)))

# join in species attributes to model assessments
best_model_assessments <- left_join(rbind(bbs_species, rls_species) %>% filter(species_name %in% unique(best_model_assessments$species_name)), 
                                    best_model_assessments)

# convert values to high and across multiple aspects
best_model_assessments <- best_model_assessments %>% 
  mutate(frequency_group = ifelse(frequency_perc >= 0.75, 'high-freq', 
                                  ifelse(frequency_perc <= 0.25, 'low-freq', NA)), 
         mean_abundance_group = ifelse(mean_abundance_perc >= 0.75, 'high-abun', 
                                       ifelse(mean_abundance_perc <= 0.25, 'low-abun', NA)), 
         sampling_group = ifelse(ecdf(n_abundance)(n_abundance) >= 0.75, 'high-data', 
                                       ifelse(ecdf(n_abundance)(n_abundance) <= 0.25, 'low-data', NA))) %>% 
  na.omit()
  

# create data for plotting
trait_data <- best_model_assessments %>% 
  group_by(dataset, cross_validation, frequency_group, mean_abundance_group, sampling_group) %>% 
  do(accuracy_median = median(.$accuracy, na.rm = T), 
     acc_upr    = quantile(.$accuracy, 0.75, na.rm = T), 
     acc_lwr    = quantile(.$accuracy, 0.25, na.rm = T), 
     
     discrimination_median = median(.$discrimination, na.rm = T), 
     dis_upr    = quantile(.$discrimination, 0.75, na.rm = T), 
     dis_lwr    = quantile(.$discrimination, 0.25, na.rm = T), 
     
     precision_median = median(.$precision, na.rm = T), 
     pre_upr    = quantile(.$precision, 0.75, na.rm = T), 
     pre_lwr    = quantile(.$precision, 0.25, na.rm = T), 
     
     all_median = median(.$aggregated_evaluation_metrics, na.rm = T), 
     all_upr    = quantile(.$aggregated_evaluation_metrics, 0.75, na.rm = T), 
     all_lwr    = quantile(.$aggregated_evaluation_metrics, 0.25, na.rm = T)) %>% 
  unnest() %>% 
  pivot_longer(., cols = c(accuracy_median, discrimination_median, precision_median, all_median), names_to = 'median', values_to = 'median_value') %>% 
  pivot_longer(., cols = c(acc_lwr, dis_lwr, pre_lwr, all_lwr), names_to = 'lwr', values_to = 'lwr_value') %>% 
  pivot_longer(., cols = c(acc_upr, dis_upr, pre_upr, all_upr), names_to = 'upr', values_to = 'upr_value') %>% 
  mutate(evaluation_group = gsub('_median', '', .$median))
  
trait_data$evaluation_group = factor(trait_data$evaluation_group, levels = c('accuracy', 'precision', 'discrimination', 'all'))
trait_data$facet_factor     = factor(paste(trait_data$frequency_group, trait_data$evaluation_group), 
                                     levels = 
                                       unique(paste(trait_data$frequency_group, trait_data$evaluation_group))[c(1,5,2,6,3,7,4,8)])

pdf('figures/species-performance-figures/abun-data-freq-effects/basic.pdf', width = 6, height = 6)
ggplot(data = trait_data %>% 
         filter(cross_validation == 'basic') %>% 
         filter(evaluation_group != 'all')) + 
  geom_rect(aes(fill = sampling_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.5) + 
  geom_point(aes(x = mean_abundance_group, y = median_value, group = dataset, colour = dataset), position=position_dodge(width=0.5)) + 
  geom_line(aes(x = mean_abundance_group, y = median_value, group = dataset, colour = dataset), position=position_dodge(width=0.5)) + 
  facet_grid(evaluation_group ~  sampling_group + frequency_group) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.05), 
        panel.grid.major.x = element_blank()) + 
  scale_colour_manual(values = colours[c(4,1)]) + 
  scale_fill_manual(values = c('grey90', 'white')) + 
  xlab(NULL) + 
  ylab('relative model rank') + 
  guides(fill = FALSE)
dev.off()





# fit models for best species and produce interaction plots ---- 

# get the best models and their assessment values
best_model_assessments <- left_join(best_models , 
                                    all_assessments_relative %>% mutate(cross_validation = .$cross_validation_2))

# check for NAs any species assessments if evaluations could not be assessed 
lapply(best_model_assessments, function(x) table(is.na(x)))

# join in species attributes to model assessments
best_model_assessments <- left_join(rbind(bbs_species, rls_species) %>% filter(species_name %in% unique(best_model_assessments$species_name)), 
                                    best_model_assessments)

# create ranked sampling column
best_model_assessments$sampling_n_perc <- ecdf(best_model_assessments$Evaluation_number)(best_model_assessments$Evaluation_number)

# create plots and model outputs 
best_model_assessments %>% 
  group_by(dataset, cross_validation) %>% 
  nest() %>% 
  .[1,] %>% 
  mutate(test = purrr::map(data, ~interaction_trait_plots(., 
                                 directory = 'figures/species-performance-figures/trait-interaction-models/', 
                                 name      = paste(dataset, cross_validation, sep = '-'), 
                                 dataset   = dataset, 
                                 width_marginal     = 3, 
                                 height_marginal    = 3)))

best_model_assessments %>% 
  group_by(dataset, cross_validation) %>% 
  nest() %>% 
  .[2,] %>% 
  mutate(test = purrr::map(data, ~interaction_trait_plots(., 
                                                          directory = 'figures/species-performance-figures/trait-interaction-models/', 
                                                          name      = paste(dataset, cross_validation, sep = '-'), 
                                                          width_marginal     = 3, 
                                                          height_marginal    = 3, 
                                                          dataset   = dataset)))

best_model_assessments %>% 
  group_by(dataset, cross_validation) %>% 
  nest() %>% 
  .[3,] %>% 
  mutate(test = purrr::map(data, ~interaction_trait_plots(., 
                                                          directory = 'figures/species-performance-figures/trait-interaction-models/', 
                                                          name      = paste(dataset, cross_validation, sep = '-'), 
                                                          width_marginal     = 3, 
                                                          height_marginal    = 3, 
                                                          dataset   = dataset)))

best_model_assessments %>% 
  group_by(dataset, cross_validation) %>% 
  nest() %>% 
  .[4,] %>% 
  mutate(test = purrr::map(data, ~interaction_trait_plots(., 
                                                          directory = 'figures/species-performance-figures/trait-interaction-models/', 
                                                          name      = paste(dataset, cross_validation, sep = '-'), 
                                                          width_marginal     = 3, 
                                                          height_marginal    = 3, 
                                                          dataset   = dataset)))



# testing creating marginal effect plots for significant terms 
best_model_assessments %>% 
   group_by(dataset, cross_validation) %>% 
   nest() %>% .$data %>% .[[2]] -> plot_data
 





