# script to produce figures that evaluate model performance across different modelling proceedures

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'grid')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# source functions ----

source('scripts/figures/functions/model-performance-functions.R')

colours = colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(4)

colours_V2 = c('orange', 'purple', 'light pink')


levels = c('glm', 'gam', 'gbm', 'rf')

metrics = c('Amae_rel_mean', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')

targets = list(#Armse    = c(0, -20, 20), 
               #Amae     = c(0, -20, 20), 
               Amae_rel_mean = c(1,0,3),
               Dintercept  = c(0, -5, 20), 
               Dslope      = c(1,  0, 2), 
               Dpearson    = c(1,  0, 1), 
               Dspearman   = c(1,  0, 1), 
               #Psd         = c(0,  0, 100), 
               Pdispersion = c(1,  0, 20)#, 
               #Pr2         = c(1,  0, 1)
               )

# load in evaluation data ----

all_assessments <- lapply(list.files('results/model_assessment_all_R2/validation', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# select only the columns to be used later 
all_assessments <- all_assessments %>% select(-family_grouped_simple, -family_grouped, -family, 
                           -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  select(dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name,
         metrics, Evaluation_number, Evaluation_message) %>% 
  
  # change abundance_response
  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))

# create plots of % of species per model type ----

# summarise models and species

length(unique(all_assessments$plot_level))
length(unique(all_assessments$species_name))

# create a full set of models and species

full_model_set <- all_assessments %>% 
  group_by(dataset, fitted_model, abundance_response) %>% 
  do(expand.grid(plot_level = .$plot_level %>% unique(), species_name = .$species_name %>% unique()))
  
full_model_set <- rbind(full_model_set %>% mutate(cross_validation_2 = 'cv'), full_model_set %>% mutate(cross_validation_2 = 'basic'))

full_model_join <- full_join(full_model_set, na.omit(all_assessments) %>% 
                               select(plot_level, species_name, dataset, fitted_model, abundance_response, cross_validation_2) %>% 
                               unique() %>% 
                               mutate(model_present = 1)) # mutate(model_present = 1))
 
# estimate proportions
model_summaries <- full_model_join %>% 
  ungroup() %>% 
  group_by(dataset, cross_validation_2, plot_level, fitted_model, abundance_response) %>% 
  do(models_fit_proportion = sum(.$model_present, na.rm = T) / length(unique(.$species_name))) %>% 
  unnest()

model_summaries$plot_level <- gsub('_', '-', gsub('gam.|rf.|glm.|gbm.', '', model_summaries$plot_level))

dir.create('figures-R2/model-performance-figures/model_counts', recursive = T)
pdf('figures-R2/model-performance-figures/model_counts/model_counts.pdf', width = 20, height = 8)
ggplot(model_summaries) + 
  geom_bar(aes(x = plot_level, y = 1-models_fit_proportion, fill = fitted_model, alpha = models_fit_proportion), 
           stat = 'identity', col = 'black') +
  facet_grid(dataset ~ cross_validation_2 + fitted_model, scales = 'free_x') + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5), 
        aspect.ratio = 1) + 
  scale_fill_manual(values = colours) + 
  ylim(c(0,1)) + 
  ylab('proportion of species with unsuccessful models') + 
  xlab('model type') + 
  scale_alpha_continuous(range = c(1, 0.2))
dev.off()

# estimate median values for accuracy, discrimination and precision for best and worst overall framework ----

framework_arrange <- all_assessments %>% 
  group_by(plot_level, cross_validation) %>% 
  do(Amae_median = median(.$Amae_rel_mean, na.rm = T), 
     Dspearman_median = median(.$Dspearman, na.rm = T), 
     Pdispersion_median = median(.$Pdispersion, na.rm = T)) %>% 
  unnest() %>% 
  split(., .$cross_validation)

framework_arrange_selection <- bind_rows(lapply(framework_arrange, function(x){
  mix_max_framework <- pivot_longer(x, Amae_median:Dspearman_median) %>% 
    group_by(name) %>% 
    do(worst_framework = .$plot_level[which.min(.$value)], 
       worst_value     = signif(.$value[which.min(.$value)],2),
       best_framework = .$plot_level[which.max(.$value)], 
       best_value     = signif(.$value[which.max(.$value)], 2)) %>% 
    unnest()
  mix_max_framework_2 <- mix_max_framework
  mix_max_framework_2[1, c(4,5)] <- mix_max_framework[1, c(2,3)]
  mix_max_framework_2[1, c(2,3)] <- mix_max_framework[1, c(4,5)]
  mix_max_framework_2
  }
  ), .id = 'dataset')

framework_arrange_dataset <- all_assessments %>% 
  mutate(cross_validation = sub('.*\\_', '', .$cross_validation)) %>% 
  group_by(cross_validation, fitted_model) %>% 
  do(Amae_median = median(.$Amae_rel_mean, na.rm = T), 
     Dspearman_median = median(.$Dspearman, na.rm = T)) %>% 
  unnest() %>% 
  split(., .$cross_validation)

framework_arrange_selection_dataset <- bind_rows(lapply(framework_arrange_dataset, function(x){
  mix_max_framework <- pivot_longer(x, Amae_median:Dspearman_median) %>% 
    group_by(name) %>% 
    do(worst_framework = .$fitted_model[which.min(.$value)], 
       worst_value     = signif(.$value[which.min(.$value)],2),
       best_framework = .$fitted_model[which.max(.$value)], 
       best_value     = signif(.$value[which.max(.$value)],2)) %>% 
    unnest()
  mix_max_framework_2 <- mix_max_framework
  mix_max_framework_2[1, c(4,5)] <- mix_max_framework[1, c(2,3)]
  mix_max_framework_2[1, c(2,3)] <- mix_max_framework[1, c(4,5)]
  mix_max_framework_2
}
), .id = 'dataset')


# save xlsx for summaries
writexl::write_xlsx(list(framework_arrange_selection, framework_arrange_selection_dataset), 
                    path = 'figures-R2/model-performance-figures/best_worst_model_comparisons.xlsx')

# what percentage of model frameworks produce 'un-acceptable' predictions ----

all_assessments %>% 
  mutate(cross_validation = sub('.*\\_', '', .$cross_validation)) %>% 
  group_by(cross_validation) %>% 
  do(n_models = length(unique(paste0(.$plot_level, .$species_name))),
     n_Dspearman_0.75 = sum(.$Dspearman > 0.75, na.rm = T), 
     n_Dspearman_0.5 = sum(.$Dspearman > 0.5, na.rm = T)) %>% 
  mutate(perc_Dspearman_0.75 = (n_Dspearman_0.75/n_models)*100, 
         perc_Dspearman_0.5 = (n_Dspearman_0.5/n_models)*100, ) %>% 
  unnest()


# plot of all metrics for each individual model type ----

nested_assessments <- all_assessments %>% 
  group_by(dataset, cross_validation_2) %>% 
  nest()

# plot_data <- nested_assessments$data[[3]] # for testing functions

# apply across both validation types
lapply(1:nrow(nested_assessments), 
       function(x){
         all_model_plots_v2(plot_data = nested_assessments$data[[x]], 
                            outlier_quantile = 0.1,
                            metrics = metrics,
                            targets = targets,
                            directory = paste0('figures-R2/model-performance-figures/all_model_all_metric_combined/'), 
                            name =  paste0(nested_assessments$dataset[x], '_', nested_assessments$cross_validation_2[x]), 
                            height = 10, 
                            width = 14)})


# plot aggregated models within metrics types ----

nested_assessments <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

# plot_data = nested_assessments$data[[1]]
lapply(1:nrow(nested_assessments),
       function(x){
         combined_assessment_metrics(plot_data = nested_assessments$data[[x]], 
                                     response = 'all',
                                     levels = levels, 
                                     metrics = metrics,
                                     targets = targets,
                                     colours = colours,
                                     directory  = 'figures-R2/model-performance-figures/all_model_boxplots/', 
                                     name = nested_assessments$cross_validation_2[x], 
                                     width = 8, 
                                     height = 10)})

# values of aggregated models within metrics types
options(scipen=999)

# make summary tables across metrics of values included in the point-line plots
summary_table <- all_assessments %>% 
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(gathered_data = purrr::map(data, 
                                    # gather by metrics and create aggregations across species and metrics
                                    ~gather(., key = metric, value = value,  metrics) %>% 
                                      mutate(value = ifelse(is.infinite(value), NA, value)) %>%
                                      group_by(fitted_model,
                                               abundance_response,
                                               metric) %>% 
                                      do(metric_median = round(signif(median(.$value, na.rm = T), 2), 2), 
                                         metric_sd    = signif(quantile(.$value, 0.75, na.rm = T) - quantile(.$value, 0.25, na.rm = T), 2)) %>% 
                                      unnest(c(metric_median, metric_sd)) %>% 
                                      mutate(metric_sd = ifelse(.$metric_sd < 0.00001, 0, .$metric_sd)))) %>%
  unnest(c(gathered_data)) %>% 
  mutate(metric_final = paste0(metric_median, ' (±', metric_sd, ')')) %>% 
  select(-metric_median, -metric_sd) %>% 
  dplyr::select(-data) %>% 
  pivot_wider(., names_from = metric, values_from = c(metric_final)) %>% 
  data.frame %>% 
  rename(., data = dataset, cv = cross_validation_2, response = abundance_response, model = fitted_model)

summary_table_full <- all_assessments %>% 
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(gathered_data = purrr::map(data, 
                                    # gather by metrics and create aggregations across species and metrics
                                    ~gather(., key = metric, value = value,  metrics) %>% 
                                      mutate(value = ifelse(is.infinite(value), NA, value)) %>% 
                                      group_by(fitted_model,
                                               abundance_response,
                                               plot_level, 
                                               metric) %>% 
                                      do(metric_median = round(signif(median(.$value, na.rm = T), 2), 2), 
                                         metric_sd    = signif(quantile(.$value, 0.75, na.rm = T) - quantile(.$value, 0.25, na.rm = T), 2)) %>% 
                                      unnest(c(metric_median, metric_sd)) %>% 
                                    mutate(metric_sd = ifelse(.$metric_sd < 0.00001, 0, .$metric_sd)))) %>% 
  unnest(c(gathered_data)) %>% 
  mutate(metric_final = paste0(metric_median, ' (±', metric_sd, ')')) %>% 
  select(-metric_median, -metric_sd) %>% 
  dplyr::select(-data) %>% 
  pivot_wider(., names_from = metric, values_from = c(metric_final)) %>% 
  data.frame %>% 
  rename(., data = dataset, cv = cross_validation_2, response = abundance_response, model = fitted_model)

# write to spreadsheet for reference
writexl::write_xlsx(list(full = summary_table_full, agg = summary_table), path = 'figures-R2/model-performance-figures/summary_table.xlsx')


# plots of rescaled values comparing between models ----

all_assessments_relative <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

# checking aggregate metrics function
plot_data = all_assessments_relative$data[[1]] %>% 
  # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
  group_by(dataset, species_name) %>% 
  nest() %>% 
  .$data %>% 
  .[[1]]

# checking plotting function
plot_data = all_assessments_relative$data[[1]] %>% 
  # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
  group_by(dataset, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                  metrics = metrics))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  mutate(dataset = all_assessments_relative$data[[1]]$dataset, 
         species_name = all_assessments_relative$data[[1]]$species_name)

# aggreagte basic models
plot_all_aggregated(all_assessments_relative$data[[1]] %>% 
                      # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                      group_by(dataset, species_name) %>% 
                      nest() %>% 
                      mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                      metrics = metrics))) %>% 
                      .$metric_aggregation %>% 
                      do.call(rbind, .) %>% 
                      mutate(dataset = all_assessments_relative$data[[1]]$dataset, 
                             species_name = all_assessments_relative$data[[1]]$species_name), 
                    directory = 'figures-R2/model-performance-figures/all_model_rescaled', 
                    colours = colours,
                    name = 'basic', 
                    levels = c('glm', 'gam', 'gbm', 'rf'))

# aggregate oob_cv models
plot_all_aggregated(all_assessments_relative$data[[2]] %>% 
                      # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                      group_by(dataset, species_name) %>% 
                      nest() %>% 
                      mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                      metrics = metrics))) %>% 
                      .$metric_aggregation %>% 
                      do.call(rbind, .) %>% 
                      mutate(dataset = all_assessments_relative$data[[2]]$dataset, 
                             species_name = all_assessments_relative$data[[2]]$species_name), 
                    directory = 'figures-R2/model-performance-figures/all_model_rescaled', 
                    colours = colours,
                    name = 'cv', 
                    levels = c('glm', 'gam', 'gbm', 'rf'))


#  plots of rescaled values comparing between models: model type only ----


# checking plotting function
plot_data = all_assessments_relative$data[[1]] %>% 
  # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
  group_by(dataset, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                  metrics = metrics))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  mutate(dataset            = all_assessments_relative$data[[1]]$dataset, 
         species_name       = all_assessments_relative$data[[1]]$species_name, 
         abundance_response = all_assessments_relative$data[[1]]$abundance_response)


plot_all_aggregated_by_model(all_assessments_relative$data[[1]] %>% 
                              # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                              group_by(dataset, species_name) %>% 
                              nest() %>% 
                              mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                              metrics = metrics))) %>% 
                              .$metric_aggregation %>% 
                              do.call(rbind, .) %>% 
                              mutate(dataset = all_assessments_relative$data[[1]]$dataset, 
                                     species_name = all_assessments_relative$data[[1]]$species_name), 
                            directory = 'figures-R2/model-performance-figures/all_model_rescaled_by_models', 
                            colours = colours,
                            name = 'basic', 
                            levels = c('glm', 'gam', 'gbm', 'rf'))


plot_all_aggregated_by_model(all_assessments_relative$data[[2]] %>% 
                               # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                               group_by(dataset, species_name) %>% 
                               nest() %>% 
                               mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                               metrics = metrics))) %>% 
                               .$metric_aggregation %>% 
                               do.call(rbind, .) %>% 
                               mutate(dataset = all_assessments_relative$data[[2]]$dataset, 
                                      species_name = all_assessments_relative$data[[2]]$species_name), 
                             directory = 'figures-R2/model-performance-figures/all_model_rescaled_by_models', 
                             colours = colours,
                             name = 'cv', 
                             levels = c('glm', 'gam', 'gbm', 'rf'))


#  plots of rescaled values comparing between models: abundance response only ----

# checking plotting function
plot_data = all_assessments_relative$data[[1]] %>% 
  # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
  group_by(dataset, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                  metrics = metrics))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  mutate(dataset            = all_assessments_relative$data[[1]]$dataset, 
         species_name       = all_assessments_relative$data[[1]]$species_name, 
         abundance_response = all_assessments_relative$data[[1]]$abundance_response)


plot_all_aggregated_by_model(all_assessments_relative$data[[1]] %>% 
                               # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                               group_by(dataset, species_name) %>% 
                               nest() %>% 
                               mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                               metrics = metrics))) %>% 
                               .$metric_aggregation %>% 
                               do.call(rbind, .) %>% 
                               mutate(dataset = all_assessments_relative$data[[1]]$dataset, 
                                      species_name = all_assessments_relative$data[[1]]$species_name), 
                             directory = 'figures-R2/model-performance-figures/all_model_rescaled_by_abundance', 
                             colours = colours,
                             name = 'basic', 
                             levels = c('glm', 'gam', 'gbm', 'rf'))


plot_all_aggregated_by_model(all_assessments_relative$data[[2]] %>% 
                               # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                               group_by(dataset, species_name) %>% 
                               nest() %>% 
                               mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                               metrics = metrics))) %>% 
                               .$metric_aggregation %>% 
                               do.call(rbind, .) %>% 
                               mutate(dataset = all_assessments_relative$data[[2]]$dataset, 
                                      species_name = all_assessments_relative$data[[2]]$species_name), 
                             directory = 'figures-R2/model-performance-figures/all_model_rescaled_by_abundance', 
                             colours = colours,
                             name = 'cv', 
                             levels = c('glm', 'gam', 'gbm', 'rf'))



# plots of ranked models ----

nested_assessments <- split(all_assessments, f = all_assessments$cross_validation)

lapply(1:length(nested_assessments), 
       function(x){
         rank_plots(plot_data = nested_assessments[[x]], 
                    levels = levels, 
                    metrics = metrics,
                    targets = targets,
                    colours = colours,
                    directory  = 'figures-R2/model-performance-figures/rank_plots/', 
                    name = unique(paste0(nested_assessments[[x]]$dataset, '_', nested_assessments[[x]]$cross_validation_2)), 
                    width = 12, 
                    height = 10)})


# distribution of the values in models selected as 'best' ----

# find the best model for a species
best_models <- all_assessments %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae_rel_mean', 'Dintercept', 'Dslope', 
                                                                              'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
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

# produce relative rankings in the assessment metrics
all_assessments_relative <- all_assessments %>% 
  group_by(dataset, cross_validation_2, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, 
                                         ~aggregate_metrics(., 
                                                            metrics = c('Amae_rel_mean', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data) %>% 
  ungroup()

# filter assessment metrics by best model and species combinations
best_model_assessments <- left_join(best_models , 
                                    all_assessments_relative %>% mutate(cross_validation = .$cross_validation_2))

# summaries for optimal models
best_model_assessments %>% filter(cross_validation == 'basic') %>% .$fitted_model %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic')) *100
# gam      gbm      glm       rf 
# 16.75291 26.19664 16.10608 40.94437 

best_model_assessments %>% filter(cross_validation == 'basic') %>% .$abundance_response %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic'))
# abun       abun-occ   abun-occ-2stage 
# 0.2238034  0.5937904  0.1824062 

best_model_assessments %>% filter(cross_validation == 'cv') %>% .$fitted_model %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic')) *100
# gam      gbm      glm       rf 
# 19.79301 21.41009 31.50065 26.84347 

best_model_assessments %>% filter(cross_validation == 'cv') %>% .$abundance_response %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic'))
# abun        abun-occ   abun-occ-2stage 
# 0.1979301   0.5769728  0.2205692 


# create barplots of best models for discrimination ----

best_model_barplot_function(all_assessments, 
                            performance_metric = 'accuracy',
                            save_dir = 'figures-R2/model-performance-figures/best-model-bar-plots/best-model-barplots-accuracy.png') 

best_model_barplot_function(all_assessments, 
                            performance_metric = 'discrimination',
                            save_dir = 'figures-R2/model-performance-figures/best-model-bar-plots/best-model-barplots-discrimination.png') 

best_model_barplot_function(all_assessments, 
                            performance_metric = 'precision',
                            save_dir = 'figures-R2/model-performance-figures/best-model-bar-plots/best-model-barplots-precision.png') 

best_model_barplot_function(all_assessments, 
                            performance_metric = 'all',
                            save_dir = 'figures-R2/model-performance-figures/best-model-bar-plots/best-model-barplots-all.png') 

# produce histograms of model performance for best models ----

# input data to function
library(viridis)
best_model_assessments %>% 
  filter(cross_validation == 'basic') %>% 
  spp_best_assessment_metrics(., 
                              metrics = metrics, 
                              targets = targets, 
                              directory = 'figures-R2/model-performance-figures/best-model-histograms', 
                              name = 'basic', 
                              width = 6, 
                              height = 8, 
                              colours = c(viridis(10, option = 3)[5], viridis(10, option = 7)[5]))

best_model_assessments %>% 
  filter(cross_validation == 'cv') %>% 
  spp_best_assessment_metrics(., 
                              metrics = metrics, 
                              targets = targets, 
                              directory = 'figures-R2/model-performance-figures/best-model-histograms', 
                              name = 'cv', 
                              width = 6, height = 8, 
                              colours = c(viridis(10, option = 3)[5], viridis(10, option = 7)[5]))


# updated verion of code

# test the new function
spp_best_assessment_metrics_V2

performance_data = rbind(
  
  best_model_assessments %>% 
  filter(cross_validation == 'basic') %>% 
  pivot_longer(., cols = Amae_rel_mean:Pdispersion) %>% 
  group_by(dataset, cross_validation) %>% 
  nest() %>% 
  mutate(outlier_lwr = purrr::map(data, ~ boxplot.stats(.$value)$stats[1]), 
         outlier_upr = purrr::map(data, ~ boxplot.stats(.$value)$stats[5])) %>% 
  unnest(c(outlier_lwr, outlier_upr)) %>% unnest(data) %>% 
  filter(value < outlier_upr, value > outlier_lwr), 
  
  best_model_assessments %>% 
    filter(cross_validation == 'cv') %>% 
    pivot_longer(., cols = Amae_rel_mean:Pdispersion) %>% 
    group_by(dataset, cross_validation) %>% 
    nest() %>% 
    mutate(outlier_lwr = purrr::map(data, ~ boxplot.stats(.$value)$stats[1]), 
           outlier_upr = purrr::map(data, ~ boxplot.stats(.$value)$stats[5])) %>% 
    unnest(c(outlier_lwr, outlier_upr)) %>% unnest(data) %>% 
    filter(value < outlier_upr, value > outlier_lwr)
  
)

#[1] "Amae_rel_mean" "Dintercept"    "Dpearson"      "Dslope"        "Dspearman"     "Pdispersion"  

performance_data$name2 <- factor(performance_data$name, labels = c('Accuracy: \n proportional MAE',
                                                     'Discrimination: \n regression slope', 
                                                     'Discrimination: \n regression intercept', 
                                                    'Discrimination: \n pearson correlation', 
                                                    'Discrimination: \n spearmans rank', 
                                                    'Precision: \n dispersion'), 
                          levels = c("Amae_rel_mean","Dslope","Dintercept","Dpearson","Dspearman","Pdispersion"))

yintercept <- ifelse(performance_data$name == 'Dintercept', 0, 1)
alpha <- ifelse(performance_data$cross_validation_2 == 'cv', 0.9, 1)


dir.create('figures-R2/model-performance-figures/best-model-performance-values/', recursive = T)
png(filename = 'figures-R2/model-performance-figures/best-model-performance-values/performance-values.png', res = 300,
    width = 1500, height = 1000)

ggplot(data = performance_data) + 
  geom_boxplot(aes(y = value, x = dataset, fill = dataset, alpha = cross_validation_2, 
                   group = paste0(cross_validation_2, dataset)), outlier.shape = NA) + 
  geom_hline(aes(yintercept = yintercept), lty = 2) + 
  facet_wrap(~name2) + 
  xlab(NULL) +
  ylab('peformance metric value') +
  theme_classic() + 
  theme(aspect.ratio = 0.75, 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = 'bottom') +
  scale_fill_manual(values = c(viridis(10, option = 3)[5], viridis(10, option = 7)[5]),
                    name = '') +
  scale_alpha_manual(values = c(1, 0.5)) + 
  guides(alpha = FALSE)

dev.off()

# summary table of best models across datasets ----

basic_best_summary <- best_model_assessments %>% 
  filter(cross_validation == 'basic') %>% 
  group_by(dataset) %>% 
  select(dataset, Amae_rel_mean:Pdispersion) %>% 
  pivot_longer(., Amae_rel_mean:Pdispersion) %>% 
  filter(name %in% metrics) %>% 
  group_by(dataset, name) %>% 
  do(median = median(.$value, na.rm = T), 
     IQR.25 = quantile(.$value, 0.25, na.rm = T), 
     IQR.75 = quantile(.$value, 0.75, na.rm = T)) %>% 
  unnest()

cv_best_summary <- best_model_assessments %>% 
  filter(cross_validation == 'cv', ) %>%
  group_by(dataset) %>% 
  select(dataset, Amae_rel_mean:Pdispersion) %>% 
  pivot_longer(., Amae_rel_mean:Pdispersion) %>% 
  filter(name %in% metrics) %>% 
  group_by(dataset, name) %>% 
  do(median = median(.$value, na.rm = T), 
     IQR.25 = quantile(.$value, 0.25, na.rm = T), 
     IQR.75 = quantile(.$value, 0.75, na.rm = T)) %>% 
  unnest()

writexl::write_xlsx(list(basic_best_summary, cv_best_summary), 'figures-R2/model-performance-figures/summary_table_aggregated.xlsx')

# plots of rescaled values comparing between models including 'best' ----

# find the best model for a species
best_models <- all_assessments %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae_rel_mean', 'Dintercept', 'Dslope', 
                                                                              'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  na.omit(.) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$aggregated_evaluation_metrics)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model) %>% 
  mutate(dataset = plyr::revalue(.$dataset, c(bbs_cv = 'bbs', rls_cv = 'rls')))

# filter assessment metrics by best model and species combinations
best_model_assessments <- left_join(best_models , 
                                    all_assessments %>% mutate(cross_validation = .$cross_validation_2))

# rbind together the best models with the rest of the model assessments
best_model_assessments$fitted_model = 'optimal'
all_assessments_best <- bind_rows(best_model_assessments, all_assessments)


# produce relative rankings in the assessment metrics
all_assessments_best_relative <- all_assessments_best %>% 
  group_by(dataset, cross_validation_2, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, 
                                         ~aggregate_metrics(., 
                                                            metrics = c('Amae_rel_mean', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data) %>% 
  ungroup()

# group for plotting
all_assessments_best_relative_plot <- all_assessments_best_relative %>% 
  group_by(cross_validation_2) %>% 
  nest()

plot_data = all_assessments_relative$data[[1]]

# aggreagte basic models
plot_all_aggregated(all_assessments_best_relative_plot$data[[1]], 
                    directory = 'figures-R2/model-performance-figures/all_model_rescaled_with_optimal', 
                    colours = c(colours, 'black'),
                    name = 'basic', 
                    levels = c('glm', 'gam', 'gbm', 'rf', 'optimal'))

# aggregate oob_cv models
plot_all_aggregated(all_assessments_best_relative_plot$data[[2]], 
                    directory = 'figures-R2/model-performance-figures/all_model_rescaled_with_optimal', 
                    colours = c(colours, 'black'),
                    name = 'cv', 
                    levels = c('glm', 'gam', 'gbm', 'rf', 'optimal'))



# correlation plots amongst best models ----

# run the previous best model code

# plots
best_model_assessments %>% 
  filter(cross_validation == 'basic', dataset == 'rls') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures-R2/model-performance-figures/correlation-metrics', 
                    name = 'rls_basic', 
                    width = 10, height = 5)

# plots
best_model_assessments %>% 
  filter(cross_validation == 'cv', dataset == 'rls') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures-R2/model-performance-figures/correlation-metrics', 
                    name = 'rls_cv', 
                    width = 10, height = 5)

# plots
best_model_assessments %>% 
  filter(cross_validation == 'basic', dataset == 'bbs') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures-R2/model-performance-figures/correlation-metrics', 
                    name = 'bbs_basic', 
                    width = 10, height = 5)
# plots
best_model_assessments %>% 
  filter(cross_validation == 'cv', dataset == 'bbs') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures-R2/model-performance-figures/correlation-metrics', 
                    name = 'bbs_cv', 
                    width = 10, height = 5)
