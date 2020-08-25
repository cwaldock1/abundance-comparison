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

levels = c('glm', 'gam', 'gbm', 'rf')

metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')

targets = list(#Armse    = c(0, -20, 20), 
               Amae     = c(0, -20, 20), 
               Dintercept  = c(0, -5, 20), 
               Dslope      = c(1,  0, 2), 
               Dpearson    = c(1,  0, 1), 
               Dspearman   = c(1,  0, 1), 
               #Psd         = c(0,  0, 100), 
               Pdispersion = c(1,  0, 20)#, 
               #Pr2         = c(1,  0, 1)
               )

# load in evaluation data ----

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

# create plots of % of species per model type ----

# summarise models and species

length(unique(all_assessments$plot_level))
length(unique(all_assessments$species_name))

# create a full set of models and species

full_model_set <- all_assessments %>% 
  group_by(dataset, fitted_model, abundance_response) %>% 
  do(expand.grid(plot_level = .$plot_level %>% unique(), species_name = .$species_name %>% unique()))
  
full_model_set <- rbind(full_model_set %>% mutate(cross_validation_2 = 'cv'), full_model_set %>% mutate(cross_validation_2 = 'basic'))

full_model_join <- full_join(full_model_set, all_assessments %>% 
                               select(plot_level, species_name, dataset, fitted_model, abundance_response, cross_validation_2) %>% 
                               unique() %>% 
                               mutate(model_present = 1))

# estimate proportions
model_summaries <- full_model_join %>% 
  ungroup() %>% 
  group_by(dataset, cross_validation_2, plot_level, fitted_model, abundance_response) %>% 
  do(models_fit_proportion = sum(.$model_present, na.rm = T) / length(unique(.$species_name))) %>% 
  unnest()

model_summaries$plot_level <- gsub('_', '-', gsub('gam.|rf.|glm.|gbm.', '', model_summaries$plot_level))

dir.create('figures/model-performance-figures/model_counts')
pdf('figures/model-performance-figures/model_counts/model_counts.pdf', width = 20, height = 8)
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
                            directory = paste0('figures/model-performance-figures/all_model_all_metric_combined/'), 
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
                                     directory  = 'figures/model-performance-figures/all_model_boxplots/', 
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
writexl::write_xlsx(list(full = summary_table_full, agg = summary_table), path = 'figures/model-performance-figures/summary_table.xlsx')


# plots of rescaled values comparing between models ----

all_assessments_relative <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

plot_data = all_assessments_relative$data[[1]]

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
                    directory = 'figures/model-performance-figures/all_model_rescaled', 
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
                    directory = 'figures/model-performance-figures/all_model_rescaled', 
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
                    directory  = 'figures/model-performance-figures/rank_plots/', 
                    name = unique(paste0(nested_assessments[[x]]$dataset, '_', nested_assessments[[x]]$cross_validation_2)), 
                    width = 12, 
                    height = 10)})


# distribution of the values in models selected as 'best' ----

# find the best model for a species
best_models <- all_assessments %>% 
  select(-Armse, -Psd) %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 
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
                                                            metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data) %>% 
  ungroup()

# filter assessment metrics by best model and species combinations
best_model_assessments <- left_join(best_models , 
                                    all_assessments_relative %>% mutate(cross_validation = .$cross_validation_2))

# summaries for optimal models
best_model_assessments %>% filter(cross_validation == 'basic') %>% .$fitted_model %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic')) *100
#gam        gbm        glm         rf 
#1.5161503  8.2399473  0.9228741 89.3210283 

best_model_assessments %>% filter(cross_validation == 'basic') %>% .$abundance_response %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic'))
#abun        abun-occ abun-occ-2stage 
#0.2175346       0.6222808       0.1601846 

best_model_assessments %>% filter(cross_validation == 'cv') %>% .$fitted_model %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic')) *100
#gam      gbm      glm       rf 
#17.93013 20.17139 30.65260 30.58668 

best_model_assessments %>% filter(cross_validation == 'cv') %>% .$abundance_response %>% table / nrow(best_model_assessments %>% filter(cross_validation == 'basic'))
#abun        abun-occ abun-occ-2stage 
#0.2043507       0.5418589       0.2471984 

# input data to function
library(viridis)
best_model_assessments %>% 
  filter(cross_validation == 'basic') %>% 
  spp_best_assessment_metrics(., 
                              metrics = metrics, 
                              targets = targets, 
                              directory = 'figures/model-performance-figures/best-model-histograms', 
                              name = 'basic', 
                              width = 6, 
                              height = 8, 
                              colours = c(viridis(10, option = 3)[5], viridis(10, option = 7)[5]))

best_model_assessments %>% 
  filter(cross_validation == 'cv') %>% 
  spp_best_assessment_metrics(., 
                              metrics = metrics, 
                              targets = targets, 
                              directory = 'figures/model-performance-figures/best-model-histograms', 
                              name = 'cv', 
                              width = 6, height = 8, 
                              colours = c(viridis(10, option = 3)[5], viridis(10, option = 7)[5]))

# summary table of best models across datasets ----

basic_best_summary <- best_model_assessments %>% 
  filter(cross_validation == 'basic') %>% 
  group_by(dataset) %>% 
  select(dataset, Armse:Pr2) %>% 
  pivot_longer(., Armse:Pr2) %>% 
  filter(name %in% metrics) %>% 
  group_by(dataset, name) %>% 
  do(median = median(.$value, na.rm = T), 
     IQR.25 = quantile(.$value, 0.25, na.rm = T), 
     IQR.75 = quantile(.$value, 0.75, na.rm = T)) %>% 
  unnest()

cv_best_summary <- best_model_assessments %>% 
  filter(cross_validation == 'cv', ) %>%
  group_by(dataset) %>% 
  select(dataset, Armse:Pr2) %>% 
  pivot_longer(., Armse:Pr2) %>% 
  filter(name %in% metrics) %>% 
  group_by(dataset, name) %>% 
  do(median = median(.$value, na.rm = T), 
     IQR.25 = quantile(.$value, 0.25, na.rm = T), 
     IQR.75 = quantile(.$value, 0.75, na.rm = T)) %>% 
  unnest()

writexl::write_xlsx(list(basic_best_summary, cv_best_summary), 'figures/model-performance-figures/summary_table_aggregated.xlsx')

# plots of rescaled values comparing between models including 'best' ----

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
                                                            metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
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
                    directory = 'figures/model-performance-figures/all_model_rescaled_with_optimal', 
                    colours = c(colours, 'black'),
                    name = 'basic', 
                    levels = c('glm', 'gam', 'gbm', 'rf', 'optimal'))

# aggregate oob_cv models
plot_all_aggregated(all_assessments_best_relative_plot$data[[2]], 
                    directory = 'figures/model-performance-figures/all_model_rescaled_with_optimal', 
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
                    directory = 'figures/model-performance-figures/correlation-metrics', 
                    name = 'rls_basic', 
                    width = 10, height = 5)

# plots
best_model_assessments %>% 
  filter(cross_validation == 'cv', dataset == 'rls') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures/model-performance-figures/correlation-metrics', 
                    name = 'rls_cv', 
                    width = 10, height = 5)

# plots
best_model_assessments %>% 
  filter(cross_validation == 'basic', dataset == 'bbs') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures/model-performance-figures/correlation-metrics', 
                    name = 'bbs_basic', 
                    width = 10, height = 5)
# plots
best_model_assessments %>% 
  filter(cross_validation == 'cv', dataset == 'bbs') %>% 
  correlation_plots(plot_data = ., 
                    metrics = metrics, 
                    targets = targets, 
                    directory = 'figures/model-performance-figures/correlation-metrics', 
                    name = 'bbs_cv', 
                    width = 10, height = 5)


