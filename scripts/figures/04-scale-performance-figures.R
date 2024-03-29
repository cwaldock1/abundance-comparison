# script to produce figures that evaluate model performance across different modelling proceedures

# load packages ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'grid', 'data.table')
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
               Pdispersion = c(1,  0, 20), 
               Pr2         = c(1,  0, 1))

# load in evaluation data ----

# load in assessments across all scales
all_assessments <- lapply(list.files('results/model_assessment_scale', full.names = T, recursive = T), readRDS)

# apply names to the data
names(all_assessments) <- str_split(list.files('results/model_assessment_scale', full.names = T, recursive = T), '/', simplify = T)[,3]

# turn into data.frame
all_assessments <- bind_rows(all_assessments, .id = 'spatial_scale')

# select only the columns to be used later 
all_assessments <- all_assessments %>% select(-family_grouped_simple, -family_grouped, -family, 
                           -transformation, -n_absence, -n_boot_absence, -abundance_response_simple) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  select(spatial_scale, dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name,
         Armse, Amae, Dintercept, Dslope, Dpearson, Dspearman, Psd, Pdispersion, Pr2, Evaluation_number, Evaluation_message) %>% 
  
  # change abundance_response
  mutate(abundance_response = plyr::revalue(.$abundance_response, c(abunocc = "abun-occ", abunocc_2stage = "abun-occ-2stage")))

# plots of rescaled values comparing between models ----

all_assessments_relative <- all_assessments %>% 
  group_by(cross_validation_2) %>% 
  nest()

# for testing
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
plot_all_aggregated_spatial_scale(all_assessments_relative$data[[1]] %>% 
                      # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                      group_by(dataset, species_name) %>% 
                      nest() %>% 
                      mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                      metrics = metrics))) %>% 
                      .$metric_aggregation %>% 
                      do.call(rbind, .) %>% 
                      mutate(dataset = all_assessments_relative$data[[1]]$dataset, 
                             species_name = all_assessments_relative$data[[1]]$species_name), 
                    directory = 'figures/scale-performance-figures/all_model_rescaled', 
                    name = 'basic', 
                    levels = c('glm', 'gam', 'gbm', 'rf'))

# aggregate oob_cv models
plot_all_aggregated_spatial_scale(all_assessments_relative$data[[2]] %>% 
                      # the species name grouping here is essential otherwise the rankings occur across all species and model combinations
                      group_by(dataset, species_name) %>% 
                      nest() %>% 
                      mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(.,
                                                                                      metrics = metrics))) %>% 
                      .$metric_aggregation %>% 
                      do.call(rbind, .) %>% 
                      mutate(dataset = all_assessments_relative$data[[2]]$dataset, 
                             species_name = all_assessments_relative$data[[2]]$species_name), 
                    directory = 'figures/scale-performance-figures/all_model_rescaled', 
                    name = 'cv', 
                    levels = c('glm', 'gam', 'gbm', 'rf'))


# distribution of the values in models selected as 'best' ----

# find the best model for a species at the finest spatial scale
best_scale_0.1 <- all_assessments %>% 
 filter(spatial_scale == 'spatial_scale_0.1') %>% 
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
 group_by(species_name, spatial_scale, cross_validation) %>% 
 do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
 unnest(cols = c('best_model')) %>% 
 select(-spatial_scale) 

# pre-able to filter the best models at a 0.1 scale from the total assessments of scale
all_best_scale_0.1 <- data.table(best_scale_0.1 %>% mutate(plot_level = best_model) %>% select(-best_model))

# read in overall best model object
overall_best_models <- readRDS('results/overall_best_models.RDS')
overall_best_models <- data.table(overall_best_models)
all_assessments_dt <- data.table(all_assessments)

# join using data.table
best_models_scale <- left_join(all_best_scale_0.1, all_assessments_dt)

# edit names for labels 
best_models_scale <- best_models_scale %>% mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
       cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
       plot_level = .$best_model) %>% 
  mutate(dataset = plyr::revalue(.$dataset, c(bbs_cv = 'bbs', rls_cv = 'rls')))

# save and read in best_models_scale object
saveRDS(best_models_scale, file = 'results/scale_best_models.RDS')
best_models_scale <- readRDS(file = 'results/scale_best_models.RDS')

# input data to function

plot_data = best_models_scale # for debugging function

# conver to data.frame
best_models_scale <- data.frame(best_models_scale)

best_models_scale %>% 
  spp_best_assessment_metrics_scale(plot_data = ., 
                                    metrics = metrics,
                                    targets = targets, 
                                    spatial_filter = 11,
                                    outlier_quantile = 0,
                                    directory = 'figures/scale-performance-figures/best-model-barplots', 
                                    name = 'combined_plot', 
                                    width = 8, 
                                    height = 4)

# barplots of the most discriminatory model ----

# find the best model for a species
best_models_each_scale <- all_assessments %>% 
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
  group_by(species_name, spatial_scale, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model) %>% 
  mutate(dataset = plyr::revalue(.$dataset, c(bbs_cv = 'bbs', rls_cv = 'rls')))

# produce relative rankings in the assessment metrics: takes 5 minutes or so
all_assessments_relative <- all_assessments %>% 
  group_by(dataset, cross_validation_2, spatial_scale, species_name) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, 
                                         ~aggregate_metrics(., 
                                                            metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion')))) %>% 
  unnest(metric_aggregation) %>% 
  dplyr::select(-data) %>% 
  ungroup()

# filter assessment metrics by best model and species combinations
best_model_assessments <- left_join(best_models_each_scale , 
                                    all_assessments_relative %>% mutate(cross_validation = .$cross_validation_2))

# saveRDS(best_model_assessments, file = 'results/scale_compiled.RDS')

best_model_assessments <- readRDS(file = 'results/scale_compiled.RDS')

# percentages for summaries in paper
best_model_assessments %>% filter(cross_validation == 'basic',spatial_scale == 'spatial_scale_10') %>% .$fitted_model %>% table /
  nrow(best_model_assessments %>% filter(cross_validation == 'basic',spatial_scale == 'spatial_scale_10')) *100

# gam       gbm       glm        rf 
# 5.000000  8.904110  3.561644 82.534247 

best_model_assessments %>% filter(cross_validation == 'cv',spatial_scale == 'spatial_scale_10') %>% .$fitted_model %>% table /
  nrow(best_model_assessments %>% filter(cross_validation == 'cv',spatial_scale == 'spatial_scale_10')) *100

# gam      gbm      glm       rf 
# 18.00450 22.65566 29.33233 30.00750 


# calculate proportions
best_model_props <- best_model_assessments %>% 
  group_by(dataset, cross_validation, spatial_scale) %>% 
  do(count_best = data.frame(table(.$fitted_model))) %>% 
  unnest(count_best) %>% 
  rename(fitted_model = Var1, count = Freq) %>% 
  group_by(dataset, cross_validation, spatial_scale) %>% 
  mutate(proportion_best = count / sum(count))

# edit levels
best_model_props$spatial_scale = gsub('spatial_scale_', '', best_model_props$spatial_scale)
best_model_props <- best_model_props %>% ungroup %>% mutate(spatial_scale = as.numeric(spatial_scale)) %>%  filter(spatial_scale < 11)
best_model_props$spatial_scale = factor(best_model_props$spatial_scale, levels = sort(unique(best_model_props$spatial_scale)))
best_model_props$fitted_model = factor(best_model_props$fitted_model, levels = levels)
best_model_props$cross_validation = factor(best_model_props$cross_validation, labels = c('within-sample', 'out-of-sample'))
best_model_props$dataset = factor(best_model_props$dataset, labels = c('breeding-bird survey', 'reef-life survey'))

# create directory and plots
dir.create('figures/scale-performance-figures/best-model-counts')
pdf('figures/scale-performance-figures/best-model-counts/best_model_counts.pdf', width = 7, height = 7)
grid.arrange(
ggplot(data = best_model_props %>% filter(dataset == 'breeding-bird survey')) + 
  geom_bar(aes(x = fitted_model, y = proportion_best, group = spatial_scale, fill = as.numeric(as.character(spatial_scale))), position = 'dodge', stat = 'identity') + 
  scale_fill_viridis_c(option = 3, end = 0.9, 
                       guide = guide_colourbar(direction = 'horizontal', barheight = 0.5, label.vjust = 0)) + 
  theme_bw() + 
  facet_grid(cross_validation ~ dataset) + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        strip.background = element_blank(), 
        legend.position = 'bottom', 
        legend.title = element_blank(),
        legend.background = element_blank(), 
        strip.text.x = element_blank(), 
        strip.text.y = element_text(size = 12)) + 
  xlab(NULL) + 
  ylab('proportion best fitted model within spatial scale'),

ggplot(data = best_model_props %>% filter(dataset == 'reef-life survey')) + 
  geom_bar(aes(x = fitted_model, y = proportion_best, group = spatial_scale, fill = as.numeric(as.character(spatial_scale))), position = 'dodge', stat = 'identity') + 
  scale_fill_viridis_c(option = 7, end = 0.9, 
                       guide = guide_colourbar(direction = 'horizontal', barheight = 0.5, label.vjust = 0)) + 
  theme_bw() + 
  facet_grid(cross_validation ~ dataset) + 
  theme(panel.grid = element_blank(),
        aspect.ratio = 1, 
        strip.background = element_blank(), 
        legend.position = 'bottom', 
        legend.title = element_blank(),
        legend.background = element_blank(), 
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 12)) + 
  xlab(NULL) + 
  ylab('proportion best fitted model within spatial scale'), 
nrow = 1)
dev.off()


# create table of scale by model performance values for best models ----

best_models_scale$spatial_scale = gsub('spatial_scale_', '', best_models_scale$spatial_scale)
best_models_scale <- best_models_scale %>% ungroup %>% mutate(spatial_scale = as.numeric(spatial_scale)) %>%  filter(spatial_scale < 11)
best_models_scale <- na.omit(best_models_scale)

writexl::write_xlsx(
list(best_models_scale %>% 
  select(-Armse, -Psd) %>% 
  filter(spatial_scale %in% c(0.1, 10)) %>% 
  group_by(cross_validation, spatial_scale) %>% 
  summarise_at(vars(c(Amae:Pdispersion)), list(median, sd)) %>% 
  ungroup() %>% 
  pivot_longer(data = .,cols = Amae_fn1:Pdispersion_fn2, names_to = 'metric', values_to = 'value') %>% 
  mutate(measure = ifelse(grepl('fn1', metric), 'median', 'sd'), 
         metric = gsub('_fn1|_fn2', '', metric)) %>% 
  pivot_wider(., names_from = measure, values_from = value) %>% 
  mutate(spatial_scale = as.numeric(gsub('spatial_scale_', '', spatial_scale))) %>% 
  .[order(.$cross_validation, .$spatial_scale),] %>% 
  pivot_wider(., names_from = c(spatial_scale), values_from = c(median, sd)) %>% 
  ungroup() %>% 
  mutate_at(vars(c(median_0.1:sd_10)), signif, digits = 2),

  best_models_scale %>% 
  select(-Armse, -Psd) %>% 
  filter(spatial_scale %in% c(0.1, 10)) %>% 
  group_by(dataset, cross_validation, spatial_scale) %>% 
  summarise_at(vars(c(Amae:Pdispersion)), list(median, sd)) %>% 
  ungroup() %>% 
  pivot_longer(data = .,cols = Amae_fn1:Pdispersion_fn2, names_to = 'metric', values_to = 'value') %>% 
  mutate(measure = ifelse(grepl('fn1', metric), 'median', 'sd'), 
         metric = gsub('_fn1|_fn2', '', metric)) %>% 
  pivot_wider(., names_from = measure, values_from = value) %>% 
  mutate(spatial_scale = as.numeric(gsub('spatial_scale_', '', spatial_scale))) %>% 
  .[order(.$cross_validation, .$spatial_scale),] %>% 
  pivot_wider(., names_from = c(spatial_scale), values_from = c(median, sd)) %>% 
  ungroup() %>% 
  mutate_at(vars(c(median_0.1:sd_10)), signif, digits = 2)),
path = 'figures/scale-performance-figures/metric_summary_table.xlsx')
  
  
# estimate percentage change between 0.1° and 10° for manuscript text
writexl::write_xlsx(
list(best_models_scale %>% 
       filter(spatial_scale %in% c(0.1, 10)) %>% 
  select(-Armse, -Psd) %>% 
  group_by(dataset, cross_validation, spatial_scale) %>% 
  summarise_at(vars(c(Amae:Pdispersion)), list(mean)) %>% 
  ungroup() %>% 
  pivot_longer(data = .,cols = Amae:Pdispersion, names_to = 'metric', values_to = 'value') %>% 
  pivot_wider(., names_from = spatial_scale, values_from = value) %>% 
  mutate(baseline = `0.1`) %>% 
  mutate(`0.1` = ((`0.1` - baseline) / baseline)*100, 
         `10` = ((`10` - baseline) / baseline)*100) %>% 
  group_by(cross_validation, metric) %>% 
  summarise_at(vars(c(`0.1`:baseline)), list(mean)) %>% 
  select(-baseline) %>% 
  pivot_longer(data = ., cols = `0.1`:`10`, names_to = 'spatial_scale', values_to = 'value') %>% 
  mutate(spatial_scale = as.numeric(gsub('spatial_scale_', '', spatial_scale))) %>% 
  .[order(.$cross_validation, .$spatial_scale),] %>% 
  pivot_wider(., names_from = 'spatial_scale', values_from = 'value') %>% 
  mutate_at(vars(c(`0.1`:`10`)), signif, digits = 2), 
  
  best_models_scale %>% 
    filter(spatial_scale %in% c(0.1, 10)) %>% 
    select(-Armse, -Psd) %>% 
    group_by(dataset, cross_validation, spatial_scale) %>% 
    summarise_at(vars(c(Amae:Pdispersion)), list(mean)) %>% 
    ungroup() %>% 
    pivot_longer(data = .,cols = Amae:Pdispersion, names_to = 'metric', values_to = 'value') %>% 
    pivot_wider(., names_from = spatial_scale, values_from = value) %>% 
    mutate(baseline = `0.1`) %>% 
    mutate(`0.1` = ((`0.1` - baseline) / baseline)*100, 
           `10` = ((`10` - baseline) / baseline)*100) %>% 
    group_by(dataset, cross_validation, metric) %>% 
    summarise_at(vars(c(`0.1`:baseline)), list(mean)) %>% 
    select(-baseline) %>% 
    pivot_longer(data = ., cols = `0.1`:`10`, names_to = 'spatial_scale', values_to = 'value') %>% 
    mutate(spatial_scale = as.numeric(gsub('spatial_scale_', '', spatial_scale))) %>% 
    .[order(.$cross_validation, .$spatial_scale),] %>% 
    pivot_wider(., names_from = 'spatial_scale', values_from = 'value') %>% 
    mutate_at(vars(c(`0.1`:`10`)), signif, digits = 2)),
path = 'figures/scale-performance-figures/metric_percentage_table.xlsx')
















