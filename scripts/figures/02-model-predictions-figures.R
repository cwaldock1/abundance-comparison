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

# theme for histograms
theme_hist <- function(){
  theme(legend.position = 'none',
        aspect.ratio = 0.5, 
        panel.grid = element_blank(),  
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
        strip.background = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line())}

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

# saveRDS(best_models, file = 'results/overall_best_models.RDS')

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



# get the slope and intercept once the model is rescaled and back transformed ----

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

# outputs
all_evals_agg <- list()
summary_evals <- list()
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
           validation_observed_mean,   validation_predict_mean, 
           validation_locations)
  
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
       validation_predict_mean = .$validation_predict_mean[[1]][.$validation_observed_mean[[1]] > 0], 
       validation_locations = .$validation_locations[[1]][.$validation_observed_mean[[1]] > 0,]) %>% 
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
  
  # apply the function abundance_assessment_metrics to a rescaled set of data, and compare to the non-rescaled set of data. 
  evaluations <-
    left_join(
      validation_data %>% 
        rowwise() %>% 
        do(metrics = abundance_assessment_metrics(predictions   = rescale_01(log10(.$validation_predict_mean+1)), 
                                                  observations  = rescale_01(log10(.$validation_observed_mean+1)), 
                                                  locations = .$validation_locations,
                                                  scale = NULL)) %>% 
        unnest(metrics) %>% 
        rename_with(., .fn = ~paste0(., '_rescaled'), .cols = colnames(.)) %>% 
        bind_cols(validation_data[,1:6], .), 
      
      validation_data %>% 
        rowwise() %>% 
        do(metrics = abundance_assessment_metrics(predictions   = .$validation_predict_mean, 
                                                  observations  = .$validation_observed_mean, 
                                                  locations = .$validation_locations,
                                                  scale = NULL)) %>% 
        unnest(metrics) %>% 
        bind_cols(validation_data[,1:6], .)
      )
  
  # get the summaries for all the metrics and save in a list
  summary_evals[[i]] <- evaluations %>% 
    select(dataset, cross_validation, Armse_rescaled:Pr2_rescaled, Armse:Pr2) %>% 
    pivot_longer(., cols = Armse_rescaled:Pr2) %>% 
    group_by(name) %>% 
    do(dataset = unique(.$dataset), 
       cross_validation = unique(.$cross_validation),
       mean_value    = median(.$value, na.rm = T), 
       sd_value       = sd(.$value, na.rm = T)) %>% 
    unnest()
  
  # select columns for plots
  evaluations <- evaluations %>% 
    select(dataset, cross_validation, 
           Dslope_rescaled, Dslope, 
           Dintercept_rescaled, Dintercept)
  
  # aggregate for points and lines
  all_evals_agg[[i]] <- evaluations %>% 
    pivot_longer(., cols = Dslope_rescaled:Dintercept) %>% 
    group_by(name) %>% 
    do(mean_value    = median(.$value, na.rm = T), 
       sd_value       = sd(.$value, na.rm = T)) %>% 
    unnest()

  # create plots for dslope and dintercept
  theme_set(theme_grey())
  require(patchwork)
  dir <- 'figures/model-prediction-figures/validation_best_models_aggregated_intercept_slope'
  dir <- paste0(dir, '/', model_type[[i]])
  dir.create(dir, recursive = T)
  col <- ifelse(unique(evaluations$dataset) == 'bbs', viridis::viridis(10, option = 3)[5], viridis::viridis(10, option = 7)[5])
  
  # combined figure for slope ----
  d_slope_xy <- 
    ggplot(data = evaluations) + 
    geom_point(aes(x = Dslope, y = Dslope_rescaled)) + 
    geom_hline(aes(yintercept = 1)) + 
    geom_vline(aes(xintercept = 1)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20)) + 
    ylab('rescaled slope') + 
    xlab('original slope')
  
  d_slope_righthist <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dslope_rescaled), alpha = 1, fill = 'gray75', col = 'transparent') +
    geom_point(data = all_evals_agg[[i]] %>% filter(name == 'Dslope_rescaled'), 
               aes(x = mean_value, y = 50), col = 'black', size = 4) +
    geom_linerange(data = all_evals_agg[[i]] %>% filter(name == 'Dslope_rescaled'), 
               aes(xmin = mean_value - sd_value, xmax = mean_value + sd_value, y = 50), col = 'black', 
               size = 1.1) +
    theme_bw() +  theme_hist() +
    theme(aspect.ratio = 5, 
          axis.text.x = element_text(angle = 270, vjust = 0.5), 
          axis.text.y = element_blank(), 
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")) + 
    xlab(NULL) + ylab(NULL) + coord_flip() + 
    xlim(pivot_longer(evaluations, Dslope:Dslope_rescaled) %>% .$value %>% range(., na.rm = T))
  
    d_slope_tophist <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dslope), alpha = 1, fill = 'gray75', col = 'transparent') + 
    geom_point(data = all_evals_agg[[i]] %>% filter(name == 'Dslope'), 
               aes(x = mean_value, y = 50), col = 'black', size = 4) +
    geom_linerange(data = all_evals_agg[[i]] %>% filter(name == 'Dslope'), 
                   aes(xmin = mean_value - sd_value, xmax = mean_value + sd_value, y = 50), col = 'black', 
                   size = 1.1) +
    theme_bw() + theme_hist() + theme(aspect.ratio = 0.2, 
          axis.text.x = element_blank(), 
          plot.margin = unit(c(0,0,0,0), "cm"), 
          axis.line.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    scale_colour_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
    scale_fill_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
    xlab(NULL) + 
    ylab(NULL) + 
    xlim(pivot_longer(evaluations, Dslope:Dslope_rescaled) %>% .$value %>% range(., na.rm = T))
  
  png(paste0(dir, '/', 'slope.png'), width = 1750, height = 1750, res = 300)
  print(d_slope_tophist + plot_spacer() + d_slope_xy + d_slope_righthist + 
    plot_layout(ncol = 2, nrow = 2, widths = c(1, 0.3), heights = c(0.3, 1)))
  dev.off()
  
  # histograms for slope ----
  
  hist_slope_rescaled <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dslope_rescaled), alpha = 1, fill = col, col = 'transparent') +
    geom_vline(aes(xintercept = 1)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          aspect.ratio = 3, 
          axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0), 
          axis.text.y = element_blank(),
          axis.title.y  = element_blank(),
          axis.line = element_line(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"), 
          plot.title = element_text(colour = 'gray20', size = 13, hjust = 0.5, vjust= -5)) + 
    coord_flip() + 
    xlim(pivot_longer(evaluations, Dslope:Dslope_rescaled) %>% .$value %>% range(., na.rm = T)  + c(-0.5, 0.5))+
    ggtitle(all_evals_agg[[i]] %>% filter(name == 'Dslope_rescaled') %>% .$mean_value %>% round(., 2)) + 
    ylab(NULL)
  
  hist_slope <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dslope), alpha = 1, fill = col, col = 'transparent') +
    geom_vline(aes(xintercept = 1)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          aspect.ratio = 3, 
          axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0), 
          axis.line = element_line(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"), 
          plot.title = element_text(colour = 'gray20', size = 13, hjust = 0.5, vjust= -5)) + 
    coord_flip() + 
    xlim(pivot_longer(evaluations, Dslope:Dslope_rescaled) %>% .$value %>% range(., na.rm = T) + c(-0.5, 0.5)) +
    ggtitle(all_evals_agg[[i]] %>% filter(name == 'Dslope') %>% .$mean_value %>% round(., 2)) + 
  ylab(NULL)
  
  png(paste0(dir, '/', 'slope_hist.png'), width = 700, height = 700, res = 300)
  print(hist_slope + hist_slope_rescaled + 
          plot_layout(ncol = 2, nrow = 1, widths = c(1, 1), heights = c(1, 1)))
  dev.off()
  
  
  # combined intercept plots ----
  
  
  d_intercept_xy <- 
    ggplot(data = evaluations) + 
    geom_point(aes(x = Dintercept, y = Dintercept_rescaled)) + 
    geom_hline(aes(yintercept = 0)) + 
    geom_vline(aes(xintercept = 0)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20)) + 
    ylab('rescaled intercept') + 
    xlab('original intercept') + 
    ylim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T))
  
  d_intercept_righthist <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dintercept_rescaled), alpha = 1, fill = 'gray75', col = 'transparent', bins = 100) +
    geom_point(data = all_evals_agg[[i]] %>% filter(name == 'Dintercept_rescaled'), 
               aes(x = mean_value, y = 50), col = 'black', size = 4) +
    geom_linerange(data = all_evals_agg[[i]] %>% filter(name == 'Dintercept_rescaled'), 
                   aes(xmin = mean_value - sd_value, xmax = mean_value + sd_value, y = 50), col = 'black', 
                   size = 1.1) +
    theme_bw() +  theme_hist() +
    theme(aspect.ratio = 5, 
          axis.text.x = element_text(angle = 270, vjust = 0.5), 
          axis.text.y = element_blank(), 
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")) + 
    xlab(NULL) + ylab(NULL) + coord_flip() + 
    xlim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T))
  
  d_intercept_tophist <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dintercept), alpha = 1, fill = 'gray75', col = 'transparent', bins = 100) + 
    geom_point(data = all_evals_agg[[i]] %>% filter(name == 'Dintercept'), 
               aes(x = mean_value, y = 50), col = 'black', size = 4) +
    geom_linerange(data = all_evals_agg[[i]] %>% filter(name == 'Dintercept'), 
                   aes(xmin = mean_value - sd_value, xmax = mean_value + sd_value, y = 50), col = 'black', 
                   size = 1.1) +
    theme_bw() + theme_hist() + theme(aspect.ratio = 0.2, 
                                      axis.text.x = element_blank(), 
                                      plot.margin = unit(c(0,0,0,0), "cm"), 
                                      axis.line.x = element_blank(), 
                                      axis.ticks.x = element_blank()) +
    scale_colour_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
    scale_fill_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
    xlab(NULL) + 
    ylab(NULL) + 
    xlim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T))
  
  
  
  png(paste0(dir, '/', 'intercept.png'), width = 1750, height = 1750, res = 300)
  print(d_intercept_tophist + plot_spacer() + d_intercept_xy + d_intercept_righthist + 
    plot_layout(ncol = 2, nrow = 2, widths = c(1, 0.3), heights = c(0.3, 1)))
  dev.off()
  
  # simple histograms for intercept ----
  
  hist_intercept_rescaled <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dintercept_rescaled), alpha = 1, fill = col, col = 'transparent') +
    geom_vline(aes(xintercept = 1)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          aspect.ratio = 3, 
          axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0), 
          axis.text.y = element_blank(),
          axis.title.y  = element_blank(),
          axis.line = element_line(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"), 
          plot.title = element_text(colour = 'gray20', size = 13, hjust = 0.5, vjust= -5)) + 
    coord_flip() + 
    xlim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T))+
    ggtitle(all_evals_agg[[i]] %>% filter(name == 'Dintercept_rescaled') %>% .$mean_value %>% round(., 2)) +
    xlim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T) + c(-0.5, 0.5)) +
    ylab(NULL)
  
  hist_intercept <- ggplot(evaluations) + 
    geom_histogram(aes(x = Dintercept), alpha = 1, fill = col, col = 'transparent') +
    geom_vline(aes(xintercept = 1)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          aspect.ratio = 3, 
          axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0), 
          axis.line = element_line(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"), 
          plot.title = element_text(colour = 'gray20', size = 13, hjust = 0.5, vjust= -5)) + 
    coord_flip() + 
    xlim(pivot_longer(evaluations, Dintercept:Dintercept_rescaled) %>% .$value %>% range(., na.rm = T) + c(-0.5, 0.5)) +
    ggtitle(all_evals_agg[[i]] %>% filter(name == 'Dintercept') %>% .$mean_value %>% round(., 2)) + 
    ylab(NULL)
  
  png(paste0(dir, '/', 'intercept_hist.png'), width = 700, height = 700, res = 300)
  print(hist_intercept + hist_intercept_rescaled + 
          plot_layout(ncol = 2, nrow = 1, widths = c(1, 1), heights = c(1, 1)))
  dev.off()
  
}

writexl::write_xlsx(summary_evals, path = 'figures/model-prediction-figures/validation_best_models_aggregated/aggregate_summaries.xlsx')

