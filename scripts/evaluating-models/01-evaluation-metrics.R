# Managing outputs of the abundance models 

# Set working directory if needed ----
#setwd('/Users/cwaldock/Dropbox/ETH REEF FUTURES/abundance-comparison')

# Load packages ----
lib_vect <- c('tidyverse', 'summarytools', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# Custom function ---- 

# managing data groupings for plots
add_family_plot_column <- function(plot_data){
  
  # edit transformations
  plot_data$transformation[is.na(plot_data$transformation)] <- ''
  plot_data$transformation <- ifelse(plot_data$transformation == '', '', paste0('-', plot_data$transformation))
  
  # edit family groups
  plot_data$family <- recode(plot_data$family, discrete = "multinomial", discrete_multinomial = "multinomial", 
                             continuous_gaussian = 'gaussian',continuous = 'gaussian', continuous_poisson = 'poisson', 
                             zip = 'poisson')
  plot_data$family[is.na(plot_data$family)] <- 'gaussian'
  
  # edit zero-inflations
  plot_data$zi <- ifelse(plot_data$zi, '-ZI', '')
  plot_data$zi[is.na(plot_data$zi)] <- ''
  
  plot_data$family_plot <- paste(plot_data$family, plot_data$zi, plot_data$transformation, sep = '')
  
  return(plot_data)
  
}


# for editing x axis
theme_remove_x <- function(){theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank(), 
                                   plot.margin = unit(c(0,0,0,0), units = 'cm'))}


# Load in all data ---- 
all_files <- list.files('results_hpc', recursive = T, full.names = T) 
predictions <- list()
for(i in 1:length(all_files)){
  load(all_files[i])
  extracted_predictions$abundance_response <- gsub('predictions_', '', str_split(all_files[i], '/')[[1]][3])
  predictions[[i]] <- extracted_predictions
}
all_files_bind <- do.call(rbind, predictions)

# Load in fish abundance groups ---- 
load("data/rls_abun_modelling_data_v2.RData")
abundance_key <- abundance_key %>% rename(., species_name = TAXONOMIC_NAME) %>% as_tibble

# Function for assessment metrics and plots ----
abundance_assessment_metrics <- function(predictions, observations, locations){

  # observations
  
  # Linear model between values
  lm_test           <- lm(observations ~ predictions)
  sum_lm            <- summary(lm_test)
  cor.test_pearson  <- cor.test(observations, predictions, method = 'pearson')
  cor.test_spearman <- cor.test(observations, predictions, method = 'spearman')
  
  # summaries from linear model
  residual_standard_error <- lm_test$sigma
  intercept <- coef(lm_test)[1]
  slope     <- coef(lm_test)[2]
  
  # Accuracy (A-overall)
  # root mean squared error between average predicted abundance 
  # across all sites and observed abundance at each site
  Armse <- sqrt(mean((observations - predictions)^2, na.rm = T)) # root mean squared error gives more weight to big deviations
  Amae  <- mean((observations - predictions), na.rm = T)  # mean absolute error weights all errors the same - if positive the value of observations is > the values of predictions
  
  # Discrimination
  Dintercept <- coef(lm_test)[1]
  Dslope     <- coef(lm_test)[2]
  Dpearson   <- cor.test_pearson$estimate
  Dspearman  <- cor.test_spearman$estimate
  
  # Precision
  Psd         <- sd(predictions) 
  Pdispersion <- sd(predictions) / sd(observations)
  Pr2         <- summary(lm_test)$r.squared
  
  metric_summary <- data.frame(Armse = Armse, 
                               Amae  = Amae, 
                               Dintercept = Dintercept, 
                               Dslope = Dslope, 
                               Dpearson = Dpearson, 
                               Dspearman = Dspearman, 
                               Psd = Psd, 
                               Pdispersion = Pdispersion, 
                               Pr2 = Pr2)
  
  return(metric_summary)
  
}

# function for seqential plots
Metrics <- c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2')

# aim, minimum, maximum
Target  <- list(Armse    = c(0, -10, 10), 
                Amae     = c(0, -10, 10), 
                Dintercept  = c(0, -5, 5), 
                Dslope      = c(1,  0, 2), 
                Dpearson    = c(1,  0, 1), 
                Dspearman   = c(1,  0, 1), 
                Psd         = c(0,  0, 10), 
                Pdispersion = c(1,  0, 2), 
                Pr2         = c(1,  0, 1))

Levels <- c('glm', 'gam', 'rf', 'gbm')

# Function that  
write_plots <- function(plot_data, directory){
  for(j in 1:length(Metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), Levels)}
    print(j)
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      select(colnames(.)[-which(colnames(.)%in%Metrics)], Metrics[j]) %>% 
      mutate(Metrics = as.numeric(.[,Metrics[j]][[1]]))
    
    t_v <- Target[[j]]
    
    metric_plot_data$Metrics[!is.finite(metric_plot_data$Metrics)] <- NA
    metric_plot_data$Metrics[which(metric_plot_data$Metrics > t_v[3])] <- t_v[3]
    metric_plot_data$Metrics[which(metric_plot_data$Metrics < t_v[2])] <- t_v[2]
    metric_plot_data <- na.omit(metric_plot_data)
    
    # set lower ylim
    y_lower <- min(c(min(metric_plot_data$Metrics,na.rm=T), t_v[2]))-0.1
    y_upper <- max(c(max(metric_plot_data$Metrics,na.rm=T), t_v[3]))+0.1
    
    plots <- ggplot(data = metric_plot_data) + 
      geom_boxplot(aes(x = abundance_class, y = Metrics)) + 
      geom_hline(aes(yintercept=t_v[1]), col = 'red') + 
      ylab(Metrics[j]) + 
      ylim(y_lower, y_upper) + 
      facet_grid(fitted_model ~ family_plot, scales = 'free', drop = T) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1), 
            aspect.ratio = 1)
    
    dir.create(directory, recursive = T) 
    pdf(file = paste0(directory, '/', Metrics[j],'.pdf'), width = 14, height = 10)
    print(plots)
    dev.off()
    
  }
}
# Estimate evaluation criteria across all data ----
final_assessment_data <- all_files_bind %>% 
  # estimate evaluation criteria across all models
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$verification_observed_mean, 
                                            .$verification_predict_mean)) %>% 
  unnest(metrics) %>% 
  # bind in model identifiers
  cbind(all_files_bind[,1:10], .) %>% 
  cbind(all_files_bind[,'abundance_response'], .) %>% 
  # join in abundance information of different groups
  left_join(abundance_key, .) %>% 
  # convert to tibble
  as_tibble() %>% 
  # changes family values to ensure consistency
  add_family_plot_column

# create plot labels ----

unique(final_assessment_data$abundance_response) # 
unique(final_assessment_data$fitted_model) # 
unique(final_assessment_data$only_abundance) # 
unique(final_assessment_data$family) # 
unique(final_assessment_data$family_plot) # 
unique(final_assessment_data$transformation) # 
unique(final_assessment_data$zi) # 

# recode values
final_assessment_data$abundance_response <- recode(final_assessment_data$abundance_response, 
                                                   abun = "a", 
                                                   abunocc_2stage = "a2", 
                                                   abunocc = "ao")

final_assessment_data$family_plot2 <- recode(final_assessment_data$family_plot, 
                                             `gaussian-log`      = "g.l", 
                                             `multinomial-log`   = "mn.l", 
                                             `gaussian-log10`    = 'g.l10',
                                             `multinomial-log10` = "mn.l10", 
                                             poisson  = 'p', 
                                             nbinom   = 'nb', 
                                             tweedie  = 'tw', 
                                             gaussian = 'g', 
                                             `poisson-ZI` = 'p.zi', 
                                             `nbinom-ZI`  = 'nb.zi', 
                                             `tweedie-ZI` = 'tw.zi')

final_assessment_data$plot_levels <- paste(final_assessment_data$fitted_model, 
      final_assessment_data$abundance_response, 
      final_assessment_data$family_plot2, sep = '.')


                           

# PLOT: assessment metric ordered by median values  ----

ggplot(data = final_assessment_data) +
  geom_boxplot(aes(x = paste(fitted_model, only_abundance, family_plot, transformation, zi), 
                   y = Pr2))


  
# Figures for abun-occ ---- 

# verification abundance occupancy 
metrics_VER_abunocc <- abun_occ_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$verification_observed_mean, 
                                            .$verification_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun_occ_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

# validation abundance occupancy
metrics_VAL_abunocc <- abun_occ_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$validation_observed_mean, 
                                            .$validation_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun_occ_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

write_plots(plot_data = metrics_VER_abunocc, directory = 'figures/verification_abunocc')

write_plots(metrics_VAL_abunocc, 'figures/valdation_abunocc')




# Figures for abun ---- 

# verification abundance occupancy 
metrics_VER_abun <- abun_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$verification_observed_mean, 
                                            .$verification_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

# validation abundance occupancy
metrics_VAL_abun <- abun_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$validation_observed_mean, 
                                            .$validation_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

write_plots(metrics_VER_abun, 'figures/verification_abun')

write_plots(metrics_VAL_abun, 'figures/valdation_abun')





# Figures for abun-2stage ---- 

# verification abundance occupancy 
metrics_VER_abun2stage <- abun2stage_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$verification_observed_mean, 
                                            .$verification_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun2stage_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

# validation abundance occupancy
metrics_VAL_abun2stage <- abun2stage_models %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$validation_observed_mean, 
                                            .$validation_predict_mean)) %>% 
  unnest(metrics) %>% 
  cbind(abun2stage_models[,1:10], .) %>% 
  left_join(abundance_key, .) %>% 
  as_tibble() %>% 
  add_family_plot_column

write_plots(metrics_VER_abun, 'figures/verification_abun2stage')

write_plots(metrics_VAL_abun, 'figures/valdation_abun2stage')





# Combine all into verification - groups by model and abundance response ----
metrics_VER_abunocc$model_response    <- 'AO'
metrics_VER_abun$model_response       <- 'A'
metrics_VER_abun2stage$model_response <- 'AO-2stage'

combined_models <- rbind(metrics_VER_abunocc, metrics_VER_abun, metrics_VER_abun2stage)

combined_models_mean <- combined_models %>% 
  gather(., key = metric, value = value,  Aoverall:Pr2) %>% 
  filter(value < 1000) %>% 
  group_by(species_name, 
           model_response, 
           fitted_model, 
           metric) %>% 
  do(metric_mean = median(.$value, na.rm = T)) %>% 
  unnest()

combined_models_mean <- na.omit(do.call(data.frame,lapply(combined_models_mean, function(x) replace(x, is.infinite(x),NA))))

plot_single <- list()
for(i in 1:length(Metrics)){
  
  print(Metrics[i])
  
metric_plot_data <- combined_models_mean %>% filter(metric %in% Metrics[i])

t_v2 <- Target[[i]]

metric_plot_data$metric_mean[which(metric_plot_data$metric_mean > t_v2[3])] <- t_v2[3]
metric_plot_data$metric_mean[which(metric_plot_data$metric_mean < t_v2[2])] <- t_v2[2]
metric_plot_data <- na.omit(metric_plot_data)

# set lower ylim
y_lower <- min(c(min(metric_plot_data$metric_mean,na.rm=T), t_v2[2]))-0.1
y_upper <- max(c(max(metric_plot_data$metric_mean,na.rm=T), t_v2[3]))+0.1

plot_single[[i]] <- ggplot(data = metric_plot_data) + 
  geom_boxplot(aes(x = paste(fitted_model, model_response, sep = '-'), 
                   y = metric_mean, 
                   group = paste(fitted_model, model_response, sep = '-'),
                   fill = fitted_model, 
                   size = 0.75), 
               outlier.shape = NA) + 
  geom_point(aes(x = fitted_model, 
                 y = metric_mean, 
                 col = fitted_model, 
                 fill = fitted_model),
             position = position_nudge(x=0.1), 
             stat = "sum", 
             show.legend = F, pch = 19, alpha = 0.5) + 
  facet_wrap(~metric, scales = 'free_y') + 
  theme_classic() + 
  xlab(NULL) + 
  ylab(Metrics[i]) + 
  ylim(y_lower, y_upper) +
  scale_size_continuous(range = c(0.5, 2)) + 
  scale_fill_manual(values = brewer.pal(4, 'Dark2')) + 
  scale_colour_manual(values = brewer.pal(4, 'Dark2')) + 
  theme(aspect.ratio = 0.75, legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}

pdf(file = 'figures/combined_plots/VER_mean_assessment_species_distributions.pdf', width = 9, height = 9)
grid.arrange(plot_single[[1]] + geom_hline(aes(yintercept = Target[[1]][1]), col = 'red'),
             plot_single[[2]]+ geom_hline(aes(yintercept = Target[[2]][1]), col = 'red'),
             plot_single[[3]]+ geom_hline(aes(yintercept = Target[[3]][1]), col = 'red'),
             plot_single[[4]]+ geom_hline(aes(yintercept = Target[[4]][1]), col = 'red'),
             plot_single[[5]]+ geom_hline(aes(yintercept = Target[[5]][1]), col = 'red'),
             plot_single[[6]]+ geom_hline(aes(yintercept = Target[[6]][1]), col = 'red'),
             plot_single[[7]]+ geom_hline(aes(yintercept = Target[[7]][1]), col = 'red'),
             plot_single[[8]]+ geom_hline(aes(yintercept = Target[[8]][1]), col = 'red'))
dev.off()             




# Combine all into verification - groups by model and abundance response ----
metrics_VAL_abunocc$model_response    <- 'AO'
metrics_VAL_abun$model_response       <- 'A'
metrics_VAL_abun2stage$model_response <- 'AO-2stage'

combined_models <- rbind(metrics_VAL_abunocc, metrics_VAL_abun, metrics_VAL_abun2stage)

combined_models_mean <- combined_models %>% 
  gather(., key = metric, value = value,  Aoverall:Pr2) %>% 
  filter(value < 1000) %>% 
  group_by(species_name, 
           model_response, 
           fitted_model, 
           metric) %>% 
  do(metric_mean = median(.$value, na.rm = T)) %>% 
  unnest()

combined_models_mean <- na.omit(do.call(data.frame,lapply(combined_models_mean, function(x) replace(x, is.infinite(x),NA))))

plot_single <- list()
for(i in 1:length(Metrics)){
  
  print(Metrics[i])
  
  metric_plot_data <- combined_models_mean %>% filter(metric %in% Metrics[i])
  
  t_v2 <- Target[[i]]
  
  metric_plot_data$metric_mean[which(metric_plot_data$metric_mean > t_v2[3])] <- t_v2[3]
  metric_plot_data$metric_mean[which(metric_plot_data$metric_mean < t_v2[2])] <- t_v2[2]
  metric_plot_data <- na.omit(metric_plot_data)
  
  # set lower ylim
  y_lower <- min(c(min(metric_plot_data$metric_mean,na.rm=T), t_v2[2]))-0.1
  y_upper <- max(c(max(metric_plot_data$metric_mean,na.rm=T), t_v2[3]))+0.1
  
  plot_single[[i]] <- ggplot(data = metric_plot_data) + 
    geom_boxplot(aes(x = paste(fitted_model, model_response, sep = '-'), 
                     y = metric_mean, 
                     group = paste(fitted_model, model_response, sep = '-'),
                     fill = fitted_model, 
                     size = 0.75), 
                 outlier.shape = NA) + 
    geom_point(aes(x = fitted_model, 
                   y = metric_mean, 
                   col = fitted_model, 
                   fill = fitted_model),
               position = position_nudge(x=0.1), 
               stat = "sum", 
               show.legend = F, pch = 19, alpha = 0.5) + 
    facet_wrap(~metric, scales = 'free_y') + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(Metrics[i]) + 
    ylim(y_lower, y_upper) +
    scale_size_continuous(range = c(0.5, 2)) + 
    scale_fill_manual(values = brewer.pal(4, 'Dark2')) + 
    scale_colour_manual(values = brewer.pal(4, 'Dark2')) + 
    theme(aspect.ratio = 0.75, legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
}

pdf(file = 'figures/combined_plots/VAL_mean_assessment_species_distributions.pdf', width = 9, height = 9)
grid.arrange(plot_single[[1]] + geom_hline(aes(yintercept = Target[[1]][1]), col = 'red'),
             plot_single[[2]]+ geom_hline(aes(yintercept = Target[[2]][1]), col = 'red'),
             plot_single[[3]]+ geom_hline(aes(yintercept = Target[[3]][1]), col = 'red'),
             plot_single[[4]]+ geom_hline(aes(yintercept = Target[[4]][1]), col = 'red'),
             plot_single[[5]]+ geom_hline(aes(yintercept = Target[[5]][1]), col = 'red'),
             plot_single[[6]]+ geom_hline(aes(yintercept = Target[[6]][1]), col = 'red'),
             plot_single[[7]]+ geom_hline(aes(yintercept = Target[[7]][1]), col = 'red'),
             plot_single[[8]]+ geom_hline(aes(yintercept = Target[[8]][1]), col = 'red'))
dev.off()             



