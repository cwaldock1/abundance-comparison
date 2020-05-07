# Managing outputs of the abundance models 

# Set working directory if needed ----
#setwd('/Users/cwaldock/Dropbox/ETH REEF FUTURES/abundance-comparison')

# Load packages ----
lib_vect <- c('tidyverse', 'summarytools', 'cowplot', 'RColorBrewer', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# Custom functions ---- 

# all functions for evaluating outputs
source('scripts/evaluating-models/functions/evaluation_functions.R')

# for editing x axis
theme_remove_x <- function(){theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank(), 
                                   plot.margin = unit(c(0,0,0,0), units = 'cm'))}


# Load in all data ---- 

# get folders of interest
result_folders <- list.dirs('results/predictions', recursive = F)

for(folder in 1:length(result_folders)){

# here in the future run iteratively for each folder of interest

all_files <- list.files('results/predictions', recursive = T, full.names = T)

# remove suitability files

all_files <- all_files[-grep('suitability', all_files)]

# get folder of interest based on loop

all_files <- all_files[grep(result_folders[folder], all_files)]

# in the future, run this over each modelling subset for scalability

bind_files <- bind_results(all_files)

# Clean data levels ----

# some of the encodings output from the models aren't well match so match values here across modelling frameworks using the clean_levels functions

clean_files <- clean_levels(bind_files)

# Calculate assessment metrics ----

model_assessment <- clean_files %>% 
  rowwise() %>% 
  do(metrics = abundance_assessment_metrics(.$verification_observed_mean, 
                                            .$verification_predict_mean)) %>% 
  unnest(metrics) %>% 
  bind_cols(clean_files[,1:13], .) %>% 
  mutate(cross_validation = gsub('results/predictions/','', result_folders[folder]))

# Attached abundance information into metrics ----

model_assessment <- left_join(model_assessment, 
                              readRDS(paste0('data/', if(model_assessment$dataset == 'bbs'){'bbs_species_properties.RDS'}else{'rls_species_properties.RDS'})) %>% 
                                rename(., species_name = TAXONOMIC_NAME))

dir.create('results/model_assessment/', recursive = T)
  
saveRDS(model_assessment, file = paste0('results/model_assessment/', gsub('results/predictions/','', result_folders[folder]), '.rds'))

}


#### END OF SCRIPT 

# Load in fish abundance groups ---- 
load("data/rls_abun_modelling_data_v2.RData")

abundance_key <- abundance_key %>% rename(., species_name = TAXONOMIC_NAME) %>% as_tibble

# function for seqential plots
metrics <- c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2')

# aim, minimum, maximum
target  <- list(Armse    = c(0, -10, 10), 
                Amae     = c(0, -10, 10), 
                Dintercept  = c(0, -5, 5), 
                Dslope      = c(1,  0, 2), 
                Dpearson    = c(1,  0, 1), 
                Dspearman   = c(1,  0, 1), 
                Psd         = c(0,  0, 10), 
                Pdispersion = c(1,  0, 2), 
                Pr2         = c(1,  0, 1))

levels <- c('glm', 'gam', 'rf', 'gbm')

# Estimate evaluation criteria across all data ----

final_assessment_data <- bind_files %>% 
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



