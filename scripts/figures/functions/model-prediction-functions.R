
# function to rescale between 0 and 1 ----

rescale_01 = function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm=T))}

# function for hex density plots of observed vs. predicted abundance ----

observed_predicted_plot <- function(input_data, 
                                    rescale = T, 
                                    model_level, # options are fitted_model or plot_level
                                    directory, 
                                    name, 
                                    width, 
                                    height, 
                                    upper_limit = 3,
                                    levels = c('glm', 'gam', 'gbm', 'rf'), 
                                    nbins = 100
                                    ){
  
  require(tidyverse)
  
  input_data$fitted_model <- factor(as.factor(input_data$fitted_model), levels)
  
  input_data %>% 
    do(verification_predict_mean = input_data$verification_predict_mean)

  # fitted model level
  if(model_level == 'fitted_model'){
  if(rescale == T){
  # fitted_model level analysis
  model_outputs <- input_data %>% 
    group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
    do(observed  = rescale_01(log10(unlist(.[,grepl('observed', names(.))]))),
       predicted = rescale_01(log10(replace(unlist(.[,grepl('predict', names(.))]), unlist(.[,grepl('predict', names(.))]) < 1, 1)))) %>% 
    ungroup() %>% 
    group_by(dataset, fitted_model) %>% 
    do(observed = unlist(.[,grepl('observed', names(.))]), 
       predicted = unlist(.[,grepl('predict', names(.))])) %>% 
    unnest(c(observed, predicted))
  }else{
    model_outputs <- input_data %>% 
      group_by(dataset, fitted_model, species_name) %>% 
      do(observed = log10(unlist(.[,grepl('observed', names(.))])), 
         predicted = log10(replace(unlist(.[,grepl('predict', names(.))]), unlist(.[,grepl('predict', names(.))]) < 1, 1))) %>% 
      unnest(c(observed, predicted))
  }}
  
  
  if(model_level == 'plot_level'){
  # plot_level model level
  if(rescale == T){
    # fitted_model level analysis
    model_outputs <- input_data %>% 
      group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
      do(observed  = rescale_01(log10(unlist(.[,grepl('observed', names(.))]))),
         predicted = rescale_01(log10(replace(unlist(.[,grepl('predict', names(.))]), unlist(.[,grepl('predict', names(.))]) < 1, 1)))) %>% 
      ungroup() %>% 
      group_by(dataset, fitted_model, plot_level) %>% 
      do(observed = unlist(.[,grepl('observed', names(.))]), 
         predicted = unlist(.[,grepl('predict', names(.))])) %>% 
      unnest(c(observed, predicted))
  }else{
    model_outputs <- input_data %>% 
      group_by(data_set, fitted_model, plot_level, species_name) %>% 
      do(observed = log10(unlist(.[,grepl('observed', names(.))])), 
         predicted = log10(replace(unlist(.[,grepl('predict', names(.))]), unlist(.[,grepl('predict', names(.))]) < 1, 1))) %>% 
      unnest(c(observed, predicted))
  }}
  
  model_outputs$predicted[model_outputs$predicted > as.numeric(quantile(model_outputs$predicted, 0.99, na.rm =T))] <- quantile(model_outputs$predicted, 0.99, na.rm =T)
  model_outputs$predicted[model_outputs$predicted < as.numeric(quantile(model_outputs$predicted, 0.01, na.rm =T))] <- quantile(model_outputs$predicted, 0.01, na.rm =T)
  
  model_outputs <- model_outputs %>% 
    mutate(observed = plyr::round_any(.$observed, 0.1)) %>% 
    group_by(dataset, fitted_model, species_name, observed) %>% 
    do(predicted = mean(.$predicted, na.rm = T)) %>% 
    unnest()
  
  # create basic plot
  base_plot <- ggplot(data = model_outputs) + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
    geom_abline() + 
    viridis::scale_fill_viridis(option = 3, begin = 0.2, end = 0.9,
                                limits = c(NA, upper_limit), 
                                na.value = viridis::viridis(10,option = 3, end = 0.9)[10], 
                                breaks = log10(c(1, 10^0.5, 10, 10^1.5, 100)), 
                                labels = signif(c(1, 10^0.5, 10, 10^1.5, 100), 1), 
                                name = 'no. obs.') +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
          aspect.ratio = 1) 
  
  # add faceting
  if(model_level == 'fitted_model'){facet_plot <- base_plot + 
    facet_wrap(~fitted_model, scales='free') + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
    geom_abline() + 
    ggtitle(label = unique(paste(input_data$dataset, ifelse(input_data$cross_validation == 'basic', 'interpolation', 'extrapolation'), sep = ' ')))}
  
  if(model_level == 'plot_level'){facet_plot <- base_plot + 
    facet_wrap(fitted_model~plot_level) + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
    geom_abline() + 
    ggtitle(label = unique(paste(input_data$dataset, ifelse(input_data$cross_validation == 'basic', 'interpolation', 'extrapolation'), input_data$abundance_response, sep = ' ')))}
  
  if(rescale == T){facet_plot <- facet_plot + scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + scale_y_continuous(breaks =  c(0, 0.2, 0.4, 0.6, 0.8, 1))}
  
  # save in nested directories
  dir.create(directory, recursive = T)
  png(file = paste0(directory, '/', name,'.png'), width = width, height = height, res = 200)
  print(facet_plot)
  dev.off()
  
}
