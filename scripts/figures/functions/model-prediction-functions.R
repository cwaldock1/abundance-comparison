
# function to rescale between 0 and 1 ----

rescale_01 = function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm=T))}

# function for hex density plots of observed vs. predicted abundance ----

observed_predicted_plot <- function(input_data, 
                                    rescale = T, 
                                    model_level, # options are fitted_model or plot_level
                                    directory, 
                                    name, 
                                    width, 
                                    height
                                    ){
  
  require(tidyverse)

  # fitted model level
  if(model_level == 'fitted_model'){
  if(rescale == T){
  # fitted_model level analysis
  model_outputs <- input_data %>% 
    group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
    do(observed  = rescale_01(log10(unlist(.[,grepl('observed', names(.))]))),
       predicted = rescale_01(log10(unlist(.[,grepl('predict', names(.))])))) %>% 
    ungroup() %>% 
    group_by(dataset, fitted_model) %>% 
    do(observed = unlist(.[,grepl('observed', names(.))]), 
       predicted = unlist(.[,grepl('predict', names(.))])) %>% 
    unnest(c(observed, predicted))
  }else{
    model_outputs <- input_data %>% 
      group_by(dataset, fitted_model) %>% 
      do(observed = unlist(.[,grepl('observed', names(.))]), 
         predicted = unlist(.[,grepl('predict', names(.))])) %>% 
      unnest(c(observed, predicted))
  }}
  
  if(model_level == 'plot_level'){
  # plot_level model level
  if(rescale == T){
    # fitted_model level analysis
    model_outputs <- input_data %>% 
      group_by(dataset, cross_validation, fitted_model, abundance_response, plot_level, species_name) %>% 
      do(observed  = rescale_01(log10(unlist(.[,grepl('observed', names(.))]))),
         predicted = rescale_01(log10(unlist(.[,grepl('predict', names(.))])))) %>% 
      ungroup() %>% 
      group_by(dataset, fitted_model, plot_level) %>% 
      do(observed = unlist(.[,grepl('observed', names(.))]), 
         predicted = unlist(.[,grepl('predict', names(.))])) %>% 
      unnest(c(observed, predicted))
  }else{
    model_outputs <- input_data %>% 
      group_by(data_set, fitted_model, plot_level) %>% 
      do(observed = unlist(.[,grepl('observed', names(.))]), 
         predicted = unlist(.[,grepl('predict', names(.))])) %>% 
      unnest(c(observed, predicted))
  }}
  
  # create basic plot
  base_plot <- ggplot(data = model_outputs) + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = 200) + 
    geom_abline() + 
    viridis::scale_fill_viridis(option = 3, end = 0.9, begin = 0.3) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill = 'grey90', colour = 'grey90')) + 
    ggtitle(label = unique(paste(input_data$dataset, ifelse(input_data$cross_validation == 'basic', 'interpolation', 'extrapolation'), input_data$abundance_response, sep = ' ')))
  
  # add faceting
  if(model_level == 'fitted_model'){facet_plot <- base_plot + 
    facet_wrap(~fitted_model, scales='free')}
  
  if(model_level == 'plot_level'){facet_plot <- base_plot + 
    facet_wrap(fitted_model~plot_level, scales='free') + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = 100)}
  
  # save in nested directories
  dir.create(directory, recursive = T)
    pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
  print(facet_plot)
  dev.off()
}
