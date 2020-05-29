
# function to rescale between 0 and 1 ----

rescale_01 = function(x){
  
  (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm=T))
  
  }

# function to unnest large data from as a data.table ----

# https://www.johannesbgruber.eu/post/a-faster-unnest/#fn1

unnest_dt2 <- function(tbl, ...) {
  
  tbl <- as.data.table(tbl)
  
  col <- ensyms(...)
  
  clnms <- syms(setdiff(colnames(tbl), as.character(col)))
  
  tbl <- as.data.table(tbl)
  
  tbl <- eval(
    expr(tbl[, lapply(.SD, unlist), by = list(!!!clnms), .SDcols = as.character(col)])
  )
  
  colnames(tbl) <- c(as.character(clnms), as.character(col))
  
  tbl
}

# function for hex density plots of observed vs. predicted abundance ----

observed_predicted_plot <- function(input_data, 
                                    rescale = T, 
                                    model_level, # options are fitted_model or plot_level
                                    directory, 
                                    name, 
                                    width, 
                                    height, 
                                    upper_limit = 1.5,
                                    nbins = 20,
                                    levels = c('glm', 'gam', 'gbm', 'rf'),
                                    option = 3
                                    ){
  
  require(tidyverse)
  
  input_data$fitted_model <- factor(as.factor(input_data$fitted_model), levels)
  
  model_outputs <- input_data
  
  # remove turn values that are <0 to NAs
  # rescale predictions
  model_outputs$predicted <- pbapply::pblapply(model_outputs$validation_predict_mean, function(x) return(ifelse(x < 0, NA, x)))
  model_outputs$predicted <- pbapply::pblapply(model_outputs$predicted, function(x) return(ifelse(x < 1, 1, x)))
  model_outputs$predicted <- pbapply::pblapply(model_outputs$predicted, function(x) rescale_01(log10(x)))
  # rescale observations
  model_outputs$observed <- pbapply::pblapply(model_outputs$validation_observed_mean, function(x) rescale_01(log10(x)))

  # unnest using data.table because dplyr is slow for this much data
  model_outputs <- as_data_frame(unnest_dt2(model_outputs %>% select(-validation_observed_mean, -validation_predict_mean), predicted, observed))
  
  # create a transformations label
  model_outputs$transformation <- gsub('gam.|gbm.|glm.|rf.|abun.|ao-2.|cc.', '', model_outputs$plot_level)
  
  # truncate the data to 99th percentiles
  model_outputs$predicted[model_outputs$predicted > as.numeric(quantile(model_outputs$predicted, 0.99, na.rm =T))] <- quantile(model_outputs$predicted, 0.99, na.rm =T)
  model_outputs$predicted[model_outputs$predicted < as.numeric(quantile(model_outputs$predicted, 0.01, na.rm =T))] <- quantile(model_outputs$predicted, 0.01, na.rm =T)
  
  # aggregate at a species level
  model_outputs <- model_outputs %>% 
    mutate(observed = plyr::round_any(.$observed, 1/nbins)) %>% 
    group_by(dataset, fitted_model, plot_level, species_name, abundance_response, transformation, observed) %>% 
    do(predicted = mean(.$predicted, na.rm = T))
  
  # unnest quickly
  model_outputs <- as_data_frame(unnest_dt2(model_outputs, predicted))
  
  
  # plots for aggregated to one panel ----
  # save non-facet wrapped plot for overall image
  if(model_level == 'aggregated'){
    
  base_plot <- ggplot(data = model_outputs) + 
    geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
    viridis::scale_fill_viridis(option = option, begin = 0.2, end = 0.9,
                                limits = c(NA, upper_limit), 
                                na.value = viridis::viridis(10,option = option, end = 0.9)[10], 
                                breaks = log10(c(1, 10^0.5, 10, 10^1.5, 100, 300, 1000)), 
                                labels = signif(c(1, 10^0.5, 10, 10^1.5, 100, 300, 1000), 1), 
                                name = 'no. obs.') +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
          aspect.ratio = 1) + 
    stat_density2d(aes(x = observed, y = predicted), col = 'black', size = 1, n = 200, alpha = 0.75) + 
    geom_abline(size = 1) 
  
  #if(density == T){base_plot + stat_density2d(aes(x = observed, y = predicted), col = 'black', size = 1, n = 200)}
    
    dir.create(directory, recursive = T)
    png(file = paste0(directory, '/', 'legend','.png'), width = 1000, height = 1000, res = 400)
    print(plot(get_legend(base_plot)))
    dev.off()
    
    png(file = paste0(directory, '/', name,'.png'), width = width, height = height, res = 400)
    print(base_plot)
    dev.off()
    
    return('plot base')
    
  }

  
  # plots for all plot levels ----
  

  plot_levels_plot <- list()
  for(i in 1:length(levels)){
    
     # create basic plot
     base_plot <- ggplot(data = model_outputs %>% filter(fitted_model == levels[i])) + 
       geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
       geom_abline() + 
       viridis::scale_fill_viridis(option = 3, begin = 0.2, end = 0.9,
                                   limits = c(NA, upper_limit), 
                                   na.value = viridis::viridis(10,option = 3, end = 0.9)[10], 
                                   breaks = log10(c(1, 10^0.5, 10, 10^1.5, 100, 300, 1000)), 
                                   labels = signif(c(1, 10^0.5, 10, 10^1.5, 100, 300, 1000), 1), 
                                   name = 'no. obs.') +
       theme_bw() + 
       theme(panel.grid = element_blank(), 
             strip.background = element_rect(fill = 'grey90', colour = 'grey90'), 
             aspect.ratio = 1)
     
     # add faceting to plot level
     facet_plot <- base_plot + 
       facet_grid(fitted_model + abundance_response ~ transformation) + 
       geom_hex(aes(x = observed, y = predicted, fill = stat(log10(count))), bins = nbins) + 
       geom_abline() + 
       scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks =  c(0, 0.5, 1))
     
     if(i == 1){
       dir.create(directory, recursive = T)
       png(file = paste0(directory, '/', 'legend','.png'), width = 1000, height = 1000, res = 400)
       print(plot(get_legend(facet_plot)))
       dev.off()
     }
     
     plot_levels_plot[[i]] <- facet_plot + theme(legend.position = 'none')
  
  }
  
  # plot pngs
  png(file = paste0(directory, '/', name,'.png'), width = width, height = height, res = 400)
  all_plots <- do.call("grid.arrange", c(plot_levels_plot, ncol = 1))
  print(all_plots)
  dev.off()
  
  
}


#  ----


