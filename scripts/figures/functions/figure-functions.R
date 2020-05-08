# functions to produce different figures

# write_plots: creates plots for each metric -----

all_model_plots <- function(plot_data, # object after running abundance_assesment_metrics and add_family_plot_column (to be renamed)
                        metrics = c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2'), 
                        targets = list(Armse    = c(0, -10, 10), 
                                       Amae     = c(0, -10, 10), 
                                       Dintercept  = c(0, -5, 5), 
                                       Dslope      = c(1,  0, 2), 
                                       Dpearson    = c(1,  0, 1), 
                                       Dspearman   = c(1,  0, 1), 
                                       Psd         = c(0,  0, 10), 
                                       Pdispersion = c(1,  0, 2), 
                                       Pr2         = c(1,  0, 1)), 
                        levels = c('glm', 'gam', 'gbm', 'rf'),
                        directory){
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      dplyr::select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) %>% 
      mutate(metrics = as.numeric(data.frame(plot_data)[,metrics[j]]))
    
    t_v <- targets[[j]]
    
    metric_plot_data$metrics[!is.finite(metric_plot_data$metrics)] <- NA
    metric_plot_data$metrics[which(metric_plot_data$metrics > t_v[3])] <- t_v[3]
    metric_plot_data$metrics[which(metric_plot_data$metrics < t_v[2])] <- t_v[2]
    metric_plot_data <- na.omit(metric_plot_data)
    
    # set lower ylim
    y_lower <- min(c(min(metric_plot_data$metrics,na.rm=T), t_v[2]))-0.1
    y_upper <- max(c(max(metric_plot_data$metrics,na.rm=T), t_v[3]))+0.1
  
    # reorder metric levels
    metric_plot_data$mean_abundance_perc_fact <- as.numeric(as.factor(metric_plot_data$mean_abundance_perc))
    metric_plot_data$frequency_perc_fact      <- as.numeric(as.factor(metric_plot_data$frequency_perc))
    
    # check ordering
    #ggplot(metric_plot_data) + 
    #  geom_point(aes(x = mean_abundance_perc_fact, y = mean_abundance_perc))
    cols <- colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(2)
    plots <-  ggplot() + 
      geom_line(data = metric_plot_data, 
                aes(x = mean_abundance_perc_fact, 
                    y = metrics), alpha = 0.25, col = cols[1], lwd = 0.1) + 
      geom_line(data = metric_plot_data, 
                aes(x = frequency_perc_fact, 
                    y = metrics), alpha = 0.25, col = 'red', lwd = 0.1) + 
      stat_smooth(data = metric_plot_data, 
                  aes(x = mean_abundance_perc_fact, 
                      y = metrics), 
                  col = cols[1], size =  0.5, se = F) + 
      stat_smooth(data = metric_plot_data, 
                  aes(x = frequency_perc_fact, 
                      y = metrics), 
                  col = 'red', size = 0.5, se = F) + 
      geom_hline(data = data.frame(t_v=t_v), aes(yintercept = t_v[1]), lty=2, col='black') + 
      ylab(metrics[j]) + 
      xlab('percentile rank \n (abundance = blue; frequency = red)') +
      facet_grid(fitted_model ~ family_grouped_simple, drop = T) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1), 
            aspect.ratio = 1, 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank())
    
    dir.create(directory, recursive = T) 
    pdf(file = paste0(directory, '/', metrics[j],'.pdf'), width = 14, height = 10)
    print(plots)
    dev.off()
    
  }
}


# function to produce plots with outputs combined across model types ----


combined_assessment_metrics <- function(plot_data, # a wide dataframe of all assessment calculations.
                                        response = 'all',
                                        directory, 
                                        name,
                                        width, 
                                        height,
                                        metrics = c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2'),
                                        targets = list(Armse    = c(0, -10, 10),
                                                       Amae     = c(0, -10, 10), 
                                                       Dintercept  = c(0, -5, 5), 
                                                       Dslope      = c(1,  0, 2), 
                                                       Dpearson    = c(1,  0, 1), 
                                                       Dspearman   = c(1,  0, 1), 
                                                       Psd         = c(0,  0, 10), 
                                                       Pdispersion = c(1,  0, 2), 
                                                       Pr2         = c(1,  0, 1)), 
                                        colours = brewer.pal(4, 'Dark2'), 
                                        levels = c('glm', 'gam', 'gbm', 'rf')){
  
  # gather assessments together
  
  plot_data <- plot_data %>% 
    gather(., key = metric, value = value,  Armse:Pr2) %>% 
    group_by(species_name,
             dataset,
             cross_validation,
             abundance_response,
             fitted_model,
             metric) %>% 
    do(metric_median = median(.$value, na.rm = T)) %>% 
    unnest(c(metric_median))
  
  # clean data
  
  plot_data <- na.omit(do.call(data.frame,lapply(plot_data, function(x) replace(x, is.infinite(x),NA))))
  
  # loop through and plot each metric independently
  
  metric_plots <- list()
  for(i in 1:length(metrics)){
    
    if(i == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    
    
    metric_plot_data <- plot_data %>% filter(metric %in% metrics[i])
    
    t_v2 <- targets[[i]]
    
    if(sum(metric_plot_data$metric_median > t_v2[3]) != 0){
    metric_plot_data$metric_median[which(metric_plot_data$metric_median > t_v2[3])] <- t_v2[3]}
    if(sum(metric_plot_data$metric_median < t_v2[2]) != 0){
    metric_plot_data$metric_median[which(metric_plot_data$metric_median < t_v2[2])] <- t_v2[2]}
    metric_plot_data <- na.omit(metric_plot_data)
    
    # set lower ylim
    if(min(metric_plot_data$metric_median,na.rm=T) > t_v2[2]){y_lower = min(metric_plot_data$metric_median,na.rm=T)-0.1}else{y_lower <- t_v2[2]}
    if(max(metric_plot_data$metric_median,na.rm=T) < t_v2[3]){y_upper = max(metric_plot_data$metric_median,na.rm=T)-0.1}else{y_upper <- t_v2[3]}
    
    if(y_lower > t_v2[1]){y_lower <- t_v2[1]}
    if(y_upper < t_v2[1]){y_upper <- t_v2[1]}
    
    p <- ggplot(data = metric_plot_data)
    
    if(response == 'all'){
      p + geom_boxplot(aes(x = abundance_response, 
                           y = metric_median, 
                           group = paste0(abundance_response, '  ', fitted_model),
                           fill = fitted_model), 
                       outlier.shape = NA, lwd = 0.5, fatten = 2) -> p2 
    }
    
    if(response == 'abundance_response'){
      p + geom_boxplot(aes(x = abundance_response, 
                           y = metric_median, 
                           group = abundance_response,
                           fill = abundance_response), 
                       outlier.shape = NA, lwd = 0.5, fatten = 2) -> p2}
    
    if(response == 'fitted_model'){
      p + geom_boxplot(aes(x = fitted_model, 
                           y = metric_median, 
                           group = fitted_model,
                           fill = fitted_model), 
                       outlier.shape = NA, lwd = 0.5, fatten = 2) -> p2}
    
      p2 + geom_hline(data = data.frame(t_v2=t_v2), aes(yintercept = t_v2[1]), lty=2, col='red') + 
      facet_wrap(~dataset) + 
      theme_classic() + 
      xlab(NULL) + 
      ylab(metrics[i]) + 
      ylim(y_lower, y_upper) +
      scale_size_continuous(range = c(0.5, 2)) + 
      scale_fill_manual(values = colours) + 
      scale_colour_manual(values = colours) + 
      theme(aspect.ratio = 0.75, legend.position = 'none', 
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + 
      ggtitle(unique(metric_plot_data$metric)) -> metric_plots[[i]]
    
    if(i == length(metrics)){
      metric_plots[[i+1]] <- get_legend(metric_plots[[i]] + 
                                          theme(legend.position = 'top', 
                                                legend.title = element_blank())) 
        
    }
  
  }
  
  n <- length(metric_plots)
  nCol <- floor(sqrt(n))
  
  dir.create(directory, recursive = T) 
  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
  all_plots <- do.call("grid.arrange", c(metric_plots, ncol=nCol))
  print(all_plots)
  dev.off()
  
}

# plots for producing a common scale of assessment criteria ----
function(plot_data,
         metrics = c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Psd', 'Pdispersion', 'Pr2'),
         targets = list(Armse    = c(0, -10, 10),
                        Amae     = c(0, -10, 10), 
                        Dintercept  = c(0, -5, 5), 
                        Dslope      = c(1,  0, 2), 
                        Dpearson    = c(1,  0, 1), 
                        Dspearman   = c(1,  0, 1), 
                        Psd         = c(0,  0, 10), 
                        Pdispersion = c(1,  0, 2), 
                        Pr2         = c(1,  0, 1)), 
         levels = c('glm', 'gam', 'gbm', 'rf'), 
         colours = brewer.pal(4, 'Dark2')
         ){
  
  require(pcaMethods)
  
  plot_data = nested_assessments$data[[1]]# remove
  
  # check relationships between metric groupings
  Accuracy  <- plot_data[c('Armse', 'Amae')]
  Discrim   <- plot_data[c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman')]
  Precision <- plot_data[c('Psd', 'Pdispersion', 'Pr2')]
  
  fix_outliers <- function(x, level = 0.01){   
    x[which(x > quantile(x, 1-level, na.rm = T))] <- quantile(x, 1-level, na.rm = T)
    x[which(x < quantile(x, level, na.rm = T))] <- quantile(x, level, na.rm = T)
    x[!is.finite(x)] <- NA
    x
  }

  
  Accuracy_1 <- data.frame(sapply(Accuracy,  fix_outliers))
  #pairs(Accuracy_1)
  Discrim_1 <- data.frame(sapply(Discrim,  fix_outliers))
  #pairs(Discrim_1)
  Precision_1 <- data.frame(sapply(Precision, function(x) fix_outliers(x, 0.025)))
  #pairs(Precision_1)
  
  # rescale values
  rescale_01 <- function(x){(x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
  Accuracy_2 <- data.frame(sapply(Accuracy_1, rescale_01))
  #pairs(Accuracy_2)
  Discrim_1$Dslope <- 1-Discrim_1$Dslope # target is 1 without boundary
  Discrim_2 <- data.frame(sapply(Discrim_1, rescale_01))
  #pairs(Discrim_2)
  Precision_1$Pdispersion <- 1-Precision_1$Pdispersion # target is 1 without boundary
  Precision_2 <- data.frame(sapply(Precision_1, rescale_01))
  #pairs(Precision_2)

  # invert range where necessary so that high numbers is always good
  invert_range <- function(x){( max(x, na.rm = T) + min(x, na.rm = T) ) - x}
  
  Accuracy_2$Armse <- invert_range(Accuracy_2$Armse)
  Accuracy_2$Amae <- invert_range(Accuracy_2$Amae)
  
  Discrim_2$Dintercept <- invert_range(Discrim_2$Dintercept)
  Discrim_2$Dslope <- invert_range(Discrim_2$Dslope)

  Precision_2$Psd <- invert_range(Precision_2$Psd)
  Precision_2$Pdispersion <- invert_range(Precision_2$Pdispersion)
  
  # aggregate over metric types 
  accuracy_metrics = names(Accuracy_2)
  discrimination_metrics = names(Discrim_2)
  precision_metrics = names(Precision_2)
  
  plot_data[accuracy_metrics] <- Accuracy_2
  plot_data[discrimination_metrics] <- Discrim_2
  plot_data[precision_metrics] <- Precision_2
  
  plot_data$accuracy       <- rowMeans(plot_data[accuracy_metrics], na.rm = T)
  plot_data$discrimination <- rowMeans(plot_data[discrimination_metrics], na.rm = T)
  plot_data$precision      <- rowMeans(plot_data[precision_metrics], na.rm = T)
  
  metric_plot_data <- plot_data %>%     
    gather(., key = metric, value = value,  c(accuracy, discrimination, precision)) %>% 
    group_by(species_name,
             dataset,
             cross_validation,
             abundance_response,
             fitted_model,
             metric) %>% 
    do(metric_median = median(.$value, na.rm = T)) %>% 
    unnest(c(metric_median))
  
  metric_plot_data$fitted_model <- factor(as.factor(metric_plot_data$fitted_model), levels)
  
  ggplot(data = metric_plot_data) + 
    geom_boxplot(aes(x = abundance_response, 
                         y = metric_median, 
                         group = paste0(abundance_response, '  ', fitted_model),
                         fill = fitted_model), 
                     outlier.shape = NA, lwd = 0.5, fatten = 2) + 
    facet_wrap(~dataset + metric) + 
    theme_classic() + 
    xlab(NULL) + 
    scale_size_continuous(range = c(0.5, 2)) + 
    scale_fill_manual(values = colours) + 
    scale_colour_manual(values = colours) + 
    theme(aspect.ratio = 0.75, legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + 
    ylab('relative performance')

}











