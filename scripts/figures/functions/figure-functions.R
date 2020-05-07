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
                        levels = c('glm', 'gam', 'rf', 'gbm'),
                        directory){
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) %>% 
      mutate(metrics = as.numeric(.[,metrics[j]]))
    
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
    
    plots <-  ggplot() + 
      geom_line(data = metric_plot_data, 
                aes(x = mean_abundance_perc_fact, 
                    y = metrics), alpha = 0.5) + 
      geom_line(data = metric_plot_data, 
                aes(x = frequency_perc_fact, 
                    y = metrics), alpha = 0.5, col = 'red') + 
      stat_smooth(data = metric_plot_data, 
                  aes(x = mean_abundance_perc_fact, 
                      y = metrics), 
                  col = 'black', size =  0.5, se = F) + 
      stat_smooth(data = metric_plot_data, 
                  aes(x = frequency_perc_fact, 
                      y = metrics), 
                  col = 'red', size = 0.5, se = F) + 
      ylab(metrics[j]) + 
      xlab('percentile rank \n (abundance = black; frequency = red)') +
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
