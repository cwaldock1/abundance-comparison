# functions to produce different figures
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  require(ggplot2)
  require(gridExtra)
  require(grid)
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}
# function to fix outliers in plotting ----
fix_outliers <- function(x, level = 0.01, lower = T, upper = T){   
  if(upper == T){x[which(x > quantile(x, 1-level, na.rm = T))] <- quantile(x, 1-level, na.rm = T)}
  if(lower == T){x[which(x < quantile(x, level, na.rm = T))]   <- quantile(x, level, na.rm = T)}
  x[!is.finite(x)] <- NA
  x
}

# all_model_plots: creates plots for each metric -----

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
    
    if(length(is.finite(metric_plot_data$metrics))!=0){ metric_plot_data$metrics[!is.finite(metric_plot_data$metrics)] <- NA}
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


# all_model_plots_v2: creates large plot across all models ---- 

all_model_plots_v2 <- function(plot_data, # object after running abundance_assesment_metrics and add_family_plot_column (to be renamed)
                               height = 10, 
                               width = 8,
                               outlier_quantile = 0.05,
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
                            directory, 
                            name){
  
  plots_all <- list()
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      dplyr::select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) %>% 
      mutate(metrics = as.numeric(data.frame(plot_data)[,metrics[j]]))
    
    # set boundaries
    lower = ifelse(metrics[j] %in% c('Dintercept'), T, F)
    upper = ifelse(metrics[j] %in% c('Dpearson', 'Dspearman', 'Pr2', 'Dintercept'), F, T)
    t_v = targets[j][[1]]
    decreasing = ifelse(metrics[j] %in% c('Armse', 'Amae', 'Dintercept', 'Psd'), T, F)
    
    # fix extreme values based on boundaries
    metric_plot_data$metrics <- fix_outliers(metric_plot_data$metrics, outlier_quantile, lower = lower, upper = upper)
    
    # calculate median value
    metric_plot_data <- metric_plot_data %>% 
      group_by(plot_level) %>% 
      nest() %>% 
      mutate(median_metric = purrr::map(data, ~median(.$metrics, na.rm = T)), 
             upr_metric    = purrr::map(data, ~as.numeric(quantile(.$metrics, 0.75, na.rm = T))),
             lwr_metric    = purrr::map(data, ~as.numeric(quantile(.$metrics, 0.25, na.rm = T)))) %>% 
      unnest()
    metric_plot_data$plot_level <- factor(metric_plot_data$plot_level, levels = unique(metric_plot_data$plot_level[order(metric_plot_data$median_metric, decreasing=decreasing)]))
    
    # remove unneccessary data
    metric_plot_data <- metric_plot_data %>% dplyr::select(plot_level, median_metric, lwr_metric, upr_metric, fitted_model) %>% unique()
    
    # produce baseplots 
      p <- ggplot(data = metric_plot_data) + 
      geom_point(aes(x = plot_level, y = median_metric, fill = fitted_model, colour = fitted_model)) + 
      geom_linerange(aes(x = plot_level, ymin = lwr_metric, ymax = upr_metric, fill = fitted_model, colour = fitted_model)) + 
      geom_hline(data = data.frame(t_v=t_v), aes(yintercept = t_v[1]), lty=2, col='black') + 
      ylab(metrics[j]) + 
      theme_bw() + 
      theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1), 
            aspect.ratio = 0.5, 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.position = 'none', 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.x = element_blank(), 
            plot.margin = unit(c(1,10,1,1), "lines")) + 
      scale_fill_manual(values = colours) + 
      scale_colour_manual(values = colours) 
      
      # produce annotations
      y_range <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
      y_values <- seq(from = min(y_range, na.rm  = T), to = max(y_range, na.rm  = T), length.out = 10)
      y_values <- y_values[10:5]#ifelse(rep(metrics[j] %in% c('Armse', 'Amae', 'Dintercept', 'Psd', 'Pdispersion'),5), y_values[1:5], rev(y_values[5:10]))
      labels <- as.character(rev(levels(metric_plot_data$plot_level)))[1:5]
      for(i in 1:5){
      p <- p + annotation_custom(grob = textGrob(label = gsub('_','',gsub('.','-',labels[i], fixed=T),fixed=T), 
                                                 hjust = 0, 
                                                 gp = gpar(cex = 0.8)),
        ymin = y_values[i],      # Vertical position of the textGrob
        ymax = y_values[i],
        xmin = length(unique(metric_plot_data$plot_level)) + 2,         # Note: The grobs are positioned outside the plot area
        xmax =  length(unique(metric_plot_data$plot_level)) + 2)
      }
      
      # remove panel so that annotations can be viewed
    gt <- ggplot_gtable(ggplot_build(p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    plots_all[[j]] <- gt
    
    # grab the legend
    if(j == length(metrics)){
      plots_all[[j+1]] <- get_legend(ggplot(data = metric_plot_data) + 
                                          geom_point(aes(x = NA, y = NA, fill = fitted_model, colour = fitted_model)) +
                                          geom_line(aes(x = NA, y = NA, fill = fitted_model, colour = fitted_model)) +
                                          scale_fill_manual(values = colours) + 
                                          scale_colour_manual(values = colours) +
                                          theme(legend.title = element_blank(), 
                                                legend.justification = 'centre'))
      
    }
    
    
  }
  
  # combine plots

  dir.create(directory, recursive = T) 

  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
  
  all_plots <- do.call("grid.arrange", c(plots_all, ncol=2))
  print(all_plots)
  dev.off()
  
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
    
    #  if(sum(metric_plot_data$metric_median > t_v2[3]) != 0){
    #  metric_plot_data$metric_median[which(metric_plot_data$metric_median > t_v2[3])] <- t_v2[3]}
    #  if(sum(metric_plot_data$metric_median < t_v2[2]) != 0){
    #  metric_plot_data$metric_median[which(metric_plot_data$metric_median < t_v2[2])] <- t_v2[2]}
    #  metric_plot_data <- na.omit(metric_plot_data)
    #  
    #  # set lower ylim
    #  if(min(metric_plot_data$metric_median,na.rm=T) > t_v2[2]){y_lower = min(metric_plot_data$metric_median,na.rm=T)-0.1}else{y_lower <- t_v2[2]}
    #  if(max(metric_plot_data$metric_median,na.rm=T) < t_v2[3]){y_upper = max(metric_plot_data$metric_median,na.rm=T)-0.1}else{y_upper <- t_v2[3]}
    #  
    #  if(y_lower > t_v2[1]){y_lower <- t_v2[1]}
    #  if(y_upper < t_v2[1]){y_upper <- t_v2[1]}
    
    
    if(response == 'all'){
      
      # aggregate to median and quantiles
      metric_plot_data_all <- metric_plot_data %>% 
        group_by(dataset,
                 cross_validation,
                 abundance_response,
                 fitted_model) %>% 
        nest() %>% 
        mutate(species_median = purrr::map(data, ~median(.$metric_median, na.rm = T)), 
               species_upr    = purrr::map(data, ~as.numeric(quantile(.$metric_median, 0.75, na.rm = T))),
               species_lwr    = purrr::map(data, ~as.numeric(quantile(.$metric_median, 0.25, na.rm = T)))) %>% 
        unnest()
      
      # remove unneccessary data
      metric_plot_data_all <- metric_plot_data_all %>% dplyr::select(species_median, 
                                                                     species_upr, 
                                                                     species_lwr, 
                                                                     dataset,
                                                                     cross_validation,
                                                                     abundance_response,
                                                                     fitted_model) %>% 
        unique()
      
      
      p <- ggplot(data = metric_plot_data_all)
      
      p + geom_point(aes(x = abundance_response, 
                           y = species_median, 
                           group = paste0(abundance_response, '  ', fitted_model),
                           col = fitted_model), 
                     position=position_dodge(width=0.75), size = 3) + 
        geom_linerange(aes(x = abundance_response, 
                           ymin = species_lwr,
                           ymax = species_upr,
                           group = paste0(abundance_response, '  ', fitted_model),
                           col = fitted_model), 
                       position=position_dodge(width=0.75)) -> p2 
    }
    
    p2 + geom_hline(data = data.frame(t_v2=t_v2), aes(yintercept = t_v2[1]), lty=2, col='black') + 
      facet_wrap(~dataset) + 
      theme_classic() + 
      xlab(NULL) + 
      ylab(metrics[i]) + 
      #ylim(y_lower, y_upper) +
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

# function for producing a common scale of assessment criteria ----

aggregate_metrics <- function(plot_data,
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
         levels = c('glm', 'gam', 'gbm', 'rf')
         ){
  
  # get metrics
  metrics_data <- na.omit(plot_data[metrics])
  
  # match targets to boundarise to express as abosulute difference from 0 (where 0 is perfect)
  metrics_data[,which(names(metrics_data) %in% c('Dslope', 'Pdispersion'))] <- abs(metrics_data[,which(names(metrics_data) %in% c('Dslope', 'Pdispersion'))] - 1)
  metrics_data[,which(names(metrics_data) %in% c('Dintercept'))] <- abs(metrics_data[,which(names(metrics_data) %in% c('Dintercept'))])
  
  # rank order and rescale
  rescale_01 <- function(x){(x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
  metrics_data <- data.frame(sapply(metrics_data, function(x) rescale_01(rank(x))))
  
  # range inversion
  invert_range <- function(x){(max(x, na.rm = T) + min(x, na.rm = T) ) - x}
  metrics_data[which(names(metrics_data) %in% c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Psd', 'Pdispersion'))] <- 
    sapply(metrics_data[which(names(metrics_data) %in% c('Armse', 'Amae', 'Dintercept', 'Dslope', 'Psd', 'Pdispersion'))], function(x) invert_range(x))
  
  # create columns for plotting
  if(length(which(names(metrics_data) %in% c('Armse', 'Amae'))) == 1){
    accuracy <- metrics_data[,which(names(metrics_data) %in% c('Armse', 'Amae'))]}else{
    accuracy <- rowSums(metrics_data[,which(names(metrics_data) %in% c('Armse', 'Amae'))])/length(which(names(metrics_data) %in% c('Armse', 'Amae')))}
  
  if(length(which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))) == 1){
    discrimination <- metrics_data[,which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))]}else{
    discrimination <- rowSums(metrics_data[,which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))])/length(which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman')))}
  
  if(length(which(names(metrics_data) %in% c("Pdispersion", "Pr2"))) == 1){precision <- metrics_data[,which(names(metrics_data) %in% c("Pdispersion", "Pr2"))]}else{
    precision <- rowSums(metrics_data[,which(names(metrics_data) %in% c("Pdispersion", "Pr2"))])/length(which(names(metrics_data) %in% c("Pdispersion", "Pr2")))}
  
  # aggregate the evaluation metrics
  aggregated_evaluation_metrics <- rowSums(cbind(accuracy, discrimination, precision))/3
  
  # index for ANY NAs
  NA_index <- rowSums(sapply(plot_data[metrics], is.na))>0
  
  plot_data$accuracy                      <- NA
  plot_data$discrimination                <- NA
  plot_data$precision                     <- NA
  plot_data$aggregated_evaluation_metrics <- NA
  
  if(sum(NA_index)!=0){
  plot_data$accuracy[-which(NA_index)]                      <- accuracy
  plot_data$discrimination[-which(NA_index)]                <- discrimination
  plot_data$precision[-which(NA_index)]                     <- precision
  plot_data$aggregated_evaluation_metrics[-which(NA_index)] <- aggregated_evaluation_metrics
  }else{
    plot_data$accuracy                      <- accuracy
    plot_data$discrimination                <- discrimination
    plot_data$precision                     <- precision
    plot_data$aggregated_evaluation_metrics <- aggregated_evaluation_metrics
  }
  
  return(plot_data)
  
}

# function to plot aggregated metrics ----

plot_all_aggregated <- function(plot_data,
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
         levels = c('glm', 'gam', 'gbm', 'rf')){
  
  metric_plot_data <- plot_data %>%     
    gather(., key = metric, value = value,  c(accuracy, discrimination, precision)) %>% 
    group_by(species_name,
             dataset,
             cross_validation,
             abundance_response,
             fitted_model,
             metric) %>% 
    # take the mean across metrics because of rescaling
    do(metric_mean = mean(.$value, na.rm = T)) %>% 
    unnest(c(metric_mean)) %>% 
    group_by(dataset,
             cross_validation,
             abundance_response,
             fitted_model,
             metric) %>% 
    do(species_median = median(.$metric_mean, na.rm = T), 
       species_upr    = quantile(.$metric_mean, 0.95, na.rm = T), 
       species_lwr    = quantile(.$metric_mean, 0.05, na.rm = T)) %>% 
    unnest(c(species_median, species_upr, species_lwr))
    
  
  metric_plot_data$fitted_model <- factor(as.factor(metric_plot_data$fitted_model), levels)
  
  plot_1 <- ggplot(data = metric_plot_data) + 
    geom_point(aes(x = abundance_response, 
                   y = species_median, 
                   group = paste0(abundance_response, '  ', fitted_model),
                   col = fitted_model),
               size = 3,
               position=position_dodge(width=0.75)) + 
    geom_linerange(aes(x = abundance_response, 
                       ymin = species_lwr,
                       ymax = species_upr,
                       group = paste0(abundance_response, '  ', fitted_model),
                       col = fitted_model), 
                   position=position_dodge(width=0.75)) +
    #geom_boxplot(aes(x = abundance_response, 
    #                     y = metric_mean, 
    #                     group = paste0(abundance_response, '  ', fitted_model),
    #                     fill = fitted_model), 
    #                 outlier.shape = NA, lwd = 0.5, fatten = 2) + 
    facet_wrap(~dataset + metric) + 
    theme_classic() + 
    xlab(NULL) + 
    scale_size_continuous(range = c(0.5, 2)) + 
    scale_fill_manual('model', values = colours) + 
    scale_colour_manual('model', values = colours) + 
    theme(aspect.ratio = 0.75, 
          legend.position = 'bottom', 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + 
    ylab('relative performance') + 
    ylim(0,1)

dir.create(directory, recursive = T) 
pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
print(plot_1)
dev.off()

}

# function to plot aggregated metrics across scales ----

plot_all_aggregated_spatial_scale <- function(plot_data,
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
                                levels = c('glm', 'gam', 'gbm', 'rf')){
  
  plot_data_DT <- data.table(plot_data %>% 
                               gather(., key = metric, value = value,  c(accuracy, discrimination, precision)) %>% 
                               filter(Evaluation_message == 'none') %>% 
                               group_by(dataset, spatial_scale, abundance_response, fitted_model, metric, species_name))
  
  plot_data_DT[,metric_mean := mean(value, na.rm = T),
               by = c('dataset',
                      'spatial_scale',
                      'abundance_response',
                      'fitted_model',
                      'metric',
                      'species_name')]
  
  plot_data_DT <- unique(plot_data_DT[,c('dataset',
                                         'spatial_scale',
                                         'abundance_response',
                                         'fitted_model',
                                         'metric',
                                         'species_name', 
                                         'metric_mean'),])
  
  key <- c('dataset',
           'spatial_scale',
           'abundance_response',
           'fitted_model',
           'metric')
  plot_data_DT[,species_median := median(metric_mean, na.rm = T), by = key]
  plot_data_DT[,species_upr := quantile(metric_mean, 0.95, na.rm = T),by = key]
  plot_data_DT[,species_lwr := quantile(metric_mean, 0.05, na.rm = T),by = key]
  
  metric_plot_data <- plot_data_DT %>% 
    select(spatial_scale,
           metric,
           dataset,
           abundance_response,
           fitted_model, 
           species_median, 
           species_upr, 
           species_lwr) %>% 
    unique() %>% 
    as_data_frame()
  
  metric_plot_data$fitted_model <- factor(as.factor(metric_plot_data$fitted_model), levels)
  
  metric_plot_data$spatial_scale <- as.numeric(gsub('spatial_scale_', '', metric_plot_data$spatial_scale))
  
  
  bbs_plot <- ggplot(data = metric_plot_data %>% filter(dataset == 'bbs')) + 
    geom_point(aes(x = spatial_scale, 
                   y = species_median, 
                   col = fitted_model), 
               size = 2) +
    geom_line(aes(x = spatial_scale, 
                   y = species_median, 
                   col = fitted_model), 
              size = 1.5)  + 
    facet_wrap(~ metric + abundance_response) + 
    theme_classic() + 
    xlab('spatial scale (degrees)') + 
    scale_fill_manual('model', values = colours) + 
    scale_colour_manual('model', values = colours) + 
    theme(aspect.ratio = 0.75, 
          legend.position = 'bottom') +
    ylab('relative performance') + 
    ylim(0,1)
  
  rls_plot <- ggplot(data = metric_plot_data %>% filter(dataset == 'rls')) + 
    geom_point(aes(x = spatial_scale, 
                   y = species_median, 
                   col = fitted_model), 
               size = 2) +
    geom_line(aes(x = spatial_scale, 
                  y = species_median, 
                  col = fitted_model), 
              size = 1.5)  + 
    facet_wrap(~ metric + abundance_response) + 
    theme_classic() + 
    xlab('spatial scale (degrees)') + 
    scale_fill_manual('model', values = colours) + 
    scale_colour_manual('model', values = colours) + 
    theme(aspect.ratio = 0.75, 
          legend.position = 'bottom') +
    ylab('relative performance') + 
    ylim(0,1)
  
  
  dir.create(directory, recursive = T) 
  pdf(file = paste0(directory, '/', 'bbs', '_', name,'.pdf'), width = width, height = height)
  print(bbs_plot)
  dev.off()
  pdf(file = paste0(directory, '/', 'rls', '_', name,'.pdf'), width = width, height = height)
  print(rls_plot)
  dev.off()
  
  
}


# function to plot model rank counts ----

rank_plots <- function(plot_data,
                       directory,
                       name,
                       width, 
                       height,
                       colours,
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
                       levels = c('glm', 'gam', 'gbm', 'rf')){
  
  # create function to apply
  get_models <- function(input){
    
    value <- ifelse(unique(input$metrics) %in% c('Armse', 'Amae', 'Dintercept', 'Psd', 'Pdispersion'), 
                        which.min(input$value), 
                        which.max(input$value))
    
    return_model <- input$plot_level[value]
    
    return(return_model)
    
  }
  
  
  model_count <- plot_data %>% 
    gather(., key = 'metrics', value = 'value', metrics) %>% 
    group_by(dataset, cross_validation_2, species_name, metrics) %>%
    do(best_model = get_models(.)) %>% 
    unnest(cols = c(best_model)) %>% 
    na.omit()
  
  model_count$fitted_model <- sub('\\..*', '', model_count$best_model)
  model_count$fitted_model <- factor(as.factor(model_count$fitted_model), levels)
  model_count$best_model <- gsub('rf.|glm.|gam.|gbm.','' , model_count$best_model)
  model_count$best_model <- gsub('abun-occ_2stage', 'abun-occ-2stage', gsub('abunocc', 'abun-occ', model_count$best_model))
  model_count$abundance_response <- sub('\\..*', '', model_count$best_model)
  model_count$abundance_response <- plyr::revalue(model_count$abundance_response, c(`abun-occ` = 'ao', `abun` = 'a', `abun-occ-2stage` = 'ao-2'))
  model_count$best_model <- gsub('abun-occ-2stage.|abun-occ.|abun.', '', model_count$best_model)
  
  
  #detach("package:cowplot", unload = TRUE)
  plot <- ggplot(data = model_count) + 
    geom_bar(aes(x=gsub('_','',gsub('.','-',best_model, fixed = T)), fill = fitted_model), col = 'black', lwd = 0.5) + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust= 0.5, size = 8), 
          panel.border = element_rect(linetype = 1), 
          panel.background = element_rect(fill = 'grey95'), 
          strip.background = element_blank(), 
          aspect.ratio=1, 
          strip.text.x=element_text(margin=margin(t=1, r=5, b=5, l=5), size = 10), 
          strip.text.y=element_text(margin=margin(t=1, r=1, b=1, l=1), size = 8)) + 
    facet_grid(metrics~fitted_model + abundance_response, scale = 'free') + 
    xlab('') + 
    ylab('model count') + 
    scale_fill_manual(values = colours, name = '') + 
    scale_colour_manual(values = colours, name = '')
  
  dir.create(directory, recursive = T) 
  
  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
  print(plot)
  dev.off()
  
}


# function to plot best-models distribution of metrics ---- 

spp_best_assessment_metrics <- function(plot_data, # object after running abundance_assesment_metrics and add_family_plot_column (to be renamed)
                                        height = 10, 
                                        width = 8,
                                        outlier_quantile = 0.05,
                                        metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2'), 
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
                                        colours = c('gray20', 'gray80'),
                                        directory, 
                                        name){
  
  plots_all <- list()
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      dplyr::select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) %>% 
      mutate(metrics = as.numeric(data.frame(plot_data)[,metrics[j]]))
    
    # set boundaries
    lower = ifelse(metrics[j] %in% c('Dintercept'), T, F)
    upper = ifelse(metrics[j] %in% c('Dpearson', 'Dspearman', 'Pr2', 'Dintercept'), F, T)
    t_v = targets[j][[1]]
    decreasing = ifelse(metrics[j] %in% c('Armse', 'Amae', 'Dintercept', 'Psd'), T, F)
    
    # fix extreme values based on boundaries
    metric_plot_data$metrics <- fix_outliers(metric_plot_data$metrics, outlier_quantile, lower = lower, upper = upper)
    
    # make plots
    require(ggridges)
    plots_all[[j]] <- ggplot(data = metric_plot_data) + 
      geom_density_ridges(aes(y = dataset, x = metrics, fill = dataset, col = dataset), alpha = 0.5, 
                          jittered_points = T, 
                          position = 'raincloud', 
                          point_alpha = 0.2, 
                          point_stroke = 0, 
                          point_size = 2) + 
      theme_bw() + 
      theme(aspect.ratio = 0.15, 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.line.x = element_line(),
            legend.position = 'none', 
            panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size = 10)) + 
      ggtitle(metrics[j]) + 
      ylab(NULL) + 
      scale_fill_manual(values = colours) + 
      scale_colour_manual(values = colours) + 
      xlim(c(min(metric_plot_data$metrics), max(metric_plot_data$metrics)))
    
    labels <- metric_plot_data %>% 
      group_by(dataset) %>% 
      do(median_value = signif(median(.$metrics, na.rm = T),2)) %>% 
      unnest() %>% .$median_value
    
    # create annotations with median values
    for(i in 1:length(labels)){
      plots_all[[j]] <- plots_all[[j]] +
        annotation_custom(grob = textGrob(label = labels[i], 
                                                 hjust = -0.2, 
                                                 gp = gpar(cex = 1, 
                                                           col = ifelse(i == 1, colours[1], colours[2])),                                          ),
                                 ymin = ifelse(i==1, 1, 2),      # Vertical position of the textGrob
                                 ymax = ifelse(i==1, 1, 2),
                                 xmin = max(metric_plot_data$metrics),         # Note: The grobs are positioned outside the plot area
                                 xmax = max(metric_plot_data$metrics))
    }
    # remove panel so that annotations can be viewed
    gt <- ggplot_gtable(ggplot_build(plots_all[[j]]))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    plots_all[[j]] <- gt
    
    
    
  }
  
  # combine plots
  
  dir.create(directory, recursive = T) 
  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height, useDingbats = F)
  all_plots <- do.call("grid.arrange", c(plots_all, ncol=1))
  print(all_plots)
  dev.off()
  
}

# function to plot best-models distribution of metrics across scales ---- 

spp_best_assessment_metrics_scale <- function(plot_data, # object after running abundance_assesment_metrics and add_family_plot_column (to be renamed)
                                              selected_dataset, 
                                              height = 10, 
                                        width = 8,
                                        outlier_quantile = 0.05,
                                        metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2'), 
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
                                        colours = c('gray20', 'gray80'),
                                        directory, 
                                        name){
  
  plots_all <- list()
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      dplyr::select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) 
    
    # set boundaries
    lower = ifelse(metrics[j] %in% c('Dintercept'), T, F)
    upper = ifelse(metrics[j] %in% c('Dpearson', 'Dspearman', 'Pr2', 'Dintercept'), F, T)
    t_v = targets[j][[1]]
    decreasing = ifelse(metrics[j] %in% c('Armse', 'Amae', 'Dintercept', 'Psd'), T, F)
    
    # remove outliers
    metric_plot_data$metrics <- fix_outliers(as.numeric(metric_plot_data[metrics[j]][[1]]), outlier_quantile, lower = lower, upper = upper)
    
    # metric_plot_data <- metric_plot_data %>%  
    #   group_by(spatial_scale, cross_validation, dataset) %>% 
    #   mutate(metrics_median = median(metrics), 
    #          metrics_upr = quantile(metrics, 0.95), 
    #          metrics_lwr = quantile(metrics, 0.05)) %>% 
    #   select(spatial_scale, cross_validation, dataset, metrics_median, metrics_upr, metrics_lwr) %>% 
    #   unique()
    
    # convert spatial scale to numeric
    metric_plot_data$spatial_scale <- gsub('spatial_scale_', '', metric_plot_data$spatial_scale)
    
    metric_plot_data$spatial_scale <- factor(metric_plot_data$spatial_scale, levels = c(0.1, 1, 5, 10, 20, 35, 50))
    
    # make plots
    require(viridis)
   plots_all[[j]] <- 
      ggplot(data = metric_plot_data) + 
      geom_boxplot(aes(x = spatial_scale, y = metrics, fill = dataset, group = paste0(dataset, spatial_scale)), lwd = 0.5) +
     geom_hline(aes(yintercept = t_v[1]), lty = 2, colour = 'gray50') +
      theme_bw() + 
        theme(aspect.ratio = 0.15, 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.line.x = element_line(),
            legend.position = 'none', 
            panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            #axis.text.y = element_blank(), 
            #axis.ticks.y = element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(size = 10)) + 
      #ggtitle(metrics[j]) + 
      ylab(NULL) + 
      scale_fill_manual(values = colours)
    
  }
  
  # combine plots
  
  dir.create(directory, recursive = T) 
  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height, useDingbats = F)
  all_plots <- do.call("grid.arrange", c(plots_all, ncol=1))
  print(all_plots)
  dev.off()
  
  pdf(file = paste0(directory, '/', 'legend','.pdf'), width = width, height = height, useDingbats = F)
  print(plot(get_legend(ggplot(data = metric_plot_data) + 
                    geom_boxplot(aes(x = spatial_scale, y = metrics, fill = dataset, group = paste0(dataset, spatial_scale)), lwd = 0.5) + 
                    scale_fill_manual(values = colours))))
  dev.off()
  
}

