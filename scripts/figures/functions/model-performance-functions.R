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
    decreasing = ifelse(metrics[j] %in% c('Armse', 'Amae', 'Dslope', 'Psd', 'Pdispersion'), T, F)
    
    # fix extreme values based on boundaries
    metric_plot_data$metrics <- fix_outliers(metric_plot_data$metrics, outlier_quantile, lower = lower, upper = upper)
    
    # calculate median value
    metric_plot_data <- metric_plot_data %>% 
      group_by(plot_level) %>% 
      nest() %>% 
      mutate(median_metric = purrr::map(data, ~median(.$metrics, na.rm = T))) %>% 
      unnest()
    metric_plot_data$plot_level <- factor(metric_plot_data$plot_level, levels = unique(metric_plot_data$plot_level[order(metric_plot_data$median_metric, decreasing=decreasing)]))
    
    # produce baseplots 
      p <- ggplot(data = metric_plot_data) + 
      geom_boxplot(aes(x = plot_level, y = metrics, fill = fitted_model, colour = fitted_model), alpha = 0.5) + 
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
      
      # produce annotatins
      y_values <- seq(from = min(metric_plot_data$metrics, na.rm  = T), to = max(metric_plot_data$metrics, na.rm  = T), length.out = 10)
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
                                          geom_boxplot(aes(x = NA, y = NA, fill = fitted_model, colour = fitted_model), alpha = 0.5) +
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
  metrics_data <- plot_data[metrics]
  
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
  if(length(which(names(metrics_data) %in% c('Armse', 'Amae'))) == 1){plot_data$accuracy <- metrics_data[,which(names(metrics_data) %in% c('Armse', 'Amae'))]}else{
    plot_data$accuracy <- rowSums(metrics_data[,which(names(metrics_data) %in% c('Armse', 'Amae'))])/ncol(metrics_data)}
  
  if(length(which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))) == 1){plot_data$discrimination <- metrics_data[,which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))]}else{
    plot_data$discrimination <- rowSums(metrics_data[,which(names(metrics_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))])/ncol(metrics_data)}
  
  if(length(which(names(metrics_data) %in% c("Pdispersion", "Pr2"))) == 1){plot_data$precision <- metrics_data[,which(names(metrics_data) %in% c("Pdispersion", "Pr2"))]}else{
    plot_data$precision <- rowSums(metrics_data[,which(names(metrics_data) %in% c("Pdispersion", "Pr2"))])/ncol(metrics_data)}
  
  # aggregate the evaluation metrics
  plot_data$aggregated_evaluation_metrics <- rowSums(metrics_data)/ncol(metrics_data)
  
  ## check relationships between metric groupings
  #Accuracy  <- plot_data[which(names(plot_data) %in% c('Armse', 'Amae'))]
  #Discrim   <- plot_data[which(names(plot_data) %in% c('Dintercept', 'Dslope', 'Dpearson', 'Dspearman'))]
  #Precision <- plot_data[which(names(plot_data) %in% c('Psd', 'Pdispersion', 'Pr2'))]
  
  # Accuracy_1 <- data.frame(sapply(Accuracy,  function(x) fix_outliers(x, 0.025)))
  # #pairs(Accuracy_1)
  # Discrim_1 <- data.frame(sapply(1:ncol(Discrim),  function(x) fix_outliers(Discrim[[x]], 0.025, lower = (names(Discrim) %in% 'Dintercept')[x], upper = (names(Discrim) %in% 'Dintercept')[x])))
  # #pairs(Discrim_1)
  # Precision_1 <- data.frame(sapply(Precision, function(x) fix_outliers(x, 0.025)))
  # #pairs(Precision_1)
  # 
  # # rescale values
  # rescale_01 <- function(x){(x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}
  # Accuracy_2 <- data.frame(sapply(Accuracy_1, rescale_01))
  # #pairs(Accuracy_2)
  # 
  # if(c('Dslope') %in% names(plot_data)){Discrim_1$Dslope <- 1-Discrim_1$Dslope} # target is 1 without boundary
  # Discrim_2 <- data.frame(sapply(Discrim_1, rescale_01))
  # #pairs(Discrim_2)
  # 
  # if(c('Pdispersion') %in% names(plot_data)){Precision_1$Pdispersion <- 1-Precision_1$Pdispersion} # target is 1 without boundary
  # Precision_2 <- data.frame(sapply(Precision_1, rescale_01))
  # #pairs(Precision_2)
  # 
  # # invert range where necessary so that high numbers is always good
  # invert_range <- function(x){( max(x, na.rm = T) + min(x, na.rm = T) ) - x}
  # 
  # if(c('Armse') %in% names(plot_data)){Accuracy_2$Armse <- invert_range(Accuracy_2$Armse)}
  # if(c('Amae') %in% names(plot_data)){Accuracy_2$Amae <- invert_range(Accuracy_2$Amae)}
  # 
  # # if(c('Dintercept') %in% names(plot_data)){Discrim_2$Dintercept <- invert_range(Discrim_2$Dintercept)}
  # if(c('Dslope') %in% names(plot_data)){Discrim_2$Dslope <- invert_range(Discrim_2$Dslope)}
  # 
  # if(c('Psd') %in% names(plot_data)){Precision_2$Psd <- invert_range(Precision_2$Psd)}
  # if(c('Pdispersion') %in% names(plot_data)){Precision_2$Pdispersion <- invert_range(Precision_2$Pdispersion)}
  # 
  # # check directions when writing function
  # #pairs(cbind(Accuracy_2, Discrim_2, Precision_2))
  # 
  # # aggregate over metric types 
  # accuracy_metrics = names(Accuracy_2)
  # discrimination_metrics = names(Discrim_2)
  # precision_metrics = names(Precision_2)
  # 
  # plot_data[accuracy_metrics] <- Accuracy_2
  # plot_data[discrimination_metrics] <- Discrim_2
  # plot_data[precision_metrics] <- Precision_2
  # 
  # plot_data$accuracy       <- rowMeans(plot_data[accuracy_metrics], na.rm = T)
  # plot_data$discrimination <- rowMeans(plot_data[discrimination_metrics], na.rm = T)
  # plot_data$precision      <- rowMeans(plot_data[precision_metrics], na.rm = T)
  
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
    unnest(c(metric_mean))
  
  metric_plot_data$fitted_model <- factor(as.factor(metric_plot_data$fitted_model), levels)
  
  plot_1 <- ggplot(data = metric_plot_data) + 
    geom_boxplot(aes(x = abundance_response, 
                         y = metric_mean, 
                         group = paste0(abundance_response, '  ', fitted_model),
                         fill = fitted_model), 
                     outlier.shape = NA, lwd = 0.5, fatten = 2) + 
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
    unnest(cols = c(best_model))
  
  model_count$fitted_model <- sub('\\..*', '', model_count$best_model)
  model_count$fitted_model <- factor(as.factor(model_count$fitted_model), levels)
  
  #detach("package:cowplot", unload = TRUE)
  plot <- ggplot(data = model_count) + 
    geom_bar(aes(x=gsub('_','',gsub('.','-',best_model, fixed = T)), fill = fitted_model), col = 'black', lwd = 0.5) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), 
          panel.border = element_rect(linetype = 1), 
          panel.background = element_rect(fill = 'grey95'), 
          strip.background = element_blank()) + 
    facet_grid(metrics~fitted_model, scale = 'free') + 
    xlab('') + 
    ylab('model count') + 
    scale_fill_manual(values = colours, name = '') + 
    scale_colour_manual(values = colours, name = '')
  
  dir.create(directory, recursive = T) 
  
  pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
  print(plot)
  dev.off()
  
}


# function to plot comparison between basic and cross-validations for each metric ---- 

plot_cv_comparison <- function(plot_data, 
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
  
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    
    fix_outliers <- function(x, level = 0.01){   
      x[which(x > quantile(x, 1-level, na.rm = T))] <- quantile(x, 1-level, na.rm = T)
      x[which(x < quantile(x, level, na.rm = T))] <- quantile(x, level, na.rm = T)
      x[!is.finite(x)] <- NA
      x
    }
    
  # set up data for plotting in loop across all metrics
  metric_plot_data <- plot_data %>% 
    group_by(fitted_model, family_grouped, species_name, plot_level) %>% 
    # subset to when both species are present in a model
    do(both_present = nrow(.)) %>% unnest() %>% filter(both_present == 2) %>% 
    left_join(., plot_data) %>% 
    dplyr::select(colnames(.)[-which(colnames(.) %in% metrics)], metrics[j]) %>% 
    mutate(metrics = as.numeric(data.frame(.)[,metrics[j]])) #%>% 
    #mutate(cross_validation_factor = ifelse(grepl('oob', .$cross_validation), 'oob', 'non_oob')) %>% 
    #select(species_name, plot_level, cross_validation_factor, metrics) %>% 
    #ungroup() %>% 
    #pivot_wider(., names_from = cross_validation_factor, values_from =  metrics) %>% 
    #unnest() %>% 
    #mutate(RR_metric = log(abs(fix_outliers(.$non_oob, 0.1))) / log(abs(fix_outliers(.$oob, 0.1)))) %>% 
    #left_join(., plot_data %>% dplyr::select(plot_level, fitted_model, abundance_response_simple)) 
  
    
  plots <- ggplot(data = metric_plot_data) +
    geom_boxplot(aes(x = plot_level, y = metrics, group = paste0(plot_level, cross_validation), fill = cross_validation)) + 
    facet_wrap(~fitted_model, scales= 'free') + 
    theme(aspect.ratio = 0.75, 
          legend.position = 'right', 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.y = element_blank()) + 
    scale_fill_manual('', values = colours) + 
    ylab(metrics[[j]])
  
  dir.create(directory, recursive = T) 
  pdf(file = paste0(directory, '/', metrics[j],'.pdf'), width = 14, height = 10)
  print(plots)
  dev.off()
  
  }
  
  }

