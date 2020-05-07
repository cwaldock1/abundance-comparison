# functions for evaluating metrics


# bind_results: binds together results folders of interest ----

bind_results <- function(file_list){
  
  predictions <- list()
  for(i in 1:length(file_list)){
    load(file_list[i])
    extracted_predictions$abundance_response <- gsub('predictions_', '', str_split(file_list[i], '/')[[1]][3])
    predictions[[i]] <- extracted_predictions
  }
  all_files_bind <- do.call(rbind, predictions)
  return(all_files_bind)
  
}


# clean data managing data groupings for plots ----

clean_levels <- function(data_input){
  
  require(tidyverse)
  
  # edit transformations
  data_input$transformation[is.na(data_input$transformation)] <- ''
  data_input$transformation <- ifelse(data_input$transformation == '', '', paste0('-', data_input$transformation))
  
  # edit family groups
  data_input$family <- recode(data_input$family, 
                             discrete = "multinomial", 
                             discrete_multinomial = "multinomial", 
                             continuous_gaussian = 'gaussian',
                             continuous = 'gaussian', 
                             continuous_poisson = 'poisson', 
                             zip = 'poisson')
  data_input$family[is.na(data_input$family)] <- 'gaussian'
  
  # edit zero-inflations
  data_input$zi <- ifelse(data_input$zi, '-ZI', '')
  data_input$zi[is.na(data_input$zi)] <- ''
  
  data_input$family_grouped <- paste(data_input$family, data_input$zi, data_input$transformation, sep = '')
  
  data_input$family_grouped_simple <- recode(data_input$family_grouped, 
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
  
  # recode abundance response
  data_input$abundance_response_simple <- recode(data_input$abundance_response,
                                                            abun = "a", 
                                                            abunocc_2stage = "a2", 
                                                            abunocc = "ao")
  
  # produce levels for plots
  data_input$plot_level <- paste(data_input$fitted_model, 
                                 data_input$abundance_response, 
                                 data_input$family_grouped_simple, sep = '.')
  
  data_input <- data_input %>% 
    select(-c(transformation, zi, only_abundance)) %>% 
    mutate(n_absence = as.character(.$n_absence), 
           n_boot_absence = as.character(.$n_boot_absence), 
           n_abundance = as.character(.$n_abundance))
  
  data_input <- data_input[names(sort(sapply(data_input, class)))]
    
  return(data_input)
  
}


# abundance_assessment_metrics: calculates assessment metrics ----

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

# write_plots: creates plots for each metric -----

# 
write_plots <- function(plot_data, # object after running abundance_assesment_metrics and add_family_plot_column (to be renamed)
                        metrics, 
                        targets, 
                        levels,
                        directory){
  for(j in 1:length(metrics)){
    
    # set up levels
    if(j == 1){plot_data$fitted_model <- factor(as.factor(plot_data$fitted_model), levels)}
    print(j)
    # set up data for plotting in loop across all metrics
    metric_plot_data <- plot_data %>% 
      select(colnames(.)[-which(colnames(.)%in%metrics)], metrics[j]) %>% 
      mutate(metrics = as.numeric(.[,metrics[j]][[1]]))
    
    t_v <- targets[[j]]
    
    metric_plot_data$metrics[!is.finite(metric_plot_data$metrics)] <- NA
    metric_plot_data$metrics[which(metric_plot_data$metrics > t_v[3])] <- t_v[3]
    metric_plot_data$metrics[which(metric_plot_data$metrics < t_v[2])] <- t_v[2]
    metric_plot_data <- na.omit(metric_plot_data)
    
    # set lower ylim
    y_lower <- min(c(min(metric_plot_data$metrics,na.rm=T), t_v[2]))-0.1
    y_upper <- max(c(max(metric_plot_data$metrics,na.rm=T), t_v[3]))+0.1
    
    plots <- ggplot(data = metric_plot_data) + 
      geom_boxplot(aes(x = abundance_class, y = metrics)) + 
      geom_hline(aes(yintercept=t_v[1]), col = 'red') + 
      ylab(metrics[j]) + 
      ylim(y_lower, y_upper) + 
      facet_grid(fitted_model ~ family_plot, scales = 'free', drop = T) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1), 
            aspect.ratio = 1)
    
    dir.create(directory, recursive = T) 
    pdf(file = paste0(directory, '/', metrics[j],'.pdf'), width = 14, height = 10)
    print(plots)
    dev.off()
    
  }
}
