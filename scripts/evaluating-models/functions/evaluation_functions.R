# functions for evaluating metrics


# bind_results: binds together results folders of interest ----

bind_results <- function(file_list){
  
  predictions <- list()
  for(i in 1:length(file_list)){
    if( i %% 10 == 0 ) cat(paste("iteration", i, "complete\n"))
    load(file_list[i])
    extracted_predictions$abundance_response <- gsub('predictions_', '', str_split(file_list[i], '/')[[1]][length(str_split(file_list[i], '/')[[1]])-1])
    predictions[[i]] <- extracted_predictions
  }
  all_files_bind <- do.call(rbind, predictions)
  return(all_files_bind)
  
}

bind_results_save <- function(file_list, directory, name){
  
  predictions <- list()
  for(i in 1:length(file_list)){
    if( i %% 50 == 0 ) cat(paste("iteration", signif((i/length(file_list))*100, 3), "% complete\n"))
    load(file_list[i])
    extracted_predictions$abundance_response <- gsub('predictions_', '', str_split(file_list[i], '/')[[1]][length(str_split(file_list[i], '/')[[1]])-1])
    predictions[[i]] <- extracted_predictions
  }
  
  all_files_bind <- do.call(rbind, predictions)
  
  dir.create(directory, recursive = T)
  
  saveRDS(all_files_bind, file = paste0(directory, '/', name, '.rds'))
  
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
    dplyr::select(-c(zi, only_abundance)) %>% 
    mutate(n_absence = as.character(.$n_absence), 
           n_boot_absence = as.character(.$n_boot_absence), 
           n_abundance = as.character(.$n_abundance))
  
  data_input <- data_input[names(sort(sapply(data_input, class)))]
    
  return(data_input)
  
}


# abundance_assessment_metrics: calculates assessment metrics ----

# predictions <- clean_files %>% rowwise() %>% .$validation_predict_mean %>% .[[1]]
# observations <- clean_files %>% rowwise() %>% .$validation_observed_mean %>% .[[1]]
# locations <- clean_files %>% rowwise() %>% .$validation_locations %>% .[[1]]

abundance_assessment_metrics <- function(predictions, observations, locations, scale = NULL){
  
  # check lengths are the same
  if(length(observations) != length(predictions) & nrow(locations) != length(observations)){return(data.frame(Armse = NA, 
                                                                                                              Amae  = NA, 
                                                                                                              Dintercept = NA, 
                                                                                                              Dslope = NA, 
                                                                                                              Dpearson = NA, 
                                                                                                              Dspearman = NA, 
                                                                                                              Psd = NA, 
                                                                                                              Pdispersion = NA, 
                                                                                                              Pr2 = NA, 
                                                                                                              Evaluation_number = NA, 
                                                                                                              Scale = if(is.null(scale)){0}else{scale}, 
                                                                                                              Evaluation_message = 'warning: length of observations and predictions does not match locations'))}
  
  # keep only observations that are abundances in the first input
  to_keep <- observations>0
  observations <- observations[to_keep]
  predictions  <- predictions[to_keep]
  locations    <- locations[to_keep,]
  
  # create dataframe of locations predictions and observations
  if(!is.null(scale)){
    
    # create dataframe
    ddata <- data.frame(observations, predictions, locations)
    
    # aggregate by scale
    rescaled_observations <- ddata %>% 
      mutate(SiteLatitude = plyr::round_any(.$SiteLatitude, scale), 
             SiteLongitude = plyr::round_any(.$SiteLongitude, scale)) %>% 
      group_by(SiteLatitude, SiteLongitude) %>% 
      do(observations = mean(.$observations), 
         predictions = mean(.$predictions)) %>% 
      unnest(c('observations', 'predictions'))
    
    #if(nrow(ddata) == nrow(rescaled_observations) | nrow(rescaled_observations) == 1){return(data.frame(Armse = NA, 
    #                                                                 Amae  = NA, 
    #                                                                 Dintercept = NA, 
    #                                                                 Dslope = NA, 
    #                                                                 Dpearson = NA, 
    #                                                                 Dspearman = NA, 
    #                                                                 Psd = NA, 
    #                                                                 Pdispersion = NA, 
    #                                                                 Pr2 = NA, 
    #                                                                 Evaluation_number = NA, 
    #                                                                 Scale = if(is.null(scale)){0}else{scale}, 
    #                                                                 Evaluation_message = 'warning: no aggregation at this scale'))}
    
    observations = rescaled_observations$observations
    predictions  = rescaled_observations$predictions
    
    }
  
  # Linear model between values
  lm_test           <- tryCatch(lm(predictions ~ observations), error = function(e) NA)
  
  # if the lm test is NA then some metrics cannot be calculated
  if(!is.na(lm_test)){
    sum_lm            <- summary(lm_test)
    coef_lm           <- coef(lm_test)
    # summaries from linear model
    residual_standard_error <- lm_test$sigma
    Dintercept <- coef_lm[1]
    Dslope     <- coef_lm[2]
    Pr2         <- sum_lm$r.squared
  }else{
    Dintercept <- NA
    Dslope     <- NA
    Pr2        <- NA
  }
  
  tryCatch(cor.test(observations, predictions, method = 'pearson'), error = function(e) NA)
  cor.test_pearson  <-   tryCatch(cor.test(observations, predictions, method = 'pearson'), error = function(e) NA)
  cor.test_spearman <-   tryCatch(cor.test(observations, predictions, method = 'spearman'), error = function(e) NA)
  sd_predictions    <- tryCatch(sd(predictions), error = function(e) NA)
  

  # Accuracy (A-overall)
  # root mean squared error between average predicted abundance 
  # across all sites and observed abundance at each site
  Armse <- sqrt(mean((predictions - observations)^2, na.rm = T)) # root mean squared error gives more weight to big deviations
  Amae  <- mean(abs((predictions  - observations)), na.rm = T)    # mean absolute error weights all errors the same - if positive the value of observations is > the values of predictions
  
  # Discrimination
  if(is.na(cor.test_pearson)){Dpearson <- NA}else{Dpearson <- cor.test_pearson$estimate}
  if(is.na(cor.test_spearman)){Dspearman <- NA}else{Dspearman  <- cor.test_spearman$estimate}
  
  # Precision
  if(is.na(sd_predictions)){Psd <- NA}else{Psd <- sd_predictions}
  if(is.na(sd_predictions)){Pdispersion <- NA}else{Pdispersion <- sd_predictions / sd(observations)}
  
  metric_summary <- data.frame(Armse = Armse, 
                               Amae  = Amae, 
                               Dintercept = Dintercept, 
                               Dslope = Dslope, 
                               Dpearson = Dpearson, 
                               Dspearman = Dspearman, 
                               Psd = Psd, 
                               Pdispersion = Pdispersion, 
                               Pr2 = Pr2, 
                               Evaluation_number = length(observations), 
                               Scale = if(is.null(scale)){0}else{scale}, 
                               Evaluation_message = 'none')
  
  return(metric_summary)
  
}

