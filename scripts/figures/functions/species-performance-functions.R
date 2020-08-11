


# produce plots from outputs of linear models ----

#input_data <- bbs_basic
#directory  <- paste0('figures/species-performance-figures/', unique(input_data$cross_validation_2), '/')
#species_traits <- c('mean_abundance', 'frequency')

plot_trait_models <- function(input_data, 
                              directory, 
                              name,
                              species_traits =  c('mean_abundance', 'frequency'), 
                              models = c('accuracy', 'discrimination', 'precision', 'aggregated_evaluation_metrics')
                              ){
  
  # create output directory
  dir.create(directory, recursive = T)
  
  # save pairs plots between variables
  png(filename = paste0(directory, 'pairs_plot_models.png'), height = 1000, width = 1000, res = 100)
  pairs(input_data[models])
  dev.off()
  
  png(filename = paste0(directory, 'pairs_plot_traits.png'), height = 1000, width = 1000, res = 100)
  pairs(input_data[species_traits])
  dev.off()
  
  
  plot_output <- rep(list(list()), length(models))
  for(j in 1:length(models)){
    
    input_data$response <- as.numeric(input_data[models[j]][[1]])
    
  for(i in 1:length(species_traits)){
    
  input_data$trait <- as.numeric(input_data[species_traits[i]][[1]])
  
  # create a model between metrics and species traits
  output         <- lm(response ~ trait*plot_level, data = input_data)
  summary_output <- summary(output)
  coefficients   <- summary_output$coefficients
  
  # estimate coefficients of interaction
  interaction_coefficients <- coefficients[grep(':plot_level', rownames(coefficients)),1]
  interaction_pvalue       <- coefficients[grep(':plot_level', rownames(coefficients)),4]
  baseline_slope           <- coefficients[grep('trait', rownames(coefficients), fixed = T)[1],1]
  
  interaction_coefficients <- baseline_slope + interaction_coefficients # estiamte coefficient values
  
  # turn into dataset 
  df_names <- sub('.*\\:plot_level', '', names(interaction_coefficients), fixed = F)
  output <- data.frame(model = sub('\\..*','', df_names), 
             abundance_response = stringi::stri_split_fixed(df_names, '.', simplify = T)[,2], 
             distribution       = stringi::stri_split_fixed(df_names, '.', simplify = T)[,3],
             transformation     = stringi::stri_split_fixed(df_names, '.', simplify = T)[,4],
             full_name          = gsub('.','-', df_names, fixed = T),
             interaction_coefficients = as.numeric(interaction_coefficients), 
             interaction_pvalue       = as.numeric(interaction_pvalue))
  
  # assign factor levels for colours
  output$model <- factor(as.factor(output$model), levels)
  
  # create plot for a given run
  plot <- ggplot(output) + 
    geom_hline(aes(yintercept = 0)) + 
    geom_bar(aes(x = reorder(full_name, -interaction_coefficients), y = interaction_coefficients, fill = model),
             col = 'black',
             stat = 'identity') + 
    geom_bar(data = output[which(interaction_pvalue > 0.05),], 
             aes(x = reorder(full_name, -interaction_coefficients), y = interaction_coefficients),
             col = 'black',
             stat = 'identity') + 
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), 
          axis.title.x = element_blank(), 
          aspect.ratio = 0.2, 
          title = element_text(size = 8), 
          legend.position = 'none') + 
    scale_fill_manual(values = colours) + 
    ylab('coefficient') + 
    ggtitle(label = paste(gsub('_', ' ',models[j]), gsub('_', ' ',species_traits[i]), sep = ' ~ '))
  
  plot_output[[j]][[i]] <- plot
  
  }
  }
  
  # combine all plots
  
  # print to pdf
  pdf(paste0(directory,'/', name, '.pdf'), width = 15, height = 15)
  all_plots <- do.call("grid.arrange", c(flatten(plot_output), ncol = 2))
  print(all_plots)
  dev.off()
  
}


# create linear models for best species and plot as three way interaction plots ----


interaction_trait_plots <- function(plot_data, 
                                    directory, 
                                    name, 
                                    width = 10, 
                                    height = 10, 
                                    width_marginal = 3, 
                                    height_marginal = 10,
                                    dataset){

require(MuMIn)

plot_data$abundance    <- as.numeric(scale(log10(plot_data$mean_abundance), scale = T, center = T))
plot_data$frequency    <- as.numeric(scale(log10(plot_data$frequency), scale = T, center = T))
plot_data$observations <- as.numeric(scale(log10(plot_data$n_per_species), scale = T, center = T))

# fit linear model
model_accuracy       <- lm(accuracy       ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')
model_discrimination <- lm(discrimination ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')
model_precision      <- lm(precision      ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')

# perform model selection
dredge_accuracy       <- dredge(model_accuracy)
dredge_discrimination <- dredge(model_discrimination)
dredge_precision      <- dredge(model_precision)

# get raw models
models_accuracy       <- get.models(model.sel(dredge_accuracy, fit = T), subset = delta < 2)
models_discrimination <- get.models(model.sel(dredge_discrimination, fit = T), subset = delta < 2)
models_precision      <- get.models(model.sel(dredge_precision, fit = T), subset = delta < 2)

# get r2
r2_accuracy       <- mean(unlist(lapply(models_accuracy, function(x) summary(x)$r.squared)))
r2_discrimination <- mean(unlist(lapply(models_discrimination, function(x) summary(x)$r.squared)))
r2_precision      <- mean(unlist(lapply(models_precision, function(x) summary(x)$r.squared)))

# get model averages
if(length(models_accuracy) == 1){
  avg_accuracy <- models_accuracy
  accuracy_coef <- data.frame(coef(summary(avg_accuracy[[1]])))
  accuracy_coef$evaluation <- 'accuracy'
  accuracy_coef$term <- rownames(accuracy_coef)
  accuracy_coef$r2   <- summary(avg_accuracy[[1]])$r.squared
  accuracy_coef <- accuracy_coef %>% dplyr::rename(test = 't.value', `p-value` = 'Pr...t..' )
}else{
  avg_accuracy <- model.avg(dredge(model_accuracy), subset = delta < 2, fit = T)
  accuracy_coef <- data.frame(summary(avg_accuracy)$coefmat.full)[,-3]
  accuracy_coef$evaluation <- 'accuracy'
  accuracy_coef$term <- rownames(accuracy_coef)
  accuracy_coef$r2   <- r2_accuracy
  accuracy_coef <- accuracy_coef %>% dplyr::rename(test = 'z.value', `p-value` = 'Pr...z..' ) 
}


if(length(models_discrimination) == 1){
  avg_discrimination <- models_discrimination
  discrimination_coef <- data.frame(coef(summary(avg_discrimination[[1]])))
  discrimination_coef$evaluation <- 'discrimination'
  discrimination_coef$term <- rownames(discrimination_coef)
  discrimination_coef$r2   <- summary(avg_discrimination[[1]])$r.squared
  discrimination_coef <- discrimination_coef %>% dplyr::rename(test = 't.value', `p-value` = 'Pr...t..' )
}else{
  avg_discrimination <- model.avg(dredge(model_discrimination), subset = delta < 2, fit = T)
  discrimination_coef <- data.frame(summary(avg_discrimination)$coefmat.full)[,-3]
  discrimination_coef$evaluation <- 'discrimination'
  discrimination_coef$term <- rownames(discrimination_coef)
  discrimination_coef$r2   <- r2_discrimination
  discrimination_coef <- discrimination_coef %>% dplyr::rename(test = 'z.value', `p-value` = 'Pr...z..' ) 
}

if(length(models_precision) == 1){
  avg_precision <- models_precision
  precision_coef <- data.frame(coef(summary(avg_precision[[1]])))
  precision_coef$evaluation <- 'precision'
  precision_coef$term <- rownames(precision_coef)
  precision_coef$r2   <- summary(avg_precision[[1]])$r.squared
  precision_coef <- precision_coef %>% dplyr::rename(test = 't.value', `p-value` = 'Pr...t..' )
}else{
  avg_precision <- model.avg(dredge(model_precision), subset = delta < 2, fit = T)
  precision_coef <- data.frame(summary(avg_precision)$coefmat.full)[,-3]
  precision_coef$evaluation <- 'precision'
  precision_coef$term <- rownames(precision_coef)
  precision_coef$r2   <- r2_precision
  precision_coef <- precision_coef %>% dplyr::rename(test = 'z.value', `p-value` = 'Pr...z..' ) 
}

model_table <- rbind(accuracy_coef, discrimination_coef, precision_coef) %>% 
  filter(term != '(Intercept)') %>% 
  dplyr::rename(coef = 'Estimate', se = 'Std..Error', `r-squared` = r2) %>% 
  select(evaluation, term, coef, se, `test`, `p-value`, `r-squared`) %>% 
  mutate(coef = signif(coef, 2), 
         se   = signif(se, 2), 
         `z-value`   = round(`test`, 2), 
         `p-value`   = ifelse(round(`p-value`, 3) == 0, 0.001, round(`p-value`, 3)), 
         `r-squared` = signif(`r-squared`, 2))

dir.create(directory, recursive = T)
writexl::write_xlsx(model_table, path = paste0(directory, '/', name,'.xlsx'))

# make prediction frames
prediction_frame <- expand.grid(abundance = seq(from = min(plot_data$abundance), to =  max(plot_data$abundance), length.out = 50), 
                                frequency = seq(from = min(plot_data$frequency), to =  max(plot_data$frequency), length.out = 50),
                                observations = seq(from = min(plot_data$observations), to =  max(plot_data$observations), length.out = 50))

# predict for each model 
if(length(models_accuracy) == 1){accuracy_preds <- predict(avg_accuracy[[1]], prediction_frame, se.fit = T)}else{accuracy_preds <- predict(avg_accuracy, prediction_frame, se.fit = T)}
if(length(models_discrimination) == 1){discrimination_preds <- predict(avg_discrimination[[1]], prediction_frame, se.fit = T)}else{discrimination_preds <- predict(avg_discrimination, prediction_frame, se.fit = T)}
if(length(models_precision) == 1){precision_preds <- predict(avg_precision[[1]], prediction_frame, se.fit = T)}else{precision_preds <- predict(avg_precision, prediction_frame, se.fit = T)}

prediction_frame$accuracy       <- accuracy_preds$fit
prediction_frame$discrimination <- discrimination_preds$fit
prediction_frame$precision      <- precision_preds$fit

prediction_frame$accuracy_upr       <- accuracy_preds$fit + accuracy_preds$se.fit 
prediction_frame$discrimination_upr <- discrimination_preds$fit + discrimination_preds$se.fit 
prediction_frame$precision_upr      <- precision_preds$fit + precision_preds$se.fit 

prediction_frame$accuracy_lwr       <- accuracy_preds$fit - accuracy_preds$se.fit 
prediction_frame$discrimination_lwr <- discrimination_preds$fit - discrimination_preds$se.fit 
prediction_frame$precision_lwr      <- precision_preds$fit - precision_preds$se.fit 

theme_interaction <- function(){theme(aspect.ratio = 1, 
                                     panel.background = element_blank(), 
                                     panel.border = element_blank(), 
                                     panel.grid = element_blank(), 
                                     legend.position = 'none', 
                                     axis.ticks = element_blank())}

# labels <- c('accuracy', 'discrimination', 'precision')
# panel <- list()
# for(i in 1:length(labels)){
#   
#   legend <- ggplot(data = prediction_frame) + 
#     geom_point(aes(x = prediction_frame[,labels[i]], y = prediction_frame[,labels[i]], fill = prediction_frame[,labels[i]])) + 
#     theme_bw() +
#     scale_fill_viridis_c(option = colour, name =  labels[i]) + 
#     theme(legend.justification = 'centre')
#   
#   panel[[i]] <- grid.arrange(
#   ggplot(data = prediction_frame) + 
#     geom_tile(aes(x = abundance, y = frequency, fill = prediction_frame[,labels[i]])) + 
#     theme_bw() +
#     theme_interaction() + 
#     scale_fill_viridis_c(option = colour) + 
#     xlab('mean abundance') + ylab('frequency of occurrence'),
#   
#   ggplot(data = prediction_frame) + 
#     geom_tile(aes(x = abundance, y = observations, fill = prediction_frame[,labels[i]])) + 
#     theme_bw() +
#     theme_interaction() + 
#     scale_fill_viridis_c(option = colour) + 
#     xlab('mean abundance') + ylab('number of observations'),
#   
#   ggplot(data = prediction_frame) + 
#     geom_tile(aes(x = frequency, y = observations, fill = prediction_frame[,labels[i]])) + 
#     theme_bw() +
#     theme_interaction() + 
#     scale_fill_viridis_c(option = colour) + 
#     xlab('frequency of occurrence') + ylab('number of observations'),
#   
#   get_legend(legend),
#   nrow = 1)
# }
# 
# pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
# print(do.call('grid.arrange', panel))
# dev.off()

# plot only the significant marginal effects or interactions ----

sig_term <- model_table %>% filter(`p-value` < 0.05)

input_data <- list()
plot_list <- list()

colour = ifelse(dataset == 'bbs', 3, 7)

plot_list <- lapply(1:nrow(sig_term), function(i){
  
  # first check if there is an interaction...
  interaction <- grepl(':', sig_term$term[i], fixed = T)
  if(interaction == T){
    
    terms <- str_split(sig_term$term[i], ':', simplify = F)[[1]]
    filter_values   <- unique(prediction_frame[,terms[2]])[c(1, 25, 50)]
    input_data[[i]] <- prediction_frame[,c(terms, 
                                           sig_term$evaluation[i], 
                                           paste0(sig_term$evaluation[i], '_upr'),
                                           paste0(sig_term$evaluation[i], '_lwr'))]
    input_data[[i]] <- input_data[[i]][which(input_data[[i]][,2] %in% filter_values),]
    input_data[[i]] <- distinct(input_data[[i]], across(contains(terms)), .keep_all = TRUE)
  
    plot_list <- ggplot() + 
    geom_ribbon(data = input_data[[i]], aes(x = input_data[[i]][,terms[1]], 
                                            ymax=input_data[[i]][, paste0(sig_term$evaluation[i], '_upr')], 
                                            ymin=input_data[[i]][, paste0(sig_term$evaluation[i], '_lwr')], 
                                            group = input_data[[i]][,terms[2]], 
                                            fill = input_data[[i]][,terms[2]]), alpha = 0.5) + 
    geom_line(data = input_data[[i]], aes(x = input_data[[i]][,terms[1]], 
                                          y = input_data[[i]][,sig_term$evaluation[i]],
                                          group = input_data[[i]][,terms[2]])) + 
    scale_fill_gradientn(colours = viridis(3, option = colour), 
                         name = terms[2]) + 
    xlab(terms[1]) + 
    ylab(sig_term$evaluation[i]) + 
      theme(aspect.ratio = 1)
    
    # what to do for three way interaction
    if(length(terms) == 3){
      
      lty_values <- unique(prediction_frame[,terms[3]])[c(1, 25, 50)]
      
      input_data[[i]] <- prediction_frame[,c(terms, 
                                             sig_term$evaluation[i], 
                                             paste0(sig_term$evaluation[i], '_upr'),
                                             paste0(sig_term$evaluation[i], '_lwr'))]
      term_2 <- input_data[[i]][which(input_data[[i]][,2] %in% filter_values),]
      term_2 <- distinct(term_2, across(contains(terms[1:2])), .keep_all = TRUE)
      term_3 <- input_data[[i]][which(input_data[[i]][,3] %in% lty_values),]
      term_3 <- distinct(term_3, across(contains(terms[c(1,3)])), .keep_all = TRUE)
      
      term_2_plot <- ggplot() + 
        geom_ribbon(data = term_2, aes(x = term_2[,terms[1]], 
                                       ymax=term_2[, paste0(sig_term$evaluation[i], '_upr')], 
                                       ymin=term_2[, paste0(sig_term$evaluation[i], '_lwr')], 
                                       group = term_2[,terms[2]], 
                                       fill = term_2[,terms[2]]), alpha = 0.5) + 
        geom_line(data = term_2, aes(x = term_2[,terms[1]], 
                                     y = term_2[,sig_term$evaluation[i]],
                                     group = term_2[,terms[2]])) + 
        scale_fill_gradientn(colours = viridis(3, option = colour, begin = 0.2, end = 0.8), 
                             name = terms[2]) + 
        xlab(terms[1]) + 
        ylab(sig_term$evaluation[i]) + 
        theme(aspect.ratio = 1)
      
      term_3_plot <- ggplot() + 
        geom_ribbon(data = term_3, aes(x = term_3[,terms[1]], 
                                       ymax=term_3[, paste0(sig_term$evaluation[i], '_upr')], 
                                       ymin=term_3[, paste0(sig_term$evaluation[i], '_lwr')], 
                                       group = term_3[,terms[3]], 
                                       fill = term_3[,terms[3]]), alpha = 0.5) + 
        geom_line(data = term_3, aes(x = term_3[,terms[1]], 
                                     y = term_3[,sig_term$evaluation[i]],
                                     group = term_3[,terms[3]])) + 
        scale_fill_gradientn(colours = viridis(3, option = colour, begin = 0.2, end = 0.8), 
                             name = terms[3]) + 
        xlab(terms[1]) + 
        ylab(sig_term$evaluation[i]) + 
        theme(aspect.ratio = 1)
      
      dir.create(paste0(directory, '_marginal'), recursive = T)
      pdf(file = paste0(paste0(directory, '_marginal'), '/', name,'_3-way', '.pdf'), width = 8, height = 4)
      print(patchwork::wrap_plots(list(term_2_plot, term_3_plot)))
      dev.off()
      
      plot_list <- NULL
      
      return(plot_list)
      
      }}else{
  
  input_data[[i]] <- prediction_frame[,c(sig_term$term[i], 
                                         sig_term$evaluation[i], 
                                         paste0(sig_term$evaluation[i], '_upr'),
                                         paste0(sig_term$evaluation[i], '_lwr'))]
  input_data[[i]] <- distinct(input_data[[i]], across(contains(sig_term$term[i])), .keep_all = TRUE)
  
  plot_list <- ggplot() + 
    geom_ribbon(data = input_data[[i]], aes(x = input_data[[i]][,sig_term$term[i]], 
                    ymax=input_data[[i]][, paste0(sig_term$evaluation[i], '_upr')], 
                    ymin=input_data[[i]][, paste0(sig_term$evaluation[i], '_lwr')]), 
                fill = 'gray50', alpha = 0.5) + 
    geom_line(data = input_data[[i]], aes(x = input_data[[i]][,sig_term$term[i]], 
                  y = input_data[[i]][,sig_term$evaluation[i]])) + 
    xlab(sig_term$term[i]) + 
    ylab(sig_term$evaluation[i]) + 
    theme(aspect.ratio = 1)
  
  }
  
  return(plot_list)
  
})

plot_list <- purrr::compact(plot_list)

# create directory for marginal plots
dir.create(paste0(directory, '_marginal'), recursive = T)
pdf(file = paste0(paste0(directory, '_marginal'), '/', name,'.pdf'), width = width_marginal, height = height_marginal)
print(patchwork::wrap_plots(plot_list))
dev.off()

}






