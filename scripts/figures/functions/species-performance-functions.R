


# produce plots from outputs of linear models ----

input_data <- bbs_basic
directory  <- paste0('figures/species-performance-figures/', unique(input_data$cross_validation_2), '/')
species_traits <- c('mean_abundance', 'frequency')

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
                                    option){

require(MuMIn)

plot_data$abundance <- scale(log10(plot_data$mean_abundance), scale = T, center = T)
plot_data$frequency      <- scale(log10(plot_data$frequency), scale = T, center = T)
plot_data$observations     <- scale(log10(plot_data$Evaluation_number), scale = T, center = T)

# fit linear model
model_accuracy       <- lm(accuracy       ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')
model_discrimination <- lm(discrimination ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')
model_precision      <- lm(precision      ~ abundance * frequency * observations, data = plot_data, na.action = 'na.fail')

# perform model selection
avg_accuracy       <- model.avg(dredge(model_accuracy), subset = delta < 2, fit = T)
avg_discrimination <- model.avg(dredge(model_discrimination), subset = delta < 2, fit = T)
avg_precision      <- model.avg(dredge(model_precision), subset = delta < 2, fit = T)

# create table of coefficients and save
accuracy_coef <- data.frame(summary(avg_accuracy)$coefmat.full)
accuracy_coef$evaluation <- 'accuracy'
accuracy_coef$term <- rownames(accuracy_coef)
discrimination_coef <- data.frame(summary(avg_discrimination)$coefmat.full)
discrimination_coef$evaluation <- 'discrimination'
discrimination_coef$term <- rownames(discrimination_coef)
precision_coef <- data.frame(summary(avg_precision)$coefmat.full)
precision_coef$evaluation <- 'precision'
precision_coef$term <- rownames(precision_coef)

model_table <- rbind(accuracy_coef, discrimination_coef, precision_coef) %>% 
  select(-Adjusted.SE) %>% 
  filter(term != '(Intercept)') %>% 
  rename(coef = 'Estimate', se = 'Std..Error', `z-value` = z.value, `p-value` = Pr...z..) %>% 
  select(evaluation, term, coef, se, `z-value`, `p-value`) %>% 
  mutate(coef = signif(coef, 2), 
         se   = signif(se, 2), 
         `z-value` = round(`z-value`, 2), 
         `p-value` = round(`p-value`, 3))

dir.create(directory, recursive = T)
writexl::write_xlsx(model_table, path = paste0(directory, '/', name,'.xlsx'))

# make prediction frames
prediction_frame <- expand.grid(abundance = seq(from = min(plot_data$abundance), to =  max(plot_data$abundance), length.out = 50), 
                                frequency = seq(from = min(plot_data$frequency), to =  max(plot_data$frequency), length.out = 50),
                                observations = seq(from = min(plot_data$observations), to =  max(plot_data$observations), length.out = 50))

# predict for each model 
prediction_frame$accuracy       <- predict(avg_accuracy, prediction_frame)
prediction_frame$discrimination <- predict(avg_discrimination, prediction_frame)
prediction_frame$precision      <- predict(avg_precision, prediction_frame)

theme_interaction <- function(){theme(aspect.ratio = 1, 
                                     panel.background = element_blank(), 
                                     panel.border = element_blank(), 
                                     panel.grid = element_blank(), 
                                     legend.position = 'none', 
                                     axis.ticks = element_blank())}

labels <- c('accuracy', 'discrimination', 'precision')
panel <- list()
for(i in 1:length(labels)){
  
  legend <- ggplot(data = prediction_frame) + 
    geom_point(aes(x = prediction_frame[,labels[i]], y = prediction_frame[,labels[i]], fill = prediction_frame[,labels[i]])) + 
    theme_bw() +
    scale_fill_viridis_c(option = option, name =  labels[i]) + 
    theme(legend.justification = 'centre')
  
  panel[[i]] <- grid.arrange(
  ggplot(data = prediction_frame) + 
    geom_tile(aes(x = abundance, y = frequency, fill = prediction_frame[,labels[i]])) + 
    theme_bw() +
    theme_interaction() + 
    scale_fill_viridis_c(option = option) + 
    xlab('mean abundance') + ylab('frequency of occurrence'),
  
  ggplot(data = prediction_frame) + 
    geom_tile(aes(x = abundance, y = observations, fill = prediction_frame[,labels[i]])) + 
    theme_bw() +
    theme_interaction() + 
    scale_fill_viridis_c(option = option) + 
    xlab('mean abundance') + ylab('number of observations'),
  
  ggplot(data = prediction_frame) + 
    geom_tile(aes(x = frequency, y = observations, fill = prediction_frame[,labels[i]])) + 
    theme_bw() +
    theme_interaction() + 
    scale_fill_viridis_c(option = option) + 
    xlab('frequency of occurrence') + ylab('number of observations'),
  
  get_legend(legend),
  nrow = 1)
}

pdf(file = paste0(directory, '/', name,'.pdf'), width = width, height = height)
print(do.call('grid.arrange', panel))
dev.off()

}






