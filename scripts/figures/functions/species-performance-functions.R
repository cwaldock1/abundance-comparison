


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
