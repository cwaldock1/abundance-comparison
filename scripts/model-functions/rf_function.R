
# Function for fitting random forest abundance models 

# load in abundance data
# load("data/rls_abun_modelling_data.RData")
# abundance = rls_abun_fitting

# load in covariates
# load("data/rls_covariates.RData")
# covariates = rls_xy[c('SiteCode', 'SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]
# covariates[,4] <- as.numeric(scale(log(abs(covariates[,4])), center = T))
# get species names
# species_name <- as.character(names(abundance)[4])

# discrete
# discrete <- T

# transformation
# transformation = 'log'

rf_function <- function(abundance = abundance, 
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         discrete = NA,       # option is T or F 
                         species_name = species_name, 
                         base_dir = 'results/modelfits'){
  
  require(tidyverse)
  require(randomForest)
  
  # filter to focal species
  abundance <- abundance[c('SiteCode', 'SiteLongitude', 'SiteLatitude', species_name)]
  names(abundance)[4] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
  }
    if(transformation == 'log10'){
      abundance$abundance <- log10(abundance$abundance+1)
    }}
  
  # convert to discrete values based on groupings as in Howard, C., Stephens, P. A., Pearce-Higgins, J. W., Gregory, R. D. & Willis, S. G. Improving species distribution models: the value of data on abundance. Methods Ecol. Evol. 5, 506–513 (2014).
  if(discrete == T){
    abundance$abundance <- round(abundance$abundance)     # round logged abundances
    abundance$abundance[abundance$abundance > 6] <- 6     # truncate abundances
    abundance$abundance <- as.factor(abundance$abundance) # turn into factors for random forests
  }
  
  # rename and get general names for covariates for generalism
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[4+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[3+i] <- paste0('cov', i)}
  covNames_new <- names(covariates) # randomForests take matrix which can be subset with this object 
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  
  # fit the random forests
  # ?randomForest
  model_fit <- randomForest(x = abundance[covNames_new], 
                            y = abundance$abundance, 
                            ntree = 1000, 
                            importance=FALSE)
  #importance(model_fit)
  #rp <- response.plot2(models = c('model_fit'), 
  #               Data = data.frame(abundance[covNames_new]), 
  #               show.variables = covNames_new, 
  #               fixed.var.metric = 'mean', 
  #               plot = F, 
  #               use.formal.names = T)
  
  #par(mfrow = c(3,3))
  #for(i in 1:length(rp)){
  #  plot(rp[[i]][,2] ~ rp[[i]][,1], type = 'l')
  #}
  
  #ggplot(data = data.frame(rp), aes(x = expl.val, y = pred.val, lty = pred.name)) + 
  #  geom_line() + 
  #  rp.gg.theme + 
  #  facet_grid(~expl.name, scales = 'free_x')
  
  # save output into appropriate folder system
  # transformation/model/family/zi/
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/',
                      'rf' , '/', 
                      if(discrete == T){'discrete/'}else{'continuous/'})
  model_path <- paste0(base_dir, '/', model_dir)
  dir.create(model_path, recursive = T)
  save(model_fit, file = paste0(model_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)
  
}
  
