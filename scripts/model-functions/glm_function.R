# function to fit glms 

# buildmer for model fitting and stepwise model selection of glmmTMB
# remotes::install_github("cvoeten/buildmer"); https://github.com/cvoeten/buildmer

# load in abundance data
#load("data/rls_abun_modelling_data_v2.RData")
#abundance = rls_abun_fitting[[1]]

# load in covariates
#load("data/rls_covariates.RData")
#covariates = rls_xy[c('SiteLongitude', 'SiteLatitude',
#                      'Depth_GEBCO', 
#                      'robPCA_1', 'robPCA_2', 'robPCA_3', 'robPCA_4', 'robPCA_5', 'robPCA_6')]

# get species names
#species_name <- unique(abundance$TAXONOMIC_NAME)

# model
#model <- 'glm'

# transformation
#transformation = NA

# family
#family = 'poisson'
#family = 'tweedie'

# function to fit glms
glm_function <- function(abundance = abundance, 
                         covariates = covariates, 
                         transformation = NA, # option is NA, log, log10
                         model = NA, # option is lm or glm
                         family = NA, # option is NA, poisson, or nbinom, or tweedie
                         zi = F,     # option is F or T
                         species_name = species_name, 
                         base_dir = 'results/modelfits'){
  require(tidyverse)
  require(glmmTMB)
  require(buildmer)
  
  if(model == 'lm' & is.na(transformation)){stop('attempting to fit a lm to raw abundance data')}
  
  # filter to focal species
  abundance <- abundance[c('SiteLongitude', 'SiteLatitude', 'Num')]
  names(abundance)[3] <- 'abundance'
  
  # join together abundance and covariates dataframes 
  abundance <- left_join(abundance, covariates)
  
  # apply appropriate transformation
  if(is.na(transformation)){NULL}else{ if(transformation == 'log'){
    abundance$abundance <- log(abundance$abundance+1)
  }
  if(transformation == 'log10'){
    abundance$abundance <- log10(abundance$abundance+1)
  }}

  # response variable name
  response <- 'abundance~'
  
  # rename covariates
  covNames_org <- names(covariates)
  covNames_org <- covNames_org[-which(covNames_org %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  for(i in 1:length(covNames_org)){names(abundance)[3+i] <- paste0('cov', i)}
  for(i in 1:length(covNames_org)){names(covariates)[2+i] <- paste0('cov', i)}
  
  # create formula with new covariate names
  covNames_new <- names(covariates)
  covNames_new <- covNames_new[-which(covNames_new %in% c('SiteCode', 'SiteLongitude', 'SiteLatitude'))]
  covNames_new_2 <- paste0('I(', covNames_new, '^2',')')
  covNames_combined <- paste0(c(covNames_new, covNames_new_2), collapse = '+')
  model_formula <- as.formula(paste0(response, covNames_combined))
  
  # build formula structure using buildmer
  model_tab <- tabulate.formula(model_formula)
  model_tab[-1, 'block'] <- rep(covNames_new, 2)
    
  if(model == 'lm'){
  
  # fit model with backwards stepwise regression
  model_fit <- buildglmmTMB(formula = model_tab,
                            data = abundance, 
                            dep = 'abundance') 
  
  # remove quadratic terms that are non-significant and refit model
  coef_table       <- coef(summary(model_fit))$cond               # get coefficient table
  coef_names       <- rownames(coef(summary(model_fit))$cond)[-1]
  coef_names_ns    <- coef_names[which(coef_table[-1,4] > 0.05)]  # get coefficients to remove
  
  # if there are non-significant I(cov^2) terms
  if(sum(grepl('I(', coef_names_ns, fixed = T)) > 0){
    model_terms <- coef_names[-which(coef_names %in% coef_names_ns[grep('I(', coef_names_ns, fixed = T)])] # get non-significant quadraties from table
    model_tab_2 <- model_tab[which(model_tab$term %in% model_terms),]    # subset model_tab by all significant terms
    model_tab <- rbind(model_tab[1,], model_tab_2) # put back intercept
    
  # refit model
  model_fit <- buildglmmTMB(formula = model_tab,
                            data = abundance, 
                            dep = 'abundance')

  }else{print('all terms are significant')}
  
  }
  
  if(model == 'glm'){
    
    
    # fit model with backwards stepwise regression
    if(!family %in% c('poisson', 'nbinom', 'tweedie')){stop('family for glm must be one of poisson or nbinom')}
    
    if(family == 'poisson'){model_fit <- buildglmmTMB(formula = model_tab,
                              data = abundance, 
                              dep = 'abundance', 
                              family = 'poisson')}
    
    if(family == 'nbinom'){model_fit <- buildglmmTMB(formula = model_tab,
                                                      data = abundance, 
                                                      dep = 'abundance', 
                                                      family = nbinom1(link = 'log'))}
    
    if(family == 'tweedie'){model_fit <- buildglmmTMB(formula = model_tab,
                                                     data = abundance, 
                                                     dep = 'abundance', 
                                                     family = tweedie(link = 'log'))}
    
    # remove quadratic terms that are non-significant and refit model
    coef_table       <- coef(summary(model_fit))$cond               # get coefficient table
    coef_names       <- rownames(coef(summary(model_fit))$cond)[-1]
    coef_names_ns    <- coef_names[which(coef_table[-1,4] > 0.05)]  # get coefficients to remove
    
    # if there are non-significant I(cov^2) terms
    if(sum(grepl('I(', coef_names_ns, fixed = T)) > 0){
      
    model_terms <- coef_names[-which(coef_names %in% coef_names_ns[grep('I(', coef_names_ns, fixed = T)])] # get non-significant quadraties from table
    model_tab_2 <- model_tab[which(model_tab$term %in% model_terms),]    # subset model_tab by all significant terms
    model_tab <- rbind(model_tab[1,], model_tab_2) # put back intercept
      
      # refit model
      if(family == 'poisson'){model_fit <- buildglmmTMB(formula = model_tab,
                                                        data = abundance, 
                                                        dep = 'abundance', 
                                                        family = 'poisson')}
      
      if(family == 'nbinom'){model_fit <- buildglmmTMB(formula = model_tab,
                                                       data = abundance, 
                                                       dep = 'abundance', 
                                                       family = nbinom1(link = 'log'))}
      
      if(family == 'tweedie'){model_fit <- buildglmmTMB(formula = model_tab,
                                                        data = abundance, 
                                                        dep = 'abundance', 
                                                        family = tweedie(link = 'log'))}
      
    }else{print('All terms are significant')}
    
  }
  
  
  
  # for zero-inflation take the selected model with no zero-inflation and remodel covariate set
  if(zi == T){
    
    # get the zi model structure
    model_structure <- formula(model_fit@model)
    zi_formula <- formula(paste('~', paste(attr(terms(model_structure), 'term.labels'), collapse = '+')))
    
    if(family == 'poisson'){model_fit <- glmmTMB(formula = model_structure,
                                                 zi = zi_formula, 
                                                 data = abundance, 
                                                 family = 'poisson')}
    
    if(family == 'nbinom'){model_fit <- glmmTMB(formula = model_structure,
                                                zi = zi_formula,
                                                data = abundance,                   
                                                family = nbinom1(link = 'log'))}
    
    if(family == 'tweedie'){model_fit <- glmmTMB(formula = model_structure,
                                                zi = zi_formula,
                                                data = abundance,                   
                                                family = tweedie(link = 'log'))}
    
    
  }
  
  # save output into appropriate folder system
  # transformation/model/family/zi/
  model_dir <- paste0(if(is.na(transformation)){'raw'}else{transformation}, '/', 
                      model , '/', 
                      if(is.na(family)){'gaussian'}else{family}, '/', 
                      if(zi==F){'no_ZI'}else{'ZI'}
                      )
  model_path <- paste0(base_dir, '/', model_dir)
  dir.create(model_path, recursive = T)
  save(model_fit, file = paste0(model_path, '/', gsub(' ', '_', species_name), '.RData'), recursive = T)

}



