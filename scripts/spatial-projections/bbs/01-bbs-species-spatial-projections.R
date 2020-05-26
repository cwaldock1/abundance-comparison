# script to identify the best performing model and produce a spatial projection from that model 

# load libraries ----
lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'raster')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load functions ---- 

# for aggregate_metrics
source('scripts/figures/functions/model-performance-functions.R')

# select best fitted model for each model type based on a concensus metrics ----

# read and clean assessment data as in the script 01-model-performance-figures.R

all_assessments <- lapply(list.files('results/model_assessment', full.names = T), readRDS)

all_assessments <- do.call(rbind, all_assessments)

# select only the columns to be used later 
all_assessments <- all_assessments %>% dplyr::select(-family_grouped_simple, -family_grouped, -family, 
                                              -transformation, -n_absence, -n_boot_absence, -abundance_response_simple, 
                                              -mean_abundance, -frequency, -mean_abundance_perc, -frequency_perc, -n_abundance) %>% 
  
  # create new cross validation level
  mutate(cross_validation_2 = gsub('rls_', '', gsub('bbs_', '', .$cross_validation))) %>% 
  
  # select final columns for this script
  dplyr::select(dataset, cross_validation, cross_validation_2, fitted_model, abundance_response, plot_level, species_name,
         Armse, Amae, Dintercept, Dslope, Dpearson, Dspearman, Psd, Pdispersion, Pr2)

# overall this code finds the best fitted model for each species within a type of model and cross-validation and dataset combination. Use this object to left join to the true predictions and create plots
# comparing each modelling frameworks overall predictions vs. observations. 
best_models <- all_assessments %>% 
  dplyr::select(-Armse, -Psd) %>% 
  # estimate the relative metric performance within a cross validation and dataset
  group_by(cross_validation_2, dataset) %>% 
  nest() %>% 
  mutate(metric_aggregation = purrr::map(data, ~aggregate_metrics(., 
                                                                  metrics = c('Amae', 'Dintercept', 'Dslope', 'Dpearson', 'Dspearman', 'Pdispersion', 'Pr2')))) %>% 
  .$metric_aggregation %>% 
  do.call(rbind, .) %>% 
  # find the best fitting model for each species within each fitted_model
  group_by(species_name, cross_validation) %>% 
  do(best_model = .$plot_level[which.max(.$discrimination)]) %>% 
  unnest(cols = c('best_model')) %>% 
  mutate(dataset = gsub('_basic|_oob_cv','', .$cross_validation),
         cross_validation = gsub('bbs_|rls_', '', .$cross_validation), 
         plot_level = .$best_model)

# select only bbs and oob_cv to get the best out the bag validation
best_models_spp <- best_models %>% filter(cross_validation == 'oob_cv', dataset == 'bbs')



# read in spatial datasets as a large matrix ----

spatial_datasets <- list.files('data/bbs_rasters', pattern = 'gri', full.names = T)

test_raster <- lapply(spatial_datasets, function(x){
  
  ras <- raster(x)
  crs(ras) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
  ras_points <- data.frame(rasterToPoints(ras))
  #ras_points$x <- signif(ras_points$x, digits = 4)
  #ras_points$y <- signif(ras_points$y, digits = 4)
  return(ras_points)

})

# ensure x and y are identical (previous code seems to have added 0.00001 to some of the xy values in conversions)
test_raster <- lapply(test_raster, function(x){
  x$x <- test_raster[[1]]$x 
  x$y <- test_raster[[1]]$y
  return(x)
})

# apply together
test_raster_all <- Reduce(left_join, test_raster)

# check for NAs 
lapply(test_raster_all, function(x) table(is.na(x)))

# select only covariates used in modelling
load('data/bbs_covariates.RData')
covariates = bbs_xy[c('SiteLongitude', 
                      'SiteLatitude',
                      "robPCA_1", 
                      "robPCA_2", 
                      "robPCA_3", 
                      "human_pop",
                      'primary_forest', 
                      'Elevation_GEBCO')]

# grab relevant columns and rename
covariates_usa <- test_raster_all[c('x', 
                                    'y',
                                    "merra.pc1",
                                    "merra.pc2", 
                                    "merra.pc3", 
                                    "human_pop",
                                    'primary_forest', 
                                    'Elevation.relative.to.sea.level')]

covNames_org <- names(covariates_usa)
covNames_org <- covNames_org[-which(covNames_org %in% c('x', 'y'))]
for(i in 1:length(covNames_org)){names(covariates_usa)[2+i] <- paste0('cov', i)}
covNames_new <- names(covariates_usa) # randomForests take matrix which can be subset with this object 
covNames_new <- covNames_new[-which(covNames_new %in% c('x', 'y'))]

# read in a random model and see predictions
load('/Volumes/Simulation/conor/abundance-comparison/results/bbs_cv_all/suitability/Acanthis hyemalis.RData')
load('/Volumes/Simulation/conor/abundance-comparison/results/bbs_cv_all/model_abunocc_2stage/log10/lm/gaussian/no_ZI/Acanthis_hyemalis.RData')
library(randomForest)
library(gamm4)
library(buildmer)
test_prediction <- as.numeric(predict(model_fit[[1]]@model, covariates_usa))
hist(test_prediction)

test_raster <- rasterFromXYZ(cbind(covariates_usa, test_prediction))
plot(test_raster[[7]])


