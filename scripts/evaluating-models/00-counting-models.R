#### Script to check modelling outputs and whether models fit for each class of model

# Set working directory if needed ----

setwd('/Volumes/Simulation/conor/abundance-comparison')

# Load packages ----

lib_vect <- c('stringi', 'tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)


# list files and create summaries ----

# all the information is stored in the predictions on which types of models are used etc 
# this needs to be summarised across all folders

# function to recursively list files
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

# get all the relevant prediction folders
files <- unlist(lapply(list.files('results', full.names = T)[1:4], function(x) list.dirs.depth.n(x, n = 2)))
prediction_files <- files[grepl('predictions', files)]
  
# perform loop to list and combine all the files for predictions
all_predictions <- list()
for(i in 1:length(prediction_files)){
  print(i/length(prediction_files))
  # list all the files
  sub_file_names <- list.files(prediction_files[i], recursive = T, full.names = T)
  
  # loop over all files and extract values
  sub_predictions <- list()
  for(j in 1:length(sub_file_names)){
  load(sub_file_names[j])
  sub_predictions[[j]] <- extracted_predictions
  }
  
  # create master list with all extracted predictions in each subfolder
  all_predictions[[i]] <- do.call(rbind, sub_predictions)
  names(all_predictions)[i] <-  prediction_files[i]
}


# convert list of dataframes into large dataframe
all_predictions_df <- bind_rows(all_predictions, .id = 'folder')

# manage columns for summaries ----
# managing data groupings for plots
add_family_plot_column <- function(plot_data){
  
  # edit transformations
  plot_data$transformation[is.na(plot_data$transformation)] <- ''
  plot_data$transformation <- ifelse(plot_data$transformation == '', '', paste0('-', plot_data$transformation))
  
  # edit family groups
  plot_data$family <- recode(plot_data$family, discrete = "multinomial", discrete_multinomial = "multinomial", 
                             continuous_gaussian = 'gaussian',continuous = 'gaussian', continuous_poisson = 'poisson', 
                             zip = 'poisson')
  plot_data$family[is.na(plot_data$family)] <- 'gaussian'
  
  # edit zero-inflations
  plot_data$zi <- ifelse(plot_data$zi, '-ZI', '')
  plot_data$zi[is.na(plot_data$zi)] <- ''
  
  plot_data$family_plot <- paste(plot_data$family, plot_data$zi, plot_data$transformation, sep = '')
  
  return(plot_data)
  
}

# manage famile columsn into standard categories
final_predictions <- add_family_plot_column(all_predictions_df) %>% 
  select(folder, dataset, species_name, fitted_model, only_abundance, family, transformation, zi)

# check out summaries of numbers of species ---- 

# summary of numbers of species across each model
summary_predictions <- final_predictions %>% 
  group_by(folder, dataset, fitted_model, only_abundance, family, transformation, zi) %>% 
  do(n_species = length(.$species_name)) %>% 
  unnest() %>% 
  data.frame()

# which models have fewer that 50 species? 
summary_predictions %>% filter(n_species !=50)

# check in the gbm
summary_predictions %>% filter(n_species != 50) %>% filter(fitted_model == 'gbm', family == 'multinomial')

# which species are not fitted by gbms? can this be debugged further or is it more fundamental 
sp_gbm <- final_predictions %>% 
  filter(fitted_model == 'gbm', family == 'multinomial') %>% 
  .$species_name %>% 
  unique

# species for which GBMs did no fit. 
unique(final_predictions$species_name)[unique(final_predictions$species_name) %in% sp_gbm] # try and debug manually from here for m
which(unique(final_predictions$species_name) %in% sp_gbm) # try and debug manually from here for m
i = 1


# summaries of overall modelling 
table(final_predictions$folder, final_predictions$dataset, final_predictions$fitted_model)

# create a batch file for re-fitting models that don't run on first instance ----

# function to find missing species across folders for a particular dataset 
find_missing_mofos <- function(prediction_folder, all_species){ 
  
  test <-  list.files(prediction_folder)
  test_location <- stri_locate_all(test, fixed = '_')
  located_string <- lapply(test_location, function(x){ 
    n <- nrow(x)
    return(as.numeric(x[n-1,1]))})
  
  species_in_folder <- do.call(c, lapply(1:length(test), function(x){
    gsub('.RData', '', str_sub(test[x], located_string[[x]]+1, -1))
  }))
  
  model_in_folder <- do.call(c, lapply(1:length(test), function(x){
    gsub('.RData', '', str_sub(test[x], 0, located_string[[x]]-1))
  }))
  
  test_2 <- tibble(model_in_folder, species_in_folder) %>% 
    group_by(model_in_folder) %>% 
    do(missing_species = list(.$species_in_folder))
  
  all_missing <- tibble(model_in_folder, species_in_folder) %>% 
    group_by(model_in_folder) %>% 
    do(missing_species = all_species[-grep(paste(all_species, collapse = '|'), .$species_in_folder)]) 
  
  return(unique(unlist(all_missing$missing_species)))
  
  }

# read in bbs data
load("data/bbs_abun_modelling_data.RData")

# find missing species for each folder and list all together as a full set of missing species
bbs_predictions <- prediction_files[grep('bbs', prediction_files)]

list_bbs_missing <- lapply(bbs_predictions, function(x) 
  find_missing_mofos(x, 
                     all_species = gsub(' ', '_', unique(abundance_key$TAXONOMIC_NAME))))

bbs_missing <- unique(do.call(c, list_bbs_missing))


# create batch file for missing species
i <- which(gsub(' ', '_',abundance_key$TAXONOMIC_NAME) %in% bbs_missing)

Nodes <- c('@RunAsMultiple,@NodeSet_16')

file_out <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-models%bbs%bbs_config_basic.R',
                        ' ', 
                        i)

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/fitting-models/batchscripts/batch_bbs_basic_debug.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)

fileConn<-file("scripts/fitting-models/batchscripts/batch_bbs_CV_debug.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


# read in bbs data
load("data/rls_abun_modelling_data_V2.RData")

# find missing species for each folder and list all together as a full set of missing species
rls_predictions <- prediction_files[grep('rls', prediction_files)]

list_rls_missing <- lapply(rls_predictions, function(x) 
  find_missing_mofos(x, 
                     all_species = gsub(' ', '_', unique(abundance_key$TAXONOMIC_NAME))))

rls_missing <- unique(do.call(c, list_rls_missing))


# create batch file for missing species
i <- which(gsub(' ', '_',abundance_key$TAXONOMIC_NAME) %in% rls_missing)

Nodes <- c('@RunAsMultiple,@NodeSet_16')

file_out <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-models%rls%rls_config_basic.R',
                   ' ', 
                   i)

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/fitting-models/batchscripts/batch_rls_basic_debug.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)

fileConn<-file("scripts/fitting-models/batchscripts/batch_rls_CV_debug.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)



# alternative approach to seeing number of files in folders ----

fils <- list.files("results/rls/model_abunocc_2stage", pattern=".RData", full.names = TRUE, recursive = TRUE)

data_frame(dir = dirname(fils)) %>% 
  count(dir) %>% 
  data.frame() %>% 
  filter(n < 50)
  
# find i
load("data/rls_abun_modelling_data_v2.RData")
taxa<-c()
for(i in 1:length(rls_abun_fitting)){
taxa[i]<-paste0(gsub(' ', '_', unique(rls_abun_fitting[[i]]$TAXONOMIC_NAME)), '.RData')
}

# missing taxa for a specific model
which(!taxa %in% list.files('results/rls/model_abunocc/log/lm/gaussian/no_ZI'))
which(!taxa %in% list.files('results/rls_sst/model_abunocc/log/brt/discrete/multinomial'))


