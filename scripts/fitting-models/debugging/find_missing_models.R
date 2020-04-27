# script to find missing species and create an individual batch script for the species which did not work on the first run. 
# there are sometimes sporadic errors on the server, so when you re-run the species it works, and it works locally but not first time on the server. 

# set the working directory to the server ----
setwd('/Volumes/Simulation/conor/abundance-comparison')

# load packages ----

lib_vect <- c('stringi', 'tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# load functions ----

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


# get predictions ----

# get all the relevant prediction folders

files <- unlist(lapply(list.files('results', full.names = T)[1:4], function(x) list.dirs.depth.n(x, n = 2)))
prediction_files <- files[grepl('predictions', files)]


# find missing bbs species and create batch file ----

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



# find missing rls species and create batch file ----

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
