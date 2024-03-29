# R script to pre-process abundance data from raw RLS surveys provided by Rick Stuart-Smith and Graham Edgar on 08/08/2019
# The output of this script is a data object that has 
# 1. species by xy matrix
# 2. species by prevalence and mean abundance matrix and class
options("scipen"=5, "digits"=5)

# load packages ----
lib_vect <- c('tidyverse', 'summarytools', 'rgdal')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)


# read in rls data ----
rls_raw  <- read.csv('raw-data/RLS-spatial-fish-data-extract.csv')
rls_meta <- read.csv('raw-data/RLS-site-subset-metadata.csv')
rls_meta$SiteCode <- as.character(rls_meta$SiteCode)
meow     <- readOGR(dsn = "raw-data/MEOW", layer = "meow_ecos")


# data structure ----

names(rls_raw)
#SurveyDate, full RLS data have been subset to the survey nearest to 2012.
#Biomass has a size correction
#SPECIES_NAME is the field we’d use for richness calculations, as it includes names of undescribed or unidentified species, 
#TAXONOMIC_NAME has these rolled up to genus or whatever

head(rls_raw)

# convert all to appropriate classes
rls_raw$SiteCode       <- as.character(rls_raw$SiteCode)
rls_raw$SiteLatitude   <- as.numeric(as.character(rls_raw$SiteLatitude))
rls_raw$SiteLongitude  <- as.numeric(as.character(rls_raw$SiteLongitude))
rls_raw$SurveyID       <- as.character(rls_raw$SurveyID)
rls_raw$SurveyDate     <- as.character(rls_raw$SurveyDate)
rls_raw$Depth          <- as.numeric(as.character(rls_raw$Depth))
rls_raw$Num            <- as.numeric(as.character(rls_raw$Num))
rls_raw$Sizeclass      <- as.character(rls_raw$Sizeclass)
rls_raw$Biomass        <- as.numeric(as.character(rls_raw$Biomass))
rls_raw$SPECIES_NAME   <- as.character(rls_raw$SPECIES_NAME)
rls_raw$TAXONOMIC_NAME <- as.character(rls_raw$TAXONOMIC_NAME)

view(dfSummary(rls_raw[-3]))

# find the ecoregion for each point and input to data ----

coordiantes <- data.frame(x = rls_meta$SiteLongitude, 
                          y = rls_meta$SiteLatitude)
coordiantes <- SpatialPoints(coordiantes, proj4string = meow@proj4string)
Ecoregions <- over(coordiantes, meow)
rls_meta$Ecoregion <- Ecoregions$ECOREGION

# aggregate surveys to a site level ----

# are there multiple latitudes and longitudes within a site: 
# safe to aggregate over sites but must average depths too
rls_raw %>% 
  select(SiteCode, SurveyID, SiteLatitude, SiteLongitude, Depth) %>% 
  unique() %>% 
  group_by(SiteCode) %>% 
  do(n_latlong = length(unique(paste0(.$SiteLatitude, .$SiteLongitude))), 
     n_depth   = length(unique(.$Depth))) %>% 
  unnest() %>% 
  ungroup() %>% 
  .[,2:3] %>% 
  lapply(., range)

# subset to sites in australia with bounding box of From Lat/Lon -50, 110 to Lat/Lon 3, 165 using rls_meta$SiteCode
# aggregate abundance and biomass over site for each taxonomic names
rls_sum <- rls_raw %>% 
  
  # filter out sizes not in the autralian subset of sites
  filter(SiteCode %in% unique(rls_meta$SiteCode)) %>% 
  select(SiteCode, SiteLatitude, SiteLongitude, Sizeclass, TAXONOMIC_NAME, SurveyID,
         Num, Biomass) %>% 
  
  # sum up abundance and biomass size classes within a survey
  group_by(SurveyID, TAXONOMIC_NAME) %>% 
  nest() %>% 
  mutate(Num_sum = purrr::map(data, ~round(sum(.$Num, na.rm = T))), 
         Biomass_sum = purrr::map(data, ~sum(.$Biomass, na.rm = T))) %>% 
  ungroup() %>% 
  unnest(c(Num_sum, Biomass_sum)) %>% 
  unnest(data) %>% 
  select(-Num, -Biomass, -Sizeclass) %>% 
  dplyr::rename(Num = Num_sum, Biomass = Biomass_sum) %>% 
  unique() %>% 

  # average abundance and biomass over surveys i.e., gets rid of yearly variation
  group_by(SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME) %>% 
  nest() %>% 
  mutate(Num_mean = purrr::map(data, ~round(mean(.$Num, na.rm = T))), 
         Biomass_mean = purrr::map(data, ~mean(.$Biomass, na.rm = T))) %>% 
  unnest(c(Num_mean, Biomass_mean)) %>% 
  unnest(data) %>% 
  select(-Num, -Biomass, -SurveyID) %>% 
  dplyr::rename(Num = Num_mean, Biomass = Biomass_mean) %>% 
  unique()

saveRDS(rls_sum, 'data/scrap_rls_sum_object.rds')

# calculate number of species available for analysis 
length(unique(rls_sum$TAXONOMIC_NAME))

rls_sum <- readRDS('data/scrap_rls_sum_object.rds')

# filter to ensure all species have 50 presences per species 
species_50 <- rls_sum %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(n_per_species = length(unique(.$SiteCode))) %>% 
  unnest(n_per_species) %>% 
  filter(n_per_species > 50) %>% 
  .$TAXONOMIC_NAME

# estimating properties of species abundances ----

# average frequency per species across the ecoregion it is present within
freq <- rls_sum %>% 
  reshape2::dcast(SiteCode ~ TAXONOMIC_NAME, 
                  value.var = 'Num', 
                  fun.aggregate = function(x) mean(x, na.rm = T), 
                  fill = 0) %>% 
  gather(key = 'TAXONOMIC_NAME', value = 'Num', -SiteCode) %>% 
  left_join(., rls_meta %>% select(SiteCode, Ecoregion) %>% unique()) %>% 
  group_by(TAXONOMIC_NAME, Ecoregion) %>% 
  nest() %>% 
  mutate(remove_values = purrr::map(data, ~sum(.$Num, na.rm = T))) %>% 
  ungroup() %>% 
  filter(remove_values != 0) %>% 
  unnest(remove_values) %>% 
  unnest(data) %>% 
  group_by(TAXONOMIC_NAME, Ecoregion) %>% 
  do(frequency = sum(.$Num!=0)/length(.$Num)) %>% 
  ungroup() %>% 
  unnest(frequency) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(frequency = mean(.$frequency, na.rm = T)) %>% 
  unnest(frequency) %>% 
  ungroup()
hist(freq$frequency)

# average abundance when present
abun <- rls_sum %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(mean_abundance = mean(.$Num, na.rm = T)) %>% 
  unnest(mean_abundance)
hist(log(abun$mean_abundance))

# count number of observations per species 
n_per_species <- rls_sum %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(n_per_species = length(unique(.$SiteCode))) %>% 
  unnest(n_per_species)

# combine into single dataframe
species_properties <- left_join(left_join(abun,freq), n_per_species)

# remove species with less that 50 records
species_properties <- species_properties %>% filter(TAXONOMIC_NAME %in% species_50)

# remove unknown species
species_properties <- species_properties[-which(grepl('spp.', species_properties$TAXONOMIC_NAME, fixed = T)),]
species_properties <- species_properties[-which(grepl('sp.', species_properties$TAXONOMIC_NAME, fixed = T)),]

# estimate percentiles and subset
species_properties$mean_abundance_perc <- ecdf(species_properties$mean_abundance)(species_properties$mean_abundance)
species_properties$frequency_perc      <- ecdf(species_properties$frequency)(species_properties$frequency)

# save species object
saveRDS(species_properties, 'data/rls_species_properties.RDS')


# high abundance high frequency 
species_properties %>% 
  filter(mean_abundance_perc > 0.7, mean_abundance_perc < 0.9, 
         frequency_perc > 0.7, frequency_perc < 0.9) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = T),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> ha_hf
  
# high abundance low frequency 
species_properties %>% filter(mean_abundance_perc > 0.7, mean_abundance_perc < 0.9, 
                              frequency_perc > 0.1, frequency_perc < 0.3) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = T),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> ha_lf

# low abundance high frequency 
species_properties %>% filter(mean_abundance_perc > 0.1, mean_abundance_perc < 0.3, 
                              frequency_perc > 0.7, frequency_perc < 0.9) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = F),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> la_hf

# low abundance low frequency 
species_properties %>% filter(mean_abundance_perc > 0.05, mean_abundance_perc < 0.3, 
                              frequency_perc > 0.05, frequency_perc < 0.3) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = F),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> la_lf

# average abundance and average frequency
species_properties %>% filter(mean_abundance_perc > 0.4, mean_abundance_perc < 0.6, 
                              frequency_perc > 0.4, frequency_perc < 0.6) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = T),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> aa_af


# filter to species that we are interested in (defined above)
# rls_sum <- rls_sum %>% filter(TAXONOMIC_NAME %in% c(aa_af, ha_hf, ha_lf, la_hf, la_lf))

# ALL SPECIES create species-specific datasets including 0s from an object with all available sites ----

source('scripts/data-processing/functions/get_buffered_absences.R')

# obtain object with all sites 
rls_sites <- rls_meta %>% select(SiteCode, SiteLatitude, SiteLongitude) %>% unique()

# obtain only numerical abundance from our focal species
rls_site_abun <- rls_sum %>% select(SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME, Num) %>% ungroup()

# filter to suitable species
rls_site_abun <- rls_site_abun %>% filter(TAXONOMIC_NAME %in% unique(species_properties$TAXONOMIC_NAME))

# create a list object that contains absences and presences
rls_sum_absences <- lapply(1:length(unique(rls_site_abun$TAXONOMIC_NAME)), FUN = function(x){
  
  print(x)
  
  input <- rls_site_abun %>% 
    ungroup() %>% 
    dplyr::filter(TAXONOMIC_NAME == unique(rls_site_abun$TAXONOMIC_NAME)[x])
  
  # buffer round the absences 
  buffered_absence <- get_buffered_absences(presences = input, 
                                            sites = rls_sites, 
                                            x_name = 'SiteLongitude',
                                            y_name = 'SiteLatitude',
                                            sp_name = 'TAXONOMIC_NAME')
  
  # I should take 80%/20% of BOTH the abundance and absences (rather than a subsample of the whole dataset)
  all_absences  <- buffered_absence %>% filter(Num == 0)  # absences
  all_presences <- buffered_absence %>% filter(Num != 0)  # presences
  absences_sample_fitting  <- sample(1:nrow(all_absences), round(nrow(all_absences)*0.8), replace = F)
  presences_sample_fitting <- sample(1:nrow(all_presences), round(nrow(all_presences)*0.8), replace = F)
  absences_sample_validation <- setdiff(1:nrow(all_absences), absences_sample_fitting)
  presence_sample_validation <- setdiff(1:nrow(all_presences), absences_sample_fitting)
  
  all_absences[absences_sample_fitting,]     # absences - fitting
  all_presences[presences_sample_fitting,]   # presences - fitting
  all_absences[absences_sample_validation,]  # absences - validation
  all_presences[presence_sample_validation,] # presences - validation
  
  fitting    <- rbind(all_absences[absences_sample_fitting,],    all_presences[presences_sample_fitting,])
  validation <- rbind(all_absences[absences_sample_validation,], all_presences[presence_sample_validation,])
  
  fitting$set    <- 'fitting'
  validation$set <- 'validation'
  
  list_outputs <- list(fitting = fitting, validation = validation)
  
  # create output directory
  dir.create('data/rls_all_basic/')
  
  # save RDS
  saveRDS(list_outputs, paste0('data/rls_all_basic/', 
                               as.character(gsub(' ','_',unique(buffered_absence$TAXONOMIC_NAME))),
                               '.RDS'))
  
  
})

# Data summaries for all species ----

all_rls_data <- do.call(rbind, lapply(list.files('data/rls_all_basic', full.names = T), function(x) readRDS(x)$fitting))

nrow(all_rls_data)
length(unique(all_rls_data$SiteCode))
length(unique(all_rls_data$TAXONOMIC_NAME))


# create species-specific datasets including 0s from an object with all available sites ----

source('scripts/data-processing/get_buffered_absences.R')

# obtain object with all sites 
rls_sites <- rls_meta %>% select(SiteCode, SiteLatitude, SiteLongitude) %>% unique()

# obtain only numerical abundance from our focal species
rls_site_abun <- rls_sum %>% select(SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME, Num) %>% 
  ungroup() %>% 
  filter(TAXONOMIC_NAME %in% c(aa_af, ha_hf, ha_lf, la_hf, la_lf))

# create a list object that contains absences and presences
rls_sum_absences <- lapply(1:length(unique(rls_site_abun$TAXONOMIC_NAME)), FUN = function(x){
                          
                           print(x)
  
                           get_buffered_absences(presences = rls_site_abun %>% 
                                                   filter(TAXONOMIC_NAME == unique(rls_site_abun$TAXONOMIC_NAME)[x]), 
                                                 sites = rls_sites, 
                                                 x_name = 'SiteLongitude',
                                                 y_name = 'SiteLatitude',
                                                 sp_name = 'TAXONOMIC_NAME')
                           
                           })


# create fitting and validation sets ----

# I should take 80%/20% of BOTH the abundance and absences (rather than a subsample of the whole dataset)
set.seed(123)

rls_model_data <- lapply(rls_sum_absences, function(x){
  
  all_absences  <- x %>% filter(Num == 0)  # absences
  all_presences <- x %>% filter(Num != 0)  # presences
  absences_sample_fitting  <- sample(1:nrow(all_absences), round(nrow(all_absences)*0.8), replace = F)
  presences_sample_fitting <- sample(1:nrow(all_presences), round(nrow(all_presences)*0.8), replace = F)
  absences_sample_validation <- setdiff(1:nrow(all_absences), absences_sample_fitting)
  presence_sample_validation <- setdiff(1:nrow(all_presences), absences_sample_fitting)
  
  all_absences[absences_sample_fitting,]     # absences - fitting
  all_presences[presences_sample_fitting,]   # presences - fitting
  all_absences[absences_sample_validation,]  # absences - validation
  all_presences[presence_sample_validation,] # presences - validation
  
  fitting    <- rbind(all_absences[absences_sample_fitting,],    all_presences[presences_sample_fitting,])
  validation <- rbind(all_absences[absences_sample_validation,], all_presences[presence_sample_validation,])
  
  fitting$set    <- 'fitting'
  validation$set <- 'validation'
  
  list_outputs <- list(list(fitting = fitting, validation = validation))
  names(list_outputs) <- unique(x$TAXONOMIC_NAME)
  
  return(list_outputs)
  
})

# create outputs to save - this is the random validation
rls_abun_fitting    <- lapply(rls_model_data, function(x) x[[1]][[1]])
rls_abun_validation <- lapply(rls_model_data, function(x) x[[1]][[2]])
rls_abun_list       <- rls_model_data
rls_abun            <- rbind(do.call(rbind, rls_abun_fitting), do.call(rbind, rls_abun_validation))
  
# create fitting and validation sets
# fitting_set <- sample(1:nrow(rls_abun), round(0.8*nrow(rls_abun)), replace = F)
# fitting_set
# rls_abun_fitting    <- rls_abun[sort(fitting_set), ]
# rls_abun_validation <- rls_abun[-sort(fitting_set), ]
# sum(rls_abun_fitting$SiteCode %in% rls_abun_validation$SiteCode) # double check for non-overlap.

# create a key for the common and rare classes for each species
abundance_key <- left_join(data.frame(TAXONOMIC_NAME  = c(aa_af, ha_hf, ha_lf, la_hf, la_lf),
           abundance_class = ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% aa_af, 'average',
                                    ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_hf, 'h_abun h_freq', 
                                           ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_lf, 'h_abun l_freq', 
                                                  ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% la_hf, 'l_abun h_freq', 
                                                         'l_abun l_freq'))))), 
           species_properties)

dir.create('data/rls_50_basic/')
save(rls_abun, # whole data object as datafrane
     rls_abun_list, # data object as list
     rls_abun_fitting, # listed fitting data
     rls_abun_validation, # list validation data
     abundance_key, file = 'data/rls_50_basic/rls_abun_modelling_data.RData')

# summarise for table ----
abundance_key %>% 
  group_by(abundance_class) %>% 
  do(mean_frequency = mean(.$frequency),
            sd_frequency   = sd(.$frequency), 
            mean_abundance = mean(.$mean_abundance),
            sd_abundance   = sd(.$mean_abundance)) %>% 
  unnest()


#  abundance_class mean_frequency sd_frequency mean_abundance sd_abundance
#  <fct>                    <dbl>        <dbl>          <dbl>        <dbl>
#  1 average                  0.239       0.0194          11.2         2.79 
#  2 h_abun h_freq            0.379       0.0364          74.6        12.4  
#  3 h_abun l_freq            0.143       0.0150          37.7        14.9  
#  4 l_abun h_freq            0.364       0.0267           4.18        0.378
#  5 l_abun l_freq            0.135       0.0195           3.77        0.253


left_join(rls_abun, abundance_key) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(mean_obs = length(unique(paste(.$SiteLatitude, .$SiteLongitude, .$SiteCode)))) %>% 
  unnest() %>% 
  .$mean_obs %>% mean

left_join(rls_abun, abundance_key) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  nest() %>% 
  mutate(mean_obs = purrr::map(data, ~length(unique(paste(.$SiteLatitude, .$SiteLongitude, .$SiteCode))))) %>% 
  unnest() %>% 
  group_by(abundance_class) %>% 
  do(mean_ac = mean(.$mean_obs), 
     sd_ac   = sd(.$mean_obs)) %>% 
  unnest() %>% data.frame 
  

