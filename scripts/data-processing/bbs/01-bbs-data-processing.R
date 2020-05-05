# Aim: Process data from the USA Breeding Bird Survey available at https://www.pwrc.usgs.gov/BBS/RawData/Choose-Method.cfm
# Author: Conor Waldock

# packages ----

lib_vect <- c('tidyverse', 'plyr', 'data.table')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# read in bbs raw data ----


# list the raw directory files

bbs_files <- list.files('/Users/cwaldock/Dropbox/ETH_REEF_FUTURES/abundance-comparison/raw-data/breeding-bird-survey-usa/50-StopData/1997ToPresent_SurveyWide', full.names = T)

# create abundance column across all stops

bbs_all <- lapply(bbs_files, function(x){
  
  # read csv
  
  bbs <- read.csv(x)
  
  # sum abundance across all stops
  
  abundance <- pbapply::pbapply(bbs[paste0('Stop', 1:50)], 1, sum)
  
  # remove stops
  
  bbs <- bbs[-which(colnames(bbs) %in% paste0('Stop', 1:50))]
  
  # add abundance to data
  
  bbs$Num <- as.numeric(as.character(abundance))
  
  return(bbs)
  
})

bbs_all_2 <- do.call(plyr::rbind.fill, bbs_all)

head(bbs_all_2)

# fix state code

bbs_all_2$StateNum[is.na(bbs_all_2$StateNum)] <- bbs_all_2$statenum[is.na(bbs_all_2$StateNum)]

# filter for quality ---- 

# read in routes and survey information (i.e., weather) data

route   <- read.csv('raw-data/breeding-bird-survey-usa/50-StopData/routes.csv')
weather <- read.csv('raw-data/breeding-bird-survey-usa/50-StopData/Weather.csv')

# read in species names and modify

sp_list <- read.delim('raw-data/breeding-bird-survey-usa/50-StopData/SpeciesList.txt')
sp_list <- sp_list %>% select(AOU, Genus, Species) %>% mutate(TAXONOMIC_NAME = paste(.$Genus, .$Species))

# filter by run type = 1

bbs_filter <- left_join( 
  
  weather %>% 
    filter(RunType == 1) %>% 
    dplyr::select(RouteDataID, RunType, Year, Month, Day) %>% 
    unique(),
  
  bbs_all_2
  
  ) %>% 
  
  # remove unwanted columns 
  
  select(RouteDataID, CountryNum, StateNum, Route, RunType, AOU, Num, Year, Month, Day) %>% 
  
  na.omit() %>% 
  
  # join by country, state and route numbers
  
  left_join(
    ., 
    
    route %>% select(CountryNum, StateNum, Route, Latitude, Longitude)
    
  ) %>% 
  
  select(RouteDataID, CountryNum, StateNum, Route, RunType, Latitude, Longitude, Year, Month, Day, AOU, Num)

# create a route identifier that is the unique latitude and longitude for all routes and equivalent to a site in the RLS data
bbs_filter$SiteCode <- with(bbs_filter, paste(CountryNum, StateNum, Route, Latitude, Longitude, sep = '_'))
# 4660 sites

# aggregate by year
bbs_filter_2 <- as_tibble(setDT(bbs_filter)[ , .(Num_mean = round(mean(Num, na.rm = T))), by = c('SiteCode', 'AOU')]) %>% 
  dplyr::rename(Num = Num_mean) %>% 
  unique() %>% 
  left_join(., bbs_filter %>% dplyr::select(CountryNum, StateNum, SiteCode, Longitude, Latitude) %>% unique())

# clean species data ----

# clean species data
bbs_filter_2 <- left_join(bbs_filter_2, sp_list %>% dplyr::select(AOU, Genus, Species)) 

to_remove <- c(which(grepl('spp.', bbs_filter_2$Species, fixed = T)), 
  which(grepl('sp.', bbs_filter_2$Species, fixed = T)), 
  which(grepl('/', bbs_filter_2$Species, fixed = T)), 
  which(grepl(' x ', bbs_filter_2$Species, fixed = T)), 
  which(grepl(' or ', bbs_filter_2$Species, fixed = T)), 
  which(grepl('blue', bbs_filter_2$Species, fixed = T)))

bbs_filter_2 <- bbs_filter_2[-to_remove,]

# assign species names
bbs_filter_2$TAXONOMIC_NAME <- paste(bbs_filter_2$Genus, 
                                     sub('\\ .*', '', unique(bbs_filter_2$Species[grep(' ', bbs_filter_2$Species)])), # removes sub-species
                                     sep = ' ')

# filter to ensure all species have 50 presences per species 
species_50 <- bbs_filter_2 %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(n_per_species = length(unique(.$SiteCode))) %>% 
  unnest(n_per_species) %>% 
  filter(n_per_species > 50) %>% 
  .$TAXONOMIC_NAME

# filter species to only those with 50 records
bbs_filter_2 <- bbs_filter_2 %>% filter(TAXONOMIC_NAME %in% species_50)

# filter by the specific species I want to investigate ----

# spread data for estimating frequency 

bbs_freq <- bbs_filter_2 %>% 
  # fill with 0s
  reshape2::dcast(SiteCode ~ TAXONOMIC_NAME, 
                  value.var = 'Num', 
                  fun.aggregate = function(x) mean(x, na.rm = T), 
                  fill = 0)

# apply frequency estimates

bbs_freq2 <- bbs_freq %>% 
  
  # gather back together
  gather(key = 'TAXONOMIC_NAME', value = 'Num', -SiteCode) %>% 
  
  # join in state information
  left_join(., bbs_filter_2 %>% dplyr::select(SiteCode, StateNum) %>% unique()) %>% 
  
  # group and average frequency within a state gives a better idea of occupancy rates
  group_by(TAXONOMIC_NAME, StateNum) %>% 
  
  # estimate states where the species aren't present and remove these species-state combinations
  nest() %>% 
  mutate(remove_values = purrr::map(data, ~sum(.$Num, na.rm = T))) %>% 
  filter(remove_values != 0) %>% 
  unnest(remove_values) %>% 

  # estimate average frequency within a states
  group_by(TAXONOMIC_NAME, StateNum) %>% 
  mutate(frequency = purrr::map(data, ~sum(.$Num!=0)/length(.$Num))) %>% 
  unnest(frequency) %>% 
  dplyr::select(-data) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(frequency = mean(.$frequency, na.rm = T)) %>% 
  unnest(frequency)
  
hist(bbs_freq2$frequency)


# average abundance when present

bbs_abun <- bbs_filter_2 %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(mean_abundance = mean(.$Num, na.rm = T)) %>% 
  unnest(mean_abundance) %>% 
  mutate(TAXONOMIC_NAME = as.character(.$TAXONOMIC_NAME))

hist(log(bbs_abun$mean_abundance))

# combine into single dataframe
species_properties <- left_join(bbs_abun, bbs_freq2)

# remove species with mean abundance < 3
species_properties <- species_properties %>% filter(mean_abundance >=3) %>% filter(TAXONOMIC_NAME %in% species_50)

# estimate percentiles and subset
species_properties$mean_abundance_perc <- ecdf(species_properties$mean_abundance)(species_properties$mean_abundance)
species_properties$frequency_perc      <- ecdf(species_properties$frequency)(species_properties$frequency)

# save species properties object
saveRDS(species_properties, 'data/bbs_species_properties.RDS')

# high abundance high frequency 
species_properties %>% 
  filter(mean_abundance_perc > 0.7, mean_abundance_perc < 0.9, 
         frequency_perc > 0.7, frequency_perc < 0.9) %>% 
  .[order(.$mean_abundance_perc,.$frequency_perc, decreasing = T),] %>% 
  .$TAXONOMIC_NAME %>% as.character %>% .[1:10] -> ha_hf

# high abundance low frequency 
species_properties %>% filter(mean_abundance_perc > 0.7, mean_abundance_perc < 0.90,  # increased to get 50 species
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


# # filter to species that we are interested in (defined above)

# ALL SPECIES create species-specific datasets including 0s from an object with all available sites ----

source('scripts/data-processing/functions/get_buffered_absences.R')

# obtain object with all sites 
bbs_sites <- bbs_filter_2 %>% select(SiteCode, Latitude, Longitude) %>% unique()  %>% dplyr::rename(SiteLatitude = Latitude, SiteLongitude = Longitude)

# obtain only numerical abundance from our focal species
bbs_site_abun <- bbs_filter_2 %>% select(SiteCode, Latitude, Longitude, TAXONOMIC_NAME, Num)

# change names to match across datasets
bbs_site_abun <- bbs_site_abun %>% dplyr::rename(SiteLatitude = Latitude, SiteLongitude = Longitude) %>% ungroup()

# create a list object that contains absences and presences
bbs_sum_absences <- lapply(1:length(unique(bbs_site_abun$TAXONOMIC_NAME)), FUN = function(x){
  
  print(x)
  
  # buffer round the absences 
  buffered_absence <- get_buffered_absences(presences = 
                          bbs_site_abun %>% 
                          filter(TAXONOMIC_NAME == unique(bbs_site_abun$TAXONOMIC_NAME)[x]), 
                        sites = bbs_sites, 
                        x_name = 'SiteLongitude',
                        y_name = 'SiteLatitude',
                        sp_name = 'TAXONOMIC_NAME', 
                        coord_ref = '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0 +units=m'
                        )
  
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
  dir.create('data/bbs_all_basic/')

  # save RDS
  saveRDS(list_outputs, paste0('data/bbs_all_basic/', 
                               as.character(gsub(' ','_',unique(buffered_absence$TAXONOMIC_NAME))),
                               '.RDS'))

  
})

# 50 SPECIES create species-specific datasets including 0s from an object with all available sites ----
 
# filter to the species of interest
bbs_filter_spp <- bbs_filter_2 %>% filter(TAXONOMIC_NAME %in% c(aa_af, ha_hf, ha_lf, la_hf, la_lf))

# obtain object with all sites 
bbs_sites <- bbs_filter_spp %>% select(SiteCode, Latitude, Longitude) %>% unique()  %>% dplyr::rename(SiteLatitude = Latitude, SiteLongitude = Longitude)

# obtain only numerical abundance from our focal species
bbs_site_abun <- bbs_filter_spp %>% select(SiteCode, Latitude, Longitude, TAXONOMIC_NAME, Num)

# change names to match across datasets
bbs_site_abun <- bbs_site_abun %>% dplyr::rename(SiteLatitude = Latitude, SiteLongitude = Longitude)

# buffer species absences
bbs_sum_absences <- lapply(1:length(unique(bbs_site_abun$TAXONOMIC_NAME)), FUN = function(x){
  
  print(x)
  
                           # buffer round the absences 
                            get_buffered_absences(presences = 
                                              bbs_site_abun %>% 
                                              filter(TAXONOMIC_NAME == unique(bbs_site_abun$TAXONOMIC_NAME)[x]), 
                                            sites = bbs_sites, 
                                            x_name = 'SiteLongitude',
                                            y_name = 'SiteLatitude',
                                            sp_name = 'TAXONOMIC_NAME', 
                                            coord_ref = '+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0 +units=m')
  }
  )


# create fitting and validation sets

# I should take 80%/20% of BOTH the abundance and absences (rather than a subsample of the whole dataset)
set.seed(123)

bbs_model_data <- lapply(bbs_sum_absences, function(x){
  
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
  #names(list_outputs) <- unique(bbs_sum_absences)
  
  return(list_outputs)
  
})


# create outputs to save for the 50 species runs - this is the random validation
bbs_abun_fitting    <- lapply(bbs_model_data, function(x) x[[1]][[1]])
bbs_abun_validation <- lapply(bbs_model_data, function(x) x[[1]][[2]])
bbs_abun_list       <- bbs_model_data
bbs_abun            <- rbind(do.call(rbind, bbs_abun_fitting), do.call(rbind, bbs_abun_validation))


# create a key for the common and rare classes for each species
abundance_key <- left_join(data.frame(TAXONOMIC_NAME  = c(aa_af, ha_hf, ha_lf, la_hf, la_lf),
                                      abundance_class = ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% aa_af, 'average',
                                                               ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_hf, 'h_abun h_freq', 
                                                                      ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_lf, 'h_abun l_freq', 
                                                                             ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% la_hf, 'l_abun h_freq', 
                                                                                    'l_abun l_freq'))))), 
                           species_properties)

dir.create('data/bbs_50_basic/')
save(bbs_abun, # whole data object as datafrane
     bbs_abun_list, # data object as list
     bbs_abun_fitting, # listed fitting data
     bbs_abun_validation, # list validation data
     abundance_key, file = 'data/bbs_50_basic/bbs_abun_modelling_data.RData')


# summarise for table ----


abundance_key %>% 
  group_by(abundance_class) %>% 
  do(mean_frequency = mean(.$frequency),
     sd_frequency   = sd(.$frequency), 
     mean_abundance = mean(.$mean_abundance),
     sd_abundance   = sd(.$mean_abundance)) %>% 
  unnest()

#   abundance_class mean_frequency sd_frequency mean_abundance sd_abundance
# <fct>                    <dbl>        <dbl>          <dbl>        <dbl>
# 1 average                 0.0795      0.00443           7.90       0.326 
# 2 h_abun h_freq           0.145       0.0327           17.8        1.41  
# 3 h_abun l_freq           0.0557      0.00538          18.4        2.19  
# 4 l_abun h_freq           0.160       0.0237            4.01       0.192 
# 5 l_abun l_freq           0.0468      0.00690           3.57       0.0881


left_join(bbs_abun %>% mutate(TAXONOMIC_NAME = as.character(.$TAXONOMIC_NAME)), abundance_key %>% mutate(TAXONOMIC_NAME = as.character(.$TAXONOMIC_NAME))) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(mean_obs = length(unique(paste(.$Latitude, .$Longitude, .$SiteCode)))) %>% 
  unnest() %>% 
  .$mean_obs %>% mean

# number of presence and absenes
# 3936.98


left_join(bbs_abun %>% mutate(TAXONOMIC_NAME = as.character(.$TAXONOMIC_NAME)), abundance_key %>% mutate(TAXONOMIC_NAME = as.character(.$TAXONOMIC_NAME))) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  nest() %>% 
  mutate(mean_obs = purrr::map(data, ~length(unique(paste(.$Latitude, .$Longitude, .$SiteCode))))) %>% 
  unnest() %>% 
  group_by(abundance_class) %>% 
  do(mean_ac = mean(.$mean_obs), 
     sd_ac   = sd(.$mean_obs)) %>% 
  unnest() %>% data.frame 

# number of presence and absences
# abundance_class  mean_ac    sd_ac
# 1         average 4251.228 194.5585
# 2   h_abun h_freq 4016.011 789.7255
# 3   h_abun l_freq 4097.136 378.8375
# 4   l_abun h_freq 4275.763 154.6182
# 5   l_abun l_freq 3613.603 812.3752

