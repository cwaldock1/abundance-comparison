# R script to pre-process abundance data from raw RLS surveys provided by Rick Stuart-Smith and Graham Edgar on 08/08/2019
# The output of this script is a data object that has 
# 1. species by xy matrix
# 2. species by prevalence and mean abundance matrix and class

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
  unnest(Num_sum, Biomass_sum) %>% 
  unnest(data) %>% 
  select(-Num, -Biomass, -Sizeclass) %>% 
  rename(Num = Num_sum, Biomass = Biomass_sum) %>% 
  unique() %>% 

  # average abundance and biomass over surveys
  group_by(SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME) %>% 
  nest() %>% 
  mutate(Num_mean = purrr::map(data, ~round(mean(.$Num, na.rm = T))), 
         Biomass_mean = purrr::map(data, ~mean(.$Biomass, na.rm = T))) %>% 
  unnest(Num_mean, Biomass_mean) %>% 
  unnest(data) %>% 
  select(-Num, -Biomass, -SurveyID) %>% 
  rename(Num = Num_mean, Biomass = Biomass_mean) %>% 
  unique()

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
  filter(remove_values != 0) %>% 
  unnest(remove_values) %>% 
  unnest(data) %>% 
  group_by(TAXONOMIC_NAME, Ecoregion) %>% 
  do(frequency = sum(.$Num!=0)/length(.$Num)) %>% 
  unnest(frequency) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(frequency = mean(.$frequency, na.rm = T)) %>% 
  unnest(frequency)
hist(freq$frequency)

# average abundance when present
abun <- rls_sum %>% 
  filter(Num > 0) %>% 
  group_by(TAXONOMIC_NAME) %>% 
  do(mean_abundance = mean(.$Num, na.rm = T)) %>% 
  unnest(mean_abundance)
hist(log(abun$mean_abundance))

# combine into single dataframe
species_properties <- left_join(abun,freq)

# remove species with mean abundance < 3
species_properties <- species_properties %>% filter(mean_abundance >=3) %>% filter(TAXONOMIC_NAME %in% species_50)

# remove unknown species
species_properties <- species_properties[-which(grepl('spp.', species_properties$TAXONOMIC_NAME, fixed = T)),]
species_properties <- species_properties[-which(grepl('sp.', species_properties$TAXONOMIC_NAME, fixed = T)),]

# estimate percentiles and subset
species_properties$mean_abundance_perc <- ecdf(species_properties$mean_abundance)(species_properties$mean_abundance)
species_properties$frequency_perc      <- ecdf(species_properties$frequency)(species_properties$frequency)

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


# create wide dataset from all sites ----

# create a species by site matrix and aggregate across size-classes
rls_subset_wide <- rls_sum %>% 
  select(SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME, Num) %>% 
  #reshape2::dcast(SiteCode ~ TAXONOMIC_NAME, 
  #                fun.aggregate = function(x) mean(x, na.rm = T), 
  #                fill = 0) %>% 
  spread(key = TAXONOMIC_NAME, value = Num, fill = 0) %>% 
  .[,which(colnames(.) %in% c('SiteCode', 'SiteLatitude', 'SiteLongitude', 
                              aa_af, ha_hf, ha_lf, la_hf, la_lf))]

# check that all species are present in x number of cells 
rls_subset_wide

rls_abun <- rls_subset_wide

# create fitting and validation sets
fitting_set <- sample(1:nrow(rls_abun), round(0.8*nrow(rls_abun)), replace = F)
fitting_set
rls_abun_fitting    <- rls_abun[sort(fitting_set), ]
rls_abun_validation <- rls_abun[-sort(fitting_set), ]
sum(rls_abun_fitting$SiteCode %in% rls_abun_validation$SiteCode) # double check for non-overlap.

# create a key for the common and rare classes for each species
abundance_key <- left_join(data.frame(TAXONOMIC_NAME  = c(aa_af, ha_hf, ha_lf, la_hf, la_lf),
           abundance_class = ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% aa_af, 'average',
                                    ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_hf, 'h_abun h_freq', 
                                           ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% ha_lf, 'h_abun l_freq', 
                                                  ifelse(c(aa_af, ha_hf, ha_lf, la_hf, la_lf) %in% la_hf, 'l_abun h_freq', 
                                                         'l_abun l_freq'))))), 
           species_properties)

save(rls_abun, rls_abun_fitting, rls_abun_validation, abundance_key, file = 'data/rls_abun_modelling_data.RData')
