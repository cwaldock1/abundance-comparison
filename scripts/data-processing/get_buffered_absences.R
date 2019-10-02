# function to extract 0s from buffers around species presences to obtain absences within a species range and range edges

# objects are created in lines 190-193 in script 01_rls-data-processing.R
#presences = rls_site_abun %>% filter(TAXONOMIC_NAME == unique(rls_site_abun$TAXONOMIC_NAME)[1])
#sites = rls_sites

get_buffered_absences <- function(presences, 
                                  sites, 
                                  buffer_size = 10){
  
  require(tidyverse)
  require(sp)
  require(rgeos)
  require(data.table)
  
  # convert all sites to a spatial object
  sites_coord <- coordinates(cbind(sites$SiteLongitude, sites$SiteLatitude))
  sites_coord_sp <- SpatialPoints(sites_coord, proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +units=m'))
  
  # covert species presences to a spatial object
  presences_coord <- coordinates(cbind(presences$SiteLongitude, presences$SiteLatitude))
  presences_coord_sp <- SpatialPoints(presences_coord, proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 +units=m'))
  
  # buffer the species presences
  buffered_presences <- gBuffer(presences_coord_sp, width = buffer_size, byid = T)
  
  # extract the buffered absences
  buffered_absences <- data.frame(sites_coord_sp[buffered_presences])

  # join together the location of buffered 0s with the site information
  buffered_absences_df <- left_join(buffered_absences, sites, 
                                    by = c('coords.x1' = 'SiteLongitude', 'coords.x2' = 'SiteLatitude')) %>% 
                             na.omit() %>% 
                             unique()
  
  # rename columns to match dataframe
  names(buffered_absences_df)[1:2] <- c('SiteLongitude', 'SiteLatitude')
  
  # create absences
  buffered_absences_df$TAXONOMIC_NAME <- unique(presences$TAXONOMIC_NAME)
  buffered_absences_df$Num <- 0
  
  # bind together abundance and absence information
  presences_v2  <- rbind(presences, buffered_absences_df)
  
  # sum over presences and absences using data.table for speed 
  presences_v2_DT <- data.table(presences_v2)
  setkey(presences_v2_DT, SiteCode, SiteLatitude, SiteLongitude, TAXONOMIC_NAME)
  presences_v2_DT_2 <- presences_v2_DT[, j = list(Num = sum(Num, na.rm = T)), 
                       by = key(presences_v2_DT)]
  
  # convert back to dataframe
  presences_v2 <- as_tibble(presences_v2_DT_2)
  
  return(presences_v2)
  
}




