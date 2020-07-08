# function to produce spatial polygons from presences ----

# dataset = 'rls'
# species_name = 'Chaetodon_lunula'
# save_dir = 'data/species_range_convex_hull'

convex_hull_creation <- function(species_name, 
                                 dataset, 
                                 save_dir){ 
  
  require(sp)
  require(rgeos)
  require(rgdal)
  require(tidyverse)
  
  # require(alphahull) option for alpha hulls, the value for a becomes arbitrary.. 
  
  # read in species level data
  species_data <- readRDS(list.files(paste0('data/', dataset, '_all_basic/'), pattern = paste0(unique(species_name), '.RDS'), full.names = T))$fitting
  
  # remove species absences
  xy_presences <- species_data %>% filter(Num > 0) %>% dplyr::select(SiteLongitude, SiteLatitude)
  
  # create convex hull
  con.hull.pos <- chull(xy_presences) # find positions of convex hull
  
  # get coordiantes of hull
  con.hull <- rbind(xy_presences[con.hull.pos,],xy_presences[con.hull.pos[1],]) # get coordinates for convex hull
  
  # convert hull to polygon
  sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(con.hull)), ID=1)))
  
  # add CRS to polygon
  crs(sp_poly) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  # buffer polygon by 2° latitude and longitude
  sp_poly_buffer <- gBuffer(sp_poly, width = 2)
  
  #plot(sp_poly_buffer); points(xy_presences)
  
  # save output into new folder in data
  dir.create(path = paste0(save_dir, '/', dataset))
  saveRDS(sp_poly_buffer, paste0(paste0(save_dir, '/', dataset), '/', species_name, '.RDS'))
  #writeOGR(sp_poly_buffer, dsn = paste0(save_dir, '/', species_name), layer="chull", driver="ESRI Shapefile")
  
}

