# script to match up rls sites with environmental rasters on a standardised grid 

# load packages ----

lib_vect <- c('tidyverse', 'rgdal', 'raster')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# load in rls sites for matching ----



