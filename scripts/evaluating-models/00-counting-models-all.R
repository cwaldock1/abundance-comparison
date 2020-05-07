
# work on the server

setwd('/Volumes/Simulation/conor/abundance-comparison')

# Load packages ----

lib_vect <- c('stringi', 'tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# Load functions ----

list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


# get lists of the relevant prediction folders ----

# get all the relevant prediction folders

files <- unlist(lapply(list.files('results', full.names = T, pattern = '_all'), function(x) list.dirs.depth.n(x, n = 2)))

prediction_files <- files[grepl('predictions', files)]




