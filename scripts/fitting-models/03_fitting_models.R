# script to fit each of the functions to species matrix

# load in packages ---- 
lib_vect <- c('tidyverse')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)




