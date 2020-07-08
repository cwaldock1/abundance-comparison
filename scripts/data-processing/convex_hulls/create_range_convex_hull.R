
# load in function ----

source('scripts/data-processing/convex_hulls/convex_hull_function.R')

# ----

# rls list all species
all_sp <- gsub('.RDS', '', list.files('results/spatial_projections/rls'))

lapply(1:length(all_sp), function(x){
  print(which(all_sp %in% all_sp[x]))
  convex_hull_creation(species_name = gsub('.RDS','', all_sp[x]),
                       dataset = 'rls',
                       save_dir = 'data/species_range_convex_hull')})

# rls list all species
all_sp <- gsub('.RDS', '', list.files('results/spatial_projections/bbs'))

lapply(1:length(all_sp), function(x){
  print(which(all_sp %in% all_sp[x]))
  convex_hull_creation(species_name = gsub('.RDS','', all_sp[x]),
                       dataset = 'bbs',
                       save_dir = 'data/species_range_convex_hull')})
