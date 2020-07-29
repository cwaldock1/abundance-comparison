# r script to combine spatial projections and compare scaled occupancy and abundance 

# load packages and functions ----

lib_vect <- c('tidyverse', 'ggplot2', 'rnaturalearth', 'gridExtra')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
#detach("package:raster", unload = TRUE)

# source functions with species_spearmans_spatial functions
source('scripts/figures/functions/spatial_projection_functions.R')

# estimate spearmans rank between abundance and occurrence ----

#spatial_dir = 'results/spatial_projections/rls'
#xy = na.omit(readRDS('data/rls_spatial_projection_data.rds'))
#dataset = 'rls'

rls_spear <- species_spearmans_spatial(spatial_dir = 'results/spatial_projections/rls', 
                                       xy = na.omit(readRDS('data/rls_spatial_projection_data.rds')), 
                                       dataset = 'rls')

bbs_spear <- species_spearmans_spatial(spatial_dir = 'results/spatial_projections/bbs', 
                                       xy = na.omit(readRDS('data/bbs_spatial_projection_data.rds')), 
                                       dataset = 'bbs')

rls_spear$dataset <- 'reef-life survey'
bbs_spear$dataset <- 'breeding-bird survey'

spearmans <- rbind(rls_spear, bbs_spear)

# subset to only species with high quality models
high_performance_species <- readRDS('results/high_performance_species.RDS')
spearmans <- spearmans %>% filter(species_name %in% gsub(' ','_',high_performance_species))

# create plot comparing spearmans ranks
spearmans_2 <- spearmans %>% 
  unique() %>% 
  group_by(dataset) %>% 
  mutate(median = median(cor), 
         upr    = quantile(cor, 0.9), 
         lwr    = quantile(cor, 0.1), 
         mean   = mean(cor), 
         sd     = sd(cor))


spearmans_plot <- ggplot(spearmans_2) + 
  geom_histogram(aes(x = cor, col = dataset, fill = dataset), alpha = 0.5) + 
  geom_point(data = spearmans_2 %>% dplyr::select(dataset, mean) %>% unique(), aes(x = mean, 
                                                                                   y = c(3, 7.5), size = 2)) +
  geom_segment(data = spearmans_2 %>% dplyr::select(dataset, mean, sd) %>% unique(), aes(x = mean-sd, xend = mean+sd, 
                                                                                   y = c(3, 7.5), yend = c(3, 7.5)))+
  theme_bw() + 
  theme(legend.position = 'none',
        aspect.ratio = 0.5, 
        panel.grid = element_blank(),  
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
        strip.background = element_blank())  + 
  facet_wrap(~dataset, scales = 'free') +
  scale_colour_manual(values = c(viridis::viridis(10, option = 3)[5], viridis::viridis(10, option = 7)[5])) + 
  scale_fill_manual(values = c(viridis::viridis(10, option = 3)[5], viridis::viridis(10, option = 7)[5])) + 
  xlab('spearmans rank between spatial projection \n of abundance and occurrence') + 
  ylab(NULL)

pdf(file = 'figures/spatial_projections/spearmans_cor.pdf', height = 5, width = 10)
spearmans_plot
dev.off()

# write spatial projection spearmans rank to file to use in the script 05_occ_rf_varimp.R
saveRDS(spearmans_2, file = 'results/spatial_species_spearmans_rank.RDS')
