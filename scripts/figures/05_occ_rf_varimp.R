# script to analyse and plot patterns in variable importance 



# libraries and set up---- 

lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'grid')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)





# read in variable importance results ----

# load in file names
vi_files <- list.files('results/variable_importance', recursive = T, full.names = T)

# read in all files
vi_read <- lapply(1:length(vi_files), function(x){
  
  z <- readRDS(vi_files[x])
  z[[1]]$dataset = str_split(vi_files[x], '/')[[1]][[3]]
  z[[2]]$dataset = str_split(vi_files[x], '/')[[1]][[3]]
  return(z)
  })

# seperate into variable important and prediction value objects
var_imp       <- do.call(rbind, lapply(vi_read, function(x) x[[1]]))
pred_occ_abun <- do.call(rbind, lapply(vi_read, function(x) x[[2]]))


# create plots of variable importance ----


agg_var_imp <- var_imp %>% 
  group_by(covariate, dataset) %>% 
  do(occ_imp_mean  = mean(.$occurrence_imp_rescaled, na.rm = T), 
     abun_imp_mean = mean(.$abundance_imp_rescaled, na.rm = T), 
     occ_imp_sd    = sd(.$occurrence_imp_rescaled, na.rm = T), 
     abun_imp_sd   = sd(.$abundance_imp_rescaled, na.rm = T)) %>% 
  unnest()

# edit variable names
agg_var_imp$covariate <- recode(as.character(agg_var_imp$covariate),
                                'Elevation_GEBCO' = 'depth/elevation', 
                                'Depth_GEBCO_transformed' = 'depth/elevation', 
                                'human_pop' = 'human', 
                                'human_pop_2015_50km' = 'human', 
                                'primary_forest' = 'forest', 
                                'robPCA_1' = 'climate PC1', 
                                'robPCA_2' = 'climate PC2', 
                                'robPCA_3' = 'climate PC3', 
                                'sst_mean' = 'sst', 
                                'wave_energy_mean' = 'wave energy', 
                                'reef_area_200km' = 'reef area')

agg_var_imp$dataset <- recode(as.character(agg_var_imp$dataset),
                              'rls' = 'reef-life survey', 
                              'bbs' = 'breeding-bird survey')
                                
# set colours
# select colours
colours = colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(length(unique(agg_var_imp$covariate)))


agg_var_imp_plot <- ggplot(data = agg_var_imp) + 
  geom_abline() + 
  geom_segment(aes(x = occ_imp_mean, xend = occ_imp_mean,
                   y = abun_imp_mean + abun_imp_sd, yend = abun_imp_mean - abun_imp_sd, col = covariate)) +
  geom_segment(aes(y = abun_imp_mean, yend = abun_imp_mean,
                   x = occ_imp_mean + occ_imp_sd, xend = occ_imp_mean - occ_imp_sd, col = covariate)) +
  geom_point(aes(x = occ_imp_mean, y = abun_imp_mean, col = covariate), size=5) + 
  ggrepel::geom_text_repel(aes(x = occ_imp_mean, y = abun_imp_mean, label = covariate, col = covariate), 
                           box.padding = 3.5, 
                           force = 7,
                           point.padding = 0,
                           alpha = 1, 
                           size = 5) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        aspect.ratio = 1, 
        panel.grid = element_blank(),  
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 14),
        strip.background = element_blank())  + 
  facet_wrap(~dataset) + 
  xlab('occurrence variable importance') + 
  ylab('abundance variable importance') + 
  scale_colour_manual(values = colours)

# create directory and save figure
dir.create(path = 'figures/variable_importance/var_imp_plots', recursive = T)
pdf(file = 'figures/variable_importance/var_imp_plots/agg_var_imp_plot_rescaled.pdf', height = 6, width = 12)
agg_var_imp_plot
dev.off()


# perform t-tests for each variable to estimate if the mean values differ across all species ----

# perform t.tests and tidy up
all_t.tests <- var_imp %>% 
  select(dataset, covariate, species_name, occurrence_imp_rescaled, abundance_imp_rescaled) %>% 
  group_by(covariate, dataset) %>% 
  do(t_test_results = t.test(y = .$abundance_imp_rescaled, x=.$occurrence_imp_rescaled)) %>% 
  broom::tidy(t_test_results) %>% 
  rename(., occurrence_estimate = estimate1, 
            abundance_estimate  = estimate2) %>% 
  arrange(., dataset, estimate)

# change covariate names to match figure
all_t.tests$covariate <- recode(as.character(all_t.tests$covariate),
                                'Elevation_GEBCO' = 'depth/elevation', 
                                'Depth_GEBCO_transformed' = 'depth/elevation', 
                                'human_pop' = 'human', 
                                'human_pop_2015_50km' = 'human', 
                                'primary_forest' = 'forest', 
                                'robPCA_1' = 'climate PC1', 
                                'robPCA_2' = 'climate PC2', 
                                'robPCA_3' = 'climate PC3', 
                                'sst_mean' = 'sst', 
                                'wave_energy_mean' = 'wave energy', 
                                'reef_area_200km' = 'reef area')


# write t.tests to csv to produce table
writexl::write_xlsx(all_t.tests, path = paste0('figures/variable_importance/var_imp_plots', '/', 'var_importance_t.tests','.xlsx'))


# plots of spearmans rank correlations at species level between variables ----

# extract species level spearmans rank correlations
var_spear <- var_imp %>% 
  select(dataset, species_name, spearmans_rank) %>%
  unique() %>% 
  group_by(dataset) %>% 
  mutate(median = median(spearmans_rank), 
         upr    = quantile(spearmans_rank, 0.9), 
         lwr    = quantile(spearmans_rank, 0.1), 
         mean   = mean(spearmans_rank), 
         sd     = sd(spearmans_rank))

var_spear$dataset <- recode(as.character(var_spear$dataset),
                              'rls' = 'reef-life survey', 
                              'bbs' = 'breeding-bird survey')


var_spear_plot <- ggplot(var_spear) + 
  geom_histogram(aes(x = spearmans_rank, fill = dataset)) + 
  #geom_point(data = var_spear %>% select(dataset, median) %>% unique(), aes(x = median, y = c(25/2, 10/2)))+
  #geom_segment(data = var_spear %>% select(dataset, median, upr, lwr) %>% unique(), aes(x = lwr, xend = upr, y = c(25/2, 10/2), yend = c(25/2, 10/2)))+
  geom_point(data = var_spear %>% select(dataset, mean) %>% unique(), aes(x = mean, y = c(25/2, 10/2), size = 2))+
  geom_segment(data = var_spear %>% select(dataset, mean, sd) %>% unique(), aes(x = mean-sd, xend = mean+sd, y = c(25/2, 10/2), yend = c(25/2, 10/2)))+
  theme_bw() + 
  theme(legend.position = 'none',
        aspect.ratio = 1, 
        panel.grid = element_blank(),  
        axis.text = element_text(size = 12),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 14),
        strip.background = element_blank())  + 
  facet_wrap(~dataset, scales = 'free') +
  scale_fill_manual(values = c(viridis::viridis(5, option = 3)[4], viridis::viridis(5, option = 7)[3])) + 
  xlab('spearmans rank between variable importance scores \n of abundance and occurrence models') + 
  ylab(NULL)
  
pdf(file = 'figures/variable_importance/var_imp_plots/var_spearmans_rank_plot.pdf', height = 5, width = 10)
var_spear_plot
dev.off()



