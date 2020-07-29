# script to analyse and plot patterns in variable importance 

# libraries and set up---- 

lib_vect <- c('tidyverse', 'cowplot', 'RColorBrewer', 'gridExtra', 'grid')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)
detach("package:raster", unload = TRUE)

# create vector of species that have good performing species occurrence models and species abundance models ----

# read in species occurrence results and obtain TSS
sp_projections <- list.files('results/spatial_projections', recursive = T, full.names = T)

# loop through read and create dataset of TSS
occ_TSS_values <- parallel::mclapply(1:length(sp_projections), mc.cores=10, function(x){
  #print(x)
  species_file <- sp_projections[x]
  TSS <- as.numeric(as.character(readRDS(species_file)$TSS_threshold$max.TSS))
  dataset <- str_split(species_file, '/')[[1]][3]
  species_name = gsub('.RDS', '', gsub('_', ' ', str_split(species_file, '/')[[1]][4]))
  return(data.frame(species_name, dataset, species_file, TSS))
  })
occ_TSS_values <- do.call(rbind, occ_TSS_values)
occ_TSS_values <- occ_TSS_values %>% filter(TSS > 0.75)
table(occ_TSS_values$dataset)
#  bbs  rls 
# 1025  489 

# read in species abundance results and obtain spearmans rank correlation for specific model type
all_assessments <- lapply(list.files('results/model_assessment_all/validation', full.names = T), readRDS)
all_assessments <- do.call(rbind, all_assessments)
abun_perform <- all_assessments %>% filter(abundance_response == 'abun', plot_level == 'rf.abun.g.l', cross_validation %in% c('bbs_basic', 'rls_basic')) %>% 
  select(species_name, Dspearman, dataset) %>% 
  filter(Dspearman > 0.5)

table(abun_perform$dataset)
# bbs  rls 
# 1009  460 

# see how many species per dataset are removed in total
full_sp_list <- full_join(abun_perform, occ_TSS_values) %>% na.omit
full_sp_list %>% filter(dataset == 'bbs') %>% .$species_name %>% unique %>% length
full_sp_list %>% filter(dataset == 'rls') %>% .$species_name %>% unique %>% length

# join together evaluations from species occurrence models and abundance models
high_performance_species <- unique(full_sp_list$species_name)

saveRDS(high_performance_species, 'results/high_performance_species.RDS')

# read in variable importance results ----

# read in high performance species
high_performance_species <- readRDS('results/high_performance_species.RDS')

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
  filter(species_name %in% high_performance_species) %>% 
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

# aggregate plot
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
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
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

# estimate overall spearmans rank values and save objects to file ----

var_imp$dataset <- recode(as.character(var_imp$dataset),
       'rls' = 'reef-life survey', 
       'bbs' = 'breeding-bird survey')

agg_var_imp %>% 
  group_by(dataset) %>% 
  do(spearmans_rank = cor.test(x = .$occ_imp_mean, y = .$abun_imp_mean, method = 'spearman')) %>% 
  broom::tidy(spearmans_rank) %>% 
  left_join(., var_imp %>% 
              group_by(dataset) %>% 
              do(spearmans_rank_all = cor.test(x = .$occurrence_imp, y = .$abundance_imp, method = 'spearman')) %>% 
              broom::tidy(spearmans_rank_all) %>% 
              unnest() %>% 
              rename(., 
                     'estimate_all' = 'estimate', 
                     'statistic_all' = 'statistic', 
                     'p.value_all'   = 'p.value') %>% select(estimate_all, statistic_all, p.value_all)) %>% 
  left_join(., 
            agg_var_imp %>% 
              group_by(dataset) %>% 
              do(r2 = summary(lm(.$occ_imp_mean ~ .$abun_imp_mean, data = .))$r.squared) %>% 
              unnest()) %>% 
  left_join(., var_imp %>% 
              group_by(dataset) %>% 
              do(r2_all = summary(lm(.$abundance_imp_rescaled ~ .$occurrence_imp_rescaled * .$covariate, data = .))$r.squared) %>% 
              unnest()) %>% 
  writexl::write_xlsx(., path = paste0('figures/variable_importance/var_imp_plots', '/', 'agg_var_importance_spearmans','.xlsx'))

# perform t-tests for each variable to estimate if the mean values differ across all species ----

# perform t.tests and tidy up
all_t.tests <- var_imp %>% 
  filter(species_name %in% high_performance_species) %>% 
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
  filter(species_name %in% high_performance_species) %>% 
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
        aspect.ratio = 0.5, 
        panel.grid = element_blank(),  
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
        strip.background = element_blank())  + 
  facet_wrap(~dataset, scales = 'free') +
  scale_fill_manual(values = c(viridis::viridis(10, option = 3)[5], viridis::viridis(10, option = 7)[5])) + 
  xlab('spearmans rank between variable importance scores \n of abundance and occurrence models') + 
  ylab(NULL)
  
pdf(file = 'figures/variable_importance/var_imp_plots/var_spearmans_rank_plot.pdf', height = 5, width = 10)
var_spear_plot
dev.off()

# write variable spearmans rank object to file
saveRDS(var_spear, file = 'results/varimp_species_spearmans_rank.RDS')


# correlation between correlations in variable importance scores and correlations in spatial patterns of occurrence and abundance ----

var_spear     <- readRDS('results/varimp_species_spearmans_rank.RDS')
spatial_spear <- readRDS('results/spatial_species_spearmans_rank.RDS')

# harmoize data so compatable
var_spear <- var_spear %>% 
  rename('var_spear' = 'spearmans_rank') %>%
  mutate(species_name = gsub(' ', '_', species_name)) %>%
  ungroup() %>% 
  select(dataset, species_name, var_spear)

spatial_spear <- spatial_spear %>% 
  rename('spatial_spear' = 'cor') %>% 
  ungroup() %>% 
  select(dataset, species_name, spatial_spear)

# combine
all_spear <- left_join(var_spear, spatial_spear)

# perform tests
spear_overall <- all_spear %>% 
  group_by(dataset) %>% 
  do(correlation = cor.test(x = .$var_spear, y = .$spatial_spear)) %>% 
  broom::tidy(correlation) 

# save model
spear_overall %>% writexl::write_xlsx(., path = paste0('figures/variable_importance/var_imp_plots', '/', 'spatial_vs_variable_correlation','.xlsx'))

# Pearson's product-moment correlation
# data:  all_spear$var_spear and all_spear$spatial_spear
# t = 7.8248, df = 1463, p-value = 9.679e-15
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1507559 0.2490811
# sample estimates:
#       cor 
# 0.2004231 

# aggregate for points
all_spear_agg <- all_spear %>% 
  group_by(dataset) %>% 
  do(mean_var_spear = mean(.$var_spear), 
     sd_var_spear = sd(.$var_spear), 
     mean_spatial_spear = mean(.$spatial_spear), 
     sd_spatial_spear = sd(.$spatial_spear)) %>% 
  unnest()

lib_vect <- c('tidyverse', 'gridExtra', 'grid')
install.lib<-lib_vect[!lib_vect %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(lib_vect,require,character=TRUE)

# manual theme function for histograms
theme_hist <- function(){
  theme(legend.position = 'none',
      aspect.ratio = 0.5, 
      panel.grid = element_blank(),  
      axis.text = element_text(size = 16),
      strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
      axis.title = element_text(size = 20),
      strip.background = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line())}

# create empty plot to fill gap
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"))

  
# BBS CORRELATION COMPARISON PLOTS ----

# bbs plot of both
bbs_plot_cors <- ggplot(all_spear %>% filter(dataset == 'breeding-bird survey'), aes(x = var_spear, y = spatial_spear)) +
  geom_density_2d() + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.75, col = 'white', n =200, size = 0.1) + 
  geom_point(alpha = 0.5, size = 2, stroke = 0) +
  stat_smooth(method = 'lm', col = 'black', se = F) + 
  scale_fill_viridis_c(option = 3, begin = 0.1, end = 0.9) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        aspect.ratio = 1, 
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        legend.position = 'none', 
        plot.title = element_text(vjust = -12, hjust = 0.1, size = 15), 
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  ylim(-1, 1.1) + 
  xlim(-1, 1.1) + 
  xlab('variable importance correlation') + 
  ylab('spatial correlation') + 
  ggtitle(label = bquote('rho' == .(signif(spear_overall[spear_overall$dataset == 'breeding-bird survey',]$estimate,2))~','~'\n'~'p <' ~ .(0.001))) 
  
bbs_var_hist <- ggplot(all_spear %>% filter(dataset == 'breeding-bird survey')) + 
  geom_histogram(aes(x = var_spear, fill = dataset), alpha = 1) + 
  geom_point(data = all_spear_agg %>% filter(dataset == 'breeding-bird survey'), aes(x = mean_var_spear, y=30), size = 4) +
  geom_segment(data = all_spear_agg %>% filter(dataset == 'breeding-bird survey'), 
               aes(x = mean_var_spear-sd_var_spear, xend = mean_var_spear+sd_var_spear, y = 30, yend = 30), size = 1.1) +
  theme_bw() + 
  theme_hist() +
  theme(aspect.ratio = 0.2, 
        axis.text.x = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
  scale_fill_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
  xlab(NULL) + 
  ylab(NULL) + 
  xlim(-1, 1.1)  + 
  scale_y_continuous(breaks = c(0, 50, 100))

  

bbs_spatial_hist <- ggplot(all_spear %>% filter(dataset == 'breeding-bird survey')) + 
  geom_histogram(aes(x = spatial_spear, fill = dataset), alpha = 1) + 
  geom_point(data = all_spear_agg %>% filter(dataset == 'breeding-bird survey'), aes(x = mean_spatial_spear, y=30), size = 4) +
  geom_segment(data = all_spear_agg %>% filter(dataset == 'breeding-bird survey'), 
               aes(x = mean_spatial_spear-sd_spatial_spear, xend = mean_spatial_spear+sd_spatial_spear, y = 30, yend = 30), size = 1.1) +
  theme_bw() + 
  theme_hist() +
  theme(aspect.ratio = 5, 
        axis.text.x = element_text(angle = 270, vjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_colour_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
  scale_fill_manual(values = c(viridis::viridis(10, option = 3, begin = 0.1, end = 0.9)[10])) + 
  xlab(NULL) + 
  ylab(NULL) + 
  coord_flip() + 
  xlim(-1, 1.1) + 
  scale_y_continuous(breaks = c(0, 50, 100))

library(patchwork)
png('figures/correlation_comparisons/bbs_combination.png', width = 1750, height = 1750, res = 300)
bbs_var_hist + plot_spacer() + bbs_plot_cors + bbs_spatial_hist + 
  plot_layout(ncol = 2, nrow = 2, widths = c(1, 0.2), heights = c(0.2, 1))
dev.off()



# RLS CORRELATION COMPARISON PLOTS ----

# rls plot of both
rls_plot_cors <- ggplot(all_spear %>% filter(dataset == 'reef-life survey'), aes(x = var_spear, y = spatial_spear)) +
  geom_density_2d() + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.75, col = 'white', n =200, size = 0.1) + 
  geom_point(alpha = 0.5, size = 2, stroke = 0) +
  stat_smooth(method = 'lm', col = 'black', se = F) + 
  scale_fill_viridis_c(option = 7, begin = 0.1, end = 0.9) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        aspect.ratio = 1, 
        axis.text = element_text(size = 16),
        strip.text.x = element_text(angle = 0, hjust = 0, size = 14),
        axis.title = element_text(size = 20),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        legend.position = 'none', 
        plot.title = element_text(vjust = -12, hjust = 0.1, size = 15), 
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  ylim(-1, 1.1) + 
  xlim(-1, 1.1) + 
  xlab('variable importance correlation') + 
  ylab('spatial correlation') + 
  ggtitle(label = bquote('rho' == .(signif(spear_overall[spear_overall$dataset == 'reef-life survey',]$estimate,2))~','~'\n'~'p <' ~ .(0.001))) 

rls_var_hist <- ggplot(all_spear %>% filter(dataset == 'reef-life survey')) + 
  geom_histogram(aes(x = var_spear, fill = dataset), alpha = 1) + 
  geom_point(data = all_spear_agg %>% filter(dataset == 'reef-life survey'), aes(x = mean_var_spear, y=10), size = 4) +
  geom_segment(data = all_spear_agg %>% filter(dataset == 'reef-life survey'), 
               aes(x = mean_var_spear-sd_var_spear, xend = mean_var_spear+sd_var_spear, y = 10, yend = 10), size = 1.1) +
  theme_bw() + 
  theme_hist() +
  theme(aspect.ratio = 0.2, 
        axis.text.x = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c(viridis::viridis(10, option = 7, begin = 0.1, end = 0.9)[5])) + 
  scale_fill_manual(values = c(viridis::viridis(10, option = 7, begin = 0.1, end = 0.9)[5])) + 
  xlab(NULL) + 
  ylab(NULL) + 
  xlim(-1, 1.1) + 
  scale_y_continuous(breaks = c(0,20,40, 60))

rls_spatial_hist <- ggplot(all_spear %>% filter(dataset == 'reef-life survey')) + 
  geom_histogram(aes(x = spatial_spear, fill = dataset), alpha = 1) + 
  geom_point(data = all_spear_agg %>% filter(dataset == 'reef-life survey'), aes(x = mean_spatial_spear, y=10), size = 4) +
  geom_segment(data = all_spear_agg %>% filter(dataset == 'reef-life survey'), 
               aes(x = mean_spatial_spear-sd_spatial_spear, xend = mean_spatial_spear+sd_spatial_spear, y = 10, yend = 10), size = 1.1) +
  theme_bw() + 
  theme_hist() +
  theme(aspect.ratio = 5, 
        axis.text.x = element_text(angle = 270, vjust = 0.5), 
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_colour_manual(values = c(viridis::viridis(10, option = 7, begin = 0.1, end = 0.9)[5])) + 
  scale_fill_manual(values = c(viridis::viridis(10, option = 7, begin = 0.1, end = 0.9)[5])) + 
  xlab(NULL) + 
  ylab(NULL) + 
  coord_flip() + 
  xlim(-1, 1.1) +
  scale_y_continuous(breaks = c(0, 20, 40, 60))

library(patchwork)
png('figures/correlation_comparisons/rls_combination.png', width = 1750, height = 1750, res = 300)
rls_var_hist + plot_spacer() + rls_plot_cors + rls_spatial_hist + 
  plot_layout(ncol = 2, nrow = 2, widths = c(1, 0.2), heights = c(0.2, 1))
dev.off()



rls_plot_cors <- ggplot(all_spear %>% filter(dataset == 'reef-life survey'), aes(x = var_spear, y = spatial_spear)) +
  geom_density_2d() + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.75, col = 'white', n =200, size = 0.1) + 
  geom_point(alpha = 0.5, size = 2, stroke = 0) +
  stat_smooth(method = 'lm', col = 'black', se = F) + 
  scale_fill_viridis_c(option = 3, begin = 0.1, end = 0.9) + 
  facet_wrap(~dataset)  +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        aspect.ratio = 1, 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = 'black', fill = 'transparent'), 
        legend.position = 'none', 
        plot.title = element_text(vjust = -15, hjust = 0.03, size = 10)) + 
  ylim(-1, 1.1) + 
  xlim(-1, 1.1) + 
  xlab('variable importance correlation') + 
  ylab('spatial correlation') + 
  ggtitle(label = bquote('rho' == .(round(spear_overall[spear_overall$dataset == 'reef-life survey',]$estimate,2))~','~'\n'~'p <' ~ .(0.001))) 

grid.arrange(bbs_plot_cors, rls_plot_cors, nrow = 1)

png(file = )






