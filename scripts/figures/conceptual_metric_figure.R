# create conceptual metric figure 

library(ggplot2)
library(cowplot) # load cowplot for clean graph

n = 100
slope = 1
intercept = 0
x = rnbinom(n, 100, 0.9)
error = round(rnorm(100, 0, 1))

lm_line = (slope*x + intercept)

int_adj_1 = 0.5
slope_adj_1 = 0.2

y = ((slope+slope_adj_1)*x + (intercept-int_adj_1)) + error

error_2 = round(rnorm(100, 0, 3))
int_adj_2 <- +7.5
slope_adj_2 <- -0.5

y_2 = ((slope+slope_adj)*x + (intercept+int_adj)) + error_2

# estimate residuals
residal_lines <- (y_2 - lm_line) + lm_line

cols <- colorRampPalette(c("#0099CC80","#9ECAE1","#58BC5D","#EEF559","#FF9933","red"), bias = 1)(5)

conc_fig <- ggplot() + 
  geom_linerange(aes(x = x, ymin = ifelse(lm_line < residal_lines, residal_lines, lm_line), 
                     ymax = ifelse(lm_line > residal_lines, residal_lines, lm_line)), alpha = 0.5) +
  #geom_point(aes(x = x, y = y), size = 5, alpha = 0.5, pch = 19, stroke = 0, col =cols[4]) +
  geom_point(aes(x = x, y = y_2), size = 7, alpha = 1, pch = 19, stroke = 0, col = cols[1]) +
  geom_abline(aes(slope = slope, intercept = intercept), col = 'black') + 
  #geom_abline(aes(slope = slope+slope_adj_1, intercept = intercept+int_adj_1), col = cols[4]) +
  geom_abline(aes(slope = slope+slope_adj_2, intercept = intercept+int_adj_2), col = '#1F407A') + 
  ylim(c(0, max(y_2) + 2)) + 
  xlim(c(0, (max(y) + 2))) + 
  xlab(NULL) + 
  ylab(NULL) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.75)
  


dir.create('figures/conceptual_metrics')
pdf('figures/conceptual_metrics/scatter_metrics.pdf', height = 5, width = 5)
conc_fig
dev.off()
