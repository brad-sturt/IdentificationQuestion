

library(tidyverse)
library(latex2exp)
library(gridExtra)
library(grid)
library(ggfittext)
library(ggh4x)

##################################################
#  Figure EC.8 in Appendix A.6
##################################################


data = read.csv("../data/appendix_A6.csv")


data = data %>%
  mutate(past_to_worst = (worst_case_pareto - max_previous_assortment) / max_previous_assortment,
         past_to_best = (best_case_pareto - max_previous_assortment) / max_previous_assortment,
         past_to_optimal = (optimal_obj_val - max_previous_assortment) / max_previous_assortment,
         past_to_expected = (expected_revenue_pareto - max_previous_assortment) / max_previous_assortment)
plot_conjoint = 
  data %>%
  select(iter,n, M,past_to_best, past_to_worst,past_to_optimal,past_to_expected) %>%
  distinct(iter, n, M,past_to_best, past_to_worst,past_to_optimal,past_to_expected) %>%
  filter(past_to_best > 0.001) %>%
  group_by(iter, n,M,past_to_optimal) %>%
  mutate(row = rank(past_to_worst)) %>% 
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  geom_col(aes(x=row,y=past_to_worst),fill='red') +
  geom_hline(aes(yintercept=past_to_optimal)) +
  geom_point(aes(x=row,y=past_to_expected)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),
                     breaks = seq(-6,12,1)/5) +
  theme( legend.position="none",
         axis.title.y=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank()) + 
  ggh4x::facet_grid2(M ~ iter, scales = "free_x",independent="x") 


png(file="../figures/conjoint_plot.png",
    width = 12, height = 6, units = 'in', res = 400)
grid.arrange(plot_conjoint)
dev.off()


