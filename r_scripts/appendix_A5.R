
library(ggplot2)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(grid)
library(ggfittext)
library(scales)

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

data = read.csv("../data/nested_R1_reverse_revenue_ordered_3.csv") 
data = data %>%
  mutate(improvement_over_best_past_assortment_best_case = best_case_revenue / max_previous_assortment,
         improvement_over_best_past_assortment_worst_case = worst_case_revenue / max_previous_assortment,
         past_to_worst = (worst_case_revenue - max_previous_assortment) / max_previous_assortment,
         past_to_best = (best_case_revenue - max_previous_assortment) / max_previous_assortment,
         )

figure_reverse_rev_order = data %>%
  select(iter, rev_ordered,n, past_to_best, past_to_worst) %>%
  mutate(past_to_best=round(past_to_best,3),
         past_to_worst=round(past_to_worst,3)) %>%
  distinct(iter, rev_ordered, n, past_to_best, past_to_worst) %>%
  group_by(iter, rev_ordered, n) %>%
  mutate(row = rank(past_to_worst)) %>%
  filter(rev_ordered == "false") %>%
  filter(past_to_best > 0.001) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  geom_col(aes(x=row,y=past_to_worst,fill='red')) +
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
  facet_grid(. ~ iter, scales = "free_x") 

png(file="../figures/nested_reverse_revenue_ordered.png",
    width = 10, height = 2, units = 'in', res = 1000)
grid.arrange(figure_reverse_rev_order)
dev.off()

figure_rev_order = data %>%
  select(iter, rev_ordered,n, past_to_best, past_to_worst) %>%
  mutate(past_to_best=round(past_to_best,3),
         past_to_worst=round(past_to_worst,3)) %>%
  distinct(iter, rev_ordered, n, past_to_best, past_to_worst) %>%
  group_by(iter, rev_ordered, n) %>%
  mutate(row = rank(past_to_worst)) %>%
  filter(rev_ordered == "true") %>%
  filter(past_to_best > 0.001) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  geom_col(aes(x=row,y=past_to_worst,fill='red')) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),
                     breaks = seq(-6,12,1)/10) +
  theme( legend.position="none",
         axis.title.y=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank()) + 
  facet_grid(. ~ iter, scales = "free_x") 

png(file="../figures/nested_revenue_ordered.png",
    width = 10, height = 2, units = 'in', res = 1000)
grid.arrange(figure_rev_order)
dev.off()

################################################
# Figure XXX in Section XXX
################################################

data = read.csv("../data/nested_R1_variety.csv") 
data = data %>%
  mutate(improvement_over_best_past_assortment_best_case = best_case_revenue / max_previous_assortment,
         improvement_over_best_past_assortment_worst_case = worst_case_revenue / max_previous_assortment,
         past_to_worst = (worst_case_revenue - max_previous_assortment) / max_previous_assortment,
         past_to_best = (best_case_revenue - max_previous_assortment) / max_previous_assortment,
  )
figure_variety = data %>%
  select(iter,n, past_to_best, past_to_worst) %>%
  mutate(past_to_best=round(past_to_best,3),
         past_to_worst=round(past_to_worst,3)) %>%
  distinct(iter, n, past_to_best, past_to_worst) %>%
  filter(past_to_best > 0.001) %>%
  group_by(iter, n) %>%
  mutate(row = rank(past_to_worst)) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  geom_col(aes(x=row,y=past_to_worst),fill='red') +
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
  facet_grid(. ~ iter, scales = "free_x") 

png(file="../figures/nested_variety.png",
    width = 10, height = 2, units = 'in', res = 800)
grid.arrange(figure_variety)
dev.off()











