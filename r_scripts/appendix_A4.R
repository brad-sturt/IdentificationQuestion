
library(ggplot2)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(grid)
library(ggfittext)



data = read.csv("../data/best_case_nested_R2.csv")


# Figure EC.6 in Appendix EC.4

figure_rev_order = 
  data %>%
  mutate(past_to_best = (best_case_revenue - max_previous_assortment) / 
           max_previous_assortment) %>%
  arrange(-past_to_best) %>%
  filter(rev_ordered == "rev_ordered") %>%
  mutate(row = row_number()) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  #geom_col(aes(x=row,y=past_to_worst,fill='red')) +
  # geom_line(aes(x=row, y=-past_to_worst),linetype="dotted") +
  #geom_point(aes(x=row,y=past_to_estimate),size=1)+
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  #coord_cartesian(ylim=c(0, 2))+ 
  #scale_x_continuous(breaks = c(0,20,40,60,80,100,120,140)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  ) 

figure_other = 
  data %>%
  mutate(past_to_best = (best_case_revenue - max_previous_assortment) / 
           max_previous_assortment) %>%
  arrange(-past_to_best) %>%
  filter(rev_ordered == "other") %>%
  mutate(row = row_number()) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  #geom_col(aes(x=row,y=past_to_worst,fill='red')) +
  # geom_line(aes(x=row, y=-past_to_worst),linetype="dotted") +
  #geom_point(aes(x=row,y=past_to_estimate),size=1)+
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  #coord_cartesian(ylim=c(0, 2))+ 
  #scale_x_continuous(breaks = c(0,20,40,60,80,100,120,140)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  ) 

figure_reverse_rev_order = 
  data %>%
  mutate(past_to_best = (best_case_revenue - max_previous_assortment) / 
           max_previous_assortment) %>%
  arrange(-past_to_best) %>%
  filter(rev_ordered == "reverse_rev_ordered") %>%
  mutate(row = row_number()) %>%
  ggplot() +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  #geom_col(aes(x=row,y=past_to_worst,fill='red')) +
  # geom_line(aes(x=row, y=-past_to_worst),linetype="dotted") +
  #geom_point(aes(x=row,y=past_to_estimate),size=1)+
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  #coord_cartesian(ylim=c(0, 2))+ 
  #scale_x_continuous(breaks = c(0,20,40,60,80,100,120,140)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  ) 

png(file="../figures/best_case_nested_reverse_revenue_ordered.png",
    width = 10, height = 2, units = 'in', res = 1000)
grid.arrange(figure_reverse_rev_order)
dev.off()

png(file="../figures/best_case_nested_revenue_ordered.png",
    width = 10, height = 2, units = 'in', res = 1000)
grid.arrange(figure_rev_order)
dev.off()

png(file="../figures/best_case_variety.png",
    width = 10, height = 2, units = 'in', res = 1000)
grid.arrange(figure_other)
dev.off()


