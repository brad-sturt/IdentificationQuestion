
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
##################################################
# Be sure that the current working directory is
# set to `IdentificationQuestion/r_script`
##################################################

################################################
# Figure XXX in Section XXX
################################################

# Load in the data for Figure 5 in Section 4.25

data = read.csv("../data/nested_R1_reverse_revenue_ordered_3.csv") 
data = data %>%
  mutate(improvement_over_best_past_assortment_best_case = best_case_revenue / max_previous_assortment,
         improvement_over_best_past_assortment_worst_case = worst_case_revenue / max_previous_assortment,
         past_to_worst = (worst_case_revenue - max_previous_assortment) / max_previous_assortment,
         past_to_best = (best_case_revenue - max_previous_assortment) / max_previous_assortment,
         ) #%>%
  #select(rev_ordered, iter, fraction, improvement_over_best_past_assortment_best_case, improvement_over_best_past_assortment_worst_case)

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

################################################
# Figure XXX in Section XXX
################################################

data_speed = read.csv("../data/nested_R1_speed.csv")
plot_speed = 
  ggplot(data_speed %>% group_by(n,M) %>% summarize(running_time = mean(speed), running_time_sd = sd(speed))) +
  geom_line(aes(x=M,y=running_time)) + 
  geom_point(aes(x=M,y=running_time)) + 
  geom_ribbon(aes(x=M,ymin = running_time - running_time_sd, ymax=running_time + running_time_sd),alpha=0.2) +
  theme_bw(base_size = 16, base_family = "Helvetica")  +
  ylab("Computation Time (Seconds)") +
  xlab("Number of Past Assortments (M)") #+ 
  #xlab("Number of Products (n)") + 
  #facet_grid(.~n)

png(file="../figures/speed_nested.png",
    width = 7, height = 5, units = 'in', res = 200)
grid.arrange(plot_speed)
dev.off()

data_speed = read.csv("../data/nested_R1_speed_2.csv")


ggplot(data_speed, aes(x=max_previous_assortments,y=worst_case_revenue)) + 
  geom_point() +
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(legend.position="none") +
  ylab(TeX("Worst-Case Expected Revenue of New Assortment")) +
  xlab("Expected Revenue of Best Previously-Offered Assortment") +
  geom_abline() +
  facet_grid(n ~ M)







