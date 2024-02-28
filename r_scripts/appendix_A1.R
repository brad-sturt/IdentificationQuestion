
library(ggplot2)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(grid)
library(ggfittext)

##################################################
# Be sure that the current working directory is
# set to `IdentificationQuestion/r_script`
##################################################

data = read.csv("../data/appendix_A1.csv")

# As a sanity check, confirm for each problem instance that the worst-case 
# expected revenue of the assortment found by estimate-then-optimize is never
# greater than the expected revenue under the best past assortment
data %>% 
  count(worst_case_revenue > max_previous_assortments + 1e-10)

# Create figure
figure = ggplot(data, aes(x=max_previous_assortments,y=worst_case_revenue)) + 
  geom_point() +
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(legend.position="none") +
  ylab("Worst-Case Expected Revenue of \nEstimate-Then-Optimize Assortment") +
  xlab("Expected Revenue of Best Previously-Offered Assortment")

png(file="../figures/revenue_ordered_assortment_prev.png",
    width = 8, height = 6, units = 'in', res = 300)
grid.arrange(figure)
dev.off()

# Obtain the subset of problem instances for which the worst-case expected
# revenue of the assortment found by estimate-then-optimize is strictly
# less than the expected revenue under the best past assortment
data_worse = data %>%
  filter(worst_case_revenue < max_previous_assortments - 1e-10) %>%
  # Add the best-case and worst-case relative percentage improvement of the 
  # expected revenue of the new assortment obtained using estimate-then-optimize 
  # over the expected revenue of the firm's best previously-offered assortment 
  mutate(past_to_worst = (worst_case_revenue - max_previous_assortments) / 
                            max_previous_assortments,
         past_to_best = (best_case_revenue - max_previous_assortments) / 
                            max_previous_assortments)

# Compute the number of problem instances in data_worse
nrow(data_worse)

# Compute the percentage of problem instances in data_worse for which the 
# best-case increase in expected revenue exceeds the worst-case increase in
# expected revenue (in absolute value)
(data_worse %>% 
  filter(abs(past_to_worst) > past_to_best) %>%
  nrow()) / nrow(data_worse) * 100

# Calculate summary statistics for data_worse
data_worse %>% summarize(mean(past_to_best),mean(past_to_worst))

# Create figure
plot_intervals = ggplot(data_worse %>% arrange(-past_to_worst) %>%
                          mutate(row = seq.int(nrow(data_worse)))) +
  geom_col(aes(x=row,y=past_to_best),fill='blue') + 
  geom_col(aes(x=row,y=past_to_worst,fill='red')) +
  geom_line(aes(x=row, y=-past_to_worst),linetype="dotted") +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"),
                     breaks = seq(-6,6,1)/10) +
  scale_x_continuous(breaks = c(0,20,40,60,80,100,120,140)) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  ) 

png(file="../figures/intervals.png",
    width = 6, height = 4, units = 'in', res = 1000)
grid.arrange(plot_intervals)
dev.off()

