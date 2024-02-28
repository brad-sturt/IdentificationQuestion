
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

################################################
# Figure EC.3 in Appendix A.2
################################################

data = read.csv("../data/appendix_A2_1.csv")
figure = ggplot(data, aes(x=max_previous_assortments,y=worst_case_revenue)) + 
  geom_point() +
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(legend.position="none") +
  ylab(TeX("Worst-Case Expected Revenue of New Assortment")) +
  xlab("Expected Revenue of Best Previously-Offered Assortment")
png(file="../figures/two_assortments.png",
    width = 8, height = 6, units = 'in', res = 300)
grid.arrange(figure)
dev.off()

################################################
# Figure EC.4 in Appendix A.2
################################################

data_speed = read.csv("../data/appendix_A2_2.csv")
plot_speed = ggplot(data_speed %>% group_by(n) %>% summarize(running_time = mean(speed), running_time_sd = sd(speed))) +
  geom_line(aes(x=n,y=running_time)) + 
  geom_point(aes(x=n,y=running_time)) + 
  geom_ribbon(aes(x=n,ymin = running_time - running_time_sd, ymax=running_time + running_time_sd),alpha=0.2) +
  theme_bw(base_size = 16, base_family = "Helvetica")  +
  ylab("Computation Time (Seconds)") +
  xlab("Number of Products (n)")

png(file="../figures/speed_two_assortments.png",
    width = 8, height = 6, units = 'in', res = 200)
grid.arrange(plot_speed)
dev.off()







