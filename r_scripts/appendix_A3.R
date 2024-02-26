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


################################################
# Figure EC.5 in Appendix A.3
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

