#'
#' @title Additional figures for the defence presentation
#' 
#' @description Plot miscellaneous figures for the thesis defence
#'
# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# load plotting theme
source("code/helper-plotting-theme.R")

# simulated correlation between biodiversity and some function

# set the number of samples
n <- 300

# obtain an environmental variable
df <- data.frame(env1 = runif(n = n, min = 0.05, max = 10))

# simulate biodiversity as a function of the environmental variable
df$SR <- rpois(n = 300, df$env1*1.5)

# simulate function as a function of both SR and env1
df$prod <- df$env1*1.2 + df$SR*1.25 + rnorm(n = 300, 0, 7)

# remove the zeros
df <- df[df$SR>0,]

# plot the results
p1 <- 
  ggplot(data = df,
       mapping = aes(x = SR, y = prod)) +
  ggbeeswarm::geom_quasirandom(shape = 1, size = 2) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.5) +
  ylab("Ecosystem process rate") +
  xlab("Biodiversity variable (e.g. species richness)") +
  theme_meta() +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

# export the figure for further modification
ggsave(filename = "figures-tables/def_fig_5.pdf", p1,
       unit = "cm", width = 10, height = 8, bg = "transparent")

