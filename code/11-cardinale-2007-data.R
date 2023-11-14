#'
#' @title Cardinale et al.'s (2007, PNAS) data
#' 
#' @description Plotting Cardinale et al.'s (2007) data on the log-response
#' rations for the NBE and transgressive overyielding
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(readr)

# load relevant scripts
source("code/helper-plotting-theme.R")

# load the data
card <- readr::read_csv("data/cardinale-2007-data.csv")

# plot the log response ratio of the NBE
p1 <- 
  ggplot(data = card, 
       mapping = aes(x = LRnet)) +
  geom_histogram(bins = 20, fill = "white", colour = "black", alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  ylab("Count") +
  xlab("Net Biodiversity Effect (LRR)") +
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
ggsave(filename = "figures-tables/def_fig_1.pdf", p1,
       unit = "cm", width = 10, height = 8, bg = "transparent")

# plot the log response ratio of transgressive overyielding
p2 <- 
  ggplot(data = card, 
       mapping = aes(x = LRtrans)) +
  geom_histogram(bins = 20, fill = "white", colour = "black", alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  ylab("Count") +
  xlab("Transgressive overyielding (LRR)") +
  theme_meta() +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
plot(p2)

# export the figure for further modification
ggsave(filename = "figures-tables/def_fig_2.pdf", p2,
       unit = "cm", width = 10, height = 8, bg="transparent")

### END
