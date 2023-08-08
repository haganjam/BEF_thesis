#'
#' @title Classic BEF experimental design figure
#' 
#' @description Hypothetical BEF relationships used to illustrate the difference
#' between two different BEF experimental designs.
#'

# load relevant libraries
library(dplyr)
library(ggplot2)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# simulate a BEF relationship
SR <- seq(1, 3, 0.05)
Prod <- 1 + (0.5*SR) - (0.085*(SR^2))

plot(Prod ~ SR)


# Experimental design 1

# add these values to a data.frame
set.seed(34265)
df1 <- 
  data.frame(id = c(1:7),
             SR = c(1, 1, 1, 2, 2, 2, 3),
             Prod = c(rnorm(n = 3, mean = Prod[SR == 1], sd = 0.1),
                      rnorm(n = 3, mean = Prod[SR == 2], sd = 0.1),
                      Prod[SR == 3]))

# convert to a tibble
df1 <- as_tibble(df1)

p1 <- 
  ggplot() +
  geom_point(data = df1, 
             mapping = aes(x = SR, y = Prod), colour = "red", size = 2) +
  geom_line(mapping = aes(x = SR, y = Prod), colour = "red") + 
  scale_x_continuous(limits = c(0.9, 3.1), breaks = c(1:3)) +
  scale_y_continuous(limits = c(1.2, 1.8)) +
  ylab("Ecosystem function (T = 1)") +
  xlab("Species richness (T = 0)") +
  theme_meta() +
  theme(axis.text.y = element_text(size = 2.5, colour = "white"),
        axis.ticks.y = element_blank())
plot(p1)

ggsave("figures-tables/fig_4a.svg", p1,
       unit = "cm", width = 6.5, height = 5.75)

# Experimental design 2

# add these values to a data.frame
set.seed(215)
df2 <- 
  data.frame(id = c(1:9),
             SR = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
             Prod = c(rnorm(n = 3, mean = Prod[SR == 1], sd = 0.1),
                      rnorm(n = 3, mean = Prod[SR == 2], sd = 0.1),
                      rnorm(n = 3, mean = Prod[SR == 3], sd = 0.1)))

# convert to a tibble
df2 <- as_tibble(df2)

p2 <- 
  ggplot() +
  geom_point(data = df2, 
             mapping = aes(x = SR, y = Prod), colour = "red", size = 2) +
  geom_line(mapping = aes(x = SR, y = Prod), colour = "red") + 
  scale_x_continuous(limits = c(0.9, 3.1), breaks = c(1:3)) +
  scale_y_continuous(limits = c(1.2, 1.9)) +
  ylab("Ecosystem function (T = 1)") +
  xlab("Species richness (T = 0)") +
  theme_meta() +
  theme(axis.text.y = element_text(size = 2.5, colour = "white"),
        axis.ticks.y = element_blank())
plot(p2)

ggsave("figures-tables/fig_4b.svg", p2,
       unit = "cm", width = 6.5, height = 5.75)

### END
