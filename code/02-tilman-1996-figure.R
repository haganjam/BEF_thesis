#'
#' @title Tilman (1996) figure
#' 
#' @description Loads data extracted from Tilman et al. (1996) and replot it
#' for inclusion in the thesis
#'

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggpubr)

# get the custom plotting theme
source("code/helper-plotting-theme.R")

# read the data
til_dat1 <- readr::read_csv(file = "data/tilman-1996-data-1.csv", col_names = FALSE)
print(til_dat1)

# rename the variables
names(til_dat1) <- c("species_richness", "plant_cover")
print(til_dat1)

# round the species richness variable to the nearest integer
til_dat1$species_richness <- round(til_dat1$species_richness, 0)
print(til_dat1)

# convert the 

# get the data into the correct format
til_dat1$stat <- rep(c("CI_upp", "mean", "CI_low"), length(unique(til_dat1$species_richness)))

til_dat1 <- 
  til_dat1 |>
  tidyr::pivot_wider(id_cols = "species_richness",
                     names_from = "stat",
                     values_from = "plant_cover")

# generate the curve values from the function fit to the data
# parameters extracted directly from the original paper
x <- seq(1, 24, 0.1)
y <- 27 + ((36.4*x)/(5.48+x))
til_dat1_curve <- data.frame(species_richness = x, mean = y)
plot(x, y)

# plot the data
p1 <- 
  ggplot() +
  geom_point(data = til_dat1,
             mapping = aes(x = species_richness, y = mean), size = 2, colour = "red") +
  geom_line(data = til_dat1_curve, 
            mapping = aes(x = species_richness, y = mean), colour = "red") +
  geom_errorbar(data = til_dat1,
                mapping = aes(x = species_richness, ymin = CI_low, ymax = CI_upp),
                width = 0, colour = "red") +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 12, 24 )) +
  annotate(geom = "text", x = 22, y = 32, label = bquote(r^2~" = "~0.18 ), size = 3) +
  xlab("Species richness treatment") +
  ylab("Plant cover (%)") +
  theme_meta()
plot(p1)  

# read the data
til_dat2 <- readr::read_csv(file = "data/tilman-1996-data-2.csv",
                            col_names = FALSE)
print(til_dat2)

# rename the variables
names(til_dat2) <- c("species_richness", "N_rooting_zone")
print(til_dat2)

# round the species richness variable to the nearest integer
til_dat2$species_richness <- round(til_dat2$species_richness, 0)
print(til_dat2)

# get the data into the correct format
til_dat2$stat <- rep(c("CI_upp", "mean", "CI_low"), length(unique(til_dat2$species_richness)))

til_dat2 <- 
  til_dat2 |>
  tidyr::pivot_wider(id_cols = "species_richness",
                     names_from = "stat",
                     values_from = "N_rooting_zone")

# generate the curve values from the function fit to the data
x <- seq(1, 24, 0.1)
y <- 0.17 + (0.24*exp(-0.41*x))
til_dat2_curve <- data.frame(species_richness = x, mean = y)
plot(x, y)

# plot the data
p2 <- 
  ggplot() +
  geom_point(data = til_dat2,
             mapping = aes(x = species_richness, y = mean), colour = "red",
             size = 2) +
  geom_line(data = til_dat2_curve, 
            mapping = aes(x = species_richness, y = mean), colour = "red") +
  geom_errorbar(data = til_dat2,
                mapping = aes(x = species_richness, ymin = CI_low, ymax = CI_upp),
                width = 0, colour = "red") +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 12, 24 )) +
  annotate(geom = "text", x = 22, y = 0.37, label = bquote(r^2~" = "~0.22 ), size = 3) +
  xlab("Species richness treatment") +
  ylab(expression(NO[3]~"mg"~kg^{-1})) +
  theme_meta()

# arrange the plots
p12 <- 
  ggarrange(p1, p2, labels = c("a", "b"),
          font.label = list(size = 11, face = "plain"),
          widths = c(1, 1.1))
plot(p12)

ggsave(filename = "figures-tables/fig_3.svg", plot = p12,
       unit = "cm", width = 20, height = 7.5)

### END
