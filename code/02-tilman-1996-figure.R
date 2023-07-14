
# thesis graphs: Tilman et al. (1996)

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggpubr)

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Thesis_document")

# get the meta_theme()
source("figures/Function_plotting_theme.R")

# read the data
fig1b <- read_csv(file = "figures/Tilman_et_al_1996_data_Fig_1/Fig_1b_data.csv",
                  col_names = FALSE)
print(fig1b)

# rename the variables
names(fig1b) <- c("species_richness", "plant_cover")
print(fig1b)

# round the species richness variable to the nearest integer
fig1b$species_richness <- round(fig1b$species_richness, 0)
print(fig1b)

# get the data into the correct format
fig1b$stat <- rep(c("CI_upp", "mean", "CI_low"), length(unique(fig1b$species_richness)))

fig1b <- 
  fig1b %>%
  pivot_wider(id_cols = "species_richness",
              names_from = "stat",
              values_from = "plant_cover")

# generate the curve values from the function fit to the data
x <- seq(1, 24, 0.1)
y <- 27 + ((36.4*x)/(5.48+x))
fig1b_curve <- data.frame(species_richness = x, mean = y)
plot(x, y)

# plot the data
p1 <- 
  ggplot() +
  geom_point(data = fig1b,
             mapping = aes(x = species_richness, y = mean)) +
  geom_line(data = fig1b_curve, 
            mapping = aes(x = species_richness, y = mean)) +
  geom_errorbar(data = fig1b,
                mapping = aes(x = species_richness, ymin = CI_low, ymax = CI_upp),
                width = 0.5) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 12, 24 )) +
  annotate(geom = "text", x = 22, y = 32, label = bquote(r^2~" = "~0.18 ), size = 3) +
  xlab("Species richness treatment") +
  ylab("Plant cover (%)") +
  theme_meta()
  

# read the data
fig1c <- read_csv(file = "figures/Tilman_et_al_1996_data_Fig_1/Fig_1c_data.csv",
                  col_names = FALSE)
print(fig1c)

# rename the variables
names(fig1c) <- c("species_richness", "N_rooting_zone")
print(fig1c)

# round the species richness variable to the nearest integer
fig1c$species_richness <- round(fig1c$species_richness, 0)
print(fig1c)

# get the data into the correct format
fig1c$stat <- rep(c("CI_upp", "mean", "CI_low"), length(unique(fig1c$species_richness)))

fig1c <- 
  fig1c %>%
  pivot_wider(id_cols = "species_richness",
              names_from = "stat",
              values_from = "N_rooting_zone")

# generate the curve values from the function fit to the data
x <- seq(1, 24, 0.1)
y <- 0.17 + (0.24*exp(-0.41*x))
fig1c_curve <- data.frame(species_richness = x, mean = y)
plot(x, y)

# plot the data
p2 <- 
  ggplot() +
  geom_point(data = fig1c,
             mapping = aes(x = species_richness, y = mean)) +
  geom_line(data = fig1c_curve, 
            mapping = aes(x = species_richness, y = mean)) +
  geom_errorbar(data = fig1c,
                mapping = aes(x = species_richness, ymin = CI_low, ymax = CI_upp),
                width = 0.5) +
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


ggsave(filename = "figures/Tilman_et_al_1996_data_Fig_1/Fig_1bc1.png", plot = p12,
       unit = "cm", width = 20, height = 7.5, dpi = 400)

### END
