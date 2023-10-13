#'
#' @title Carlander (1955) figure
#' 
#' @description Loads data extracted from Carlander (1955) in Hector and Wilby
#' (2009) and plots the figure.
#'

# load relevant libraries
library(ggplot2)
library(readr)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# load the Carlander (1955) data
car_dat <- readr::read_csv(file = "data/carlander-1955-data.csv", col_names = FALSE)
print(car_dat)

# check the data
names(car_dat) <- c("species_richness_ln", "biomass_pound_acre_ln")

# round off the species richness variable to a single value
car_dat$species_richness_ln <- round(car_dat$species_richness_ln, 0)

# run a regression on these data
lm_x <- lm(log(car_dat$biomass_pound_acre_ln) ~ log(car_dat$species_richness_ln))
summary(lm_x)

# get the colours
col_func()

# plot the figure
p1 <- 
  ggplot(data = car_dat,
       mapping = aes(x = species_richness_ln, y = biomass_pound_acre_ln)) +
  geom_jitter(width = 0.01, colour = "#118176", size = 1.75) +
  geom_smooth(method = "lm", size = 0.5, 
              alpha = 0.1, colour = "#118176", fill = "#118176") +
  scale_y_continuous(breaks = c(100, 200, 500, 1000), trans = "log") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20), trans = "log") +
  ylab(expression("Biomass pound"~acre^{-1}~"(ln)")) +
  xlab("Species richness (ln)") +
  annotate("text", x = 18, y = 80,
           label = lm_eqn(lm_x), parse = TRUE) +
  theme_meta()
plot(p1)

ggsave(filename = "figures-tables/fig_1.pdf", p1,
       unit = "cm", width = 10, height = 7.5)

### END
