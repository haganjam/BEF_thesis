
# thesis graphs: Carlander (1955)

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Thesis_document")

# get the meta_theme()
source("figures/Function_plotting_theme.R")

# load the Carlander (1955) data
df <- read_csv(file = "figures/Carlander_1955_data/Fig_X_data.csv",
               col_names = FALSE)
print(df)

# check the data
names(df) <- c("species_richness_ln", "biomass_pound_acre_ln")

# round off the species richness variable to a single value
df$species_richness_ln <- round(df$species_richness_ln, 0)

# run a regression on these data
lm.x <- lm(log(df$biomass_pound_acre_ln) ~ log(df$species_richness_ln))
summary(lm.x)

# extract the relevant data
lm.x <- summary(lm.x)
r2 <- round(lm.x$r.squared, 2)
Fst <- round(lm.x$fstatistic[1], 2)
df1 <- round(lm.x$fstatistic[2], 0)
df2 <- round(lm.x$fstatistic[3], 0)
pval <- round(lm.x$coefficients[2, 4], 2)

# plot the figure
ggplot(data = df,
       mapping = aes(x = species_richness_ln, y = biomass_pound_acre_ln)) +
  geom_jitter(width = 0.01) +
  geom_smooth(method = "lm", colour = "black", size = 0.5, 
              alpha = 0.25) +
  scale_y_continuous(breaks = c(100, 200, 500, 1000), trans = "log") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20), trans = "log") +
  ylab(expression("Biomass pound"~acre^{-1}~"(ln)")) +
  xlab("Species richness (ln)") +
  annotate(geom = "text", x = 15, y = 170, label = bquote(r^2~" = "~.(r2) ), size = 3) +
  annotate(geom = "text", x = 16.55, y = 120, label = bquote(F[.(df1)~","~.(df2)]~" = "~.(Fst) ), size = 3) +
  annotate(geom = "text", x = 14.75, y = 90, label = bquote("P = "~.(pval) ), size = 3) +
  theme_meta()

ggsave(filename = "figures/Carlander_1955_data/Fig_X1.png",
       unit = "cm", width = 10, height = 7.5, dpi = 400)





