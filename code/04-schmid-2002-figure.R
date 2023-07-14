
# reproduce Figs from Schmid 2002

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Thesis_document")

# get the meta_theme()
source("figures/Function_plotting_theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# generate a hump-shaped curve
x <- seq(-0.1, 10.1, 0.1)
y <- (5.5*x) - (0.55*(x^2))
plot(x, y)

# pull into a data.frame
df <- data.frame(prod = x, div = y)
print(df)

# generate a few random points under the curve
df$div_raw <- round(runif(n = nrow(df), min = 0+0.25, max = df$div-0.25), 0)

# randomly sample 30
df[sample(x = 1:nrow(df), 40), ]$div_raw <- NA

# check the summary statistics of df
summary(df)

# plot the graph
p1 <- 
  ggplot() +
  geom_point(data = df, 
             mapping = aes(x = prod, y = div_raw)) +
  geom_line(data = df,
            mapping = aes(x = prod, y = div), 
            size = 0.75,colour = "red", linetype = "dashed") +
  ylab("Species richness") +
  xlab("Productivity") +
  scale_y_continuous(limits = c(0, 14.5)) +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white", size = 1),
        axis.text.y = element_text(colour = "white", size = 1))

# plot the sideways version

# plot two within site curves
x <- seq(-0.1, 8.35, 0.1)
dfA <- data.frame(prod = x, div = 1.3^(x)-1)

x <- seq(-0.1, 5.7, 0.1)
dfB <- data.frame(prod = x, div = 1.6^(x)-1)

# get the end poitns
dfAB <- bind_rows(dfA[dfA$div == max(dfA$div),], dfB[dfB$div == max(dfB$div),])
dfAB$point <- c("A", "B")

p2 <- 
  ggplot() +
  geom_line(data = df,
            mapping = aes(x = prod, y = div), 
            size = 0.75,colour = "red", linetype = "dashed") +
  geom_line(data = dfA,
            mapping = aes(x = prod, y = div)) +
  geom_line(data = dfB,
            mapping = aes(x = prod, y = div)) +
  geom_point(data = dfAB, 
             mapping = aes(x = prod, y = div), size = 3) +
  geom_text(data = dfAB, 
             mapping = aes(x = prod, y = div, label = point), size =4,
             nudge_y = 0.65, nudge_x = 0.5) +
  scale_y_continuous(limits = c(0, 14.5)) +
  ylab("Species richness") +
  xlab("Productivity") +
  coord_flip() +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white", size = 1),
        axis.text.y = element_text(colour = "white", size = 1))

p12 <- 
  ggarrange(p1, p2, widths = c(1.5, 1),
            labels = c("a", "b"),
            font.label = list(size = 11, face = "plain")
            )

ggsave(filename = "figures/Schmid_2002/Fig_Schmid_2002.png", plot = p12,
       unit = "cm", width = 17, height = 7.5, dpi = 400)

### END
