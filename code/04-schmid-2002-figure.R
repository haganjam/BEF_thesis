#'
#' @title Schmid (2002) figure
#' 
#' @description Generates an example of the figure that Schmid (2002) presented
#' in his summary article about the BEF relationship and how it applies to 
#' experimental and non-experimental systems
#'

# load relevant functions
source("code/helper-plotting-theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# generate a hump-shaped curve
x <- seq(-0.1, 10.1, 0.1)
y <- (5.5*x) - (0.55*(x^2))
plot(x, y)

# pull into a data.frame
sch_dat <- data.frame(prod = x, div = y)
print(sch_dat)

# generate a few random points under the curve
sch_dat$div_raw <- round(runif(n = nrow(sch_dat), min = 0+0.25, max = sch_dat$div-0.25), 0)

# randomly sample 40
sch_dat[sample(x = 1:nrow(sch_dat), 40), ]$div_raw <- NA

# check the summary statistics of df
summary(sch_dat)

# plot the graph
p1 <- 
  ggplot() +
  geom_point(data = sch_dat, 
             mapping = aes(x = prod, y = div_raw), size = 2, colour = "#118176") +
  geom_line(data = sch_dat,
            mapping = aes(x = prod, y = div), 
            size = 0.75,colour = "#118176", linetype = "dashed") +
  ylab("Species richness") +
  xlab("Productivity") +
  scale_y_continuous(limits = c(0, 14.5)) +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white", size = 1),
        axis.text.y = element_text(colour = "white", size = 1))
plot(p1)

# plot the sideways version

# plot two within site curves
x <- seq(-0.1, 8.35, 0.1)
dfA <- data.frame(prod = x, div = 1.3^(x)-1)

x <- seq(-0.1, 5.7, 0.1)
dfB <- data.frame(prod = x, div = 1.6^(x)-1)

# get the end poitns
dfAB <- dplyr::bind_rows(dfA[dfA$div == max(dfA$div),], dfB[dfB$div == max(dfB$div),])
dfAB$point <- c("A", "B")

p2 <- 
  ggplot() +
  geom_line(data = sch_dat,
            mapping = aes(x = prod, y = div), 
            size = 0.75,colour = "#118176", linetype = "dashed") +
  geom_line(data = dfA,
            mapping = aes(x = prod, y = div), colour = "#118176") +
  geom_line(data = dfB,
            mapping = aes(x = prod, y = div), colour = "#118176") +
  geom_point(data = dfAB, 
             mapping = aes(x = prod, y = div), size = 3, colour = "#118176") +
  geom_text(data = dfAB, 
             mapping = aes(x = prod, y = div, label = point), size =4,
             nudge_y = 0.65, nudge_x = 0.5, colour = "#118176") +
  scale_y_continuous(limits = c(0, 14.5)) +
  ylab("Species richness") +
  xlab("Productivity") +
  coord_flip() +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white", size = 1),
        axis.text.y = element_text(colour = "white", size = 1))
plot(p2)

p12 <- 
  ggarrange(p1, p2, widths = c(1.5, 1),
            labels = c("a", "b"),
            font.label = list(size = 11, face = "plain")
            )
plot(p12)

ggsave(filename = "figures-tables/fig_6.pdf", plot = p12,
       unit = "cm", width = 17, height = 7.5)

### END
