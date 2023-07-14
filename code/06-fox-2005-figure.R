#'
#' @title Fox (2005) partition
#' 
#' @description Generates some empirical examples of the Fox (2005) extension
#' of Loreau and Hector's (2001) partition to give readers a feel for how it works.
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# load relevant data
source("code/helper-BEF-partitions.R")
source("code/helper-plotting-theme.R")

# generate the first test data.frame
tst1 <- test_data(S = 1, M1 = 200, M2 = 500, Y1 = 0, Y2 = 500)

p1 <- 
  ggplot() +
  geom_bar(data = tst1[[2]],
           mapping = aes(x = mixture, y = yield, fill = species), 
           position = "stack", stat = "identity", width = 0.3,
           colour = "black") +
  xlab("") +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8, alpha = 0.85) +
  ylab("Ecosystem function") +
  geom_segment(data = tst1[[3]],
               mapping = aes(x = as.numeric(xstart)-0.5, xend = as.numeric(xend)-0.5,
                             y = ystart, yend = yend),
               colour = "red", size = 1 ) +
  annotate(x = 3.5, y = 525, geom = "text", label = "NBE", colour = "red") +
  theme_meta()+
  scale_y_continuous(limits = c(0, 615)) +
  theme(legend.position = "none")
plot(p1)

# calculate the partition
pt1 <- local_scale_part(data = tst1[[1]], RYe = c(0.5, 0.5), part = "fox_2005")

# generate the second test data.frame
tst2 <- test_data(S = 1, M1 = 200, M2 = 500, Y1 = (100*1.428571), Y2 = (250*1.428571))

p2 <- 
  ggplot() +
  geom_bar(data = tst2[[2]],
           mapping = aes(x = mixture, y = yield, fill = species), 
           position = "stack", stat = "identity", width = 0.3,
           colour = "black") +
  xlab("") +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8, alpha = 0.85) +
  ylab("Ecosystem function") +
  geom_segment(data = tst2[[3]],
               mapping = aes(x = as.numeric(xstart)-0.5, xend = as.numeric(xend)-0.5,
                             y = ystart, yend = yend),
               colour = "red", size = 1 ) +
  annotate(x = 3.5, y = 525, geom = "text", label = "NBE", colour = "red") +
  theme_meta()+
  scale_y_continuous(limits = c(0, 615)) +
  theme(legend.position = "none")
plot(p2)

# calculate the partition
pt2 <- local_scale_part(data = tst2[[1]], RYe = c(0.5, 0.5), part = "fox_2005")

# generate the second test data.frame
tst3 <- test_data(S = 1, M1 = 200, M2 = 500, Y1 = (100*1.3), Y2 = (250*1.7))

p3 <- 
  ggplot() +
  geom_bar(data = tst3[[2]],
           mapping = aes(x = mixture, y = yield, fill = species), 
           position = "stack", stat = "identity", width = 0.3,
           colour = "black") +
  xlab("") +
  scale_fill_viridis_d(option = "C", begin = 0, end = 0.8, alpha = 0.85) +
  ylab("Ecosystem function") +
  geom_segment(data = tst3[[3]],
               mapping = aes(x = as.numeric(xstart)-0.5, xend = as.numeric(xend)-0.5,
                             y = ystart, yend = yend),
               colour = "red", size = 1 ) +
  annotate(x = 3.5, y = 615, geom = "text", label = "NBE", colour = "red") +
  theme_meta()+
  scale_y_continuous(limits = c(0, 615)) +
  theme(legend.position = "none")
plot(p3)

# calculate the partition
pt3 <- local_scale_part(data = tst3[[1]], RYe = c(0.5, 0.5), part = "fox_2005")

# bind the partition data 
pt <- bind_rows(pt1, pt2, pt3)
row.names(pt) <- NULL
pt$sample <- c("a", "b", "c")
names(pt) <- c("sample", "NBE", "TI_CE", "TD_CE", "DOM")

pt <- 
  pt |>
  tidyr::pivot_longer(cols = c("NBE", "TI_CE", "TD_CE", "DOM"),
                      names_to = "effect",
                      values_to = "value") |>
  dplyr::filter(effect != "NBE")

p4 <- 
  ggplot(data = pt,
         mapping = aes(x = sample, y = value, fill = effect)) +
  geom_bar(stat = "identity", width = 0.3, colour = "black") +
  scale_fill_viridis_d(option = "F", begin = 0, end = 0.8, alpha = 0.85) +
  theme_meta() +
  ylab("Net biodiversity effect") +
  xlab("") +
  theme(legend.title = element_blank())
plot(p4)

# arrange these plots
p1234 <- ggpubr::ggarrange(p1, p2, p3, p4, labels = c("a", "b", "c", "d"),
                           font.label = list(size = 11, face = "plain"))
plot(p1234)

ggsave(filename = "figures-tables/fig_X.svg", plot = p1234,
       unit = "cm", width = 20, height = 16)

### END
