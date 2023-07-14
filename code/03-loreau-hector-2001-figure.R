
# Loreau and Hector (2001): Selection-complementarity partition

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Thesis_document")

# get the meta_theme()
source("figures/Function_plotting_theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# define examples from Fox (2005)
f1 <- data.frame(sample = rep(c(1, 2, 3, 4, 5, 6), each = 2),
                 species = rep(c(1, 2), 6),
                 M = rep(c(500, 250), 6),
                 Y = c(300, 100, 330, 110, 360, 120, 390, 130, 420, 140, 450, 150))

# terms used:

# dY = deviation from total expected yield in mixture (i.e. Yo - Ye)
# Yo = observed mixture yield
# Ye = expected mixture yield (based on initial proportions and monoculture performance)

# Net Biodiversity Effect = Yo - Ye (i.e. dY)

# for each mixture

# Mi = monoculture of each species
# Yoi = observed yield of each species in mixture
# Yo = observed mixture yield - sum(Yoi)

# RYei = expected relative yield of each species (1/n spp)
# RYoi = observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields

# Yei = expected yield of each species in mixture (RYe*Mi)
# Ye = expected yield of mixture - sum(Yei)

# dRY = RYoi - RYei

# N = n spp in mixture

# Poi = observed proportion of each species in mixture (i.e. RYoi/sum(RYoi))


# define function for number of unique elements in a vector
n_unique <- function(x) length(unique(x))

# define function for the raw covariance
raw_cov <- function(x, y) {
  c1 <- x - mean(x)
  c2 <- y - mean(y)
  sum( (c1*c2) )/length(c1)
}

# define the preparation function
B.ef.prep <- function(df, RYe) {
  
  if (n_unique(df$species) != length(RYe) ) {
    stop("error, expected frequencies not defined")
  }
  
  # define expected relative yields
  df$RYe <- rep(RYe, n_unique(df$sample))
  
  # define observed relative yields
  df$RYo <- (df$Y/df$M)
  
  # define the change in relative yield
  df$dRY <- (df$RYo - df$RYe)
  
  # calculate expected yield for each species
  df$Ye <- (df$M*df$RYe)
  
  # calculate observed proportion of each species in mixture (po,ijk, Isbell et al. 2018)
  df$Poi <- unlist(lapply(split(df, df$sample), 
                          function(x) {
                            x$Y/sum(x$Y) }), use.names = FALSE)
  
  # calculate change in observed proportion relative to the expectation (d.po,ijk, Isbell et al. 2018)
  df$d.Poi <- (df$Poi - df$RYe)
  
  # calculate change in observed proportion (dRYo,ijk Isbell et al. 2018)
  df$d.RYoi <- (df$RYo - df$Poi)
  
  if(TRUE %in% grepl(pattern = "time|place", names(df))) {
    
    # define extra data.frames for means at different levels
    
    # species means for each time across places (pij and Mij single bar, Isbell et al. 2018)
    sm_t <- aggregate(df[, c("d.Poi", "M") ], list(df$species, df$time), mean)
    names(sm_t) <- c("species", "time", "d.Poi.t", "M.t")
    df <- merge(df, sm_t, all.x = TRUE)
    
    # species means for each place across times (pik and Mik single bar, Isbell et al. 2018)
    sm_p <- aggregate(df[, c("d.Poi", "M") ], list(df$species, df$place), mean)
    names(sm_p) <- c("species", "place", "d.Poi.p", "M.p")
    df <- merge(df, sm_p, all.x = TRUE)
    
    # overall species mean across all times and places
    sm_s <- aggregate(df[, c("d.Poi", "M") ], list(df$species), mean)
    names(sm_s) <- c("species", "d.Poi.s", "M.s")
    df <- merge(df, sm_s, all.x = TRUE)
    
    # reorder the columns
    df <- df[order(df$time,df$place,df$species), ]
    
    df.p <- list(f = df,
                 s = sm_s,
                 p = sm_p,
                 t = sm_t)
    
    return(df.p)
    
  } else ( return(df) )
  
}

# define a function to perform the Loreau and Hector (2001) partition
hector.loreau.2001.pt <- function(adf, RY.exp) {
  
  # check the there are times or places in the data.frame
  if(TRUE %in% grepl(pattern = "time|place", names(adf))) {
    stop("error, remove time and place columns from the data.frame or perform the whole spatial-temporal partition using isbell.2018.pt() function")
  }
  
  # use the B.ef.prep function to define relevant terms
  df.p <- B.ef.prep(df = adf, RYe = RY.exp)
  
  # define the number of species
  N <- n_unique(df.p$species)
  
  # calculate the net biodiversity effect
  nbe_df <- aggregate(df.p[, c("Y", "Ye") ], list(df.p$sample), sum)
  NBE.a <- (nbe_df$Y - nbe_df$Ye)
  
  # calculate the complementarity effect (trait independent complementarity sensu Fox 2005)
  comp_df <- aggregate(df.p[, c("dRY", "M")], list(df.p$sample), mean)
  CE.a <- N*(comp_df$dRY)*(comp_df$M)
  
  # calculate the selection effect
  SE.a <- sapply(split(df.p, df.p$sample), 
                 function(x) { N*raw_cov(x$dRY,x$M) })
  
  # wrap this into a data.frame
  be_out <- data.frame(sample = unique(df.p$sample),
                       net.biodiversity.effect = NBE.a,
                       complementarity.effect = CE.a,
                       selection.effect = SE.a)
  
  return(be_out)
  
}

# test the function with the Fox (2005) example
hector.loreau.2001.pt(adf = f1, RY.exp = c(0.6, 0.4))

# set-up a function to generate relevant test data
test_data <- function(S = 1, M1, M2, Y1, Y2) {
  
  df1 <- 
    data.frame(sample = rep(S, 2),
               species = c(1, 2),
               M = c(M1, M2),
               Y = c(Y1, Y2))
  
  df2 <- 
    data.frame(species = as.character(rep(1:2, 4)),
               mixture = factor(c("Sp.1", "Sp.1", "Sp.2", "Sp.2", "Mix (Exp.)", "Mix (Exp.)", "Mix", "Mix"),
                                levels = c("Sp.1", "Sp.2", "Mix (Exp.)", "Mix")),
               yield = c(M1, 0, 0, M2, M1/2, M2/2, Y1, Y2))
  
  df3 <- 
    test1.seg <- data.frame(xstart = unique(df2$mixture),
                            xend = unique(df2$mixture),
                            ystart = c(NA, NA, NA, (M1+M2)/2),
                            yend = c(NA, NA, NA, (Y1+Y2)))
  
  return(list(df1, df2, df3))
  
}

# generate the first test data.frame
tst1 <- test_data(S = 1, M1 = 200, M2 = 500, Y1 = 0, Y2 = 500)

# check whether we can calculate observed mixture functioning based on relative yield
tst1[[1]]$RYO <- tst1[[1]]$Y/tst1[[1]]$M
tst1[[1]]$RYO*tst1[[1]]$M == tst1[[1]]$Y

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
pt1 <- hector.loreau.2001.pt(adf = tst1[[1]], RY.exp = c(0.5, 0.5))

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
pt2 <- hector.loreau.2001.pt(adf = tst2[[1]], RY.exp = c(0.5, 0.5))

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
pt3 <- hector.loreau.2001.pt(adf = tst3[[1]], RY.exp = c(0.5, 0.5))

# bind the partition data 
pt <- bind_rows(pt1, pt2, pt3)
row.names(pt) <- NULL
pt$sample <- c("a", "b", "c")
names(pt) <- c("sample", "NBE", "CE", "SE")

pt <- 
  pt %>%
  pivot_longer(cols = c("NBE", "CE", "SE"),
               names_to = "effect",
               values_to = "value") %>%
  filter(effect != "NBE")

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

# arrange these plot
p1234 <- ggarrange(p1, p2, p3, p4, labels = c("a", "b", "c", "d"),
                   font.label = list(size = 11, face = "plain"))
plot(p1234)

ggsave(filename = "figures/Loreau_Hector_2001_Fig/Fig_Y.png", plot = p1234,
       unit = "cm", width = 20, height = 16, dpi = 400)

### END
