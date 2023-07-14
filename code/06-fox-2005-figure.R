#'
#' @title local_scale_part()
#' 
#' @description Function to calculate Loreau and Hector's (2001, Nature) and
#' Fox (2005, Ecology Letters) biodiversity effect partition using species 
#' mixture and monoculture data
#' 
#' @param data data.frame in the following format
#' column 1 - sample: variable specifying the unique place-time combination
#' column 2 - species: variable specifying the species name (all times and places must have all species names present)
#' column 3 - M: monoculture functioning
#' column 4 - Y: mixture function
#' @param RYe expected relative yields for the species as a numeric vector of length = (N species) and which sums to one
#' @param part either "loreau_2001" or "fox_2005"
#' 
#' @symbol Mi - monoculture of each species
#' @symbol Yoi - observed yield of each species in mixture
#' @symbol Yo - observed mixture yield - sum(Yoi)
#' @symbol RYei - expected relative yield of each species (usually 1/n species but can be anything)
#' @symbol RYoi - observed relative yield (Yoi/Mi) i.e. measures the extent to which species i overyields
#' @symbol dRY - RYoi - RYei
#' @symbol N - n spp in mixture
#' 

# load the helper functions
source("functions/helper_functions.R")

# install and load libraries required for these functions
install_if("dplyr")
install_if("assertthat")

local_scale_part <- function(data, RYe, part = "loreau_2001") {
  
  # test if the input data is a data.frame
  test_1 <- function(x) {
    
    is.data.frame(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(call$x, " is not a data.frame")
    
  }
  
  assertthat::assert_that(test_1(x = data))
  
  # test if all the required columns are present
  test_2 <- function(x) {
    
    all(names(x) %in% c("sample", "species", "M", "Y"))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is missing one of the following columns: sample, species, M, Y")
    
  }
  
  assertthat::assert_that(test_2(x = data))
  
  # test if RYe is a number
  test_3 <- function(x) {
    
    is.vector(x) & is.numeric(x)
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(call$x, " is not a numeric vector")
    
  }
  
  assertthat::assert_that(test_3(x = RYe))
  
  # test if the length of the RYe vector equals the number of species
  test_4 <- function(x, y) {
    
    assertthat::are_equal(length(x), y)
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0("x and y do not have equal length")
    
  }
  
  assertthat::assert_that(test_4(x = RYe, y = n_unique(data[["species"]])) )
  
  # test if the RYe vector sums to 1
  test_5 <- function(x) {
    
    sum(x) > 0.99
    
  }
  
  assertthat::on_failure(test_5) <- function(call, env){
    
    paste0("RYe does not sum to one")
    
  }
  
  assertthat::assert_that(test_5(x = RYe) )
  
  # get the number of species in the data
  n_sp <- n_unique(data$species)
  
  # sort the data.frame
  df <- 
    data %>%
    arrange(sample, species)
  
  # define expected relative yields
  df$RYe <- rep(RYe, n_unique(df$sample))
  
  # define observed relative yields: Prevent dividing by zero
  df$RYo <- ifelse(df$M == 0, 0, (df$Y/df$M))
  
  # define the change in relative yield
  df$dRY <- (df$RYo - df$RYe)
  
  # calculate expected yield for each species
  df$Ye <- (df$M*df$RYe)
  
  # calculate the net biodiversity effect
  nbe_df <- aggregate(df[, c("Y", "Ye") ], list(df$sample), sum)
  NBE <- (nbe_df$Y - nbe_df$Ye)
  
  if (part == "loreau_2001") {
    
    # calculate the complementarity effect (trait independent complementarity sensu Fox 2005)
    comp_df <- aggregate(df[, c("dRY", "M")], list(df$sample), mean)
    CE <- n_sp*(comp_df$dRY)*(comp_df$M)
    
    # calculate the selection effect
    SE <- sapply(split(df, df$sample), function(x) { n_sp*raw_cov(x$dRY, x$M) })
    
    # wrap this into a tibble()
    BEF_df <- tibble(sample = unique(df$sample),
                     NBE = NBE,
                     SE = SE,
                     CE = CE
    )
    
  } else if (part == "fox_2005") {
    
    # calculate trait independent complementarity sensu Fox 2005
    comp_df <- aggregate(df[, c("dRY", "M")], list(df$sample), mean)
    TI_CE <- n_sp*(comp_df$dRY)*(comp_df$M)
    
    # calculate trait dependent complementarity
    TD_CE <- sapply(split(df, df$sample), 
                    function(x) { n_sp*raw_cov(x$M, (x$RYo - (x$RYo/sum(x$RYo)) ) ) })
    
    # calculate the dominance effect
    DOM <- sapply(split(df, df$sample), 
                  function(x) { n_sp*raw_cov(x$M, ((x$RYo/sum(x$RYo)) - x$RYe) ) })
    
    # wrap this into a data.frame
    BEF_df <- tibble(sample = unique(df$sample),
                     NBE = NBE,
                     TI_CE = TI_CE ,
                     TD_CE = TD_CE,
                     DOM = DOM
    )
    
  } else {
    
    stop("Choose appropriate option for the `part` argument")
    
  }
  
  return(BEF_df)
  
}


# define examples from Fox (2005)
f1 <- data.frame(sample = rep(c(1, 2, 3, 4, 5, 6), each = 2),
                 species = rep(c(1, 2), 6),
                 M = rep(c(500, 250), 6),
                 Y = c(300, 100, 330, 110, 360, 120, 390, 130, 420, 140, 450, 150))

# define the answers from Fox (2005
f1.ans <- data.frame(SE = c(0, 2.5, 5, 7.5, 10, 12.5),
                     CE = c(0, 37.5, 75, 112.5, 150, 187.5))

# test if the function correctly calculates the biodiversity effects sensu Loreau and Hector (2001, Nature)

# use the partition to calculate the biodiversity effects
df.test <- local_scale_part(data = f1, RYe = c(0.60, 0.40), part = "loreau_2001")
SE.test <- all(near(df.test$SE, f1.ans$SE))
CE.test <- all(near(df.test$CE, f1.ans$CE))

# check if any of the biodiversity effects were incorrectly calculated
if ( any(c(SE.test, CE.test) == FALSE) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 
  
} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")
  
}

# test if the function correctly calculates the biodiversity effects sensu Fox (2005, Ecology Letters)

# use the partition to calculate the biodiversity effects
df.test <- local_scale_part(data = f1, RYe = c(0.60, 0.40), part = "fox_2005")
CE.test <- near(df.test$TI_CE, f1.ans$CE)
SE.test <- near((df.test$TD_CE+df.test$DOM), f1.ans$SE)

# check if any of the biodiversity effects were incorrectly calculated
if ( any(c(SE.test, CE.test) == FALSE) ) { 
  
  warning("Functions do not correctly calculate biodiversity effects of the test data") 
  
} else { 
  
  message("Functions correctly calculate biodiversity effects on the test data")
  
}

### END

# load relevant libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# set working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Thesis_document")

# get the meta_theme()
source("figures/Function_plotting_theme.R")

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
  pt %>%
  pivot_longer(cols = c("NBE", "TI_CE", "TD_CE", "DOM"),
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
library(ggpubr)
p1234 <- ggarrange(p1, p2, p3, p4, labels = c("a", "b", "c", "d"),
                   font.label = list(size = 11, face = "plain"))
plot(p1234)

ggsave(filename = "figures/Fox_2005_Fig/Fig_X.png", plot = p1234,
       unit = "cm", width = 20, height = 16, dpi = 400)

### END
