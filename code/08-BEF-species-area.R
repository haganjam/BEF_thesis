#'
#' @title Scaling up BEF using the species area relationship
#' 
#' @description Illustration of Tilman's idea to use the species-area
#' relationship to scale-up the BEF relationship
#'

# load relevant libraries
library(dplyr)
library(ggplot2)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# set-up the calculation as per Tilman (1999)

# set the z-values
zvals <- seq(0.15, 0.30, 0.01)

# set the experimentally determined local diversity required
SL <- 16

# set the area of the local plot
AL <- 1

# set the area of the region
AR <- 1000000

# calculate the regional richness required
SR <- vector(length = length(zvals))
for(i in 1:length(zvals)) {
  
  SR[i] <- SL*((AR/AL)^zvals[i])
  
}

# check the results
print(SR)

# check the rounded range
range(round(SR, 0))








