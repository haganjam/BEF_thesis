#'
#' @title Cardinale et al.'s (2004, Oikos) model
#' 
#' @description Using Cardinale et al.'s (2004) model to explore some
#' key theoretical results of the scale-dependence of the BEF
#' relationship
#'

# load relevant libraries
library(dplyr)
library(ggplot2)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# set the number of species
N <- 20

# set the number of patches
P <- 100

# set species intrinsic growth rates
r <- seq(log(0.2)-0.1, log(0.2)+0.1, length.out = N)
print(exp(r))

# check that the geometric mean is 0.2
exp(mean(log(exp(r))))

# set the case
case <- 1
if(case == 1) {
  K <- rnorm(n = N, mean = 500, sd = 175)
  a <- rep(1, N*N)
} else {
  K <- 200
  a <- rep(0.2, N*N)
}

# check the K parameters
hist(K)
range(K)

# convert the a parameters into a matrix
a <- matrix(a, nrow = N, ncol = N)

# set-up the starting biomass
b0 <- 100

# run the simulation


# set the output list












