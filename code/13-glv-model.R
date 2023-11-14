#'
#' @title Generalised Lotka-Volterra competition model
#' 
#' @description Functions for integrating the generalised Lotka-Volterra model for
#' competition of N-species
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(deSolve)
library(tidyr)

# Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-4] <- 0 # prevent numerical problems
    dxdt <- (x*r)*(1 - ((A %*% x)/K))
    list(dxdt)
  })
}

# function to plot output
plot_ODE_output <- function(out, Kmax){
  out <- as_tibble(out) |> tidyr::gather(species, density, -time)
  pl <- ggplot(data = out) + 
    aes(x = time, y = density, colour = species) + 
    geom_line() +
    expand_limits(y=-0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = Kmax, linetype = "dashed", colour = "red") +
    theme_classic()
  show(pl)
  return(out)
}

# general function to integrate GLV
integrate_GLV <- function(r, A, K, x0, maxtime = 200, steptime = 1){
  
  # two-species mixture
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A, K = K)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  # make into data.frame
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
  out
}

# run a few examples

# set the number of species
S <- 5

# get intrinsic growth rates
r <- runif(n = S, min = 0.5, max = 2.5)
print(r)

# get random interaction coefficients

# get the number of pairwise interactions
npairs <- ncol(combn(x = 1:S, m = 2))

# draw pairwise interaction coefficients for all possible pairs
aij <- runif(n = npairs, min = -0.5, max = 1)

# do some matrix manipulation to make a symmetric matrix
A <- diag(S)
A[upper.tri(A, diag=FALSE)] <- aij
A[lower.tri(A)] = t(A)[lower.tri(A)]
print(A)

# set the carrying capacities
K <- runif(n = S, min = 0.5, max = 2.5)

# set starting abundances of each species
x0 <- rep(x = 0.5, times = S)

# integrate the model
pops <- integrate_GLV(r, A, K, x0, maxtime = 200)

# plot the results
plot_ODE_output(pops, Kmax = max(K))

# END
