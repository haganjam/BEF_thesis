#'
#' @title Compare selection-complementarity to stabilising-equalising niche differences
#' 
#' @description Uses the Lotka-Volterra model from Loreau (2004) to compare
#' selection-complementarity effects when the linearity assumption is met and 
#' species coexist to niche and fitness differences
#'

# load relevant scripts
source("code/helper-plotting-theme.R")
source("code/helper-BEF-partitions.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(deSolve)
library(tidyr)

# generalised lotka-volterra model functions
# https://stefanoallesina.github.io/Sao_Paulo_School/intro.html#multi-species-dynamics

# Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * (r + A %*% x)
    list(dxdt)
  })
}
# general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.5){
  
  # two-species mixture
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  # make into data.frame
  out <- as.data.frame(out)
  colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
  out
}

# make an output list
n <- 1000
list_out <- vector("list", length = n)

# loop over these replicates
for(i in 1:n) {
  
  # print the run
  print(paste0("loop ", i))
  
  # get intrinsic growth rates
  r <- round(runif(2, 0.2, 1.5), 1)
  
  # intraspecific
  intra <- runif(n = 1, 0.7, 1)
  
  # get random interaction coefficients
  A <- runif(1, 0.1, 1)
  A <- -matrix(c(intra, A, A, intra), 2, 2, byrow = TRUE)
  
  # run the mixtures
  x0 <- runif(2)
  y <- integrate_GLV(r, A, x0)
  
  # each species in monoculture
  
  # species 1
  m1 <- integrate_LG(r[1], A[1, 1], x0[1])
  
  # species 2
  m2 <- integrate_LG(r[2], A[2, 2], x0[2])
  
  # calculate niche and fitness differences (Spaak et al. 2023)
  
  # niche differences
  p <- sqrt( (-A[1, 2]*-A[2, 1])/(-A[1, 1]*-A[2, 2]) )
  
  # fitness differences
  f1 <- r[2] * sqrt( (-A[1, 2]*-A[1, 1]) )
  f2 <- r[1] * sqrt( (-A[2, 1]*-A[2, 2]) )
  
  # calculate the invasion growth rate of each species
  # https://core.ac.uk/download/pdf/83636816.pdf
  
  # species 1
  inv_r1 <- r[1] - (m2[nrow(m2),][["sp_1"]]*(-A[1, 2]))
  
  # species 2
  inv_r2 <- r[2] - (m1[nrow(m1),][["sp_1"]]*(-A[2, 1]))
  
  # create a BEF_df object
  BEF_df <- data.frame(sample = 1,
                       species = c(1, 2),
                       M = c(m1[nrow(m1),][["sp_1"]], m2[nrow(m2),][["sp_1"]]),
                       Y = c(y[nrow(y),][["sp_1"]], y[nrow(y),][["sp_2"]]))
  
  # apply the Loreau and Hector (2001) partition
  nbe <- local_scale_part(data = BEF_df, RYe = c(0.5, 0.5), part = "loreau_2001")
  
  # pull all this information into a data.frame
  df <- data.frame(run = i,
                   inv_sp1 = inv_r1,
                   inv_sp2 = inv_r2,
                   inv_coex = all(c(inv_r1, inv_r2) > 0),
                   coex = all(BEF_df$Y > 0),
                   niche_overlap = p,
                   niche_difference = 1-p,
                   fitness1 = f1,
                   fitness2 = f2,
                   fitness3 = ifelse(f1 > f2, (f1/f2), (f2/f1)),
                   nbe = nbe$NBE,
                   se = nbe$SE,
                   ce = nbe$CE,
                   ryt = sum(BEF_df$Y/BEF_df$M),
                   coex_nf = ifelse(f1 > f2, p < (f2/f1), p < (f1/f2)))
  row.names(df) <- NULL
  
  # add to the list
  list_out[[i]] <- df
  
}

# bind into a data.frame
df_out <- dplyr::bind_rows(list_out)
head(df_out)

# check if any pairs are not coexisting
any(df_out$coex == FALSE)
any(df_out$inv_coex == FALSE)
all(df_out$inv_coex == df_out$coex_nf)

# plot the results
ggplot(data = df_out,
       mapping = aes(x = niche_difference, y = fitness2/fitness1, colour = ryt)) +
  geom_point(shape = 1, alpha = 0.5) +
  theme_bw()

ggplot(data = df_out,
       mapping = aes(x = niche_difference, y = fitness2/fitness1, colour = coex_nf)) +
  geom_point(shape = 1, alpha = 0.5) +
  theme_bw()

ggplot(data = df_out,
       mapping = aes(x = niche_difference, y = fitness3, colour = coex_nf)) +
  geom_point(shape = 1, alpha = 0.5) +
  theme_bw()

# fitness differences
ggplot(data = df_out,
       mapping = aes(x = fitness, y = se, colour = inv_coex)) +
  geom_point(shape = 1, alpha = 0.5) +
  theme_bw()

# plot niche versus fitness differences
ggplot(data = df_out,
       mapping = aes(x = niche, y = fitness, colour = inv_coex)) +
  geom_point(alpha = 0.5) +
  theme_bw()


