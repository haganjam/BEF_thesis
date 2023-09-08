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
P <- 2000

# get a vector with the different diversities for each patch
divT <- rep(seq(1:N), each = P/N)

# set the number of time steps
t <- 300

# set species intrinsic growth rates
r <- seq(log(0.2)-0.1, log(0.2)+0.1, length.out = N)
r <- exp(r)
print(exp(r))

# check that the geometric mean is 0.2
exp(mean(log((r))))

# set the case
case <- 1
if(case == 1) {
  K <- rnorm(n = N, mean = 500, sd = 175)
  a <- 1
} else {
  K <- 200
  a <- 0.2
}

# check the K parameters
hist(K)
range(K)

# set-up the starting biomass
b0 <- 500

# set an output list for each patch
pbi <- vector("list", length = P)

# loop over each patch
for(k in 1:P) {
  
  # get the species list
  sp <- sample(x = 1:N, size = divT[k])
  
  # set-up a list with starting abundances
  bi <- vector("list", length = t)
  bi[[1]] <- rep(0, N)
  bi[[1]][sp] <- rep(b0/divT[k], divT[k])
  
  # run the simulation
  for(i in 2:t) {
    btn <- vector(length = N)
    for(j in 1:N) {
      btn[j] <- bi[[i-1]][j]*exp(r[j] * (1 - ((bi[[i-1]][j] + (a*sum(bi[[i-1]][-j])) )/K[j])))
    }
    
    bi[[i]] <- btn
  }
  
  # bind into a data.frame
  bi_df <- dplyr::as_tibble(do.call("rbind", bi))
  head(bi_df)
  
  # rename to species
  names(bi_df) <- paste0("sp", 1:N)
  
  # bring into the long format
  bi_df <- 
    bi_df |> 
    dplyr::mutate(time = 1:t) |> 
    tidyr::pivot_longer(cols = dplyr::starts_with("sp"),
                        names_to = "species",
                        values_to = "abundance")
  
  pbi[[k]] <- bi_df
  
}


# plot the results
ggplot(data = bi_df,
       mapping = aes(x = time, y = abundance, colour = species)) +
  geom_line()









