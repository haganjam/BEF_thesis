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

# set the number of landscapes
L <- 100

# set the number of patches within a landscape
P <- 20

# get a vector with the different diversities for each landscape
divT <- rep(seq(1:N), each = L/N)
length(divT)

# set the number of time steps
t <- 300

# set species intrinsic growth rates
r <- seq(log(0.2)-0.1, log(0.2)+0.1, length.out = N)
r <- exp(r)
print(exp(r))

# check that the geometric mean is 0.2
exp(mean(log((r))))

# set the case:
# 1 - partitioning between patches in a landscape
# 2 - partitioning within patches in a landscape
case <- 2

# set the parameters
if(case == 1) {
  
  # set the K values to vary by patch type (Thompson et al. 2021)
  Kmax <- 1000
  sd <- 1.5
  K <- sapply(1:N, function(x) {
    # define the environment
    env <- c((x+x-1):x, (x+1):P)
    # define the patch specific K
    Kmax*(exp((x-( env ))/(2*(sd)))^2)
  }  )
  
  # round to the nearest integer
  K <- lapply(K, function(x) round(x, 0) )
  
  # convert to a matrix
  K <- do.call("rbind", K)
  
  # make sure the K is non-zero
  K[K == 0] <- 0.1
  
  # set a to 1
  a <- 1
  
} else {
  
  # set the carrying capacities to be equal
  K <- matrix(rep(200, N*P), nrow = P, ncol = N)
  
  # set alpha as either 1 or 0.2
  a <- c(0.2)
  
}

# set-up the starting biomass
b0 <- 500

# loop over each landscape
lbi <- vector("list", length = L)
for(l in 1:L) {
  
  # get the species list
  sp <- sample(x = 1:N, size = divT[l])
  
  # set an output list for each patch
  pbi <- vector("list", length = P)
  
  # loop over each patch
  for(k in 1:P) {
    
    # set-up a list with starting abundances
    bi <- vector("list", length = t)
    bi[[1]] <- rep(0, N)
    bi[[1]][sp] <- rep(b0/divT[l], divT[l])
    
    # get the patch specific carrying capacities
    Kp <- K[k,]
    
    # run the simulation
    for(i in 2:t) {
      btn <- vector(length = N)
      for(j in 1:N) {
        btn[j] <- bi[[i-1]][j]*exp(r[j] * (1 - ((bi[[i-1]][j] + (a*sum(bi[[i-1]][-j])) )/Kp[j])))
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
      dplyr::mutate(landscape = l,
                    patch = k,
                    init_SR = divT[l],
                    time = 1:t) |> 
      tidyr::pivot_longer(cols = dplyr::starts_with("sp"),
                          names_to = "species",
                          values_to = "abundance")
    
    pbi[[k]] <- bi_df
    
  }
  
  # bind the output into a data.frame
  pbi_df <- dplyr::bind_rows(pbi)
  
  # write into the landscape level output
  lbi[[l]] <- pbi_df 
  
}

# bind the output into a data.frame
lbi_df <- dplyr::bind_rows(lbi)
head(lbi_df)

# patch scale summary
lbi_patch <- 
  lbi_df |> 
  dplyr::filter(time == t) |> 
  dplyr::group_by(landscape, patch) |> 
  dplyr::summarise(init_SR = first(init_SR),
                   real_SR = sum(abundance > 0),
                   total_abun = sum(abundance)) 

# plot the BEF relationship at the patch scale
ggplot(data = lbi_patch,
       mapping = aes(x = init_SR, y = total_abun)) +
  geom_point() +
  geom_smooth()

# landscape scale summary
lbi_land <- 
  lbi_df |> 
  dplyr::filter(time == t) |> 
  dplyr::group_by(landscape) |> 
  dplyr::summarise(init_SR = first(init_SR),
                   real_SR = sum(abundance > 0),
                   total_abun = sum(abundance)) 

# plot the BEF relationship at the patch scale
ggplot(data = lbi_land,
       mapping = aes(x = init_SR, y = total_abun)) +
  geom_point() +
  geom_smooth()




