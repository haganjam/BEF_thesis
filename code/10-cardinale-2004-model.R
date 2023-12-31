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
L <- 600

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

# set the case corresponding to Fig. 9

# a - patch-scale, within patch niche partitioning
# c - landscape-scale, within patch niche partitioning
# b - patch-scale, no niche partitioning
# d - landscape-scale, between patch niche partitioning

# set ac or bd
case <- "bd"

# plots for the defence or thesis
defence <- TRUE

if(case == "ac") {
  
  # how many patches?
  P <- 20
  
  # set the carrying capacities to be equal
  K <- matrix(rep(200, N*P), nrow = P, ncol = N)
  
  # set alpha to 0.2
  a <- 0.2
  
} else if(case == "bd") {
  
  # how many patches?
  P <- 20
  
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
  
  # set alpha to 1
  a <- 1
  
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
  
  # extract the final time-point
  pbi_df <- dplyr::filter(pbi_df, time == t)
  
  # write into the landscape level output
  lbi[[l]] <- pbi_df 
  
}

# bind the output into a data.frame
lbi_df <- dplyr::bind_rows(lbi)
head(lbi_df)

# patch scale summary
lbi_patch <- 
  lbi_df |> 
  dplyr::group_by(landscape, patch) |> 
  dplyr::summarise(init_SR = first(init_SR),
                   real_SR = sum(abundance > 0),
                   total_abun = sum(abundance)) |> 
  dplyr::ungroup()

# get maximum monoculture
pm <- max(lbi_patch[lbi_patch$init_SR == 1, ]$total_abun)

# power function using nls()
nlm1 <- nls(total_abun ~ a*(init_SR^b), data = lbi_patch, start = list(a = 1 , b = 2))

# extract model parameters and standard errors
nlm1_sum <- summary(nlm1)
nlm1_sum$coefficients

# make a predicted dataset
pred_patch <- data.frame(init_SR = seq(0.8, 20, 0.1))
x <- investr::predFit(nlm1, pred_patch, interval = "confidence")
pred_patch <- dplyr::bind_cols(pred_patch, x)

# plot the BEF relationship at the patch scale
y <- 
  ggplot() +
  geom_hline(yintercept = pm, linetype = "longdash") +
  geom_jitter(data = lbi_patch,
              mapping = aes(x = init_SR, y = total_abun),
              colour = "#118176", alpha = 0.15, width = 0.2) +
  geom_ribbon(data = pred_patch,
              mapping = aes(x = init_SR, ymin = lwr, ymax = upr),
              alpha = 0.05, fill = "#118176") +
  geom_line(data = pred_patch,
            mapping = aes(x = init_SR, y = fit), size = 0.5, colour = "#118176") +
  xlab("Initial species richness") +
  ylab("Community biomass") +
  theme_meta()
plot(y)

# plot for the defence
if(defence) {
  yd1 <- 
    ggplot() +
    geom_hline(yintercept = pm, linetype = "longdash") +
    geom_jitter(data = lbi_patch |> dplyr::group_by(init_SR) |> dplyr::sample_n(10),
                mapping = aes(x = init_SR, y = total_abun),
                colour = "#118176", alpha = 1, width = 0.2) +
    geom_ribbon(data = pred_patch,
                mapping = aes(x = init_SR, ymin = lwr, ymax = upr),
                alpha = 0.05, fill = "#118176") +
    geom_line(data = pred_patch,
              mapping = aes(x = init_SR, y = fit), size = 0.5, colour = "#118176") +
    xlab("Initial species richness") +
    ylab("Community biomass") +
    theme_meta() +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
}

# assign the plot object
if(case == "ac") {
  assign("a1", y)
} else {
  assign("b1", y)
}

# plot the regional scale cases

# landscape scale summary
lbi_land <- 
  lbi_df |> 
  dplyr::group_by(landscape, patch) |>
  dplyr::summarise(init_SR = first(init_SR),
                   real_SR = sum(abundance > 0),
                   total_abun = sum(abundance)) |>
  dplyr::group_by(landscape) |>
  dplyr::summarise(init_SR = first(init_SR),
                   real_SR_m = mean(real_SR),
                   real_SR_sd = mean(real_SR),
                   total_abun_m = mean(total_abun),
                   total_abun_sd = sd(total_abun)
  ) |>
  dplyr::ungroup()

# power function using nls()
nlm2 <- nls(total_abun_m ~ a*(init_SR^b), data = lbi_land, start = list(a = 1 , b = 2))

# extract model parameters and standard errors
nlm2_sum <- summary(nlm2)
nlm2_sum$coefficients

# make a predicted dataset
pred_land <- data.frame(init_SR = seq(0.8, 20, 0.1))
x <- investr::predFit(nlm2, pred_land, interval = "confidence")
pred_land <- dplyr::bind_cols(pred_land, x)

# get maximum monoculture
lm <- max(lbi_land[lbi_land$init_SR == 1, ]$total_abun_m)

# make a factor for plotting
lbi_land$init_SR_f <- factor(lbi_land$init_SR) 

# plot the BEF relationship at the landscape scale
y <- 
  ggplot() +
  geom_hline(yintercept = lm, linetype = "longdash") +
  geom_jitter(data = lbi_land,
              mapping = aes(x = init_SR, y = total_abun_m), 
              colour = "#118176", alpha = 0.25, size = 2, width = 0.2) +
  geom_ribbon(data = pred_land,
              mapping = aes(x = init_SR, ymin = lwr, ymax = upr),
              alpha = 0.05, fill = "#118176") +
  geom_line(data = pred_land,
            mapping = aes(x = init_SR, y = fit), colour = "#118176", size = 0.5) +
  xlab("Initial species richness") +
  ylab("Community biomass") +
  theme_meta()
plot(y)

# make a plot for the defence
if(defence){
  yd2 <- 
  ggplot() +
  geom_hline(yintercept = lm, linetype = "longdash") +
  geom_jitter(data = lbi_land |> dplyr::group_by(init_SR) |> dplyr::sample_n(10),
              mapping = aes(x = init_SR, y = total_abun_m), 
              colour = "#118176", alpha = 0.25, size = 2, width = 0.2) +
  geom_ribbon(data = pred_land,
              mapping = aes(x = init_SR, ymin = lwr, ymax = upr),
              alpha = 0.05, fill = "#118176") +
  geom_line(data = pred_land,
            mapping = aes(x = init_SR, y = fit), colour = "#118176", size = 0.5) +
  xlab("Initial species richness") +
  ylab("Community biomass") +
  theme_meta() +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
}

# combine the plots for the defence and export
if(defence){
  yd <- cowplot::plot_grid(yd1,yd2, align = c("hv"))
  ggsave(filename = "figures-tables/def_fig_4.pdf", yd,
         unit = "cm", width = 18, height = 7.5, bg = "transparent")
}

# assign the plot object
if(case == "ac") {
  assign("c1", y)
} else {
  assign("d1", y)
}

# sort out the axis labels
a1 <- a1 + xlab(NULL)
b1 <- b1 + xlab(NULL) + ylab(NULL)

# sort out the axis labels
c1 <- c1 
d1 <- d1 + ylab(NULL)

# arrange the plot
pq <-
  cowplot::plot_grid(a1, b1, c1, d1, align = c("hv"),
                     nrow = 2, ncol = 2, 
                     labels = c("a", "b", "c", "d"), label_size = 11,
                     label_fontface = "plain",
                     rel_widths = c(1,1.05),
                     rel_heights = c(1.05, 1))
plot(pq)

# export the figure for further modification
ggsave(filename = "figures-tables/fig_9.pdf", pq,
       unit = "cm", width = 20, height = 14)

### END
