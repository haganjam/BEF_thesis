#'
#' @title Thompson et al. (2021) model
#' 
#' @description Use Thompson et al.'s (2021) to reproduce the theoretical 
#' results obtained by Bond and Chase (2002)
#'

# load relevant libraries
library(dplyr)
library(ggplot2)
library(mcomsimr)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# set the specific parameters and just vary dispersal

# simulate using the mcomsimr models
mc1 <- 
  simulate_MC(
  patches = 100,
  species = 50,
  dispersal = 0.01,
  plot = FALSE,
  torus = TRUE,
  kernel_exp = 0.1,
  env1Scale = 2,
  timesteps = 1200,
  burn_in = 800,
  initialization = 200,
  max_r = 5,
  min_env = 0,
  max_env = 1,
  env_niche_breadth = 0.3,
  optima_spacing = "even",
  intra = 1,
  min_inter = 1,
  max_inter = 1,
  comp_scaler = 0.05,
  extirp_prob = 0
)

# get the positive time-steps
raw_df <- 
  mc1$dynamics.df |>
  dplyr::filter(time > 0) |>
  dplyr::select(time, env, patch, species, optima, env_niche_breadth, max_r, N) |>
  dplyr::arrange(time, patch, species)

# summarise the BEF quantities: Patch scale
patch_df <- 
  raw_df |>
  dplyr::group_by(time, patch) |>
  dplyr::summarise(N_total = sum(N),
                   alpha_div = sum(N > 0)) |>
  dplyr::ungroup() |>
  dplyr::filter(time == max(time))
plot(patch_df$alpha_div, patch_df$N_total)

# summarise the BEF quantities: Landscape scale
land_df <- 
  raw_df |>
  dplyr::group_by(time) |>
  dplyr::filter(N > 0) |>
  dplyr::summarise(N_total = mean(N),
                   gamma_div = length(unique(species))) |>
  dplyr::ungroup() |>
  dplyr::filter(time == max(time))
print(land_df)  

