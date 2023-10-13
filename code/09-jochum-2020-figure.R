#'
#' @title Re-plot Jochum et al.'s (2020) data
#' 
#' @description Replot Jochum et al.'s (2020) of BEF relationships
#' in realistic community compositions
#'

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# get the meta_theme()
source("code/helper-plotting-theme.R")

# list the files in this directory
file_list <- list.files("data")
print(file_list)

# extract the jochum data
file_list <- file_list[grepl(pattern = "jochum", file_list)]
print(file_list)

# read in the files into a list
df_list <- vector("list", length = length(file_list))
for (i in 1:length(file_list)) {
  
  # read the file and extract the names
  df <- read_csv(file = paste0("data/", file_list[i]))
  df_names <- names(df)
  
  # use the name to specify the function id directly
  df$function_id <- df_names[3]
  
  # rename the columns
  names(df) <- c("con_uncon", "log_diversity", "function_value", "function_id")
  
  # set the slope identity based on the dataset
  df$slope <- gsub(".csv", "\\1", file_list[i])
  
  # add a column specifying the country
  df[["country"]] <- 
    if(grepl(pattern = "USA", file_list[i])) {
      "USA"
    } else {
      "DEU"
    }
  
  # reorder the columns
  df <-
    df |>
    dplyr::select(country, slope, con_uncon, log_diversity, function_id, function_value)
  
  # write into a list
  df_list[[i]] <- df
  
}

# bind this list into a data.frame
df_joch <- bind_rows(df_list)

# check the unique functions
unique(df_joch[["function_id"]])

# extract aboveground biomass, root biomass and sol carbon
funcs <- c("sqrt_biomass", "soil_organic_C", "log_root_biomass",
           "log_biomass", "log_soil_C", "root_biomass")

# filter the dataset
df_joch <-
  df_joch |>
  dplyr::filter(function_id %in% funcs)

# plot the germany sites
df_DEU <- 
  df_joch |> 
  dplyr::filter(country == "DEU")

# function ids
funcs <- c("sqrt_biomass", "log_root_biomass", "soil_organic_C")

# y-labs
ylabs <- c(expression(sqrt(Biomass~(g~m^{-2}))),
           expression(ln(Root~biomass~(g~m^{-2}))),
           "Soil organic C (%)")

# x-labs
xlabs <- c("", "ln( Species richness )", "")

# make the unconstrained plots
uncon_plots <- vector("list", length = length(funcs))

# set the plot titles
titles <- c("Unconstrained", "", "")

# loop over each variable
for(i in 1:length(funcs)) {
  
  # filter the data
  df_plot <-
    df_DEU |> 
    dplyr::filter(con_uncon == "uncon") |> 
    dplyr::filter(function_id == funcs[i])
  
  # fit the model
  lm1 <- lm(function_value ~ log_diversity, data = df_plot)
  
  # plot the data
  p1 <- 
    ggplot(data = df_plot,
         mapping = aes(x = log_diversity, y = function_value)) +
    geom_point(size = 1.75, colour = "#C9795B") +
    geom_smooth(method = "lm", colour = "#C9795B", fill = "#C9795B", alpha = 0.05,
                size = 0.5) +
    ggtitle(titles[i]) +
    ylab(ylabs[i]) +
    xlab(xlabs[i]) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.4,
             label = lm_eqn(lm1), parse = TRUE) +
    theme_meta()
  
  uncon_plots[[i]] <- p1
  
}

# merge the plots
q1 <-
  cowplot::plot_grid(plotlist = uncon_plots, nrow = 1, ncol = 3,
                   labels = c("a", "b", "c"), label_size = 11,
                   label_fontface = "plain",
                   vjust = 4)

# make the constrained plots
con_plots <- vector("list", length = length(funcs))

# set the plot titles
titles <- c("Constrained", "", "")

# loop over each variable
for(i in 1:length(funcs)) {
  
  # filter the data
  df_plot <-
    df_DEU |> 
    dplyr::filter(con_uncon == "con") |> 
    dplyr::filter(function_id == funcs[i])
  
  # fit the model
  lm1 <- lm(function_value ~ log_diversity, data = df_plot)
  
  # plot the data
  p1 <- 
    ggplot(data = df_plot,
           mapping = aes(x = log_diversity, y = function_value)) +
    geom_point(size = 1.75, colour = "#5685C1") +
    geom_smooth(method = "lm", colour = "#5685C1", fill = "#5685C1", alpha = 0.05,
                size = 0.5) +
    ylab(ylabs[i]) +
    xlab(xlabs[i]) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.4,
             label = lm_eqn(lm1), parse = TRUE) +
    ggtitle(titles[i]) +
    theme_meta()
  
  con_plots[[i]] <- p1
  
}

# merge the plots
q2 <- 
  cowplot::plot_grid(plotlist = con_plots, nrow = 1, ncol = 3,
                   labels = c("d", "e", "f"), label_size = 11,
                   label_fontface = "plain", 
                   vjust = 4)

# combine the two sets of plots
q12 <- cowplot::plot_grid(q1, q2, nrow = 2, ncol = 1)

ggsave(filename = "figures-tables/fig_8.pdf", plot = q12,
       unit = "cm", width = 20, height = 15)

### END
