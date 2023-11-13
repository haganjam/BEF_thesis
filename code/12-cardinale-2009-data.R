#'
#' @title Qiu and Cardinale's (2021, Ecology) data
#' 
#' @description Plot the balance between selection and complementarity effects
#'
# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# load relevant scripts
source("code/helper-plotting-theme.R")

# download the raw data

# data come from Cardinale et al. (2009): Effects of biodiversity on the functioning of ecosystems: a summary of 164 experimental manipulations of species richness, Ecology
# link to the archive is: https://esapubs.org/archive/ecol/E090/060/

# load the raw data directly from the link
qiu <- readr::read_csv(url("https://esapubs.org/archive/ecol/E090/060/BEF_summary_v2_Aug2008.csv"),
                       na = c(".", "NA"),)

# check the parsing specifications
spec(qiu)

# view the data
View(qiu)

# subset out the relevant data points
names(qiu)

# subset the selection and complementarity effects
qiu_part <- 
  qiu |> 
  dplyr::select(Entry, Study, Expt, Consumer, SE, SDSE, CE, SDCE)

# get non-nas
qiu_part <- qiu_part[complete.cases(qiu_part[,c(5, 7)]),]

# remove extreme protist example
qiu_part <- 
  qiu_part |>
  dplyr::filter(SE < 1000, SE > -1000 )|>
  dplyr::filter(CE < 1000, CE > -1000 )

# plot the results
summary(qiu_part)

# get the mean and sd
msd <- 
  qiu_part |>
  dplyr::summarise(SE_m = mean(SE),
                   SE_sd = sd(SE),
                   CE_m = mean(CE),
                   CE_sd = sd(CE))

p1 <- 
  ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = qiu_part,
             mapping = aes(x = SE, y = CE), 
             shape = 1, alpha = 0.3) +
  geom_errorbar(data = qiu_part, 
                mapping = aes(x = SE, ymin = CE-SDCE, ymax = CE+SDCE),
                size = 0.2, alpha = 0.3) +
  geom_errorbarh(data = qiu_part, 
                 mapping = aes(y = CE,xmin = SE-SDSE, xmax = SE+SDSE),
                 size = 0.2, alpha = 0.3) +
  geom_errorbar(data = msd,
                mapping = aes(x = SE_m, ymin = CE_m-CE_sd, ymax = CE_m+CE_sd),
                size = 0.5, alpha = 1, colour = "red") +
  geom_errorbarh(data = msd,
                 mapping = aes(y = CE_m, xmin = SE_m-SE_sd, xmax = SE_m+SE_sd),
                 size = 0.5, alpha = 1, colour = "red") +
  geom_point(data = msd,
             mapping = aes(x = SE_m, y = CE_m), colour = "red", shape = 23, size = 3.5,
             fill = "white", stroke = 0.75) +
  scale_x_continuous(limits = c(-500, 900), breaks = seq(-500, 900, 250)) +
  scale_y_continuous(limits = c(-500, 900), breaks = seq(-500, 900, 250)) +
  ylab("Complementarity Effect") +
  xlab("Selection Effect") +
  theme_meta()
plot(p1)

# how many estimates?
nrow(qiu_part)
length(unique(qiu_part$Study))

# export the figure for further modification
ggsave(filename = "figures-tables/def_fig_3.pdf", p1,
       unit = "cm", width = 10, height = 8)




