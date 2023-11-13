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

# plot the results
ggplot(data = qiu_part)




