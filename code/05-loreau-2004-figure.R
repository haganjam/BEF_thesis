#'
#' @title Theoretical results from Loreau (2004)
#' 
#' @description Uses the Lotka-Volterra model from Loreau (2004) to reproduce
#' the basic theoretical results consistently reported in the BEF literature
#'

# load relevant scripts
source("code/helper-plotting-theme.R")

# load relevant libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# simulate parameter values of K1, K2, a12 and a21 compatible with coexistence
n <- 1000
coex <- vector("list", length = n)
for(i in 1:n) {
  
  K1 <- runif(n = 1, min = 0.5, max = 3)
  K2 <- runif(n = 1, min = 0.5, max = 3)
  a12 <- runif(n = 1, min = 0.25, max = 1.25)
  a21 <- runif(n = 1, min = 0.25, max = 1.25)
  
  cond1 <- (a21 < (K2/K1))
  cond2 <- ((K2/K1) < (1/a12))
  
  if( all(c(cond1, cond2)) ) {
    
    coex[[i]] <- data.frame(K1 = K1,
                            K2 = K2,
                            a12 = a12,
                            a21 = a21,
                            coex = TRUE)
    
  } else {
    
    coex[[i]] <- data.frame(K1 = K1,
                            K2 = K2,
                            a12 = a12,
                            a21 = a21,
                            coex = FALSE)
    
  }
  
}

# calculate equilibrium abundances
BEF_df <-
  
  lapply(coex, function(y) {
  
  r1 <- y
  
  if(r1$coex == TRUE) {
    
    # calculate the x value
    r1$x <- r1$K2/r1$K1
    
    N1 <- with(data = r1,
               (K1*(1-(x*a12)))/(1-(a12*a21)) )
    
    N2 <- with(data = r1,
               (K1*(x-a21))/(1-(a12*a21)) )
    
    # pull into a data.frame from a classic BEF experiment
    df <- data.frame(coex = TRUE,
                     SR = c(1, 1, 2),
                     EF = c(r1$K1, r1$K2, (N1+N2)),
                     TOY = (N1+N2) > max(c(r1$K1, r1$K2)),
                     RYT = ((N1/r1$K1) + (N2/r1$K2)) )
    
  } else {
    
    # pull into a data.frame from a classic BEF experiment
    df <- data.frame(coex = FALSE,
                     SR = c(1, 1, 2),
                     EF = c(r1$K1, r1$K2, max(c(r1$K1, r1$K2))),
                     TOY = FALSE,
                     RYT = 1)
    
  }
  
  return(df)
  
} )

# bind this into a data.frame
BEF_df <- dplyr::bind_rows(BEF_df, .id = "ID")
head(BEF_df)

BEF_df <- 
  BEF_df |>
  dplyr::rename(`Transgressive overyielding` = TOY)

# check for any negative relationships
ID_in <- 
  BEF_df |>
  dplyr::group_by(ID, SR) |>
  dplyr::summarise(mean_EF = mean(EF)) |>
  dplyr::summarise(diff_EF = diff(mean_EF)) |>
  dplyr::filter(diff_EF < 0) |>
  dplyr::pull(ID)

# make a plot of this result for stable coexistence
df_plot1 <- 
  BEF_df |>
  dplyr::filter(coex == TRUE) |>
  dplyr::filter(ID %in% sample(unique(ID), 30))

p1 <- 
  ggplot() +
  geom_jitter(data = df_plot1,
              mapping = aes(x = SR, y = EF, group = ID, colour = `Transgressive overyielding`), 
              width = 0.01, alpha = 9, shape = 16, size = 1.5) +
  geom_smooth(data = df_plot1,
              mapping = aes(x = SR, y = EF, group = ID, colour = `Transgressive overyielding`),
              method = "lm", se = FALSE, size = 0.25) +
  ylab("Abundance") +
  xlab("Species richness") +
  theme_meta() +
  scale_x_continuous(breaks = c(1, 2), limits = c(0.9, 2.1)) +
  scale_y_continuous(limits = c(0.3, 4.1)) +
  scale_colour_manual(values = col_func()[4:5]) +
  theme(legend.position = "top",
        legend.key=element_blank())
plot(p1)

# make a plot of this result for unstable coexistence
df_plot2 <- 
  BEF_df |>
  dplyr::filter(coex == FALSE) |>
  dplyr::filter(ID %in% sample(unique(ID), 30))

p2 <- 
  ggplot() +
  geom_jitter(data = df_plot2,
              mapping = aes(x = SR, y = EF, group = ID, colour = `Transgressive overyielding`), 
              width = 0.01, alpha = 9, shape = 16, size = 1.5) +
  geom_smooth(data = df_plot2,
              mapping = aes(x = SR, y = EF, group = ID, colour = `Transgressive overyielding`),
              method = "lm", se = FALSE, size = 0.25) +
  guides(color=guide_legend(override.aes=list(colour=NA))) +
  ylab("Abundance") +
  xlab("Species richness") +
  theme_meta() +
  scale_x_continuous(breaks = c(1, 2), limits = c(0.9, 2.1)) +
  scale_y_continuous(limits = c(0.3, 4.1)) +
  scale_colour_manual(values = col_func()[4:5]) +
  theme(legend.position = "top") +
  theme(legend.position = "top",
        legend.key=element_blank())
plot(p2)

p12 <- ggarrange(p1, p2, common.legend = TRUE, labels = c("a", "b"),
                 font.label = list(size = 11, face = "plain"))

ggsave(filename = "figures-tables/fig_7.pdf", plot = p12,
       unit = "cm", width = 19, height = 10.5)

### END
