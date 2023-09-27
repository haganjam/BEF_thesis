
# customised plotting theme
theme_meta <- 
  function(base_family = "") {
    theme(panel.background=element_rect(fill="white", colour="black", linetype="solid"),
          panel.grid = element_blank(),
          axis.text = element_text(colour="black",size=10),
          legend.key = element_rect(fill = NA))
  }

# get the equation from the model object (m) i.e. lm1
# function to do this
lm_eqn <- function(m){
  eq <- substitute(italic(r)^2~"="~r2,
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq))
}

### END