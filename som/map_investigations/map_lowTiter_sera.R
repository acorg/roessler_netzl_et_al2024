# load required packages
rm(list = ls())
library(Racmacs)
library(dplyr)
library(ggplot2)

# load functions
source("./functions/map_functions.R")


xlim_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x 
ylim_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x 


# function to quickly remove antigens
remove_ag_and_optimize <- function(map, ag){
  
  # now make map without JN.1 and BA.2.86
  wo_JN1 <- removeAntigens(map, ag)
  ag_names <- agNames(wo_JN1)
  start_react <- rep(0, length(ag_names))
  start_react[ag_names == "Alpha"] <- NA
  
  wo_JN1 <- optimizeAgReactivity(wo_JN1, fixed_ag_reactivities = start_react, reoptimize = FALSE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
  
  wo_JN1 <- optimize_and_realign_map(wo_JN1, map, 1000, 2, option_list = list(ignore_disconnected = TRUE,dim_annealing = TRUE))
  
  return(wo_JN1) 
}


# compare maps with and without XBB.1.5 conv sera
comp_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_CH11_alpha_adj.ace")
full_map <- read.acmap("data/maps/map_threshold20_all_ags.ace")
map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")


no_ch11 <- remove_ag_and_optimize(full_map, ag = "CH.1.1")


png("som/map_investigations/map_wo_lowTiter_sera.png", width = 10, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:4), ncol=2, byrow = TRUE))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.3, 0.1, 0.3))

plot(procrustesMap(comp_map, map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "A", cex = 1.5)
plot(procrustesMap(no_ch11, map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "B", cex = 1.5)

plot(procrustesMap(comp_map, map, antigens = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "C", cex = 1.5)
plot(procrustesMap(no_ch11, map, antigens = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "D", cex = 1.5)


dev.off()

