# load required packages
rm(list = ls())
library(Racmacs)
library(dplyr)
library(ggplot2)

# load functions
source("./functions/map_functions.R")


xlim_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x + 1
ylim_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x + 1


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
comp_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace")
comp_wo <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woJN1BA286_alpha_adj.ace")
map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_alpha_adj.ace")
final_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_woJN1BA286_alpha_adj.ace")


plot(comp_map, plot_stress = TRUE, optimization_number = 189)
plot(map, plot_stress = TRUE, optimization_number = 999)

wo_JN1 <- remove_ag_and_optimize(map, "JN.1")
wo_BA286 <- remove_ag_and_optimize(map, "BA.2.86")


png("som/map_investigations/map_wo_XBB15_optima.png", width = 20, height = 9, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:8), ncol=4, byrow = TRUE))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.3, 0.1, 0.3))

plot(final_map, plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "A", cex = 2)
plot(procrustesMap(wo_JN1, final_map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "B", cex = 2)
plot(procrustesMap(wo_BA286, final_map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "C", cex = 2)
plot(procrustesMap(map, final_map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "D", cex = 2)

plot(procrustesMap(comp_wo, final_map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "E", cex = 2)
plot(procrustesMap(comp_wo, final_map, sera = FALSE, optimization_number = 623), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "F", cex = 2)
plot(procrustesMap(comp_map, final_map, sera = FALSE), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "G", cex = 2)
plot(procrustesMap(comp_map, final_map, sera = FALSE, optimization_number = 189), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "H", cex = 2)

dev.off()

