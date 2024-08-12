# load required packages
rm(list = ls())
library(Racmacs)
library(dplyr)
library(ggplot2)

# load functions
source("./functions/map_functions.R")


xlim_zoom <- read.csv("./data/metadata/xlim_zoom.csv")$x 
ylim_zoom <- read.csv("./data/metadata/ylim_zoom.csv")$x 

# compare maps with and without XBB.1.5 conv sera
map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_alpha_adj.ace")
final_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_woJN1BA286_alpha_adj.ace")

srSize(map) <- srSize(map) - 4
agSize(map) <- agSize(map) - 4



ptDrawingOrder(map) <- rev(ptDrawingOrder(map))

png("som/map_investigations/map_wJN1BA286.png", width = 4, height = 5, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:1), ncol=1, byrow = TRUE))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.3, 0.1, 0.3))

plot(map, plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
dev.off()

