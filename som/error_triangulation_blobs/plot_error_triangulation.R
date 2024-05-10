rm(list = ls())
library(Racmacs)
source("functions/map_functions.R")
set.seed(100)

figure_dir <- "som/error_triangulation_blobs/"
map_dir <- "./data/maps/"

map_files <- list.files(map_dir, pattern = ".ace", full.names = FALSE)
map_files <- map_files[grepl("singleTP", map_files)]

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x


for(map_f in map_files){
  
  map <- read.acmap(file.path(map_dir, map_f))
  
  map_positioned <- remove_na_sera(map)
  
  fig_name <- gsub(".ace", ".png", map_f)
  
  png(file.path(figure_dir, paste0("error_triangulation_", fig_name)), width = 6.5, height = 3, units = 'in', res=300, pointsize = 18)
  layout(matrix(c(1:2), ncol = 2, byrow = T))
  par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))
  plot(map_positioned, show_error = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
  plot(triangulationBlobs(relaxMap(map_positioned), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
  
  dev.off()
  
}