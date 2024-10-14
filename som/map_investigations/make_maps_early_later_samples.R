# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")
source("./functions/long_map_info.R")



map_dir <- file.path("data", "maps")
figure_dir <- file.path("som", "map_investigations")
# load data
alignment_names <- read.csv("./data/metadata/map-names.csv", row.names = "Antigen")
alignment_map <- read.acmap("data/maps/alignment_map.ace")


xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x +1
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x


agNames(alignment_map) <- alignment_names[agNames(alignment_map),]

# optimise reactivity for JN.1 and Alpha. reduce alpha and JN.1 by 2 fold
map_files <- list.files(map_dir, pattern = ".ace", full.names = FALSE)
map_files <- map_files[grepl("_woXBB15conv_CH11_alpha_adj", map_files)]


# make map with newer studies
data_long <- read.csv("data/titer_data/titer_data_long.csv")

# only keep samples from later time points
days_to_remove <- c("D7", "D9", "D10")
for(map_f in map_files){
  
  og_map <- read.acmap(file.path(map_dir, map_f))
  
  sera_to_remove <- srNames(og_map)[grepl(paste(days_to_remove, collapse = "|"), srNames(og_map))]
  
  map <- removeSera(og_map, sera_to_remove)
  map <- optimize_and_realign_map(map, og_map, 1000, 2, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
  
  png(file.path(figure_dir, gsub(".ace", "_later_days.png", map_f)), width = 10, height = 3, units = 'in', res=300, pointsize = 18)
  layout(matrix(c(1:3), ncol = 3, byrow = T))
  par(oma=c(0, 0, 0, 0),  mar=c(0, 0.25, 0, 0.25))
  plot(map, show_error = FALSE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
#   text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
  plot(triangulationBlobs(relaxMap(remove_na_coords(map)), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
 #  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
  plot(procrustesMap(map, og_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
#  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "C", cex = 1.2)
  
  dev.off()
  
  
}


# only keep samples from early days
for(map_f in map_files){
  
  og_map <- read.acmap(file.path(map_dir, map_f))
  
  sera_to_remove <- srNames(og_map)[!grepl(paste(days_to_remove, collapse = "|"), srNames(og_map))]
  
  map <- removeSera(og_map, sera_to_remove)
  map <- optimize_and_realign_map(map, og_map, 1000, 2, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
  
  png(file.path(figure_dir, gsub(".ace", "_earlier_days.png", map_f)), width = 10, height = 3, units = 'in', res=300, pointsize = 18)
  layout(matrix(c(1:3), ncol = 3, byrow = T))
  par(oma=c(0, 0, 0, 0),  mar=c(0, 0.25, 0, 0.25))
  plot(map, show_error = FALSE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
 # text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
  plot(triangulationBlobs(relaxMap(remove_na_coords(map)), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
#  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
  plot(procrustesMap(map, og_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
 # text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "C", cex = 1.2)
  
  dev.off()
  
  
}
