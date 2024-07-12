rm(list = ls())
library(Racmacs)
source("./functions/map_functions.R")
set.seed(100)


figure_dir <- "som/bootstrapping/"
map_dir <- "./data/maps/"

map_files <- c("map_threshold20_all_ags_singleTP_woXBBBQ11conv_alphaJN1ba286_adjScan.ace")

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x + 1
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x



for(map_f in map_files){
  
  print(map_f)
 
  neut <- read.acmap(file.path(map_dir, map_f))
  
  neut <- remove_na_coords(neut)
  
  neutBootTA <- bootstrapMap(
    neut,
    "bayesian",
    bootstrap_repeats = 500, # was 1000
    bootstrap_ags = TRUE,
    bootstrap_sr = TRUE,
    reoptimize = TRUE,
    optimizations_per_repeat = 1000,
    ag_noise_sd = 0.7,
    titer_noise_sd = 0.7,
    options = list(ignore_disconnected = TRUE,
                   dim_annealing = TRUE)
  )
  
  print("after TA")
  map_n <- gsub(".ace", "_neutBootTA_bayesian500_1000.ace", map_f)
  save.acmap(neutBootTA, file.path(figure_dir, map_n))
  
  neutBootTA <- bootstrapMap(
    neut,
    "noisy",
    bootstrap_repeats = 500,
    bootstrap_ags = TRUE,
    bootstrap_sr = TRUE,
    reoptimize = TRUE,
    optimizations_per_repeat = 1000,
    ag_noise_sd = 0.7,
    titer_noise_sd = 0.7,
    options = list(ignore_disconnected = TRUE,
                   dim_annealing = TRUE)
  )
  
  print("after noisy")
  map_n <- gsub(".ace", "_neutBootTA_noisy500_1000.ace", map_f)
  save.acmap(neutBootTA, file.path(figure_dir, map_n))
  
  neutBootAS <- bootstrapMap(
    neut,
    "resample",
    bootstrap_repeats = 500,
    bootstrap_ags = TRUE,
    bootstrap_sr = TRUE,
    reoptimize = TRUE,
    optimizations_per_repeat = 1000,
    ag_noise_sd = 0.7,
    titer_noise_sd = 0.7,
    options = list(ignore_disconnected = TRUE,
                   dim_annealing = TRUE)
  )
  
  map_n <- gsub(".ace", "_neutBootTA_resample500_1000.ace", map_f)
  save.acmap(neutBootAS, file.path(figure_dir, map_n))
  
  print("after resample")
  
  map_f <- gsub(".ace", "", map_f)
  
  neutBootTA <- read.acmap(paste0("./som/bootstrapping/", map_f, "_neutBootTA_bayesian500_1000.ace"))
  neutBootTABlobs <- bootstrapBlobs(neutBootTA, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
  
  neutBootResample <- read.acmap(paste0("./som/bootstrapping/", map_f, "_neutBootTA_resample500_1000.ace"))
  neutBootResampleBlobs <- bootstrapBlobs(neutBootResample, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
  
  neutBootNoisy <- read.acmap(paste0("./som/bootstrapping/", map_f, "_neutBootTA_noisy500_1000.ace"))
  neutBootNoisyBlobs <- bootstrapBlobs(neutBootNoisy, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
  
  # Plot the figure
  fig_name <- paste0(map_f, "_bbootstrap_TA.png")
  png(file.path(figure_dir, fig_name), width = 9, height = 2.5, units = 'in', res=300, pointsize = 18)
  layout(matrix(c(1, 2, 3), ncol=3))
  par(oma=c(0, 0, 0, 0), mar=c(0, 0.5, 0, 0.5))
  plot(neutBootTABlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
 # text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'A', cex = 1.4)
  plot(neutBootResampleBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
#  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'B', cex = 1.4)
  plot(neutBootNoisyBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
 # text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'C', cex = 1.4)
  dev.off()
  
}



