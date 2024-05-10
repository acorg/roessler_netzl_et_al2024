library(Racmacs)
set.seed(100)

source("./functions/map_functions.R")

figure_dir <- file.path("som", "bootstrapping", "woXBBconvBQ11conv")

neut <- read.acmap('./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x 

labels <- data.frame(
  row.names = agNames(neut),
  val = LETTERS[c(1:length(agNames(neut)))]
)

png(file.path(figure_dir, "bootstrapping-antigens.png"), width = 10, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:nrow(labels)), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(ag in agNames(neut)){
  print(ag)
  
  newMap <- removeAntigens(neut, c(ag))
  newMap <- optimize_and_realign_map(newMap,
                                     neut, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                         dim_annealing = TRUE))
  
  newMap <- realignMap(newMap, neut)
  
  save.acmap(map = newMap, filename = file.path(figure_dir, paste0("wo_ag_",ag,".ace")))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = ag, cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[ag, ], cex = 1.2)
}

dev.off()

png(file.path(figure_dir, "bootstrapping-antigens.png"), width = 12, height = 14, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:nrow(labels)), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(ag in agNames(neut)){

  newMap <- read.acmap(filename = file.path(figure_dir, paste0("wo_ag_",ag,".ace")))
  srOutlineWidth(newMap) <- 1

  p <- procrustesMap(newMap, neut, sera = FALSE)

  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = ag, cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.6, ylim_no_zoom[2]-0.6, labels[ag, ], cex = 1.2)
}

dev.off()

