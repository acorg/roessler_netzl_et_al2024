library(Racmacs)
set.seed(100)

source("./functions/titer_lineplot_functions.R")
source("./functions/map_functions.R")

figure_dir <- file.path("som", "bootstrapping", "woXBBconvBQ11conv")

neut <- read.acmap('./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

target_srgs <- c('Wuhan conv.', 'alpha conv.',
                 'beta conv.',
                 'gamma conv.', 'delta conv.', 'BA.1 conv.', 'BA.2.12.1 conv.', 'BA.4 conv.', 'BA.5 conv.',
                 'XBB.1.5 conv.',
                 'Wuhan vax. (single dose)','Beta vax. (single dose)', 'Wuhan vax. (two doses)', 'XBB.1.5 vax. (two doses)')

labels <- data.frame(
  row.names = target_srgs,
  val = LETTERS[1:length(target_srgs)]
)


png(file.path(figure_dir, "bootstrapping-sera.png"), width = 12, height = 12, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:length(target_srgs)), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(srGroup in target_srgs){
    
    print(srGroup)
  
  newMap <- removeSera(neut, srNames(neut)[as.character(srGroups(neut)) == srGroup])
  newMap <- optimize_and_realign_map(newMap,
                                     neut, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                       dim_annealing = TRUE))
  
  
  newMap <- realignMap(newMap, neut)
  
  save_text <-srGroup
  
  save.acmap(map = newMap, filename = file.path(figure_dir, paste0("wo_",save_text,".ace")))
  #srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- srGroup
 
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[srGroup, ], cex = 1.2)
}

dev.off()


# when maps already exist
png(file.path(figure_dir, "bootstrapping-sera.png"), width = 12, height = 12, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:16), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(srGroup in target_srgs){
  print(srGroup)
  
  save_text <- srGroup
  newMap <- read.acmap(file.path(figure_dir, paste0("wo_",save_text,".ace")))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- capitalize(gsub(" conv.| vax.", "", srGroup))
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[srGroup, ], cex = 1.2)
}

dev.off()

