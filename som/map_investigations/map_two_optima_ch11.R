rm(list = ls())
library(Racmacs)
source("functions/map_functions.R")

old_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_alpha_adj.ace")
no_ch11 <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")
lims <- Racmacs:::mapPlotLims(old_map, sera = TRUE)
xlim_zoom <- lims$xlim +1
ylim_zoom <- lims$ylim + 1


png("som/map_investigations/ch11_optima_effect.png", 12, height =3, units = 'in', res=300, pointsize = 12)

layout(matrix(c(1:3), ncol = 3, byrow = F))
par(mar = rep(0.2, 4))

  plot(old_map, plot_stress = TRUE, optimization_number = 1, xlim = xlim_zoom, ylim = ylim_zoom)
  text(xlim_zoom[1] + 0.5, ylim_zoom[2]-0.4, "A", cex = 2, pos = 4)
  plot(procrustesMap(no_ch11, old_map, sera = FALSE), plot_stress = TRUE, optimization_number = 1, xlim = xlim_zoom, ylim = ylim_zoom)
  text(xlim_zoom[1] + 0.5, ylim_zoom[2]-0.4, "B", cex = 2, pos = 4)
  plot(procrustesMap(no_ch11, old_map, antigens = FALSE), plot_stress = TRUE, optimization_number = 1, xlim = xlim_zoom, ylim = ylim_zoom)
  text(xlim_zoom[1] + 0.5, ylim_zoom[2]-0.4, "C", cex = 2, pos = 4)
dev.off()

# do old map in 3d to see
map_3d <- optimizeMap(old_map, 3, 1000, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
agSize(map_3d) <- agSize(map_3d) - 15
srSize(map_3d) <- srSize(map_3d) - 10
Racmacs::view(map_3d)



