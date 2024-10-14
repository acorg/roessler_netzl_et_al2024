# labelled map 
rm(list = ls())
library(Racmacs)
source("functions/map_functions.R")

fig_dir_color <- file.path("figures", "labelled_map", "lighterXBB")
dir.create(fig_dir_color, showWarnings = FALSE)

move_coords <- function(map, at = 2, by = -0.5){
  agCoords(map)[,at] <- agCoords(map)[,at] + by
  srCoords(map)[,at] <- srCoords(map)[,at] + by

  return(map)
}

# move all map y-coords down by 0.5
# read in map
# map <- move_coords(read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace"), at = 1, by = 0.5)
# map_alpha_unadj <- move_coords(read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv.ace"), at = 1, by = 0.5)
# map_alpha_unadj <- move_coords(map_alpha_unadj, at = 2, by = -1)


map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")
map_alpha_unadj <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11.ace")

lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
lims_no_zoom <- Racmacs:::mapPlotLims(map, sera = TRUE)

xlim_zoom <- round(lims$xlim)
ylim_zoom <- round(lims$ylim)-1

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

write.csv(ylim_no_zoom, "./data/metadata/ylim_no_zoom.csv")
write.csv(xlim_no_zoom, "./data/metadata/xlim_no_zoom.csv")

write.csv(ylim_zoom, "./data/metadata/ylim_zoom.csv")
write.csv(xlim_zoom, "./data/metadata/xlim_zoom.csv")


map <- move_coords(map, at = 2, by = -0.3)
map <- move_coords(map, at = 1, by = 0.2)

map_alpha_unadj <- move_coords(map_alpha_unadj, at = 2, by = -0.3)
map_alpha_unadj <- move_coords(map_alpha_unadj, at = 1, by = 0.2)


# Setup plotting function
doplot <- function(map, xlims, ylims, show_labels = TRUE) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  srSize(map) <- srSize(map) - 4
  agSize(map) <- agSize(map) - 4
  # Plot the regular map
  srOutlineWidth(map) <- 0.8
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.6)
  plot(map, xlim = xlims, 
       ylim =ylims, fill.alpha = 0.9,
       plot_stress = TRUE)
  
  if(show_labels){
    
    # Plot labels
    label_adjustments <- matrix(0, numAntigens(map), 2)
    rownames(label_adjustments) <- agNames(map)
    label_adjustments["B.1.351",] <- c(0.9, 0)
    label_adjustments["P.1.1",] <- c(-0.9, 0)
    label_adjustments["B.1.1.7+E484K",] <- c(0.9, -0.5)
    label_adjustments["BA.1",] <- c(-0.4, 0.7)
    label_adjustments["BA.2",] <- c(0, -0.6)
    label_adjustments["B.1.1.7",] <- c(0.3, -0.6)
    label_adjustments["D614G",] <- c(0, -0.5)
    label_adjustments["B.1.617.2",] <- c(0,-0.6)
    label_adjustments["BA.5.3.2",] <- c(0, 0.7)
    
    labels <- agNames(map)
    names(labels) <- agNames(map)
    labels["B.1.351"] <- "beta\n(B.1.351)"
    labels["P.1.1"] <- "gamma\n(P.1.1)"
    labels["BA.1"] <- "BA.1 omicron\n(B.1.1.529+BA.1)"
    labels["BA.2"] <- "BA.2 omicron\n(B.1.1.529+BA.2)"
    labels["B.1.617.2"] <- "delta\n(B.1.617.2)"
    labels["B.1.1.7"] <- "alpha\n(B.1.1.7)"
    labels["B.1.1.7+E484K"] <- "alpha + E484K\n(B.1.1.7+E484K)"
    labels["BA.5.3.2"] <- "BA.5.3.2 omicron\n(B.1.1.529+BA.5)"
    
    label_size <- rep(1, numAntigens(map))
    names(label_size) <- agNames(map)
    
    text(
      agCoords(map) + label_adjustments,
      cex = label_size,
      label = labels,
      font = 1
    )
  }
  

}

png("figures/labelled_map/map_no_zoom.png", 7, 7, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_no_zoom, ylim_no_zoom, FALSE)
dev.off()

png("figures/labelled_map/map_zoom.png", 4, 5, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_zoom, ylim_zoom, FALSE)
dev.off()


png("figures/labelled_map/alpha_unadj_proc_map_zoom.png", 7.5, height = 4, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1:2), ncol = 2, byrow = T))
par(mar = rep(0.5, 4))
doplot(map_alpha_unadj, xlim_zoom, ylim_zoom, FALSE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "A", cex = 1.2)
doplot(procrustesMap(map, map_alpha_unadj,sera = FALSE), xlim_zoom, ylim_zoom, FALSE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "B", cex = 1.2)
dev.off()


