library(Racmacs)
library(dplyr)
library(ggplot2)

xlim_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x+2
xlim_zoom[2] <- xlim_zoom[2] + 2
ylim_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

reduce_point_size <- function(map){
  srSize(map) <- srSize(map) - 5
  agSize(map) <- agSize(map) - 8
  
  return(map)
}

# do automated reactivity adjustment, map has the lowest stress
og_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace")


ag_names <- agNames(og_map)
start_react <- agReactivityAdjustments(og_map)

adj_maps_alpha <- list()

start_react[ag_names %in% c("JN.1", "BA.2.86")] <- NA
adj_map_optim <- optimizeAgReactivity(og_map, fixed_ag_reactivities = start_react,
                                      reoptimize = TRUE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
adj_maps_alpha[["both optim"]] <- adj_map_optim

react_adj <- c(seq(-0.65, -0.8, by = -0.01), c(-0.25, -0.5, -0.685, -0.75, -1, -1.5, -2))


for(adj in react_adj){
  start_react[ag_names %in% c("JN.1", "BA.2.86")] <- adj
  adj_map_optim <- optimizeAgReactivity(og_map, fixed_ag_reactivities = start_react, reoptimize = TRUE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))

  adj_maps_alpha[[paste(adj)]] <- adj_map_optim
}


# plot here map stresses
map_stresses <- do.call(rbind, lapply(adj_maps_alpha, function(x) mapStress(x)))

reactivities <- rownames(map_stresses)
mapdf <- data.frame(map_stresses)
mapdf$reactivities <- reactivities

mapdf <- rbind(mapdf, data.frame("map_stresses" =  map_stresses_1, "reactivities" = rownames(map_stresses_1)))

mapdf %>%
  filter(reactivities != "both optim") %>%
  ggplot(aes(x = as.numeric(reactivities), y = map_stresses)) +
  geom_hline(yintercept = mapStress(og_map), color = "grey40", linetype = "dashed") +
  geom_hline(yintercept = min(map_stresses), linetype = "dashed", color = "red") +
  ylab("Map Stress") +
  xlab("Reactivity adjustmens JN.1, BA.2.86") +
  geom_point(shape = 21) +
  theme_bw() -> p

ggsave("som/map_investigations/map_stress_scan.png", p, dpi = 300, width = 3, height = 4)


png("som/map_investigations/map_lower_JN1BA286AlphaAdj_scan_optim.png", width = 9, height = 6, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:4), ncol=2, byrow = FALSE))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.3, 0.1, 0.3))

plot(reduce_point_size(adj_maps_alpha$`both optim`), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.5, ylim_zoom[2]-0.5, "A", cex = 1.2)
text(xlim_zoom[1]+5, ylim_zoom[2]-0.5, "Reactivity optimised", cex = 0.8)
plot(reduce_point_size(procrustesMap(adj_maps_alpha$`both optim`, og_map, sera = FALSE)), xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.5, ylim_zoom[2]-0.5, "C", cex = 1.2)

plot(reduce_point_size(adj_maps_alpha$`-0.685`), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.5, ylim_zoom[2]-0.5, "B", cex = 1.2)
text(xlim_zoom[1]+5, ylim_zoom[2]-0.5, "Reactivity -0.685", cex = 0.8)
plot(reduce_point_size(procrustesMap(adj_maps_alpha2$`-0.685`, og_map, sera = FALSE)), xlim = xlim_zoom, ylim = ylim_zoom)
text(xlim_zoom[1]+0.5, ylim_zoom[2]-0.5, "D", cex = 1.2)

dev.off()

save.acmap(adj_maps_alpha$`-0.685`, "data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alphaJN1ba286_adjScan.ace")
