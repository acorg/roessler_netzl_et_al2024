rm(list = ls())
library(Racmacs)
library(dplyr)
library(ggplot2)

# load functions
source("./functions/map_functions.R")


xlim_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x + 1
ylim_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x + 1

reduce_point_size <- function(map){
  srSize(map) <- srSize(map) - 5
  agSize(map) <- agSize(map) - 8
  
  return(map)
}

# do automated reactivity adjustment, map has the lowest stress
comp_map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace")
map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv.ace")


# remove XBB.1.5 conv from map
map <- removeSera(map, srNames(map)[srGroups(map) == "XBB.1.5 conv."])
map <- optimize_and_realign_map(map, comp_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                                dim_annealing = TRUE))

# adjust alpha reactivity
ag_names <- agNames(map)

start_react <- rep(0, length(ag_names))
start_react[ag_names == "Alpha"] <- NA

map_adj <- optimizeAgReactivity(map, fixed_ag_reactivities = start_react, reoptimize = FALSE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
# print(paste(map_f, agReactivityAdjustments(map)))
map_adj <- optimize_and_realign_map(map_adj, comp_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                                dim_annealing = TRUE))


# stable optimum with JN.1 and BA.2.86
plot(map_adj, plot_stress = TRUE)
plot(map_adj, plot_stress = TRUE, optimization_number = 999)
Racmacs::view(remove_na_coords(map_adj))

# Optimum of map with XBB.1.5 conv is also stable, but we do find other optimum. Ca 10 stress difference
plot(comp_map, plot_stress = TRUE, optimization_number = 219)
plot(comp_map, plot_stress = TRUE, optimization_number = 1)


remove_ag_and_optimize <- function(map, ag){
  
  # now make map without JN.1 and BA.2.86
  wo_JN.1 <- removeAntigens(map, ag)
  ag_names <- agNames(wo_JN.1)
  start_react <- rep(0, length(ag_names))
  start_react[ag_names == "Alpha"] <- NA
  
  wo_JN1 <- optimizeAgReactivity(wo_JN.1, fixed_ag_reactivities = start_react, reoptimize = FALSE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
  # print(paste(map_f, agReactivityAdjustments(map)))
  wo_JN1 <- optimize_and_realign_map(wo_JN1, map, 1000, 2, option_list = list(ignore_disconnected = TRUE,dim_annealing = TRUE))
 
  return(wo_JN1) 
}

wo_JN1 <- remove_ag_and_optimize(map, "JN.1")
wo_ba286 <- remove_ag_and_optimize(map, "BA.2.86")
wo_both <- remove_ag_and_optimize(map, c("JN.1", "BA.2.86"))

w_XBB_wo_ag <- remove_ag_and_optimize(comp_map, c("JN.1", "BA.2.86"))
plot(wo_JN1)
plot(wo_ba286)
plot(wo_both, plot_stress = TRUE)
plot(w_XBB_wo_ag, plot_stress = TRUE)
plot(w_XBB_wo_ag, plot_stress = TRUE, optimization_number = 400)


# what to show
# w XBB.1.5, w XBB.1.5 other optimum, w XBB.1.5 wo Ags, other optimum
# wo XBB.1.5, wo XBB.1.5 wo JN.1, wo XBB.1.5 wo BA.2.86, wo JN.1 both optima
png("som/map_investigations/map_wo_XBB15_optima.png", width = 14, height = 6, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:8), ncol=4, byrow = TRUE))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.3, 0.1, 0.3))

plot(reduce_point_size(comp_map), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
plot(reduce_point_size(comp_map), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 219)
plot(reduce_point_size(w_XBB_wo_ag), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 1)
plot(reduce_point_size(w_XBB_wo_ag), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 400)

plot(reduce_point_size(map_adj), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom)
plot(reduce_point_size(wo_JN1), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 1)
plot(reduce_point_size(wo_ba286), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 1)
plot(reduce_point_size(wo_both), plot_stress = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, optimization_number = 1)

dev.off()



save.acmap(wo_both, "data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj_woJN1BA286.ace")


mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 

wo_both <- apply_color(wo_both, mapColors)

Racmacs::view(apply_color(wo_both, mapColors))


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
