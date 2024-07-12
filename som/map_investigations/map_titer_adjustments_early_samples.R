# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)
library(dplyr)
# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")
figure_dir <- "som/map_investigations/"

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x +1
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x 

xlim_zoom <- read.csv("./data/metadata/xlim_zoom.csv")$x
ylim_zoom <- read.csv("./data/metadata/ylim_zoom.csv")$x
ylim_zoom[2] <- ylim_zoom[2] +1

mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 


# calculate avg fold change from D10 for samples with 2 time points
data_long <- read.csv("data/titer_data/titer_data_long.csv")

# remove non titrated samples
data_long %>%
  group_by(sr_info) %>%
  mutate(non_titrated = length(Titer_thresh20[Titer_thresh20 == "*"]),
         all_titers = length(Titer_thresh20)) %>%
  filter(non_titrated < all_titers) -> data_long

data_long %>%
  select(NHP.ID, Time.point.post.challenge, sr_info) %>%
  unique() %>%
  group_by(NHP.ID) %>%
  mutate(nr_time_points = length(Time.point.post.challenge)) %>%
  filter(nr_time_points > 1) %>%
  pull(sr_info) -> sera_two_timepoints

# add human gmt as highlighted sample

data_long %>%
  filter(sr_info %in% sera_two_timepoints) %>%
  mutate(sr_name = sr_info,
         titer = Titer_thresh20,
         titertype = ifelse(grepl("<", titer), 2, 1),
         logtiter = ifelse(titertype == 1, log2(as.numeric(titer)/10), log2(20/10)),
         nhp_sr_group = paste0(sr_group, ": ",NHP.ID)) -> data_long_sub


# calculate mean fc for >LOD only
data_long_sub %>%
  filter(titertype == 1) %>%
  group_by(NHP.ID, ag_name, sr_group) %>%
  summarize(fc = logtiter - logtiter[Time.point.post.challenge == "D10"]) %>%
  filter(fc > 0) %>%
  ungroup() %>%
  summarize(mean_fc = round(2^mean(fc))) %>%
  pull(mean_fc) -> d10


# adjust reactivity of >LOD for early samples accordingly
base_map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alphaJN1ba286_adjScan.ace")

base_table <- titerTable(base_map)
base_log <- logtiterTable(base_map)

adj_table_d10 <- base_table[,grepl("D10|D7|D9", colnames(base_table))]
adj_table_later <- base_table[,!grepl("D10|D7|D9", colnames(base_table))]
adj_table_d10[!(adj_table_d10 %in% c("*", "<20"))] <- as.character(as.numeric(adj_table_d10[!(adj_table_d10 %in% c("*", "<20"))])*d10) # FC from D10 between 4 and 8

adj_d10_map <- base_map
titerTable(adj_d10_map) <- cbind(adj_table_later, adj_table_d10)
adj_d10_map <- optimize_and_realign_map(adj_d10_map, base_map, 2000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                                   dim_annealing = TRUE))

png(file.path(figure_dir, "map_D10titersx4_sera.png"), width = 10, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:3), ncol = 3, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0, 0.25, 0, 0.25))

plot(adj_d10_map, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
#text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "A", cex = 1.2)
plot(triangulationBlobs(relaxMap(remove_na_coords(adj_d10_map)), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
plot(procrustesMap(adj_d10_map, base_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
#text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "B", cex = 1.2)

dev.off()
