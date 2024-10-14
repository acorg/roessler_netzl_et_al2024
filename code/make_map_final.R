# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)
library(dplyr)
library(tidyverse)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

map_dir <- file.path("data", "maps")

do_thresh_1 <- FALSE
# load data
mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 

alignment_names <- read.csv("./data/metadata/map-names.csv", row.names = "Antigen")
alignment_map <- read.acmap("data/maps/alignment_map.ace")


agNames(alignment_map) <- alignment_names[agNames(alignment_map),]

keep_ags <- c('Wuhan', 'Alpha', 'Beta', 'Delta', 'BA.1', 'BA.2', 'BA.5', 'XBB.1.5', 'BA.2.86', 'JN.1', "KP.3", "KP.2", "KZ.1.1.1")

# load titer data
titer_data_long <- read.csv("data/titer_data/titer_data_long.csv") %>%
  select(!X)

make_table <- function(titer_data, titer_col){
  titer_data %>%
    select(sr_info, ag_name, all_of(titer_col)) %>%
    pivot_wider(names_from = "sr_info", values_from = titer_col) %>%
    column_to_rownames("ag_name") -> tab
  
  not_titrated_sera <- apply(tab, 2, function(x){
    vals <- unique(x == "*") == TRUE
    !(FALSE %in% vals)
  } )
  
  tab <- tab[, !not_titrated_sera]
  
  return(tab)
}

# ------------------------------------------ MAKE MAPS --------------------------------------
table_thresh20 <- make_table(titer_data_long, "Titer_thresh20")

map <- make_map(table_thresh20, mapColors["Color"], alignment_map, 2000, big_ags = keep_ags, options = list(ignore_disconnected = TRUE,
                                                                                                            dim_annealing = FALSE))

save.acmap(map, "./data/maps/map_threshold20_all_ags.ace")


# Make map with no double NHP samples
titer_data_long %>%
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
  ungroup() %>%
  mutate(sampled_day = as.numeric(gsub("D", "", Time.point.post.challenge))) %>%
  filter(sampled_day <= 10) %>%
  pull(sr_info) -> sera_two_timepoints

map_files <- list.files(map_dir, pattern = ".ace", full.names = TRUE)

map_files <- map_files[(grepl("ags.ace", map_files))]

for(map_f in map_files){
  
  map <- read.acmap(map_f)
  
  map <- removeSera(map, sera_two_timepoints)
  map <- removeSera(map, srNames(map)[grepl("23-01", srNames(map))])
  
  map <- optimize_and_realign_map(removeSera(map, srNames(map)[as.character(srGroups(map)) %in% c("XBB conv.", "BQ.1.1 conv.")]),
                                  alignment_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                             dim_annealing = TRUE))
  
  save.acmap(map, gsub(".ace", "_singleTP_woXBBBQ11conv.ace", map_f))
  
}

# remove all XBB.1.5 conv sera
map_files <- list.files(map_dir, pattern = ".ace", full.names = TRUE)
map_files <- map_files[grepl("woXBBBQ11conv.ace", map_files)]
for(map_f in map_files){
  
  map <- read.acmap(map_f)
  map <- removeSera(map, srNames(map)[srGroups(map) == "XBB.1.5 conv."])
  map <- optimize_and_realign_map(map, alignment_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                             dim_annealing = TRUE))
  
  save.acmap(map, gsub(".ace", "_woXBB15conv.ace", map_f))
  
}

# remove JN.1 and BA.2.86
map_files <- list.files(map_dir, pattern = ".ace", full.names = TRUE)
map_files <- map_files[grepl("woXBBBQ11conv", map_files)]
for(map_f in map_files){
  
  map <- read.acmap(map_f)
  map <- removeAntigens(map, c("CH.1.1"))
  map <- optimize_and_realign_map(map, alignment_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                             dim_annealing = TRUE))
  
  save.acmap(map, gsub(".ace", "_CH11.ace", map_f))
  
}

# optimise reactivity for Alpha. reduce alpha by 2 fold
map_files <- list.files(map_dir, pattern = ".ace", full.names = TRUE)
map_files <- map_files[grepl("woXBBBQ11conv", map_files)]

# down here automated alpha adjust
for(map_f in map_files){
  
  map <- read.acmap(map_f)
  ag_names <- agNames(map)
  
  start_react <- rep(0, length(ag_names))
  start_react[ag_names == "Alpha"] <- NA
  
  map <- optimizeAgReactivity(map, fixed_ag_reactivities = start_react, reoptimize = FALSE, options = list(ignore_disconnected = TRUE, dim_annealing = TRUE))
  #print(paste(map_f, agReactivityAdjustments(map)))
  map <- optimize_and_realign_map(map, alignment_map, 1000, 2, option_list = list(ignore_disconnected = TRUE,
                                                                                  dim_annealing = TRUE))
  
  
  save.acmap(map, gsub(".ace", "_alpha_adj.ace", map_f))
  
}

