# labelled map 
rm(list = ls())
library(Racmacs)

figure_dir <- "figures/labelled_map"
map_dir <- "data/maps/"

move_coords <- function(map, at = 2, by = -0.5){
  agCoords(map)[,at] <- agCoords(map)[,at] + by
  srCoords(map)[,at] <- srCoords(map)[,at] + by
  
  return(map)
}

reoptimize_subset_map <- function(map, sub_ags, sub_srs){
  
  sub_align <- subsetMap(map, antigens = sub_ags, sera = sub_srs)
  sub_align <- optimizeMap(sub_align, 2, 1000)
  sub_align <- realignMap(sub_align, map)
  
  return(sub_align)
  
}

change_map_colors <- function(map, map_colors){
  
  for(x in 1:length(agNames(map))){
    
    temp_name <- agNames(map)[x]

    if(temp_name %in% rownames(map_colors)){
      agFill(map)[x] <- map_colors[temp_name, 1]
    }
        
  }
  
  return(map)
  
}

sr_group_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";",
                            row.names = "SerumGroup")

mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
mapColors[rownames(sr_group_colors), "Color"] <- sr_group_colors$Color 

target_map <- "map_threshold20_all_ags_singleTP_woXBBBQ11conv_alphaJN1ba286_adjScan.ace"
full_map <-read.acmap(file.path("./data/maps", target_map))
full_map <- change_map_colors(full_map, mapColors)

alignment_names <- read.csv("./data/metadata/map-names.csv", row.names = "Antigen")
alignment_map <- read.acmap("./data/maps/alignment_map.ace")

agNames(alignment_map) <- alignment_names[agNames(alignment_map),]

duke_map <- read.acmap("./data/maps/Wilks_et_al_map_ndsubset_no_outliers_slope_adjusted.ace")

old_names <- agNames(duke_map)
agNames(duke_map) <- alignment_names[agNames(duke_map),1]
agNames(duke_map)[agNames(duke_map) == "NA"] <- old_names[agNames(duke_map) == "NA"]

# change map colors
alignment_map <- change_map_colors(alignment_map, mapColors)
duke_map <- change_map_colors(duke_map, mapColors)

# make common maps
subset_ags <- intersect(agNames(full_map),agNames(alignment_map))


srGroups(alignment_map) <- as.character(srGroups(alignment_map))
srGroups(alignment_map) <- gsub("WT", "Wuhan", srGroups(alignment_map))
srGroups(alignment_map) <- gsub("BNT/BNT", "Wuhan vax. (two doses)", srGroups(alignment_map))
srGroups(alignment_map) <- gsub("\\/alpha\\+E484K", "", srGroups(alignment_map))


agSize(alignment_map)[match(agNames(full_map), agNames(alignment_map))[!is.na(match(agNames(full_map), agNames(alignment_map)))]] <- agSize(full_map)[match(agNames(alignment_map), agNames(full_map))[!is.na(match(agNames(alignment_map), agNames(full_map)))]]
agSize(alignment_map)[agNames(alignment_map) == "BA.5"] <- 18
agSize(alignment_map)[agNames(alignment_map) == "BQ.1.18"] <- 14

map_srs <- as.character(unique(srGroups(full_map)))
alignment_srs <- unique(srGroups(alignment_map))
subset_sr_groups <- intersect(map_srs, alignment_srs)


if(file.exists("data/maps/roessler_et_al2023_subNHP.ace")){
  sub_align <- read.acmap("data/maps/roessler_et_al2023_subNHP.ace")
  
} else {
  
  sub_align <- reoptimize_subset_map(alignment_map, subset_ags, srNames(alignment_map)[as.character(srGroups(alignment_map)) %in% subset_sr_groups])
  
  save.acmap(sub_align, "data/maps/roessler_et_al2023_subNHP.ace")
  
}


# make subset NHP map
nhp_kimpel <- reoptimize_subset_map(full_map, subset_ags, srNames(full_map)[as.character(srGroups(full_map)) %in% subset_sr_groups])

# duke map
srGroups(duke_map) <- gsub("B.1.617.2", "delta conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub("B.1.1.7", "alpha conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub("B.1.351", "beta conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub("P.1", "gamma conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub("BA.1", "BA.1 conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub("2x mRNA-1273", "Wuhan vax. (two doses)", srGroups(duke_map))
srGroups(duke_map) <- gsub("D614G", "Wuhan conv.", srGroups(duke_map))
srGroups(duke_map) <- gsub(" new", "", srGroups(duke_map))
subset_sr_groups <- intersect(map_srs, srGroups(duke_map))

agSize(duke_map)[match(agNames(full_map), agNames(duke_map))[!is.na(match(agNames(full_map), agNames(duke_map)))]] <- agSize(full_map)[match(agNames(duke_map), agNames(full_map))[!is.na(match(agNames(duke_map), agNames(full_map)))]]

subset_ags <- intersect(agNames(full_map),agNames(duke_map))

nhp_duke <- reoptimize_subset_map(full_map, subset_ags, srNames(full_map)[as.character(srGroups(full_map)) %in% subset_sr_groups])


if(file.exists("data/maps/Wilks_et_al_subNHP.ace")){
  sub_duke <- read.acmap("data/maps/Wilks_et_al_subNHP.ace")
  
} else {
  
  sub_duke <- reoptimize_subset_map(duke_map, subset_ags, srNames(duke_map)[as.character(srGroups(duke_map)) %in% subset_sr_groups])
  
  save.acmap(sub_duke, "data/maps/Wilks_et_al_subNHP.ace")
  
}

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
xlim_no_zoom[1] <- xlim_no_zoom[1] + 2 
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x
ylim_no_zoom[1] <- ylim_no_zoom[1] + 1
ylim_no_zoom[2] <- ylim_no_zoom[2] - 1

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



#ylim_no_zoom[1] <- ylim_no_zoom[1] -1
# 
# 
# png(file.path(figure_dir, "proc_map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj_wSubset.png"), width = 8, height = 9, units = 'in', res=300, pointsize = 18)
# layout(matrix(c(1:6), ncol = 2, byrow = T))
# par(mar = rep(0.5, 4))
# doplot(procrustesMap(alignment_map, full_map, scaling = TRUE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
# doplot(procrustesMap(move_coords(duke_map, at = 1, by = -3), full_map, scaling = TRUE), xlim_no_zoom-3, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]-3+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
# doplot(procrustesMap(sub_align, nhp_kimpel, scaling = TRUE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "C", cex = 1.2)
# doplot(procrustesMap(sub_duke, nhp_duke, scaling = TRUE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "D", cex = 1.2)
# doplot(procrustesMap(nhp_kimpel, full_map, scaling = FALSE, sera = FALSE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "E", cex = 1.2)
# doplot(procrustesMap(nhp_duke, full_map, scaling = FALSE, sera = FALSE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "F", cex = 1.2)
# dev.off()
# 
# png(file.path(figure_dir, "proc_map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.png"), width = 12, height = 5, units = 'in', res=300, pointsize = 18)
# layout(matrix(c(1:2), ncol = 2, byrow = T))
# par(mar = rep(0.5, 4))
# doplot(procrustesMap(alignment_map, full_map, scaling = TRUE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
# text(xlim_no_zoom[2]-0.2, ylim_no_zoom[1]+ 0.2, "Human data set 1", cex = 0.8, adj = c(1, 0), col = "grey30")
# doplot(procrustesMap(move_coords(duke_map, at = 1, by = -3), full_map, scaling = TRUE), xlim_no_zoom-3, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]-3+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
# text(xlim_no_zoom[2]-3- 0.2, ylim_no_zoom[1] + 0.2, "Human data set 2", cex = 0.8, adj = c(1, 0), col = "grey30")
# dev.off()

ylim_no_zoom[2] <-ylim_no_zoom[2] + 1
png(file.path(figure_dir, paste0("proc_", gsub(".ace", "", target_map), "_show_only_commonAgs.png")), width = 12, height = 5, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:2), ncol = 2, byrow = T))
par(mar = rep(0.5, 4))
doplot(procrustesMap(removeAntigens(alignment_map, agNames(alignment_map)[!agNames(alignment_map) %in% agNames(full_map)]), full_map, scaling = TRUE), xlim_no_zoom, ylim_no_zoom, FALSE)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
text(xlim_no_zoom[2]-0.2, ylim_no_zoom[1]+ 0.2, "Human data set 1", cex = 0.8, adj = c(1, 0), col = "grey30")
doplot(procrustesMap(move_coords(move_coords(removeAntigens(duke_map, agNames(duke_map)[!agNames(duke_map) %in% agNames(full_map)]), at = 1, by = -2), at =2, by = +2), full_map, scaling = TRUE), xlim_no_zoom-1, ylim_no_zoom+2, FALSE)
text(xlim_no_zoom[1]-1+0.4, ylim_no_zoom[2]+2-0.4, "B", cex = 1.2)
text(xlim_no_zoom[2]-1- 0.2, ylim_no_zoom[1]+2 + 0.2, "Human data set 2", cex = 0.8, adj = c(1, 0), col = "grey30")
dev.off()


# 
# png(file.path(figure_dir, "proc_map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj_unscaled.png"), width = 12, height = 5, units = 'in', res=300, pointsize = 18)
# layout(matrix(c(1:2), ncol = 2, byrow = T))
# par(mar = rep(0.5, 4))
# doplot(procrustesMap(alignment_map, full_map, scaling = FALSE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
# doplot(procrustesMap(duke_map, full_map, scaling = FALSE), xlim_no_zoom, ylim_no_zoom, FALSE)
# text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
# dev.off()
# 
# 
