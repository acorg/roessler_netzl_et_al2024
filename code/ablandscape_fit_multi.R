#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(ablandscapes)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

fit_ags <- "all" # options are "all", "sub"
sub_ags <- c('CB.1', 'BR.3', 'CH.1.1','BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF')

source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")

figure_dir <- file.path("figures", "landscapes", "gmt_landscapes")
suppressWarnings(dir.create(figure_dir, recursive = T))

# if you want to exclude some antigens in the fit
ags_to_exclude <- c("")
# Read the base map
map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alphaJN1ba286_adjScan.ace")

sr_group_gmt <- read.csv("./data/titer_data/sr_group_gmt_threshold20.csv") %>%
  select(!X)

map <- removeAntigens(map, ags_to_exclude)
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)

if(fit_ags == "sub"){
  padding <- 1
  ags_to_fit_lndscp <- agNames(map)[!(agNames(map) %in% sub_ags)]
} else {
  ags_to_fit_lndscp <- agNames(map)
}

# read the full map
map_orig <- map

sr_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";", row.names = "SerumGroup")

# set the single exposure groups
single_exposure_sr_groups <- levels(srGroups(map))

single_exposure_sr <- srNames(map_orig)[as.character(srGroups(map_orig)) %in% single_exposure_sr_groups]

# subset the map to only multi exposure sera
map_long <- long_map_info(map_orig)

map_long %>%
  select(titer, ag_name, sr_name, sr_group) -> titerdata

titerdata %>%
  group_by(
    sr_group
  ) -> titerdata

titerdata %>%
  group_map(
    get_titertable
  ) -> titertables

lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable[,ags_to_fit_lndscp],
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = map,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
    
  }
)


titertables_groups <- group_data(titerdata)

titertables_groups$sr_group <- as.character(titertables_groups$sr_group)

titerdata$gmt <- sr_group_gmt$logtiter[match(interaction(titerdata$sr_group, titerdata$ag_name), interaction(sr_group_gmt$sr_group, sr_group_gmt$ag_name))]

# Add impulses
titerdata %>%
  # manually set GMT's that are lower than that to LOD2
  mutate(gmt = ifelse(gmt < log2(0.8), log2(0.8), gmt))-> gmt_data


# angle for html page
# angle for html page
angle <- list(
  rotation = c(-1.3594,0.0060, -0.0626), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5
  # zoom = 1.1646 # higher is more zoomed out
)

lndscp_list <- list()
data3js <- base_plot_data3js(map, lndscp_fits, agNames(map), lims, agNames(map))

two_dose_groups <- grep("two doses", titertables_groups$sr_group)

two_dose_landscapes <- plot_landscapes_from_list(data3js, titertables_groups[two_dose_groups, ], lndscp_fits[two_dose_groups], map, gmt_data, highlighted_ags = agNames(map),
                                                 ag_plot_names =  c("Wuhan", "BA.1", "JN.1", "HK.3"),
                                                 ag_labelled = c("Wuhan", "BA.1", "JN.1", "HK.3"),
                                                lndscp_colors = sr_colors,
                                                 add_ag_label = TRUE)


lndscp <-r3js(
  two_dose_landscapes,
  rotation = angle$rotation,
  zoom = angle$zoom
)

lndscp_list[["two_dose_vax"]] <- lndscp
save_name <- file.path(figure_dir, paste0("two_dose_vax_", fit_ags, "_gmt_landscapes"))
plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)


# plot landscapes
for(srg in 1:length(unique(titertables_groups$sr_group))){
  
        target_rows <- srg
        lndscp_fits_t <- lndscp_fits[target_rows]
        titertables_groups_t <- titertables_groups[target_rows,]
      
        lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors)
    
        lndscp <-r3js(
          lndscp_3js,
          rotation = angle$rotation,
          zoom = angle$zoom
        )

        
        srg_n <- titertables_groups$sr_group[target_rows]
        srg_n <- gsub("/", "_", srg_n)
        srg_n <- gsub(" ", "", srg_n)
        save_name <- file.path(figure_dir, paste0(srg_n,"_", fit_ags, "_gmt_landscapes"))
        plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)
        
        lndscp_list[[srg_n]] <- lndscp
}

lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups, lndscp_fits, map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors,
  show_gmts = FALSE,
  hide_buttons = FALSE)
        
lndscp <-r3js(
          lndscp_3js,
          rotation = angle$rotation,
          zoom = angle$zoom
        )

save_name <- file.path(figure_dir, paste0(fit_ags, "_all_gmt_landscapes"))
plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)

lndscp_list[["all"]] <- lndscp
saveRDS(lndscp_list, paste0("data/landscape_fit/biv_boosts_ags_",fit_ags,".rds"))


