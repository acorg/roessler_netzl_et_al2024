rm(list = ls())
library(Racmacs)
library(tidyverse)
library(titertools)
library(stringr)
library(patchwork)

# define titerplot function
source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")
source("./functions/sr_group_color_functions.R")
source("functions/map_functions.R")

sr_group_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";",
                            row.names = "SerumGroup")


ag_order <- read.csv("data/metadata/ag_order.csv")$x

ymax <- 11.5 #12.5 with label

map <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv.ace")
map_ch11 <- read.acmap("data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv.ace")
srGroups(map) <- gsub("BA.4 conv.|BA.5 conv.", 'BA.4/5 conv.', srGroups(map))

data_long <- long_map_info(map)

data_long %>%
  pull(sr_group) %>%
  as.character() %>%
  unique() -> target_groups

target_groups <- c('Wuhan vax. (single dose)','Beta vax. (single dose)','Wuhan vax. (two doses)', 'XBB.1.5 vax. (two doses)',
                   'Wuhan conv.', 'alpha conv.','beta conv.',
                   'gamma conv.', 'delta conv.', 'BA.1 conv.', 'BA.2.12.1 conv.', 'BA.4/5 conv.') 

paste(target_groups, collapse = "', '")

target_ags <- c('Wuhan', 'Alpha', 'Delta', 'Beta', 'BA.1', 'BA.2', 'BA.5', 'XBB.1.5', 'BA.2.86', 'JN.1', "KP.3", "KP.2", "KZ.1.1.1")

#------------------- Threshold <20

if(file.exists("data/titer_data/sr_group_gmt_threshold20_unadj.csv")){
  sr_group_gmt <- read.csv("data/titer_data/sr_group_gmt_threshold20_unadj.csv") %>%
    select(!X)
  
} else {
  sr_group_gmt <- calc_titertools_sr_group_gmt(data_long, thresh = 20)
  write.csv(sr_group_gmt, "data/titer_data/sr_group_gmt_threshold20_unadj.csv")
  
}

if(file.exists("data/titer_data/fc_from_homologous_threshold20_unadj.csv")){
  fc_df <- read.csv("data/titer_data/fc_from_homologous_threshold20_unadj.csv") %>%
    select(!X)
  
} else {
  fc_df <- calc_fc_from_homologous(data_long)
  write.csv(fc_df, "data/titer_data/fc_from_homologous_threshold20_unadj.csv")
  
}


target_sr_groups_a <- target_groups[1:4]
target_sr_groups_b <- target_groups[5:length(target_groups)]

titerplot20_a <- do_titer_plot_fc_label(map, 4, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        target_sr_groups = target_sr_groups_a,
                                        ag_order = ag_order,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE,
                                        facet_levels = target_sr_groups_a) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

targets_b_main <- target_sr_groups_b
titerplot20_b <- do_titer_plot_fc_label(map, 4, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% targets_b_main),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% targets_b_main),
                                        target_sr_groups = targets_b_main,
                                        ag_order = ag_order,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE,
                                        facet_levels = targets_b_main) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# (plot_spacer() + titerplot20_b + plot_spacer() + plot_layout(widths = c(0.25, 3, 0.4))) / titerplot20_a + plot_layout(heights = c(3.5, 1)) + plot_annotation(tag_levels = 'A') -> titerplot20

titerplot20_b / titerplot20_a + plot_layout(heights = c(2.2, 1)) + plot_annotation(tag_levels = 'A') -> titerplot20

ggsave("figures/titer_lineplots/sr_group_titerlineplot_threshold20_map_unadj_noLabel_wCH11.png", titerplot20, dpi = 300, width = 12, height = 8)

ag_order_no_ch11 <- ag_order[!grepl("CH.1.1", ag_order)]

titerplot20_a <- do_titer_plot_fc_label(map, 4, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        target_sr_groups = target_sr_groups_a,
                                        ag_order = ag_order_no_ch11,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE,
                                        facet_levels = target_sr_groups_a) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

targets_b_main <- target_sr_groups_b
titerplot20_b <- do_titer_plot_fc_label(map, 4, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% targets_b_main),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% targets_b_main),
                                        target_sr_groups = targets_b_main,
                                        ag_order = ag_order_no_ch11,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE,
                                        facet_levels = targets_b_main) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# (plot_spacer() + titerplot20_b + plot_spacer() + plot_layout(widths = c(0.25, 3, 0.4))) / titerplot20_a + plot_layout(heights = c(3.5, 1)) + plot_annotation(tag_levels = 'A') -> titerplot20

titerplot20_b / titerplot20_a + plot_layout(heights = c(2.2, 1)) + plot_annotation(tag_levels = 'A') -> titerplot20

ggsave("figures/titer_lineplots/sr_group_titerlineplot_threshold20_map_unadj_noLabel.png", titerplot20, dpi = 300, width = 12, height = 8)



#--------------------------------------------- Load alpha adjusted data
data_long_ch11 <- data_long
sr_gmt_ch11 <- sr_group_gmt


map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")

data_long <- long_map_info(remove_na_coords(map)) %>%
  mutate(titer = titer_adjusted,
         logtiter = logtiter_adjusted)

data_long %>%
  select(sr_group, sr_name) %>%
  unique() %>%
  group_by(sr_group) %>%
  summarize(n = length(sr_name)) -> sr_table

sr_table_np <- long_map_info(map) %>%
  select(sr_group, sr_name) %>%
  unique() %>%
  group_by(sr_group) %>%
  summarize(n_all = length(sr_name)) %>%
  cbind(sr_table[,2]) %>%
  mutate(n_np = n_all - n,
         "Serum Group" = sr_group,
    "Sample size" = ifelse(n_np > 0, paste0(n_all, " (", n_np, ")"), n)) %>%
  select(!n_all:n_np)

sr_table_np[match(target_groups, sr_table_np$sr_group),] %>%
  filter(!is.na(sr_group)) %>%
  select(!sr_group)-> sr_table_np

write.csv(sr_table_np, "data/titer_data/sr_table.csv", row.names = FALSE)


map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")

srGroups(map) <- gsub("BA.4 conv.|BA.5 conv.", 'BA.4/5 conv.', srGroups(map))

data_long <- long_map_info(map) %>%
  mutate(titer = titer_adjusted,
         logtiter = logtiter_adjusted)


if(file.exists("data/titer_data/sr_group_gmt_threshold20.csv")){
  sr_group_gmt <- read.csv("data/titer_data/sr_group_gmt_threshold20.csv") %>%
    select(!X)
  
} else {
  sr_group_gmt <- calc_titertools_sr_group_gmt(data_long, thresh = 20)
  write.csv(sr_group_gmt, "data/titer_data/sr_group_gmt_threshold20.csv")
  
}

if(file.exists("data/titer_data/fc_from_homologous_threshold20.csv")){
  fc_df <- read.csv("data/titer_data/fc_from_homologous_threshold20.csv") %>%
    select(!X)
  
} else {
  fc_df <- calc_fc_from_homologous(data_long)
  write.csv(fc_df, "data/titer_data/fc_from_homologous_threshold20.csv")
  
}

sr_group_gmt <- rbind(sr_group_gmt,
                      sr_gmt_ch11 %>%
                        filter(ag_name == "CH.1.1"))

data_long <- rbind(data_long,
                   data_long_ch11 %>%
                     filter(ag_name == "CH.1.1"))

titerplot20_a <- do_titer_plot_fc_label(map_ch11, 4, thresh = 20, fc_label = FALSE, adj_titers = TRUE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% target_sr_groups_a),
                                        target_sr_groups = target_sr_groups_a,
                                        ag_order = ag_order,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        facet_levels = target_sr_groups_a,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

titerplot20_b <- do_titer_plot_fc_label(map_ch11, 4, thresh = 20, fc_label = FALSE, adj_titers = TRUE,
                                        sr_group_gmt_plotdata = sr_group_gmt %>%
                                          filter(sr_group %in% target_sr_groups_b),
                                        fc_df = fc_df %>%
                                          filter(sr_group %in% target_sr_groups_b),
                                        target_sr_groups = target_sr_groups_b,
                                        ag_order = ag_order,
                                        sr_group_colors = sr_group_colors,
                                        ymax = ymax,
                                        facet_levels = target_sr_groups_b,
                                        show_gmt_conf = TRUE)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

titerplot20_b / titerplot20_a + plot_layout(heights = c(2.2, 1)) + plot_annotation(tag_levels = 'A') -> titerplot20

ggsave("figures/titer_lineplots/sr_group_titerlineplot_threshold20_map_noLabel_wCH11.png", titerplot20, dpi = 300, width = 12, height = 8)



