rm(list = ls())
library(Racmacs)
library(tidyverse)
library(titertools)
library(stringr)
library(ggsci)
library(patchwork)

# define titerplot function
source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")
source("./functions/sr_group_color_functions.R")

sr_group_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";",
                            row.names = "SerumGroup")

ag_order <- read.csv("data/metadata/ag_order.csv")$x

#------------------- Threshold <1
ymax <- 8.5


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

time_points <- unique(data_long_sub$Time.point.post.challenge)

time_colors <- data.frame("Color" = pal_npg("nrc")(length(time_points)))
rownames(time_colors) <- time_points

data_long_sub$sr_group <- factor(data_long_sub$sr_group, levels = unique(data_long_sub$sr_group))

facet_levels <- unique(data_long_sub[order(data_long_sub$sr_group), "nhp_sr_group"])$nhp_sr_group

titerplot20 <- do_titer_plot_fc_label(map, 3, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                      sr_group_gmt_plotdata = NULL,
                                      fc_df = NULL,
                                      target_sr_groups = unique(data_long_sub$sr_group),
                                      ag_order = ag_order,
                                      sr_group_colors = time_colors,
                                      ymax = ymax,
                                      data_long = data_long_sub,
                                      do_ag_tiles = FALSE,
                                      show_gmt = FALSE,
                                      titertools_gmt = FALSE,
                                      facet_var = "nhp_sr_group",
                                      color_var = "Time.point.post.challenge",
                                      facet_levels = facet_levels,
                                      legend_name = "Timepoint post challenge") +
  theme(legend.position = "top")

# calculate fc from one day to next
data_long_sub %>%
  ungroup() %>%
  select(!all_of(c("X", "sr_group_time", "sr_info", "sr_name", "titer", "titertype"))) %>%
  select(!Titer:Titer_thresh20) %>%
  pivot_wider(names_from = "Time.point.post.challenge", values_from = logtiter) %>%
  mutate(d10_to_d20 = D20 - D10,
         d10_to_d21 = D21 - D10,
         d10_to_d28 = D28 - D10) %>% 
  select(!D10:D28) %>%
  pivot_longer(cols = paste0("d10_to_d", c("20", "21", "28")), names_to = "fc_from", values_to = "logtiter") %>%
  filter(!is.na(logtiter)) %>%
  mutate(nhp_sr_group = factor(nhp_sr_group, levels = facet_levels))-> fc_df


titerplot20 + 
  geom_line(data = fc_df %>% mutate(Time.point.post.challenge = "D10"), aes(color = fc_from, group = nhp_sr_group), 
            color = "grey60") + 
  scale_y_continuous(
    labels  = function(x) round(2^x*10),
    name = "Titer",
    breaks = c(0:8),
   sec.axis = sec_axis(~., name="Fold change from D10",
                       labels = function(x) 2^x,
                       breaks = c(0:8))
  ) + 
  theme(
    axis.text.y.right =  element_text(color = "grey60"),
    axis.title.y.right = element_text(color = "grey60")
  ) -> titerplot20_fc


data_long_sub %>%
  mutate(sr_group_big = sr_group,
         sr_group = paste0(sr_group, "_", Time.point.post.challenge)) -> data_long_day_sr

target_srgs <- unique(data_long_day_sr$sr_group)

sr_group_plot <- do_titer_plot_fc_label(map, 1, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                                      sr_group_gmt_plotdata = NULL,
                                      fc_df = NULL,
                                      target_sr_groups = target_srgs,
                                      ag_order = ag_order,
                                      sr_group_colors = time_colors,
                                      ymax = ymax,
                                      data_long = data_long_day_sr,
                                      do_ag_tiles = FALSE,
                                      show_gmt = TRUE,
                                      titertools_gmt = FALSE,
                                      facet_var = "sr_group_big",
                                      color_var = "sr_group",
                                      facet_levels = c("Wuhan conv.", "alpha conv.","beta conv.", "BA.1 conv."),
                                      legend_name = "Timepoint post challenge")
ggsave("som/titers_over_time/sr_group_day_titerplots.png", sr_group_plot + 
         theme(legend.position = "top"), dpi = 300, width = 8, height = 8)

# titerplot of fold changes
fc_df %>%
  mutate(sr_group_big = factor(sr_group, levels = c("Wuhan conv.", "alpha conv.","beta conv.", "BA.1 conv.")),
         fc_from = gsub("d", "D", fc_from),
         sr_group = paste0(sr_group, ": ", gsub("_", " ", fc_from)),
         sr_name = nhp_sr_group,
         titertype = 1,
         titer = 2^logtiter*10,
         Timepoint.post.challenge = "D10") -> fc_df_srg


fc_df_groups <- fc_df_srg[order(fc_df_srg$sr_group_big),] %>%
  pull(sr_group) %>%
  unique()

fc_df_colors <- data.frame(Color = rep("grey60", length(fc_df_groups)))
rownames(fc_df_colors) <- fc_df_groups

do_titer_plot_fc_label(map, 1, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                       sr_group_gmt_plotdata = NULL,
                       fc_df = NULL,
                       target_sr_groups = fc_df_groups,
                       ag_order = ag_order,
                       sr_group_colors = fc_df_colors,
                       ymax = ymax,
                       data_long = fc_df_srg,
                       do_ag_tiles = FALSE,
                       show_gmt = TRUE,
                       titertools_gmt = FALSE,
                       facet_var = "sr_group",
                       color_var = "sr_group",
                       facet_levels = fc_df_groups,
                       legend_name = "Timepoint post challenge",
                       do_below_thresh_tile = FALSE,
                       clamp_below_thresh = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(
    labels  = function(x) round(2^x, 2),
    name = "Fold change from D10",
    limits = c(-2, 5),
    breaks = c(-2:8),
  ) +
  coord_cartesian(
    ylim = shrinkrange(c(-2, 5), 0.05)
  ) -> fc_df_plot


fc_df_colors <- data.frame("Color" =  sr_group_colors[paste(c("Wuhan", "alpha", "beta", "BA.1"), "conv."), "Color"])
rownames(fc_df_colors) <- fc_df_groups

temp_plots <- list()
for(srg in fc_df_groups){
  
  do_titer_plot_fc_label(map, 1, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                         sr_group_gmt_plotdata = NULL,
                         fc_df = NULL,
                         target_sr_groups = srg,
                         ag_order = ag_order,
                         sr_group_colors = fc_df_colors,
                         ymax = ymax,
                         data_long = fc_df_srg,
                         do_ag_tiles = FALSE,
                         show_gmt = TRUE,
                         titertools_gmt = FALSE,
                         facet_var = "sr_group",
                         color_var = "sr_group",
                         facet_levels = srg,
                         legend_name = "Timepoint post challenge",
                         do_below_thresh_tile = FALSE,
                         clamp_below_thresh = FALSE,
                         order_by_gmt = TRUE) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(
      labels  = function(x) round(2^x, 2),
      name = "Fold change from D10",
      limits = c(-2, 5),
      breaks = c(-2:8),
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(-2, 5), 0.05)
    ) -> fc_df_plot
  
  temp_plots[[srg]] <- fc_df_plot
}

singles <- patchwork::wrap_plots(temp_plots, ncol = 1, axis_titles = "collect")


fc_df_srg %>%
  mutate(sr_group_big = "Idvls") -> fc_df_srg

fc_df_colors_temp <- fc_df_colors$Color
names(fc_df_colors_temp) <- rownames(fc_df_colors)

do_titer_plot_fc_label(map, 1, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                       sr_group_gmt_plotdata = NULL,
                       fc_df = NULL,
                       target_sr_groups = fc_df_groups,
                       ag_order = ag_order,
                       sr_group_colors = fc_df_colors,
                       ymax = ymax,
                       data_long = fc_df_srg,
                       do_ag_tiles = FALSE,
                       show_gmt = TRUE,
                       titertools_gmt = FALSE,
                       facet_var = "sr_group_big",
                       color_var = "sr_group",
                       facet_levels = c("Idvls"),
                       legend_name = "Serum group",
                       do_below_thresh_tile = FALSE,
                       clamp_below_thresh = FALSE,
                       order_by_gmt = TRUE,
                      use_ag_number = TRUE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(
    labels  = function(x) round(2^x, 2),
    name = "Fold change from D10",
    limits = c(-2, 5),
    breaks = c(-2:8),
  ) +
  coord_cartesian(
    ylim = shrinkrange(c(-2, 5), 0.05)
  ) + 
  facet_wrap(~sr_group_big, ncol = 1,
             labeller = labeller(sr_group_big = c("Idvls" = "Individuals", "NA" = "GMT"))) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = c(0.9, 0.9)) + 
  scale_color_manual(values = fc_df_colors_temp,
                     name = "Serum group",
                     labels = capitalize(gsub(" conv.| vax.","", names(fc_df_colors_temp)))) + 
  scale_fill_manual(values = fc_df_colors_temp,
                    name = "Serum group",
                    labels = capitalize(gsub(" conv.| vax.","", names(fc_df_colors_temp))))-> fc_df_plot

(singles + plot_layout(axis_titles = "collect")) | fc_df_plot -> fc_comb
ggsave("som/titers_over_time/fc_plots.png", fc_comb, dpi = 300, width = 10, height = 8)



titerplot20_fc /(sr_group_plot + fc_df_plot) + plot_layout(widths = c(3, 1, 1)) + plot_annotation(tag_levels = "A") -> plot_comb

ggsave("som/titers_over_time/titerplots.png", plot_comb, dpi = 300, width = 12, height = 14)


# make plot ranked by fold change per group
color_var <- "sr_group"
calc_half_thresh_sr_group_gmt(fc_df_srg, 20) %>%
  arrange(color_var, "logtiter") %>%
  group_by_at(vars(one_of(color_var))) %>% 
  mutate(rank = rank(desc(logtiter), ties.method = "first")) -> fc_rank


temp_plots <- list()
for(srg in c("Wuhan conv.", "alpha conv.","beta conv.", "BA.1 conv.")){
  
  ag_order_temp <- fc_rank[grepl(srg, fc_rank$sr_group),]
  
  ag_order_temp[order(ag_order_temp$rank),] %>%
    pull(ag_name) -> ag_order_temp
  
  do_titer_plot_fc_label(map, 1, thresh = 20, fc_label = FALSE, adj_titers = FALSE,
                         sr_group_gmt_plotdata = NULL,
                         fc_df = NULL,
                         target_sr_groups = target_srgs[grepl(srg, target_srgs)],
                         ag_order = ag_order_temp,
                         sr_group_colors = time_colors,
                         ymax = ymax,
                         data_long = data_long_day_sr,
                         do_ag_tiles = FALSE,
                         show_gmt = TRUE,
                         titertools_gmt = FALSE,
                         facet_var = "sr_group_big",
                         color_var = "sr_group",
                         do_below_thresh_tile = TRUE,
                         clamp_below_thresh = TRUE,
                         order_by_gmt = FALSE,
                         facet_levels = c("Wuhan conv.", "alpha conv.","beta conv.", "BA.1 conv."),
                         legend_name = "") + 
    theme(legend.position = c(0.8,0.8),
          legend.key = element_blank(),
          legend.background=element_blank()) -> temp_plot
  
  temp_plots[[srg]] <- temp_plot
}

singles <- patchwork::wrap_plots(temp_plots, ncol = 1, axis_titles = "collect") + plot_layout(axis_titles = "collect")

(singles | fc_comb) + plot_layout(widths = c(1, 2)) -> titer_fc_ranks

ggsave("som/titers_over_time/titer_fc_plots_ranked.png", titer_fc_ranks, dpi = 300, width = 12, height = 8)
