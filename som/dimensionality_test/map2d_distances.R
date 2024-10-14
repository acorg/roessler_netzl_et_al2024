rm(list = ls())
library(Racmacs)
library(ggplot2)
set.seed(100)

map_dir <- "./data/maps/"


mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 

ag_order <- read.csv("data/metadata/ag_order.csv")$x

# move DV.7.1 and CH.1.1 to the back
ag_order <- c(ag_order[!(grepl("CH.1.1", ag_order))])


map_files <- list.files(map_dir, pattern = ".ace", full.names = TRUE)
map_file <- map_files[grepl("_woXBB15conv_CH11_alpha_adj", map_files)]

map <- read.acmap(map_file)


# Plot a 3D map
map3D <- map


ag_coords <- agCoords(map3D)                              
sr_coords <- srCoords(map3D)

ag_dists <- dist(ag_coords)


ag_dist_df <- reshape2::melt(as.matrix(ag_dists), varnames = c("ag", "item"))


ag_dist_df %>%
  filter(value != 0) %>%
  ggplot(aes(x = ag, y = value, color = item)) + 
  geom_point()

# Rank distances
ag_dist_df %>%
  group_by(ag) %>%
  mutate(dist_rank = order(value)) -> ag_dist_df

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# mean rank
ag_dist_df %>%
  filter(value > 0) %>%
  group_by(item) %>%
  summarize(mean_dist = mean(value),
            median_dist = median(value),
            mean_rank = mean(dist_rank),
            mode_rank = Mode(dist_rank),
            median_rank = median(dist_rank)) %>% View()

# compute here serum distances from variants
sr_dist_df <- data.frame("ag" = c(),
                         "sr" = c(),
                         "value" = c(),
                         "sr_group" = c())

for(ag in rownames(ag_coords)){
  
  for(sr_num in 1:nrow(sr_coords)){
    
    temp_coords <- rbind(ag_coords[ag,], sr_coords[sr_num,])
    
    temp_df <- data.frame("ag" = ag,
                          "sr" = rownames(sr_coords)[sr_num],
                          "value" = as.numeric(dist(temp_coords)),
                          "sr_group" = as.character(srGroups(map3D))[sr_num])
    
    sr_dist_df <- rbind(sr_dist_df, temp_df)
    
  }
  
}


sr_dist_df %>%
  filter(!is.na(value)) -> sr_dist_df


sr_dist_df %>%
  ggplot(aes(x = sr_group, y = value)) + 
  geom_point() + 
  facet_wrap(~ag)


ag_colors <- mapColors$Color
names(ag_colors) <- rownames(mapColors)

target_groups <- c('Wuhan vax. (single dose)','Beta vax. (single dose)','Wuhan vax. (two doses)', 'XBB.1.5 vax. (two doses)',
                   'Wuhan conv.', 'alpha conv.','beta conv.',
                   'gamma conv.', 'delta conv.', 'BA.1 conv.', 'BA.2.12.1 conv.', 'BA.4 conv.', 'BA.5 conv.')

sr_dist_df$sr_group <- factor(sr_dist_df$sr_group, levels = target_groups)

sr_dist_df %>%
  ggplot(aes(x = ag, y = value, color = ag)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = ag_colors) +
  facet_wrap(~sr_group) + 
  theme_bw() + 
  ylim(c(0,10)) +
  scale_x_discrete(limits = ag_order,
                   name = "Variant") +
  ylab("Euclidean distance in 2D") +
  theme(strip.background.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) -> p_3d


ggsave("som/dimensionality_test/distance_from_sera_2d.png", p_3d, width = 12, height = 8, dpi = 300)

# average distance from sera per variant
sr_dist_df %>%
  group_by(ag) %>%
  summarise(mean_dist = mean(value),
            sd = sd(value),
            lower = Rmisc::CI(value)["lower"],
            upper = Rmisc::CI(value)["upper"]) -> avg_sr_distance

avg_sr_distance %>%
  ggplot(aes(x = ag, y = mean_dist, fill = ag)) + 
 # geom_point(color = "black", shape = 21) + 
  geom_pointrange(aes(ymin = lower, ymax = upper, color = ag)) +
  scale_color_manual(values = ag_colors) +
  theme_bw() + 
  scale_x_discrete(limits = ag_order,
                   name = "Variant") +
  ylab("Mean euclidean distance from sera in 2D") +
  ylim(c(0,7)) +
  theme(strip.background.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) -> p_3d


ggsave("som/dimensionality_test/mean_distance_from_sera_2d.png", p_3d, width = 6, height = 4, dpi = 300)


