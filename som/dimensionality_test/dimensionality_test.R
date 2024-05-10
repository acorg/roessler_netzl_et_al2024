rm(list = ls())
library(Racmacs)
library(ggplot2)
set.seed(100)
map <- read.acmap('./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace')


# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 1000,
  replicates_per_dimension = 100,
  options = list(dim_annealing = FALSE,
                 ignore_disconnected = TRUE)
)

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))

saveRDS(df, "./som/dimensionality_test/dimension_test_result_threshold20_wo_dim_annealing.rds")


# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 1000,
  replicates_per_dimension = 100,
  options = list(dim_annealing = TRUE,
                 ignore_disconnected = TRUE)
)

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))

saveRDS(df, "./som/dimensionality_test/dimension_test_result_threshold20.rds")


# do the plot
df_dim <- readRDS("./som/dimensionality_test/dimension_test_result_threshold20.rds")%>%
  mutate("Dimensional annealing" = TRUE)
df <- readRDS("./som/dimensionality_test/dimension_test_result_threshold20_wo_dim_annealing.rds") %>%
  mutate("Dimensional annealing" = FALSE) %>%
  rbind(.,df_dim) 
ggplot(data=df, aes(x=dimensions, y=rmse, color = `Dimensional annealing`)) +
  geom_line()+
  geom_point() +
  theme_bw() +
  ylim(c(0,2))+
  xlab('Number of dimensions') +
  ylab('Mean RMSE of detectable titers') +
  theme(strip.background = element_blank(),
        legend.position = "top") ->dp


png("./som/dimensionality_test/dimension_test_threshold20_boht.png", 4, 4, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
dp
dev.off()


# Plot a 3D map
map3D <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 3,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none",
  options = list(ignore_disconnected = TRUE,
                 dim_annealing = TRUE)
)

map3D <- applyPlotspec(map3D, map)
map3D <- realignMap(map3D, map)

p <- procrustesMap(map3D, map, sera=FALSE)
agSize(p) <- 4
srSize(p) <- 3
Racmacs::view(p)


map1D <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 1,
  number_of_optimizations = 1000,
  minimum_column_basis    = "none",
  options = list(ignore_disconnected = TRUE,
                 dim_annealing = TRUE)
)

p <- procrustesMap(map1D, map, sera=FALSE)
agSize(p) <- 4
srSize(p) <- 3
Racmacs::view(p)