library(Racmacs)
set.seed(100)
source("./functions/map_functions.R")

# load data
mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Variable', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 



sample_w_replacement <- TRUE

map_dir <- file.path("data", "maps")
figure_dir <- file.path("map_diagnostics", "subset_serum_maps")

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x +1
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x +1


n_samples <- c(1:6)


map_files <- list.files(map_dir, pattern = ".ace", full.names = FALSE)
map_files <- map_files[grepl("_woXBB15conv_woJN1BA286_alpha_adj", map_files)]

for(map_f in map_files){
  
  neut <- read.acmap(file.path(map_dir, map_f))
  
  
  sr_groups <- as.character(unique(srGroups(neut)))
  
  for(n in n_samples) {
    
    png(file.path(figure_dir, paste0(n, "_", gsub(".ace", ".png", map_f))), width = 8, height = 3, units = 'in', res=300, pointsize = 18)
    layout(matrix(c(1:10), ncol = 5, byrow = T))
    par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))
    
    for(i in 1:10){
      
      sera <- c()
    
      for(sg in sr_groups){
        
       
        sr_names <- srNames(neut)[as.character(srGroups(neut)) == sg]
        
        n_sera <- length(sr_names)
        
      
        if(sample_w_replacement){
          
          sera <- c(sera, sr_names[sample.int(n_sera, n, replace = TRUE)])
          
        } else {
          if(n_sera > n){
            sera <- c(sera, sr_names[sample.int(n_sera, n, replace = FALSE)])
          } else {
            sera <- c(sera, sr_names)
          }
        }
        
        
      }
      
      
      new_titertable <- adjustedTiterTable(neut)[,sera]
  
      colnames(new_titertable) <- paste0(colnames(new_titertable), "-", 1:ncol(new_titertable))
      sub_map <- make_map(table = new_titertable, mapColors["Color"], neut, 1000, options =  list(ignore_disconnected = TRUE,
                                                                                                  dim_annealing = TRUE))
      
      plot(procrustesMap(sub_map, neut, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
           grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
      
    }
    
    dev.off()
  }
    
    
  }
  
  
#   for(n in n_samples){
# 
#     png(paste0("som/subset_serum_maps/",n,"_sera.png"), width = 8, height = 3, units = 'in', res=300, pointsize = 18)
#     layout(matrix(c(1:10), ncol = 5, byrow = T))
#     par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))
# 
#     for(i in 1:10){
# 
#       p <- read.acmap(paste0("som/subset_serum_maps/n_serum_",n,"_run_", i, ".ace"))
# 
#       plot(procrustesMap(p, neut, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
#            grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
# 
#     }
# 
#     dev.off()
# 
#   }
# 
# 
# 
# }

