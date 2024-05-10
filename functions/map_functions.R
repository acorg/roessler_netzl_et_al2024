makeMap <- function(table, baseMap = NULL, dimensions = 2,
                    nOptimisations = 500, mcb = "none", dilution_stepsize = 1, options = list()) {
  # Take a titer table, make a map. Optionally apply plotspec and re-align
  # map to an already existing map.
  m <- acmap(titer_table = table)
  
  dilutionStepsize(m) <- dilution_stepsize
  
  m <- optimizeMap(
    map                     = m,
    number_of_dimensions    = dimensions,
    number_of_optimizations = nOptimisations,
    minimum_column_basis    = mcb,
    options = options
  )
  
  agNames(m) <- rownames(table)
  srNames(m) <- colnames(table)
  
  if(!(is.null(baseMap))) {
    
    m <- applyPlotspec(m, baseMap)
    m <- realignMap(m, baseMap)
  }
  
  ptDrawingOrder(m) <- rev(seq_len(numPoints(m)))
  
  return(m)
}

apply_color <- function(map, map_colors) {
  
  agFill(map) <- map_colors[agNames(map),1]

  srOutline(map) <- map_colors[as.character(srGroups(map)),1]

  return(map)
}

apply_style <- function(map){
  
  nag <- numAntigens(map)
  nsr <- numSera(map)
  N <- nsr + nag
  ptDrawingOrder(map) <- c(N:(nag+1), 1:nag)
  
  srOutlineWidth(map) <- 1
  srSize(map) <- 9
  agSize(map) <- 18
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.7)
  
  sublineages <- c("B.1.1.7+E484K", "BE.1.1", "BA.5.2.1", "CB.1", "BQ.1.3", "BQ.1.1",
    "BQ.1.18", "BR.3", "CH.1.1", "BF.7")
  agSize(map)[agNames(map) %in% sublineages] <- 15
  
  return (map)
  
}


# make the map
make_map <- function(table, map_colors_info, alignment_map, n_optim =1000, options = list(ignore_disconnected = TRUE),
                     big_ags = c("Wuhan")) {
  
  # make map
  map <- acmap(
    ag_names = rownames(table),
    sr_names = colnames(table),
    titer_table = table
  )
  dilutionStepsize(map) <- 0
  
  sr_groups_map <- unlist(lapply(colnames(titerTable(map)), function(x) {
    strsplit(x, "_")[[1]][1]
  }))
  
  
  #set serum groups
  srGroups(map) <- factor(
    sr_groups_map)
  
  #set antigen colours
  ag_colors <- map_colors_info[agNames(map),]
  agFill(map) <- ag_colors
  
  #set Serum colours
  
  sr_colors <- map_colors_info[as.character(srGroups(map)),]
  srOutline(map) <- sr_colors
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.7)
  
  # Set styles
  srOutlineWidth(map) <- 1
  srSize(map) <- 9
  agSize(map) <- 14
  
  agSize(map)[agNames(map) %in% big_ags] <- 18
  
  ptDrawingOrder(map) <- rev(ptDrawingOrder(map))
  
  
  map_optim <- optimizeMap(
    map,
    number_of_dimensions = 2,
    number_of_optimizations = n_optim,
    minimum_column_basis = "none",
    options = options
  )
  
  # realign map
  map_optim <- realignMap(map_optim, target_map = alignment_map)
  
  return(map_optim)
}

remove_na_sera <- function(map){
  
  na_sera <- srCoords(map)
  na_sera <- rownames(na_sera)[is.na(na_sera)[,1]]
  
  map <- removeSera(map, na_sera)
  
  return(map)
}

remove_na_antigens <- function(map){
  
  na_sera <- agCoords(map)
  na_sera <- rownames(na_sera)[is.na(na_sera)[,1]]
  
  map <- removeAntigens(map, na_sera)
  
  return(map)
}

remove_na_coords <- function(map){
  
  map <- remove_na_sera(map)
  map <- remove_na_antigens(map)
  
  return(map)
  
}

optimize_and_realign_map <- function(map, alignment_map, optimization_numbers = 1000, dim = 2, 
                                     option_list = list(ignore_disconnected = TRUE)){
  
  map <- optimizeMap(map, dim, optimization_numbers, options = option_list)
  map <- realignMap(map, alignment_map)
  
  return(map)
}

