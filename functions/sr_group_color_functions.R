# create mix color for serum groups

blend_colors <- function(colors) {
  blend <- apply(col2rgb(colors), 1, mean)/255
  return(blend)
  
}

blend_hex <- function(colors){
    blended <- blend_colors(colors)
    blended_hex <- rgb(blended[1], blended[2], blended[3])
    return(blended_hex)
}

blend_sr_group_colors <- function(sr_group, sr_group_colors, split_by = NULL) {
  
  # split up sr group
  sr_group_split <- str_split(sr_group,pattern = "_")[[1]]

 
  target_identifiers <- sr_group_split[sr_group_split!= "NA"]
 
  sr_group_colors[target_identifiers[!(target_identifiers %in% rownames(sr_group_colors))], "Color"] <- sr_group_colors[target_identifiers[1], "Color"]
  # select colors
  colors <- sr_group_colors[target_identifiers, "Color"]
  colors <- colors[!is.na(colors)]
  
  # blend them
  blended_hex <- blend_hex(colors)
 
  return(blended_hex)
}