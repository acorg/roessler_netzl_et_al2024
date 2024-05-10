capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

calc_fc_from_homologous <- function(data_long){
  
  fc_df <- data.frame(sr_group = c(),
                      ag1 = c(),
                      ag2 = c(),
                      fc = c(),
                      fc_lower = c(),
                      fc_upper = c())
  
  ag_names <- unique(data_long$ag_name)
  # calculate fold changes per sr group and ag
  # have to optimise ag selection to not do everything double
  for(sr in unique(data_long$sr_group)){
    

    ag1 <- tolower(str_split(sr, " conv| vax", n = 2)[[1]][1])
    if(ag1 == "ba.4") {
      ag1 <- "ba.5"
    }
    
    if(ag1 %in% tolower(ag_names)){
      
      for(ag2 in ag_names){
       
        titers1 <- data_long %>% filter(sr_group == sr) %>%
          filter(tolower(ag_name) == tolower(ag1)) %>%
          pull(titer_adjusted)
        
        titers2 <- data_long %>% filter(sr_group == sr) %>%
          filter(ag_name == ag2) %>%
          pull(titer_adjusted)
        
        diff <- log2diff(titers1, titers2, dilution_stepsize = 0)
        
        temp_df <- data.frame(sr_group = sr,
                              ag1 = ag1,
                              ag2 = ag2,
                              fc = diff["mean", "estimate"],
                              fc_lower = diff["mean", "lower"],
                              fc_upper = diff["mean", "upper"])
        
        fc_df <- rbind(fc_df, temp_df)
        
      }
      
    }
    
  }
  
  
  return(fc_df)
  
}

calc_fc_from_homologous_gmt <- function(data_long){
  
  fc_df <- data.frame(sr_group = c(),
                      ag1 = c(),
                      ag2 = c(),
                      fc = c())
  
  ag_names <- unique(data_long$ag_name)
 
  # calculate fold changes per sr group and ag
  # have to optimise ag selection to not do everything double
  for(sr in unique(data_long$sr_group)){
    
    ag1 <- tolower(str_split(sr, " conv", n = 2)[[1]][1])
    
    if(ag1 == "ba.4") {
      ag1 <- "ba.5"
    }
    
    if(ag1 %in% tolower(ag_names)){
      
      for(ag2 in ag_names){
       
        titers1 <- data_long %>% filter(sr_group == sr) %>%
          filter(tolower(ag_name) == tolower(ag1)) %>%
          pull(logtiter)
        
        titers2 <- data_long %>% filter(sr_group == sr) %>%
          filter(ag_name == ag2) %>%
          pull(logtiter)
        
        temp_df <- data.frame(sr_group = sr,
                              ag1 = ag1,
                              ag2 = ag2,
                              fc = titers2 - titers1)
        
        fc_df <- rbind(fc_df, temp_df)
        
      }
      
    }
    
  }
  
  
  return(fc_df)
  
}


calc_titertools_sr_group_gmt <- function(data_long, thresh){
  
  data_long %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    reframe(
      gmt = t(titertools::gmt(titer, ci_level = 0.95, dilution_stepsize = 0)[1, 1:3]),
      all_below_thresh = length(titer[titer == paste0("<", thresh)]) == length(titer)
    )  %>% 
    ungroup() %>%
    mutate(logtiter = gmt[,1],
           lower = gmt[,2],
           upper = gmt[,3],
           titer = ifelse(all_below_thresh, paste0("<", thresh), 2^logtiter*10)) %>%
    select(!gmt) -> sr_group_gmt_plotdata
  
  return(sr_group_gmt_plotdata)
}

calc_half_thresh_sr_group_gmt <- function(data_long, thresh){
  
  data_long %>%
    mutate(logtiter = ifelse(titertype == 2, log2(thresh/20), logtiter)) %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    summarize(
      lower = Rmisc::CI(na.omit(logtiter))["lower"],
      upper = Rmisc::CI(na.omit(logtiter))["upper"],
      logtiter = Rmisc::CI(na.omit(logtiter))["mean"],
      all_below_thresh = length(titer[titer == paste0("<", thresh)]) == length(titer)
    )  %>% 
    ungroup() %>%
    mutate(titer = ifelse(all_below_thresh, paste0("<", thresh), 2^logtiter*10),
           titer = ifelse(titer == "NaN", NA, titer))-> sr_group_gmt_plotdata
  
  sr_group_gmt_plotdata[is.na(sr_group_gmt_plotdata)] <- NA
  
  return(sr_group_gmt_plotdata)
}


# calculate relative value from different data frame
calculate_realtive_titer <- function(data_big, data_small, column_name = "gmt", new_col_name = "relative_titer"){
  
  data_big_gmts <- data_big[match(interaction(data_small$ag_name, data_small$sr_group), interaction(data_big$ag_name, data_big$sr_group)),column_name]
  data_big_gmts <- (2^data_big_gmts*10)$gmt
  data_small_gmts <- (2^data_small[,column_name]*10)$gmt
  
 # return(data_small_gmts/data_big_gmts)
  fold_change <- unlist(lapply(c(1:length(data_big_gmts)), function(x){
    
    if(is.na(data_small_gmts[x])){
      ""
    } else {
      paste0(round(data_small_gmts[x]/data_big_gmts[x],1))
    }
  }))
  
  data_big[new_col_name] <- ""
  data_big[match(interaction(data_small$ag_name, data_small$sr_group), interaction(data_big$ag_name, data_big$sr_group)),new_col_name] <- fold_change

  return(data_big)  
}

# gmt and fold change label 
get_fold_change_label <- function(data, ymax) {
  
  data %>%
    mutate(titer = ifelse(all_below_thresh, titer, round(as.numeric(titer))),
           fc_label = ifelse(fold_change < 0, paste0("-", round(2^abs(fold_change), 1)), round(2^fold_change, 1)),
           label = paste0(titer, "\n", fc_label),
           label = gsub("NA", "", label),
           y = ymax - 0.5) -> data
  
  return(data)
}

# gmt and fold change label 
get_gmt_label <- function(data, ymax) {
  
  data$titer <- round(2^data$logtiter*10)
  data$titer[is.na(data$titer)] <- ""
  data$y <- ymax - 0.5
  data$label <- paste0(data$titer)
  
  return(data)
}


do_titer_plot_fc_label <- function(map, facet_col = 3, thresh = 20, fc_label =TRUE, adj_titers = TRUE,
  highlight_samples = NULL,
  sr_group_gmt_plotdata = NULL,
  fc_df = NULL,
  target_sr_groups = NULL,
  ag_order = NULL,
  sr_group_colors,
  color_var = "sr_group",
  facet_var = "sr_group",
  data_long = NULL,
  titertools_gmt = TRUE,
  titertools_fc = TRUE,
  show_idvls = TRUE,
  show_gmt_conf = FALSE,
  ymax = 11.5,
  facet_levels = c('Wuhan conv.', 'alpha conv.', 'beta conv.', 'gamma conv.', 'delta conv.', 'BA.1 conv.', "BA.2.12.1 conv.", "BA.4 conv.", "BA.5 conv.", 'BQ.1.1 conv.', 'XBB conv.', 'XBB.1.5 conv.'),
  do_ag_tiles = TRUE,
  show_gmt = TRUE,
  legend_name = "Serum group",
  do_below_thresh_tile = TRUE,
  clamp_below_thresh = TRUE,
  order_by_gmt = FALSE,
  use_ag_number = FALSE) {
  
  # Set the antigen x axis order
  if(is.null(ag_order)){
    ag_order <- agNames(map)
  } 
  
  if(is.null(data_long)){
    # Get long info
    plotdata <- long_map_info(map)
  } else {
    plotdata <- data_long
  }
  
  if(is.null(target_sr_groups)){
    target_sr_groups <- as.character(unique(plotdata$sr_group))
  } else {
    plotdata <- plotdata %>%
      filter(sr_group %in% target_sr_groups)
  }
  
  if(adj_titers){
    plotdata$logtiter <- plotdata$logtiter_adjusted
    plotdata$titer <- plotdata$titer_adjusted
  }
  
  ymin <- min(c(-0.5, log2(thresh/20)-0.5))
  ymax <- ymax
  
  # Get gmts by serum group
  if(is.null(sr_group_gmt_plotdata)){
    if(titertools_gmt){
      sr_group_gmt_plotdata <- calc_titertools_sr_group_gmt(plotdata, thresh = thresh)
    } else {
      sr_group_gmt_plotdata <- calc_half_thresh_sr_group_gmt(plotdata, thresh = thresh)
    }
    
    if(order_by_gmt) {
     
      ag_order <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$logtiter, decreasing = TRUE), ] %>%
        pull(ag_name) %>%
        unique()
    }
    
  }
  
  if(fc_label) {
    if(is.null(fc_df)){
      if(titertools_fc){
        fc_df <- calc_fc_from_homologous(plotdata)
      } else {
        fc_df <- calc_fc_from_homologous_gmt(sr_group_gmt_plotdata)
      }
      
    }
    # add label
    sr_group_gmt_plotdata$fold_change <- fc_df$fc[match(interaction(sr_group_gmt_plotdata$sr_group, sr_group_gmt_plotdata$ag_name), interaction(fc_df$sr_group, fc_df$ag2))]
    sr_group_gmt_plotdata %>% 
      get_fold_change_label(.,ymax = ymax) -> sr_group_gmt_plotdata
    
  } else {
    sr_group_gmt_plotdata %>% 
      get_gmt_label(.,ymax = ymax) -> sr_group_gmt_plotdata
  }
  
  plotdata %>% 
    mutate(sr_group = factor(sr_group, levels = target_sr_groups),
           logtiter = ifelse((titertype == 2 | logtiter < log2(thresh/20)) & clamp_below_thresh, log2(thresh/20), logtiter)) -> plotdata
  
  sr_group_gmt_plotdata %>%
    mutate(sr_group = factor(sr_group, levels = target_sr_groups),
           logtiter = ifelse((all_below_thresh | logtiter < log2(thresh/20)) & clamp_below_thresh, log2(thresh/20), logtiter)) -> sr_group_gmt_plotdata
  
  if(!facet_var %in% colnames(sr_group_gmt_plotdata) & show_gmt){
    
    sr_group_gmt_plotdata %>%
      mutate(sr_group_info = str_split(sr_group, "_"),
             sr_group_big = sapply(sr_group_info, function(x) {
               x[1]
             })) -> sr_group_gmt_plotdata
    
    sr_group_gmt_plotdata[[facet_var]] <- factor(sr_group_gmt_plotdata[[facet_var]], levels = facet_levels)
  }
  
  plotdata[[facet_var]] <- factor(plotdata[[facet_var]], levels = facet_levels)
  
  if(use_ag_number){
    
    sr_group_gmt_plotdata %>%
      arrange(color_var, "logtiter") %>%
      group_by_at(vars(one_of(color_var))) %>% 
      mutate(rank = rank(desc(logtiter), ties.method = "first")) -> sr_group_gmt_plotdata
    
    sr_group_gmt_plotdata %>%
      select(all_of(c(color_var, "ag_name", "rank"))) -> ranks_by_sr_group
    
    plotdata$ag_name <- ranks_by_sr_group$rank[match(interaction(plotdata$sr_group, plotdata$ag_name), interaction(ranks_by_sr_group[[color_var]],ranks_by_sr_group$ag_name))]
    
    sr_group_gmt_plotdata %>%
      mutate(ag_name = rank) -> sr_group_gmt_plotdata
    
    ag_order <- c(1:max(sr_group_gmt_plotdata$ag_name))
    
    
  }
  
  
  sr_groups <- plotdata %>%
    pull(color_var) %>%
    unique() %>%
    as.character()
  
  color_vals <- sapply(sr_groups, function(x) blend_sr_group_colors(x, sr_group_colors))
  fill_vals <- color_vals
  
# Do the plot
  plotdata %>%
    filter(!(sr_name %in% highlight_samples)) %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_var,
        fill = color_var
      )
    ) -> gp
  
  if(show_idvls){
    gp <- gp + 
      geom_line(
        aes(group = sr_name),
        alpha = ifelse(show_gmt, 0.2, 1),
        size = 0.8
      ) + 
      geom_point(
        alpha = ifelse(show_gmt, 0.2, 1),
      ) 
  }
  
  if(!is.null(highlight_samples)){
    
    sub_samples <- plotdata[grepl(paste(highlight_samples, collapse = "|"), plotdata$sr_name),]
    
    sub_samples <- plotdata %>%
      filter(sr_name %in% highlight_samples)
    
    gp <- gp + 
      geom_line(data = sub_samples, aes(group = sr_name),
                alpha =1,
                size =1.4,
                color = "white") +
      geom_line(data = sub_samples, aes(group = sr_name),
                alpha =1,
                size = 0.7,
                color = "#00b300") +
      geom_point(data = sub_samples, aes(group = sr_name),
                 alpha =1,
                 color = "#00b300",
                 fill = "white",
                 shape = 21) 
  }
  
  if(show_gmt_conf){
    gp + 
      geom_pointrange(data = sr_group_gmt_plotdata,
                      aes(ymin = lower,
                          ymax = upper),
                      alpha = 1) -> gp
  }  
  
  if(show_gmt){
    
    gp +
      geom_line(
        data = sr_group_gmt_plotdata,
        aes(group = sr_group),
        size = 1.2,
        alpha = 1
      ) +
      geom_point(
        data = sr_group_gmt_plotdata,
        size = 2.5,
        alpha = 1
      ) +
      geom_point(
        data = sr_group_gmt_plotdata,
        size = 1.5,
        fill = "white",
        shape = 21
      ) -> gp
    
  }
    
    gp + scale_color_manual(
      values = color_vals,
      name = legend_name
    ) +
    scale_fill_manual(
      values = fill_vals,
      name = legend_name
    ) +
    scale_y_titer(
      ymin = min(c(1, log2(thresh/10))),
      ymax = ymax - 0.5
    ) + 
    scale_x_discrete(
      limits = ag_order,
      expand = expansion(add = ifelse(do_below_thresh_tile, 0, 0.5))
    ) +
  facet_wrap(
    as.formula(paste("~", facet_var)),
    ncol = facet_col,
    labeller = as_labeller(function(x) capitalize(gsub(" conv.| vax.", "", x)))
  ) + 
    labs(
      x = ifelse(use_ag_number, "Rank", "Virus variant")
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(ymin, ymax), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 14)
    ) -> gp
  
 
  # Annotate region below detectable
  if(do_below_thresh_tile){
    gp <- gp +
      annotate(
        "tile",
        x = ag_order,
        y = min(c(0, log2(thresh/20))),
        height = 1+log2(thresh/10),
        fill = "grey50",
        color = NA,
        alpha = 0.3
      )
  }
  
  if(do_ag_tiles){
    # Annotate colors for each antigen
    for (n in seq_len(numAntigens(map))) {
      gp <- gp +
        annotate(
          "tile",
          x = agNames(map)[n],
          y = ymin-0.25,
          height = 1,
          fill = agFill(map)[n],
          color = NA
        ) 
    }  
  }
  
  
  # Label fold change from homologous
  if(fc_label){
    gp <- gp +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y, label = label),
              color = "black",
              size = 2
    )
  }
  
  return(gp)
  
}



do_ag_titer_plot_fc_label <- function(map, facet_col = 3, thresh = 20, fc_label =TRUE, adj_titers = TRUE,
                                      highlight_samples = NULL,
                                      sr_group_gmt_plotdata = NULL,
                                      fc_df = NULL,
                                      target_sr_groups = NULL,
                                      ag_order = NULL,
                                      ymax = 11.5) {
  
  # Set the antigen x axis order
  # Set the antigen x axis order
  if(is.null(ag_order)){
    ag_order <- agNames(map)
  } 
  
  # Get long info
  plotdata <- long_map_info(map)
  
  if(is.null(target_sr_groups)){
    target_sr_groups <- as.character(unique(plotdata$sr_group))
  }
  
  if(adj_titers){
    plotdata$logtiter <- plotdata$logtiter_adjusted
    plotdata$titer <- plotdata$titer_adjusted
  }
  
  ymin <- min(c(-0.5, log2(thresh/20)-0.5))
  
  # Get gmts by serum group
  if(is.null(sr_group_gmt_plotdata)){
    sr_group_gmt_plotdata <- calc_titertools_sr_group_gmt(plotdata, thresh = thresh)
  }
  
  if(fc_label) {
    if(is.null(fc_df)){
      fc_df <- calc_fc_from_homologous(plotdata)
    }
    # add label
    sr_group_gmt_plotdata$fold_change <- fc_df$fc[match(interaction(sr_group_gmt_plotdata$sr_group, sr_group_gmt_plotdata$ag_name), interaction(fc_df$sr_group, fc_df$ag2))]
    sr_group_gmt_plotdata %>% 
      get_fold_change_label(.,ymax = ymax) -> sr_group_gmt_plotdata
    
  } else {
    sr_group_gmt_plotdata %>% 
      get_gmt_label(.,ymax = ymax) -> sr_group_gmt_plotdata
  }
  
  plotdata %>% 
    filter(ag_name %in% ag_order) %>%
    filter(sr_group %in% target_sr_groups) %>%
    mutate(sr_group = factor(sr_group, levels = target_sr_groups),
           logtiter = ifelse(titertype == 2 | logtiter < log2(thresh/20), log2(thresh/20), logtiter),
           ag_name = factor(ag_name, levels = ag_order)) -> plotdata
  
  sr_group_gmt_plotdata %>%
    filter(ag_name %in% ag_order) %>%
    filter(sr_group %in% target_sr_groups) %>%
    mutate(sr_group = factor(sr_group, levels = target_sr_groups),
           logtiter = ifelse(all_below_thresh | logtiter < log2(thresh/20), log2(thresh/20), logtiter),
           ag_name = factor(ag_name, levels = ag_order)) -> sr_group_gmt_plotdata
  
  
  plotdata$ag_name <- factor(plotdata$ag_name, levels = ag_order)
  sr_group_gmt_plotdata$ag_name <- factor(sr_group_gmt_plotdata$ag_name, levels = ag_order)
  

  fill_vals <- agFill(map)
  names(fill_vals) <- agNames(map)
  
  # Do the plot
  plotdata %>%
    ggplot(
      aes(
        x = sr_group,
        y = logtiter,
        color = ag_name,
        fill = sr_group
      )
    ) + 
    geom_point(
      alpha = 0.2
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 2.5,
      alpha = 1
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
     scale_color_manual(
       values = fill_vals
    ) +
    scale_fill_manual(
       values = srGroupOutline(map)
    ) +
    scale_y_titer(
      ymin = min(c(1, log2(thresh/10))),
      ymax = ymax - 0.5
    ) + 
    scale_x_discrete(
      limits = target_sr_groups,
      expand = expansion(add = 0)
    ) +
    # scale_shape_manual(
    #   values = c(
    #     "FALSE" = 21,
    #     "TRUE" = 25
    #   )
    # ) +
    # scale_size_manual(
    #   values = c(
    #     "FALSE" = 1.2,
    #     "TRUE" = 0.8
    #   )
  # ) +
  facet_wrap(
    vars(ag_name),
    ncol = facet_col
  ) + 
    labs(
      x = "Serum group"
    ) +
    coord_cartesian(
      ylim = shrinkrange(c(ymin, ymax), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none"
    ) -> gp
  
 
  # Annotate region below detectable
  gp <- gp +
    annotate(
      "tile",
      x = target_sr_groups,
      y = min(c(0, log2(thresh/20))),
      height = 1+log2(thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    )


  # Annotate colors for each antigen
  for (n in seq_along(target_sr_groups)) {
   
    gp <- gp +
      annotate(
        "tile",
        x = target_sr_groups[n],
        y = ymin-0.25,
        height = 1,
        fill = srGroupOutline(map)[target_sr_groups[n]],
        color = NA
      )
  }
   
  # Label fold change from homologous
  gp <- gp +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = sr_group, y = y, label = label),
              color = "black",
              size = 2
    )
  
  return(gp)
  
}

# Function for accounting for a sample's indivdual affect
# calculate GMT per sample, subtract serum group GMT from that to get reactivity bias
# subtract sample's reactivity bias from sample titrations
adjust_individual_effect <- function(data, dil_stepsize = 0){
  
  # calculate gmt per sample
  data %>%
    group_by(sr_name) %>%
    mutate(gmt_sample = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  # calculate gmt per serum group
  titers_variation_adjusted %>%
    group_by(sr_group) %>%
    mutate(gmt_sr_group = meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = dil_stepsize)$mean) ->titers_variation_adjusted 
  
  
  # adjust serum reactivity bias
  titers_variation_adjusted %>%
    mutate(reactivity_bias = gmt_sample - gmt_sr_group,
           logtiter = logtiter - reactivity_bias,
          # logtiter = logtiter -gmt_sample,
           titer = 2^logtiter*10) -> titers_variation_adjusted
  
  return(titers_variation_adjusted)
  
}