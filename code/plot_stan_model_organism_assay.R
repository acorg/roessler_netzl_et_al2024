# Setup workspace
rm(list = ls())
library(labbook)
library(tidyverse)
library(patchwork)
library(Racmacs)
library(ggh4x)
library(ggpubr)

source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")

col_vals <- c(hcl.colors(3, "Dark 3"), "black")
col_vals[1] <- "#d42b52"

names(col_vals) <- c("NHP", "Roessler", "Wilks", "Modelled GMT")

facet_labels <- c("NHP" = "NHP", 
                  "Roessler" = "Human data set 1", 
                  "Wilks" = "Human data set 2", 
                  "Modelled GMT" = "Modelled GMT",
                  "Raw titers"= "Raw titers",
                  "Adj. for serum reactivity" = "Adj. for serum reactivity",
                  "Adj. for assay effect" = "Adj. for assay effect",
                  "Adj. for organism reactivity" = "Adj. for organism reactivity",
                  "Adj. for serum, assay & organism reactivity" = "Adj. for serum, assay & organism reactivity",
                  "Wuhan vax. (two doses)" = "Wuhan (two doses)" , 
                  "Wuhan conv." = "Wuhan",
                  "alpha conv." = "Alpha", 
                  "beta conv." = "Beta",
                  "delta conv." = "Delta", 
                  "BA.1 conv." = "BA.1",
                  "GMT" = "GMT")


labelling_function <- function(x){
  facet_labels[x]
}
# Function to plot results split by source
plot_sr_group_results_sub <- function(mapdata, ag_order = NULL, facet_by = "facet_col", ag_means_all = FALSE,
                                      show_only_gmt = FALSE,
                                      show_only_raw_final = FALSE,
                                      wide_format = FALSE,
                                      sr_group_levels = c("Wuhan vax. (two doses)", 
                                                          "Wuhan conv.",
                                                          "alpha conv.", "beta conv.",
                                                          "delta conv.", "BA.1 conv.")) {
  
  # Apply adjustments
  mapdata$source <- factor(mapdata$source)
  
  mapdata <- mapdata %>%
    filter(sr_group %in% sr_group_levels)
  
  # Get gmt data
  ag_gmts <- mapdata |> distinct(ag_name, sr_group, ag_mean) 
  
  mapdata -> mapdata_subset
  
  # add black GMT line to all panels
  ag_gmts |> 
    filter(
      ag_name %in% filter(mapdata_subset, !is.na(logtiter))$ag_name
    ) %>%
    mutate(facet_col = "NHP") -> ag_gmts_sub
  
  ag_gmts_subset <- rbind(ag_gmts_sub, ag_gmts_sub %>%
                     mutate(facet_col = "Roessler"),
                     ag_gmts_sub %>%
                       mutate(facet_col = "Wilks"))
  
  if(ag_means_all){
    ag_gmts_subset <- rbind(ag_gmts_subset, ag_gmts_sub %>%
                              mutate(facet_col = "GMT"))
  }
  
  mapdata_subset |> 
    filter(
      ag_name %in% filter(mapdata_subset, !is.na(logtiter))$ag_name
    ) -> mapdata_subset
  
  # Set antigen order
  if(is.null(ag_order)){
    ag_order <- ag_gmts_subset |> arrange(-ag_mean) |> pluck("ag_name") 
  }
  
  # Do the plot
  mapdata_subset |>
    mutate(
      plot_type = "Raw titers"
    ) -> plotdata_raw
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - sr_effect ,
      plot_type = "Adj. for serum reactivity"
    ) -> plotdata_sr_effect
  
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - sr_effect - organism_magnitude - assay_effect,
      plot_type = "Adj. for serum, assay & organism reactivity"
    ) -> plotdata_sr_effect_magnitude
  
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - organism_magnitude,
      plot_type =  "Adj. for organism reactivity"
    ) -> plotdata_magnitude
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - assay_effect,
      plot_type =  "Adj. for assay effect"
    ) -> plotdata_assay
  
  
  plotdata <- bind_rows(
    plotdata_raw,
    plotdata_sr_effect,
    plotdata_sr_effect_magnitude,
    plotdata_magnitude,
    plotdata_assay
  ) |> 
    mutate(
      facet_col = source,
      plot_type = factor(
        plot_type,
        c(
          "Raw titers",
          "Adj. for serum reactivity",
          "Adj. for assay effect",
          "Adj. for organism reactivity",
          "Adj. for serum, assay & organism reactivity"
        )
      )
    )
  
  # calculate gmt now
  plotdata %>%
    group_by(source, plot_type, ag_name, sr_group) %>%
    mutate(logtiter = ifelse(titer == "<16" & plot_type == "Raw titers", log2(8/10), logtiter),
           logtiter = ifelse(titer == "<20" & plot_type == "Raw titers", log2(10/10), logtiter),
            lower = Rmisc::CI(na.omit(logtiter))["lower"],
           upper = Rmisc::CI(na.omit(logtiter))["upper"],
           logtiter = mean(logtiter, na.rm =TRUE)) %>%
    ungroup() %>%
    mutate(facet_col = "GMT",
           sr_name = paste0("GMT", source)) %>%
    select(source, plot_type, ag_name, sr_group, sr_name, lower, upper, logtiter, facet_col) %>%
    unique() %>%
    mutate(alpha_val = "GMT")-> plotdata_source
  
  if(facet_by == "source"){
    plotdata_source <- rbind(plotdata_source, plotdata_source %>%
                               mutate(facet_col = source)) 
  }

  plotdata <- plyr::rbind.fill(plotdata %>%
                                 mutate(alpha_val = "Organism"), plotdata_source)
 
 plotdata %>%
  mutate(alpha_val = ifelse(grepl("GMT", sr_name), "GMT", "Organsim"),
         facet_col = factor(facet_col, levels = c("NHP", "Roessler","Wilks", "GMT")),
         sr_group = factor(sr_group, levels = sr_group_levels)) -> plotdata 
 
 ag_gmts_subset <- ag_gmts_subset %>%
   filter(sr_group %in% unique(plotdata$sr_group)) %>%
   mutate(facet_col = factor(facet_col, levels = c("NHP", "Roessler","Wilks", "GMT")),
          sr_group = factor(sr_group, levels = sr_group_levels))
 
 if(show_only_gmt){
   plotdata <- plotdata %>%
     filter(plot_type %in% c("Raw titers",
                             "Adj. for serum, assay & organism reactivity")) %>%
     filter(facet_col == "GMT")
            
   ag_gmts_subset <- ag_gmts_subset %>%
     filter(facet_col == "GMT") %>%
     mutate(facet_col = as.character(facet_col))
 }
 
 if(show_only_raw_final){
   plotdata <- plotdata %>%
     filter(plot_type %in% c("Raw titers",
                             "Adj. for serum, assay & organism reactivity"))
 }
 
 plotdata %>%
   mutate(tile_height = ifelse(source == "Roessler", log2(16/10), log2(20/10))) %>%
   select(sr_group, plot_type, source, tile_height, ag_name, facet_col) %>%
   mutate(ag_name = "Wuhan",
          source = ifelse(facet_col == "GMT", "NHP", source),
          tile_height = ifelse(facet_col == "GMT",log2(20/10), tile_height)) %>%
   unique() -> tile_frame
 
 plotdata |> 
    ggplot(
      aes(
        x = ag_name,
        color = source,
        alpha = alpha_val
      )
    ) + 
    geom_ribbon(
      aes(ymin = lower,
          ymax = upper,
          group = sr_name,
          fill = source),
      color = NA,
      alpha = 0.15
    ) + 
    geom_line(
      aes(
        y = logtiter,
        group = sr_name)
    ) + 
    geom_line(
      aes(group = sr_name,
          y = logtiter),
      data = plotdata |> filter(!is.na(logtiter)),
      linetype = "dotted"
    ) + 
    scale_alpha_discrete(range = c(0.4, 1), 
                         limits = c("Organsim", "GMT"),
                         guide = "none") +
    geom_rect(
      data = tile_frame,
      xmin = -Inf,
      xmax = Inf,
      aes(ymax = tile_height), 
      ymin = -Inf,
      fill = "black",
      color = NA,
      alpha = 0.1
    ) +
    scale_x_discrete(
      limits = ag_order,
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = function(x) 2^x*10,
                       breaks = seq(-1,14,2)) +
    scale_color_manual(values = col_vals,
                       labels = facet_labels[names(col_vals)],
                       name = "Dataset", drop = FALSE) +
    scale_fill_manual(values = col_vals, 
                      labels = facet_labels[names(col_vals)],
                      name = "Dataset", drop = FALSE) + 
    coord_cartesian(
      ylim = c(-2, 14)
    ) +
    labs(
      x = "Virus variant",
      y = "Titer"
    ) + 
    titerplot_theme() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5),
          panel.spacing.x = unit(1, "lines"),
          axis.text = element_text(size = 8),
          strip.text.y= element_text(size = 9),
          strip.text.x= element_text(size = 9))-> gp
 
 if(show_only_gmt){
  
   
   gp <- gp + 
     geom_line(
       data = ag_gmts_subset,
       mapping = aes(
         y = ag_mean,
         group = "average"
       ),
       color = "black",
       alpha =1,
       linetype = "solid"
     )
   
   if(wide_format){
     gp <- gp + facet_nested(plot_type ~ sr_group,
                             drop = TRUE,
                             labeller = as_labeller(labelling_function,
                                                    label_wrap_gen(multi_line = TRUE)),
                             nest_line = element_line(color = "grey20")
     ) 
     
   } else {
     
     gp <- gp + facet_nested(sr_group ~ plot_type,
                             drop = TRUE,
                             labeller = as_labeller(labelling_function,
                                                    label_wrap_gen(multi_line = TRUE)),
                             nest_line = element_line(color = "grey20")
     ) 
   }
     
   
 } else {
   
   gp <- gp + 
     geom_line(
       data = ag_gmts_subset,
       mapping = aes(
         y = ag_mean,
         group = "average"
       ),
       color = "black",
       alpha =1
     ) +
     facet_nested(plot_type ~ sr_group + facet_col,
                drop = TRUE,
                labeller = as_labeller(labelling_function,
                                                   label_wrap_gen(multi_line = TRUE)),
                nest_line = element_line(color = "grey20")
   ) 
   
 }
  
  
  return(gp)
}

figure_dir <- file.path("figures", "titerplots_stan")
if(!exists(figure_dir)){
  dir.create(figure_dir)
}



source_names <- c("NHP_Lentivirus" =  "NHP",
                  "Human_Authentic Virus" = "Roessler",
                  "Human_Lentivirus" = "Wilks")

for(data in c( "stan_mapdata_assay_organism_sigma_per_dataset")){
  
  
  mapdata <- readRDS(file.path("data/titer_data/", paste0(data, ".rds")))
  
  mapdata$source <- source_names[mapdata$source]
  
  
  wrapped_sub <- plot_sr_group_results_sub(mapdata, ag_order = unique(mapdata$ag_name), facet_by = "facet_col", ag_means_all = TRUE,
                                           sr_group_levels = c("Wuhan vax. (two doses)",
                                                               "Wuhan conv.",
                                                               "alpha conv."))
  
  ggsave(file.path(figure_dir, paste0("facet_ag_order_titerplots_stan_gmt_",data,"_set1.png")), plot = wrapped_sub, width = 14, height = 8)
  
  
  wrapped_sub <- plot_sr_group_results_sub(mapdata, ag_order = unique(mapdata$ag_name), facet_by = "facet_col", ag_means_all = TRUE,
                                           sr_group_levels = c("beta conv.",
                                                               "delta conv.",
                                                               "BA.1 conv."))
  
  ggsave(file.path(figure_dir, paste0("facet_ag_order_titerplots_stan_gmt_",data,"_set2.png")), plot = wrapped_sub, width = 14, height = 8)
  
  
   
   wrapped_sub <- plot_sr_group_results_sub(mapdata, ag_order = unique(mapdata$ag_name), facet_by = "facet_col", ag_means_all = TRUE,
                                            show_only_gmt = TRUE,
                                            wide_format = TRUE,
                                              sr_group_levels = c("Wuhan vax. (two doses)",
                                                                  "Wuhan conv.",
                                                                  "alpha conv.",
                                                                  "beta conv.",
                                                                  "delta conv.",
                                                                  "BA.1 conv.")) + 
     theme(legend.position = "top") + 
     theme(axis.text = element_text(size = 8))
   
  
   ggsave(file.path(figure_dir, paste0("facet_ag_order_titerplots_stan_only_gmt_",data,"_Fig3.png")), wrapped_sub, width = 12, height = 6)
   
}

mapdata %>%
  mutate(source = facet_labels[source]) %>%
  filter(titertype != 0) %>%
  select(source, ag_name, sr_group, sr_name) %>%
  group_by(sr_group, ag_name) %>%
  count(source) %>%
  pivot_wider(names_from = "source", values_from = "n") %>%
  mutate(sr_group = facet_labels[sr_group])-> sr_table

sr_table[is.na(sr_table)] <- 0

write.csv(sr_table, "data/titer_data/bayesian_comp_titer_table.csv", row.names = FALSE)

#------------------------- Show sigma error per dataset
mapdata %>%
  mutate(final_titer = imputed_logtiter - sr_effect - organism_magnitude  - assay_effect) %>%
  group_by(source) %>%
  summarize(measured = sd(logtiter, na.rm = TRUE),
            "from imputed titers" = sd(final_titer, na.rm = TRUE),
            fitted = dataset_error) %>%
  unique() %>% 
  pivot_longer(cols = c("fitted", "measured", "from imputed titers"), names_to = "Data", values_to = "Sigma Error") %>%
  mutate("Dataset" = facet_labels[source]) -> sigma

# plot sigma
sigma %>%
  ggplot(aes(x = Dataset, y = `Sigma Error`, colour = Data)) +
  geom_point() + 
  theme_bw() + 
  ylim(c(0, NA)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) -> sigma_error_plot

ggsave(file.path(figure_dir, "sigma_error_plot.png"), sigma_error_plot, dpi = 300, width = 5, height = 4)

#------------------------ Calculate organism and assay magnitudes
mapdata %>%
  select(organism, organism_magnitude, assay, assay_effect) %>%
  unique() %>%
  pivot_wider(names_from = organism, values_from = organism_magnitude) %>%
  mutate("Human - NHP reactivity" = Human - NHP,
         "Lentivirus - Authentic virus effect" = assay_effect[assay == "Lentivirus"] - assay_effect) %>%
  select(all_of(c("Human - NHP reactivity", "Lentivirus - Authentic virus effect" ))) %>%
  pivot_longer(cols = c("Human - NHP reactivity", "Lentivirus - Authentic virus effect" ), names_to = "Variable", values_to = "Estimate") %>%
  filter(Estimate > 0) %>%
  mutate(Estimate = round(2^Estimate, 2))-> difference_esitmate


#------------------ Do SOM Titer correlation plot
do_correlation_plot <- function(data, x_data, y_data){
  
  ggscatter(data, x = x_data, y= y_data,
            add = "reg.line",  # Add regression line
            add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE # Add confidence interval
  ) + stat_cor() + 
    theme_bw() + 
    scale_y_continuous(labels = function(x) round(2^x*10),
                       limits = c(0,10),
                       name = paste(y_data, "GMT")) + 
    scale_x_continuous(labels = function(x) round(2^x*10),
                       limits = c(0,10),
                       name = paste(x_data, "GMT"))  -> p
  
  return(p)
  
}


mapdata %>%
  mutate(final_titer = imputed_logtiter - sr_effect - organism_magnitude  - assay_effect) %>%
  group_by(ag_name, sr_group, source) %>%
  mutate(GMT = mean(final_titer),
         source = facet_labels[source]) %>%
  select(source, ag_name, sr_group, GMT) %>%
  unique() %>%
  pivot_wider(names_from = source, values_from = GMT)-> data_final

final_nhp_wilks <- do_correlation_plot(data_final, "NHP", "Human data set 2")
final_nhp_roessler <- do_correlation_plot(data_final, "NHP", "Human data set 1")
final_human <- do_correlation_plot(data_final, "Human data set 1", "Human data set 2")

final_plot <- final_nhp_roessler + final_nhp_wilks + final_human 

mapdata %>%
  mutate(final_titer = imputed_logtiter,
         source = facet_labels[source]) %>%
  group_by(ag_name, sr_group, source) %>%
  mutate(GMT = mean(final_titer)) %>%
  select(source, ag_name, sr_group, GMT) %>%
  unique() %>%
  pivot_wider(names_from = source, values_from = GMT)-> data_final

final_nhp_wilks <- do_correlation_plot(data_final, "NHP", "Human data set 2")
final_nhp_roessler <- do_correlation_plot(data_final, "NHP", "Human data set 1")
final_human <- do_correlation_plot(data_final, "Human data set 1", "Human data set 2")

before_plot <- final_nhp_roessler + final_nhp_wilks + final_human 


comb <- before_plot / final_plot + plot_annotation(tag_levels = 'A')
ggsave(file.path(figure_dir, "gmt_corr_plot.png"), comb, dpi = 300, width = 12, height = 8)

