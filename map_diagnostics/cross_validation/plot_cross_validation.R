# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ggplot2)

source("functions/titer_lineplot_functions.R")
results <- readRDS("map_diagnostics/cross_validation/cross_validation_mapthreshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.rds")

# Read the map
map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_woXBB15conv_CH11_alpha_adj.ace")

agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}


# Set detectable results subset
detectable_results <- filter(results, measured_titer_type == 1)

# do histogram of pred - measured
mean <- round(mean(detectable_results$residual, na.rm = T),2)
sd <- round(sd(scale(detectable_results$residual, scale = F), na.rm = T), 2)

hist_diff <- ggplot(detectable_results) +
  geom_histogram(aes(x = residual), fill = "grey50", alpha = 0.8, bins = 100) +
  xlim(c(-15, 15)) +
  geom_vline(xintercept = mean, linetype = "dashed") +
  labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
  theme_bw()

ggsave(plot = hist_diff, filename = "./map_diagnostics/cross_validation/histogram_residuals.png", width = 5, height = 4, dpi = 300)

# change names
ag_order <- read.csv("data/metadata/ag_order.csv")[,3]
ag_pretty <- data.frame(
  row.names = ag_order,
  val = ag_order
)

sr_pretty <- data.frame(
  row.names = c('Wuhan vax. (single dose)','Beta vax. (single dose)','Wuhan vax. (two doses)', 'XBB.1.5 vax. (two doses)',
                'Wuhan conv.', 'alpha conv.','beta conv.',
                'gamma conv.', 'delta conv.', 'BA.1 conv.', 'BA.2.12.1 conv.', 'BA.4 conv.', 'BA.5 conv.',
                'XBB.1.5 conv.'),
  val = c('1x Wuhan vax.','1x Beta vax.','2x Wuhan vax.', '2x XBB.1.5 vax.',
          'Wuhan conv.', 'alpha conv.','beta conv.',
          'gamma conv.', 'delta conv.', 'BA.1 conv.', 'BA.2.12.1 conv.', 'BA.4 conv.', 'BA.5 conv.',
          'XBB.1.5 conv.')
)

sr_pretty$val <- capitalize(gsub(" conv.| vax.", "", sr_pretty$val))

detectable_results$ag_pretty <- factor(ag_pretty[as.character(detectable_results$ag_name),], levels = ag_pretty$val)
detectable_results$sr_pretty <- factor(sr_pretty[as.character(detectable_results$sr_group),], levels = sr_pretty$val)


    # Antigen and serum group tab
        detectable_results %>%
          ggplot(
            aes(
              x = predicted_logtiter,
              y = measured_logtiter,
              color = ag_name
            )
          ) +
          labs(x = "Predicted log2 titer",
               y = "Measured log2 titer") +
          # geom_smooth() +
          geom_point(
            alpha = 0.1
          ) +
          geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
          ) +
          scale_color_manual(
            # values = Racmacs:::srGroupOutline(map)
            values = agFillScale(map)
          ) +
          xlim(c(-10,10))+
          ylim(c(0, 10))+
          facet_grid(
            rows = vars(sr_pretty),
            cols = vars(ag_pretty),
            labeller = label_wrap_gen(multi_line = TRUE)
          ) +
          theme_bw() +
          theme(legend.position = "none",
               strip.text.x = element_text(size = 6),
               strip.text.y = element_text(size = 6))-> gp


ggsave(plot = gp, filename = "./map_diagnostics/cross_validation/scatter_pred_vs_measured.png", width = 14, height = 10, dpi = 300)
 