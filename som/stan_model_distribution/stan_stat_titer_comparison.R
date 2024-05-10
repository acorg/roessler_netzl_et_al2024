# Setup workspace
rm(list = ls())
library(labbook)
library(patchwork)
library(Racmacs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(rstatix)
source("functions/long_map_info.R")

set.seed(100)

mapdata <- readRDS("data/titer_data/stan_mapdata_assay_organism_sigma_per_dataset.rds") %>%
  mutate(modelled_titer = imputed_logtiter - sr_effect - organism_magnitude - assay_effect)

unique(mapdata$source)
datasets  <- c("NHP_Lentivirus" = "NHP", 
               "Human_Authentic Virus" = "Roessler",
               "Human_Lentivirus" = "Wilks")

# compare titers across datasets to see if they are statistically comparable
mapdata %>%
  select(ag_name, sr_group, source, modelled_titer) %>%
  mutate(dataset = datasets[source])-> data_sub

data_sub %>%
  mutate(organism = ifelse(dataset == "NHP", "NHP", "Human")) -> data_sub

data_sub %>%
  group_by(sr_group, ag_name, organism) %>%
  shapiro_test(modelled_titer) %>%
  filter(p < 0.05)

# Data not normally distributed so we will do wilcox test 
data_sub %>%
  group_by(sr_group, ag_name) %>%
  wilcox_test(modelled_titer ~ organism) %>%
  filter(p < 0.05) %>%
  add_significance() %>%
  mutate(Variable = "Modelled titer",
         Variant = ag_name,
         `Serum group` = sr_group,
         `Significance level` = p.signif) %>%
  select(`Serum group`, Variant, group1:p, `Significance level`)-> titer_diffs

write.csv(titer_diffs, file = "som/stan_model_distribution/significant_titer_differences.csv", row.names = FALSE)
