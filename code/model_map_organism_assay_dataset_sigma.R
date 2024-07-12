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

srExtraMatrix <- function (map, optimization_number = 1) 
{
  matrix(srExtra(map), nrow = numAntigens(map), ncol = numSera(map), 
         byrow = T)
}

# Function to model the map
model_map <- function(map, 
                      ag_mean_prior_mean = 6, 
                      ag_mean_prior_sigma = 20,
                      sigma_prior_alpha = 2,
                      sigma_prior_beta = 5,
                      sigma_prior_sr_effect = 10,
                      sigma_prior_dataset_bias = 0.4, 
                      sigma_prior_organism_magnitude = 2,
                      sigma_prior_assay_effect = 2,
                      fitted_dataset_magnitude = NULL) {
  
  # Fetch data
  titer_matrix <- titerTable(map)
  ag_name_matrix <- Racmacs:::agNameMatrix(map)
  sr_extra_matrix <- srExtraMatrix(map)
  sr_name_matrix <- Racmacs:::srNameMatrix(map)
  sr_group_matrix <- Racmacs:::srGroupMatrix(map)
  assay_matrix <- apply(sr_extra_matrix, 1:2, function(x){
    strsplit(x, "_")[[1]][2]
  })
  organism_matrix <- apply(sr_extra_matrix, 1:2, function(x){
    strsplit(x, "_")[[1]][1]
  })
  
  
  # Vectorise data
  titers <- as.vector(titer_matrix)
  ag_name <- as.vector(ag_name_matrix)
  sr_name <- as.vector(sr_name_matrix)
  sr_extra <- as.vector(sr_extra_matrix)
  sr_group <- as.vector(sr_group_matrix)
  assay <- as.vector(assay_matrix)
  organism <- as.vector(organism_matrix)
  
  # Remove NA values
  na_titers <- titers == "*" | titers == "."
  titers <- titers[!na_titers]
  ag_name <- ag_name[!na_titers]
  sr_name <- sr_name[!na_titers]
  sr_extra <- sr_extra[!na_titers]
  sr_group <- sr_group[!na_titers]
  assay <- assay[!na_titers]
  organism <- organism[!na_titers]
  
  # Factor variables
  ag_name <- factor(ag_name)
  sr_name <- factor(sr_name)
  sr_extra <- factor(sr_extra)
  sr_group <- factor(sr_group)
  assay <- factor(assay)
  organism <- factor(organism)
  
  # Calculate titer limits
  titer_lims <- titertools:::calc_titer_lims(titers, dilution_stepsize = dilutionStepsize(map))
  
  # Build model
  model <- cmdstanr::cmdstan_model("models/titer_comparison_organism_assay_set_sigma_dataset.stan")
  
  
  if(is.null(fitted_dataset_magnitude)){
    fitted_dataset_magnitude <- rep(0, nlevels(sr_extra))
  } 
  
  # # Get long version of the map data
  mapdata <- long_map_info(map) |> rename(source = sr_extra)
  
  #calculate sigma error
  
  # Set input data
  stan_data <- list(
    N = length(titers),
    N_ags = nlevels(ag_name),
    N_srs = nlevels(sr_name),
    N_sr_groups = nlevels(sr_group),
    N_datasets = nlevels(sr_extra),
    N_organisms = nlevels(organism) ,
    N_assays = nlevels(assay),
    upper_lims = titer_lims$max_titers,
    lower_lims = titer_lims$min_titers,
    ags = as.numeric(ag_name),
    srs = as.numeric(sr_name),
    sr_groups = as.numeric(sr_group),
    datasets = as.numeric(sr_extra),
    organisms = as.numeric(organism),
    assays = as.numeric(assay),
    ag_mean_prior_mean = ag_mean_prior_mean,
    ag_mean_prior_sigma = ag_mean_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta,
    sigma_prior_sr_effect = sigma_prior_sr_effect,
    sigma_prior_organism_magnitude = sigma_prior_organism_magnitude,
    sigma_prior_assay_effect = sigma_prior_assay_effect,
    organism_magnitude_mean = rep(0, nlevels(organism)),
    assay_effect_mean = rep(0, nlevels(assay))
  )
  
  if(is.null(fitted_dataset_magnitude)){
    fitted_dataset_magnitude <- rep(0, nlevels(sr_extra))
  } 
  
  ag_means <- matrix(NA, nlevels(ag_name), nlevels(sr_group))
  ag_level_matrix <- matrix(levels(ag_name), nlevels(ag_name), nlevels(sr_group), byrow = F)
  sr_group_level_matrix <- matrix(levels(sr_group), nlevels(ag_name), nlevels(sr_group), byrow = T)
  ag_means[] <- vapply(seq_along(ag_means), \(i){
    mean(titer_lims$log_titers[ag_name == ag_level_matrix[i] & sr_group == sr_group_level_matrix[i]])
  }, numeric(1))
  ag_means[is.na(ag_means)] <- stan_data$ag_mean_prior_mean
  
  
  stan_init <- list(
    ag_means = ag_means,
    sr_effects = rep(0, nlevels(sr_name)),
    dataset_bias = array(0, c(nlevels(ag_name), nlevels(sr_group), nlevels(sr_extra))),
    sigma_error = rep(1, nlevels(sr_extra)),
   assay_effects = rep(0, nlevels(assay)),
   organism_magnitude = rep(0, nlevels(organism))
  )

  
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    save_warmup = FALSE,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 1000,
    max_treedepth = 20,
    refresh = 500 # print update every 500 iters
  )
  
  # Optimize model
  model_optim <- model$optimize(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init
    )
  )
  
 model_optimized <- fit
  
  # Get generated quantities
  model_pars <- model_optimized$summary()[-1,]
  model_par_list <- as.list(model_pars$mean)
  names(model_par_list) <- model_pars$variable
  
  
  model_optimized$summary("imputed_logtiters") |>
    mutate(sr_name = sr_name,
           ag_name = ag_name) |> 
    rename(
      imputed_logtiter = mean
    ) |>
    select(
      sr_name,
      ag_name,
      imputed_logtiter
    ) -> imputed_logtiters
  
  # Organize results
  sr_effects <- model_optimized$summary("sr_effects") |> mutate(sr_name = levels(sr_name)) |> select(sr_name, mean) |> rename(sr_effect = mean)
  organism_magnitude <- model_optimized$summary("organism_magnitude") |> mutate(organism = levels(organism)) |> select(organism, mean) |> rename(organism_magnitude = mean)
  assay_effect <- model_optimized$summary("assay_effects") |> mutate(assay = levels(assay)) |> select(assay, mean) |> rename(assay_effect = mean)
  
  ag_gmts <- model_optimized$summary("ag_means") |>
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) |> 
    separate(
      indices,
      sep = ",",
      into = c("ag", "group")
    ) |> 
    mutate(
      ag_name = as.character(levels(ag_name)[as.numeric(ag)]),
      sr_group = as.character(levels(sr_group)[as.numeric(group)])
    ) |> 
    select(
      ag_name, sr_group, mean
    ) |> rename(
      ag_mean = mean
    )
  
  dataset_error <- model_optimized$summary("sigma_error") |> 
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) |> 
    separate(
      indices,
      sep = ",",
      into = c("dataset")
    ) |> 
    mutate(
      source = as.character(levels(sr_extra)[as.numeric(dataset)])
    ) |> 
    select(
      source, mean
    ) |> 
    rename(
      dataset_error = mean
    )
  
  
  mapdata$assay <- sapply(mapdata$source, function(x){
    strsplit(x, "_")[[1]][2]
  })
  
  mapdata$organism <- sapply(mapdata$source, function(x){
    strsplit(x, "_")[[1]][1]
  })
  
  # Bind in estimates
  mapdata |> 
    left_join(sr_effects, by = "sr_name") |> 
    left_join(organism_magnitude, by = "organism") |> 
    left_join(assay_effect, by = "assay") |> 
    left_join(dataset_error, by = c("source")) |> 
    left_join(imputed_logtiters, by = c("ag_name", "sr_name")) |> 
    left_join(ag_gmts, by = c("ag_name", "sr_group")) -> mapdata
  
 
  return(list("mapdata" = mapdata, "fit" = fit, "optim" = model_optim))
  
}

nhp_map <- read.acmap("./data/maps/map_threshold20_all_ags_singleTP_woXBBBQ11conv_alpha_adj.ace")
roessler_map <- read.acmap("./data/maps/roessler_et_al2023_subNHP.ace")
duke_map <- read.acmap("./data/maps/Wilks_et_al_subNHP.ace")

srExtra(nhp_map) <- rep("NHP_Lentivirus", length(srExtra(nhp_map)))
srExtra(roessler_map) <- rep("Human_Authentic Virus", length(srExtra(roessler_map)))
srExtra(duke_map) <- rep("Human_Lentivirus", length(srExtra(duke_map)))

common_ags <- intersect(intersect(agNames(nhp_map), agNames(roessler_map)), intersect(agNames(nhp_map), agNames(duke_map)))
common_srs <- intersect(intersect(srGroups(nhp_map), srGroups(roessler_map)), intersect(srGroups(nhp_map), srGroups(duke_map)))

merged_map <- mergeMaps(list(nhp_map, roessler_map, duke_map))
merged_map <- subsetMap(merged_map, antigens = common_ags, sera = srNames(merged_map)[as.character(srGroups(merged_map)) %in% common_srs])

# Do Model
model_both <- model_map(merged_map, ag_mean_prior_mean = 4, sigma_prior_beta = 1.5, sigma_prior_alpha = 3, sigma_prior_sr_effect = 4, sigma_prior_dataset_bias = 1e-3,
                        sigma_prior_organism_magnitude = 4, sigma_prior_assay_effect = 4) 

saveRDS(
  object = model_both$mapdata,
  "data/titer_data/stan_mapdata_assay_organism_sigma_per_dataset.rds"
)

# Print diagnose
model_both$fit$cmdstan_diagnose()

#--------------------------------------------  Stat comparison of assay and organism effect --------------- 

# Set up of variables

draws_assay <- model_both$fit$draws(variables = "assay_effects", format = "df")
draws_organism <- model_both$fit$draws(variables = "organism_magnitude", format = "df")

map <- merged_map
sr_extra_matrix <- srExtraMatrix(map)
titer_matrix <- titerTable(map)
ag_name_matrix <- Racmacs:::agNameMatrix(map)
sr_extra_matrix <- srExtraMatrix(map)
sr_name_matrix <- Racmacs:::srNameMatrix(map)
sr_group_matrix <- Racmacs:::srGroupMatrix(map)
assay_matrix <- apply(sr_extra_matrix, 1:2, function(x){
  strsplit(x, "_")[[1]][2]
})
organism_matrix <- apply(sr_extra_matrix, 1:2, function(x){
  strsplit(x, "_")[[1]][1]
})

# Vectorise data
titers <- as.vector(titer_matrix)
ag_name <- as.vector(ag_name_matrix)
sr_name <- as.vector(sr_name_matrix)
sr_extra <- as.vector(sr_extra_matrix)
sr_group <- as.vector(sr_group_matrix)
assay <- as.vector(assay_matrix)
organism <- as.vector(organism_matrix)

# Remove NA values
na_titers <- titers == "*" | titers == "."
titers <- titers[!na_titers]
ag_name <- ag_name[!na_titers]
sr_name <- sr_name[!na_titers]
sr_extra <- sr_extra[!na_titers]
sr_group <- sr_group[!na_titers]
assay <- assay[!na_titers]
organism <- organism[!na_titers]

# Factor variables
ag_name <- factor(ag_name)
sr_name <- factor(sr_name)
sr_extra <- factor(sr_extra)
sr_group <- factor(sr_group)
assay <- factor(assay)
organism <- factor(organism)

sr_extra_levels <- levels(sr_extra)
sr_group_levels <- levels(sr_group)
ag_means <- matrix(NA, nlevels(ag_name), nlevels(sr_group))

ag_level_matrix <- matrix(levels(ag_name), nlevels(ag_name), nlevels(sr_group), byrow = F)
sr_group_level_matrix <- matrix(levels(sr_group), nlevels(ag_name), nlevels(sr_group), byrow = T)

for(r in 1:nrow(ag_means)){
  for(c in 1:ncol(ag_means)){
    ag_means[r,c] <- paste0(ag_level_matrix[r, c], ", ", sr_group_level_matrix[r, c])
  }
}

sr_effects <- levels(sr_name)

datasets <- c("Human_Authentic Virus" = "Roessler", "Human_Lentivirus" = "Wilks","NHP_Lentivirus"= "NHP")
sigma_error <- datasets[levels(sr_extra)]

#----------------------------- Do the comparison ---------------------------
colnames(draws_assay)[1:2] <- levels(assay)
colnames(draws_organism)[1:2] <- levels(organism)


draws_assay %>%
  pivot_longer(cols = levels(assay), names_to = "variable", values_to = 'estimate') -> draws_long

draws_long %>%
  group_by(variable) %>%
  summarize(mean_effect = mean(estimate))

model_both$fit$summary("assay_effects")

summary_df <- model_both$fit$summary("assay_effects") %>%
  select(variable, mean, median) %>%
  pivot_longer(cols = c("mean", "median"), names_to = "Data", values_to = "estimate") %>%
  rbind(model_both$optim$summary("assay_effects") %>%
          mutate(Data = "Optimization")
  ) %>%
  mutate(variable = ifelse(grepl("1", variable), levels(assay)[1], levels(assay)[2]))

draws_long %>%
  filter(variable == "Authentic Virus") %>%
  shapiro_test(estimate)

draws_long %>%
  filter(variable == "Lentivirus") %>%
  shapiro_test(estimate)

draws_long %>%
  t_test(estimate ~ variable) %>%
  add_significance()

draws_long %>%
  ggplot(aes(x = variable, y = estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_signif(comparisons = list(c("Authentic Virus", "Lentivirus")),
              map_signif_level = TRUE,
              test = "t.test") +
  scale_y_continuous(labels = function(x) round(2^x, 3),
                     name = "Estimate",
                     breaks = seq(-10, 10, 2)) +
  geom_point(data = summary_df, aes(color = Data), shape = 1, position = position_dodge(width = 0.3)) +
  xlab("Assay effect") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) -> assay_plot


##-------------------------  for organism
draws_organism %>%
  pivot_longer(cols = levels(organism), names_to = "variable", values_to = 'estimate') -> long_organism

long_organism %>%
  filter(variable == "Human") %>%
  shapiro_test(estimate)

long_organism %>%
  filter(variable == "NHP") %>%
  shapiro_test(estimate)

long_organism %>%
  t_test(estimate ~ variable) %>%
  add_significance()


summary_df <- model_both$fit$summary("organism_magnitude") %>%
  select(variable, mean, median) %>%
  pivot_longer(cols = c("mean", "median"), names_to = "Data", values_to = "estimate") %>%
  rbind(model_both$optim$summary("organism_magnitude") %>%
          mutate(Data = "Optimization")
  ) %>%
  mutate(variable = ifelse(grepl("1", variable), levels(organism)[1], levels(organism)[2]))

long_organism %>%
  ggplot(aes(x = variable, y = estimate)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_signif(comparisons = list(c("Human", "NHP")),
              map_signif_level = TRUE, # function(p) sprintf("p = %.2g", p),
              test = "t.test") + 
  scale_y_continuous(labels = function(x) round(2^x, 3),
                     name = "Estimate",
                     breaks = seq(-10, 10, 2)) +
  geom_point(data = summary_df, aes(color = Data), shape = 1, position = position_dodge(width = 0.3)) +
  xlab("Organism effect") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))-> organism_plot

organism_plot + theme(legend.position = "none") + assay_plot + theme(legend.position = "none") -> comb

ggsave("som/stan_model_distribution/assay_organism_distribution.png", comb, width = 8, height = 4.5, dpi = 300)

#---------------  do it for serum and antigen
plot_violin_plot_distribution <- function(model_both, target_var, xlab_label,
                                          lower_break = -10, upper_break = 10,
                                          steps = 2){
  
  draws_df <- model_both$fit$draws(variables = target_var, format = "df")

  sub_df <- draws_df[, grepl(target_var, colnames(draws_df))]
  
  sub_df %>%
    pivot_longer(cols = colnames(sub_df), names_to = "variable", values_to = "estimate") -> sub_df_long
  
  sub_df_long$variable <- sapply(sub_df_long$variable, function(x) eval(parse(text = x)))
  
  summary_df <- model_both$fit$summary(target_var) %>%
    select(variable, mean, median) %>%
    pivot_longer(cols = c("mean", "median"), names_to = "Data", values_to = "estimate") %>%
    rbind(model_both$optim$summary(target_var) %>%
            mutate(Data = "Optimization")
    )
  
  summary_df$variable <- sapply(summary_df$variable, function(x) eval(parse(text = x)))
  
  if(target_var == "sr_effects"){
    sr_names <- c(1:length(sr_effects))
    names(sr_names) <- sr_effects
    
    summary_df$variable <- factor(sr_names[summary_df$variable], levels = sr_names)
    sub_df_long$variable <- factor(sr_names[sub_df_long$variable], levels = sr_names)
  } 
  
  if(target_var == "sigma_error"){
    facet_labels <- c("NHP" = "NHP", 
                      "Roessler" = "Human data set 1", 
                      "Wilks" = "Human data set 2")
    
    summary_df$variable <- facet_labels[summary_df$variable]
    sub_df_long$variable <- facet_labels[sub_df_long$variable]
  }
  
  sub_df_long %>%
    ggplot(aes(x = variable, y = estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
    scale_y_continuous(labels = function(x) round(2^x, 3),
                       breaks = seq(lower_break, upper_break, steps),
                       name = "Estimate") +
    geom_point(data = summary_df, aes(color = Data), shape = 1, position = position_dodge(width = 0.3)) +
    xlab(xlab_label) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> p
  
  return(p)
  
}


ag_plot <- plot_violin_plot_distribution(model_both, "ag_means", "Variant, Serum group") + 
  theme(legend.position = "none")
sigma_plot <- plot_violin_plot_distribution(model_both, "sigma_error", "Sigma per dataset",
                                            lower_break = 0, upper_break = 2, steps = 0.5) + 
  theme(legend.position = "none")

sr_plot <- plot_violin_plot_distribution(model_both, "sr_effects", "Serum",
                                            lower_break = -10, upper_break = 10, steps = 2) + 
  theme(legend.position = "top")

sr_plot / (ag_plot + sigma_plot + comb + plot_layout(widths = c(3, 1, 2))) + plot_annotation(tag_levels = 'A') -> ag_sr_sigma

ggsave("som/stan_model_distribution/ag_sigma_sr_effect_distribution_test.png", ag_sr_sigma, width = 20, height = 12, dpi = 300)

sr_name_df <- data.frame("Serum number" = c(1:length(sr_effects)),
                         "Serum name" = sr_effects,
                         "Serum group" = srGroups(map)[match(srNames(map), sr_effects)],
                         "Dataset" = datasets[srExtra(map)[match(srNames(map), sr_effects)]])

write.csv(sr_name_df, file = "som/stan_model_distribution/sr_name_table.csv", row.names = FALSE)
