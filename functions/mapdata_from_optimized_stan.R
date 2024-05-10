# Function to model the map
mapdata_from_optimized_stan <- function(map, 
                                   model_optimized,
                                   model_file = "models/titer_comparison_magnitude_no_bias.stan",
                      ag_mean_prior_mean = 6, # was 6
                      ag_mean_prior_sigma = 20,
                      sigma_prior_alpha = 2,
                      sigma_prior_beta = 5,
                      sigma_prior_sr_effect = 10,
                      sigma_prior_dataset_bias = 0.4, # was 0.4,
                      sigma_prior_dataset_magnitude = 2,
                      fitted_dataset_magnitude = NULL) {
  
  # Fetch data
  titer_matrix <- titerTable(map)
  ag_name_matrix <- Racmacs:::agNameMatrix(map)
  sr_extra_matrix <- srExtraMatrix(map)
  sr_name_matrix <- Racmacs:::srNameMatrix(map)
  sr_group_matrix <- Racmacs:::srGroupMatrix(map)
  
  # Vectorise data
  titers <- as.vector(titer_matrix)
  ag_name <- as.vector(ag_name_matrix)
  sr_name <- as.vector(sr_name_matrix)
  sr_extra <- as.vector(sr_extra_matrix)
  sr_group <- as.vector(sr_group_matrix)
  
  # Remove NA values
  na_titers <- titers == "*" | titers == "."
  titers <- titers[!na_titers]
  ag_name <- ag_name[!na_titers]
  sr_name <- sr_name[!na_titers]
  sr_extra <- sr_extra[!na_titers]
  sr_group <- sr_group[!na_titers]
  
  # Factor variables
  ag_name <- factor(ag_name)
  sr_name <- factor(sr_name)
  sr_extra <- factor(sr_extra)
  sr_group <- factor(sr_group)
  
  # Calculate titer limits
  titer_lims <- titertools:::calc_titer_lims(titers, dilution_stepsize = dilutionStepsize(map))
  
  # Build model
  model <- cmdstanr::cmdstan_model(model_file)
  
  
  if(is.null(fitted_dataset_magnitude)){
    fitted_dataset_magnitude <- rep(0, nlevels(sr_extra))
  } 
  
  # Set input data
  stan_data <- list(
    N = length(titers),
    N_ags = nlevels(ag_name),
    N_srs = nlevels(sr_name),
    N_sr_groups = nlevels(sr_group),
    N_datasets = nlevels(sr_extra),
    upper_lims = titer_lims$max_titers,
    lower_lims = titer_lims$min_titers,
    ags = as.numeric(ag_name),
    srs = as.numeric(sr_name),
    sr_groups = as.numeric(sr_group),
    datasets = as.numeric(sr_extra),
    ag_mean_prior_mean = ag_mean_prior_mean,
    ag_mean_prior_sigma = ag_mean_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta,
    sigma_prior_sr_effect = sigma_prior_sr_effect,
#    sigma_prior_dataset_bias =sigma_prior_dataset_bias,
    sigma_prior_dataset_magnitude = sigma_prior_dataset_magnitude, # was 10. The higher here, the more variation is allowed (penalty = value/sigma). So less tringent fitting, eg we come out here with ag_mean = prior value
    dataset_magnitude_mean = fitted_dataset_magnitude  
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
    dataset_magnitude = fitted_dataset_magnitude,
    sigma_error = matrix(1, nlevels(sr_extra), nlevels(sr_group))
  )
  

  # Get generated quantities
  model_pars <- model_optimized$summary()[-1,]
  model_par_list <- as.list(model_pars$estimate)
  names(model_par_list) <- model_pars$variable
  
  imputed_logtiters <- model$generate_quantities(
    fitted_params = do.call(posterior::draws_matrix, model_par_list),
    data = stan_data
  )$summary()

  imputed_logtiters |>
    mutate(
      sr_name = sr_name,
      ag_name = ag_name
    ) |>
    rename(
      imputed_logtiter = mean
    ) |>
    select(
      sr_name,
      ag_name,
      imputed_logtiter
    ) -> imputed_logtiters
  
  # model_optimized$summary("imputed_logtiters") |>
  #   mutate(sr_name = sr_name,
  #          ag_name = ag_name) |> 
  #   rename(
  #     imputed_logtiter = estimate
  #   ) |>
  #   select(
  #     sr_name,
  #     ag_name,
  #     imputed_logtiter
  #   ) -> imputed_logtiters
  
  # Organize results
  sr_effects <- model_optimized$summary("sr_effects") |> mutate(sr_name = levels(sr_name)) |> select(sr_name, estimate) |> rename(sr_effect = estimate)
  dataset_magnitude <- model_optimized$summary("dataset_magnitude") |> mutate(source = levels(sr_extra)) |> select(source, estimate) |> rename(dataset_magnitude = estimate)
  
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
      ag_name, sr_group, estimate
    ) |> rename(
      ag_mean = estimate
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
      into = c("dataset", "group")
    ) |> 
    mutate(
      source = as.character(levels(sr_extra)[as.numeric(dataset)]),
      sr_group = as.character(levels(sr_group)[as.numeric(group)])
    ) |> 
    select(
      source, sr_group, estimate
    ) |> 
    rename(
      dataset_error = estimate
    )
  
  # dataset_bias <- model_optimized$summary("dataset_bias") |> 
  #   mutate(
  #     indices = gsub("^.*\\[", "[", variable),
  #     indices = gsub("\\[", "", indices),
  #     indices = gsub("\\]", "", indices)
  #   ) |> 
  #   separate(
  #     indices,
  #     sep = ",",
  #     into = c("ag", "group", "dataset")
  #   ) |> 
  #   mutate(
  #     ag_name = as.character(levels(ag_name)[as.numeric(ag)]),
  #     sr_group = as.character(levels(sr_group)[as.numeric(group)]),
  #     source = as.character(levels(sr_extra)[as.numeric(dataset)])
  #   ) |> 
  #   select(
  #     source, sr_group, ag_name, estimate
  #   ) |> 
  #   rename(
  #     dataset_sr_group_bias = estimate
  #   )
  
  # Get long version of the map data
  mapdata <- long_map_info(map) |> rename(source = sr_extra)
  
  # Bind in estimates
  mapdata |> 
    left_join(sr_effects, by = "sr_name") |> 
    left_join(dataset_magnitude, by = "source") |> 
  #  left_join(dataset_bias, by = c("ag_name", "source", "sr_group")) |> 
    left_join(imputed_logtiters, by = c("ag_name", "sr_name")) |> 
    left_join(ag_gmts, by = c("ag_name", "sr_group")) -> mapdata
  
  
  return(mapdata)
  
}