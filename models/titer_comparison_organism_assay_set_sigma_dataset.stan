
functions {
  
  real normal_int_censored_likelihood(real lower_lim, real upper_lim, real mu, real sigma) {
    
    real result;
    
    if (is_inf(lower_lim) && is_inf(upper_lim)) {
      
      result = 0;
      
    } else if (lower_lim == upper_lim) {
      
   //   real sigma_force = 0.5;
   
      result = normal_lpdf(
        lower_lim | mu, sigma
      );
      
    } else if (!is_inf(lower_lim) && !is_inf(upper_lim)) {
  
      result = log_diff_exp(
        normal_lcdf(upper_lim | mu, sigma),
        normal_lcdf(lower_lim | mu, sigma)
      );
      
    } else if (!is_inf(lower_lim) && is_inf(upper_lim)) {
      
      result = normal_lccdf(
        lower_lim | mu, sigma
      );
      
    } else if (is_inf(lower_lim) && !is_inf(upper_lim)) {
      
      result = normal_lcdf(
        upper_lim | mu, sigma
      );
      
    }
    
    return result;
    
  }
  
  real normal_lub_rng(real mu, real sigma, real lower_lim, real upper_lim) {
    
    if (lower_lim == upper_lim) {
      
      return lower_lim;
      
    } else {
      
      real p_lower_lim = normal_cdf(lower_lim, mu, sigma);
      real p_upper_lim = normal_cdf(upper_lim, mu, sigma);
      real u = uniform_rng(p_lower_lim, p_upper_lim);
      real y = mu + sigma * inv_Phi(u);
      return y;
      
    }
  }
  
}

data {
  int<lower=1> N;
  int<lower=1> N_ags;
  int<lower=1> N_srs;
  int<lower=1> N_sr_groups;
  int<lower=1> N_datasets;
  int<lower=1> N_organisms;
  int<lower=1> N_assays;
  vector[N] upper_lims;
  vector[N] lower_lims;
  array[N] int ags;
  array[N] int srs;
  array[N] int sr_groups;
  array[N] int datasets;
  array[N] int organisms;
  array[N] int assays;
  real ag_mean_prior_mean;
  real<lower=0> ag_mean_prior_sigma;
  real<lower=0> sigma_prior_alpha;
  real<lower=0> sigma_prior_beta;
  //vector<lower=0>[N_datasets] sigma_error;
  real<lower=0> sigma_prior_sr_effect;
  real<lower=0> sigma_prior_organism_magnitude;
  real<lower=0> sigma_prior_assay_effect;
 // vector[N_datasets] dataset_magnitude_mean;
 vector[N_organisms] organism_magnitude_mean;
 vector[N_assays] assay_effect_mean;
}

parameters {
  matrix[N_ags, N_sr_groups] ag_means;
  vector[N_srs] sr_effects;
  vector[N_organisms] organism_magnitude;
  vector[N_assays] assay_effects;
 // vector[N_datasets] dataset_magnitude;
  vector<lower=0>[N_datasets] sigma_error;
}

model {
  
  // Calculate likelihood of priors
  sr_effects ~ normal(0, sigma_prior_sr_effect);
  
   for (i in 1:N_datasets) {
     sigma_error[i] ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
   }
   
  
  for (i in 1:N_sr_groups) {
    
 //   sigma_error[,i] ~ inv_gamma(sigma_prior_alpha, sigma_prior_beta);
    ag_means[,i] ~ normal(ag_mean_prior_mean, ag_mean_prior_sigma);
    
  }
  
  
  for (i in 1:N_organisms) {
    organism_magnitude[i] ~ normal(organism_magnitude_mean[i], sigma_prior_organism_magnitude);
  }
  
  for (i in 1:N_assays) {
    assay_effects[i] ~ normal(assay_effect_mean[i], sigma_prior_assay_effect);
  }
  
  // Work out likelihood of each titer
  for (i in 1:N) {
    
    // Set variables
    int ag = ags[i];
    int sr = srs[i];
    int sr_group = sr_groups[i];
    int dataset = datasets[i];
    int assay = assays[i];
    int organism = organisms[i];
    
    // Add in mixed effect for each antigen and serum
    real logtiter = ag_means[ag, sr_group] + 
      sr_effects[sr] + 
      organism_magnitude[organism] + 
      assay_effects[assay]; 
      
    
    // Work out titer likelihood
    target += normal_int_censored_likelihood(
      lower_lims[i],
      upper_lims[i],
      logtiter, 
      sigma_error[dataset]
    );
    
  }
  
}

generated quantities {

  vector[N] imputed_logtiters;

  // Work out likelihood of each titer
  for (i in 1:N) {

    // Set variables
    int ag = ags[i];
    int sr = srs[i];
    int sr_group = sr_groups[i];
    int dataset = datasets[i];
    int assay = assays[i];
    int organism = organisms[i];

    // Add in mixed effect for each antigen and serum
     real logtiter = ag_means[ag, sr_group] +
      sr_effects[sr] + 
      organism_magnitude[organism] + 
      assay_effects[assay]; 

      // Add in mixed effect for each antigen and serum
    imputed_logtiters[i] = normal_lub_rng(
      logtiter,
      sigma_error[dataset],
      lower_lims[i],
      upper_lims[i]
      );




  }

}