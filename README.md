# roessler_netzl_et_al2024
This repository contains the code for the NHP antigenic cartography manuscript by Annika Rössler, Antonia Netzl, et. al., 2024. Please cite the original publication if any data or code from this repository is used. 

The repository's DOI was created with Zenodo (https://doi.org/10.25495/7gxk-rd71)

Raw data can be found in the `data` directory and will be added upon publication. The code for the analyses shown in the main manuscript can be found in the `code` directory. The `alignment_map.ace` in the `data/maps` directory is the human map from Rössler, Netzl, et al. Nat Commun 14, 5224 (2023) (https://doi.org/10.1038/s41467-023-41049-4)
and the `Wilks_et_al_map_ndsubset_no_outliers_slope_adjusted.ace` is the human map from Wilks, Mühlemann, et al., Science (2023) (https://doi.org/10.1126/science.adj0070).


To obtain a titer table for antigenic map construction, execute the `excel_to_titertable.R` script. To construct maps, first the `make_map_final.R` script. Antibody landscapes are fit in the `ablandscapes_fit_multi.R`. The `3D_landscapes.Rmd` creates an html script with 3D illustrations of the landscapes. The Bayesian titer comparison is done in the `model_map_organism_assay_dataset_sigma.R` script and plots are generated in the `plot_stan_model_organism_assay.R`. The STAN model is stored in the `model` directory.

All SOM analyses can be found in the `som` directory, the Bayesian analysis SOM plots are generated in the scripts in the `code` directory. 

The `function` directory contains utility functions, the `data/metadata` directory contains metadata such as colour information.

All analyses were performed in R version 4.2.2 (2022-10-31).
R Core Team (2022). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna,
  Austria. URL https://www.R-project.org/.
  
Antigenic maps were constructed using the Racmacs package, Version 1.1.35:
Wilks S (2022). _Racmacs: R Antigenic Cartography Macros_. https://acorg.github.io/Racmacs,
  https://github.com/acorg/Racmacs.
  
Antibody landscapes were constructed using the ablandscapes package, Version 1.1.0: 
Wilks S (2021). _ablandscapes: Making Antibody landscapes Using R_. R package
  version 1.1.0, <https://github.com/acorg/ablandscapes>.

Geometric mean titers and fold changes were calculated using the titertools package, Version 0.0.0.9001:
Wilks SH (2022). _titertools: A statistical toolkit for the annalysis of censored
  titration data_. R package version 0.0.0.9001.
  
The Bayesian modelling was performed with cmdstanr: 
Gabry J, Češnovar R (2022). _cmdstanr: R Interface to 'CmdStan'_.
  https://mc-stan.org/cmdstanr/, https://discourse.mc-stan.org.
