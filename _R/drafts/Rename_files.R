files <- c("Abundances", "community_end", "community_start", "Equitability", "perf_Sp_mean", "rank_sp", "Shannon", "sp_rich", "Spearman", "Species_richness")

for (seed in Seeds) {
  
  for(n_observed_axes in nb_obs_axes){
    
    models <- c(glue::glue("Perf_know_start_10_mort_fixed_seeds_dep_abund_10_axes_seed_{seed}"),
                glue::glue("Part_know_start_10_mort_fixed_seeds_dep_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"),
                glue::glue("Part_know_IV_start_10_mort_fixed_seeds_dep_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"))
    
    for(model in models){
    
      for(file in files){
        old_name <- paste0(file, "_", model, ".RData")
        new_name <- paste0(file, ".RData")
      
        file.rename(here::here("outputs", model, old_name), here::here("outputs", model, new_name))
      }
    }
  }
}


for (seed in Seeds) {
  
  for(n_observed_axes in nb_obs_axes){
    
    old_models <- c(glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_seed_{seed}"),
                    glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"),
                    glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"))
    
    new_models <- c(glue::glue("Perf_know_start_10_mort_fixed_seeds_dep_abund_10_axes_seed_{seed}"),
                glue::glue("Part_know_start_10_mort_fixed_seeds_dep_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"),
                glue::glue("Part_know_IV_start_10_mort_fixed_seeds_dep_abund_10_axes_{n_observed_axes}_obs_seed_{seed}"))
    
    for(model in 1:length(old_models)){
      
        old_name <- old_models[model]
        new_name <- new_models[model]
        
        file.rename(here::here("outputs", old_name), here::here("outputs", new_name))
    }
  }
}