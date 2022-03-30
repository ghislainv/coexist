## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# Number of observed axes in partial models
nb_obs_axes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

# Seeds for reproducibility: it controls the environment × species parameters configuration.
# Seeds <- sample(1:10^6, 10)
# save(Seeds, file=here::here("outputs", "Seeds.RData"))
load(here::here("outputs", "Seeds.RData"))

source(file = here::here("_R", "call_libraries.R"))

for(configuration in 1:length(Seeds)){
  
  seed <- Seeds[configuration]
  
  # ========================================
  # Launch perfect knowledge model
  # ========================================
  
  perf_know <- TRUE
  IV <- FALSE
  
  source(file = here::here("_R", "Basic_parameters.R"))
  source(file = here::here("_R", "launch_model.R"))
  
  launch_model()
  
  # ========================================
  # Infer observed intraspecific variability
  # ========================================

  for(n_observed_axis in nb_obs_axes){

    if(perf_know==TRUE){
      source(here::here("_R", "infer_IV.R"))
      infer_IV(model, n_observed_axis)
    }
  }
  
  # ========================================
  # Launch partial knowledge models
  # ========================================
    
  perf_know <- FALSE
  
  for(n_observed_axis in nb_obs_axes){
    
    # Without IV
    IV <- FALSE
    source(file = here::here("_R", "Basic_parameters.R"))
    
    launch_model()
    
    #With IV
    IV <- TRUE
    source(file = here::here("_R", "Basic_parameters.R"))
    
    launch_model()
    
  }
}

# ===============
# Compare models
# ===============

source(here::here("_R", "Model_comparison.R"))
compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side)
Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)

# =========================
# End of file
# =========================