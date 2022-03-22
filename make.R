## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# Number of observed axes in partial models
nb_obs_axes <- c(1, 3, 5, 7)

# Seeds for reproducibility: it controls the environment × species parameters configuration.
# Seeds <- sample(1:10^6, 10)
# save(Seeds, file=here::here("outputs", "Seeds.RData"))
load(here::here("outputs", "Seeds.RData"))

for(configuration in 1:10){
  
  seed <- Seeds[configuration]
  
  # ========================================
  # Launch perfect knowledge model
  # ========================================
  
  perf_know <- TRUE
  IV <- FALSE
  
  source(file = here::here("Basic_parameters.R"))
  source("launch_model.R")
  
  launch_model()
  
  # ========================================
  # Infer observed intraspecific variability
  # ========================================
  
  for(n_observed_axis in nb_obs_axes){
  
    if(perf_know==TRUE){
      source(file = here::here("_R", "call_libraries.R"))
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
compare_models()

# =========================
# End of file
# =========================