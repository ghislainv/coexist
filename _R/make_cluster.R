## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================


load(here::here("Array_simulations.RData"))

source(file = here::here("_R", "call_libraries.R"))

source(file=here::here("_R", "Math_functions.R"))

source(file=here::here("_R", "Plot_functions_cluster.R"))

source(file=here::here("_R", "Generate_environment_cluster.R"))

source(file=here::here("_R", "Species_parameters_cluster.R"))

source(file = here::here("_R", "launch_script_cluster.R"))

source(here::here("_R", "infer_IV_cluster.R"))

# Number of observed axes in partial models
nb_obs_axes <- c(0:10)

# Seeds for reproducibility: it controls the environment × species parameters configuration.
# Seeds <- sample(1:10^6, 10)
# save(Seeds, file=here::here("outputs", "Seeds.RData"))
load(here::here("Seeds.RData"))

s=1
  
  # ========================================
  # Launch perfect knowledge model
  # ========================================
  
  perf_know <- TRUE
  IV <- FALSE
  
  n_observed_axes <- 0
  
  launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
  
  # ========================================
  # Infer observed intraspecific variability
  # ========================================
  
  for(n_observed_axes in nb_obs_axes){
    
    if(perf_know==TRUE){
      infer_IV(n_observed_axes)
    }
  }
  
  # ========================================
  # Launch partial knowledge models
  # ========================================
  
  perf_know <- FALSE
  
  for(n_observed_axes in nb_obs_axes){
    
    # Without IV
    IV <- FALSE
    
    launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
    
    #With IV
    IV <- TRUE
    
    launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
  }
#}

# # ===============
# # Compare models
# # ===============
# 
# source(here::here("_R", "Model_comparison.R"))
# compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes)
# Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)
# 
# # =========================
# # End of file
# # =========================