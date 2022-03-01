## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================


# Seed for reproducibility
seed <- sample(1:10^6, 1)

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

if(perf_know==TRUE){
  source(file = "./call_libraries.R")
  source("infer_IV.R")
  infer_IV(model, n_observed_axis)
}

# ========================================
# Launch partial knowledge models
# ========================================
# Without IV
perf_know <- FALSE
source(file = here::here("Basic_parameters.R"))

launch_model()

#With IV
IV <- TRUE
source(file = here::here("Basic_parameters.R"))

launch_model()

# ===============
# Compare models
# ===============

source("Model_comparison.R")
compare_models()

# =========================
# End of file
# =========================