## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

source(file = here::here("Basic_parameters.R"))

source("launch_model.R")
launch_model()

# ========================================
# Infer observed intraspecific variability
# ========================================

if(perf_know==TRUE){
  source("infer_IV.R")
  infer_IV(model, n_observed_axis)
}

# ===============
# Compare models
# ===============

source("Model_comparison.R")
compare_models()

# =========================
# End of file
# =========================