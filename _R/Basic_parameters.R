# Figure width
fig_width <- 16.6 # in cm

# Number of niche/environmental axes
n_axes <- 10

# Perfect of partial knowledge
if(perf_know==FALSE){part_know<-TRUE} else{part_know<-FALSE}

#Random optimum of species
randomOptSp<-TRUE

#Initialization of the landscape
start_full_landscape<-FALSE
start_one_ind_per_species<-FALSE
start_ten_ind_per_species<-TRUE

# Niche width (used only if randomOptSp == FALSE)
niche_width <- 0.25

# Basal mortality
theta <- 0.01

# Strength of unsuitability for mortality
b <- -0.5

# Fecundity (used only if the nb_seeds_dep_abund==TRUE)
fecundity <- 0.5

# Number of repetitions
nrep <- 2

# Number of generations
ngen <-100

############
#Model name#
############

if(perf_know==TRUE){
  mod <- "Perf_know"
}else{
  if(IV==FALSE){
    mod <- "Part_know"
  }else{
    mod <- "Part_know_IV"
  }
}

if(nb_seeds_dep_abund==TRUE){
  nb_seeds <- "seeds_dep_abund"
}else{nb_seeds <- "seeds_indep_abund"}

if(mortality_fixed==TRUE){
  mort <- "mort_fixed"
}
if(mortality_stocha==TRUE){
  if(mortality_stocha_basal==TRUE){
    mort <- "mort_stocha_basal"
  }else{mort <- "mort_stocha"}
}
if(mortality_proportion==TRUE){
    mort <- "mort_prop"
}

if(start_ten_ind_per_species==TRUE){
  start <- "start_10"
}else{
  if(start_one_ind_per_species==TRUE){
    start <- "start_1"
  }else{
    if(start_full_landscape==TRUE){
      start <- "full"
    }
  }
}

model <- glue::glue("{mod}_{start}_{mort}_{nb_seeds}_{n_axes}_axes")

if(part_know==TRUE){
  model <- glue::glue("{model}_{n_observed_axes}_obs_seed_{seed}")
  model_perf <- glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}")
}else{model <- glue::glue("{model}_seed_{seed}")}

##############
# Load files #
##############

#Load species parameters and environmental variables
if(part_know==TRUE){
  # Charge files from the model without IV
  # Use model run just before
  load(file=here::here("outputs", model_perf, glue::glue("lm_fit_{n_observed_axes}_obs_axes.RData")))
  load(here::here("outputs", model_perf, "sites.RData"))
  load(here::here("outputs", model_perf, "env.RData"))
  load(file=here::here("outputs", model_perf, "niche_optimum.RData"))
  load(file=here::here("outputs", model_perf, glue::glue("Inferred_species_parameters_{n_observed_axes}_obs_axes.RData")))
  load(file=here::here("outputs", model_perf, "env_entrelac.RData"))
}

#Number of species
if(part_know==TRUE){
  nsp <- length(unique(model.frame(lm_fit)$Species))
} else {nsp <- 20}

#Load IV
if(IV==TRUE){
  load(file=here::here("outputs", model_perf, glue::glue("V_intra_{n_observed_axes}_obs_axes.RData")))
}

#Size of the side of the environmental matrix
if(part_know==TRUE){
  nsite_side <- sqrt(nrow(sites))
} else {nsite_side <- 25}

nsite <- nsite_side^2
