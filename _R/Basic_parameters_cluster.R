# Figure width
#fig_width <- 16.6 # in cm

# Number of niche/environmental axes
n_axes <- 15

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
fec <- 0.5

# Number of repetitions
#nrep <- 10

# Number of generations
ngen <-10000

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
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_sites.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_env.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Inferred_species_parameters_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
  #load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_env_entrelac.RData")))
}

#Number of species
if(part_know==TRUE){
  nsp <- nrow(Inferred_species_parameters)
} else {nsp <- 20}

#Load IV
if(IV==TRUE){
  load(paste0(directory_writing, "/outputs/", glue::glue("V_intra_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
}

#Size of the side of the environmental matrix
if(part_know==TRUE){
  nsite_side <- sqrt(nrow(sites))
} else {nsite_side <- 25}

nsite <- nsite_side^2
