# Figure width
fig_width <- 16.6 # in cm

# Number of niche/environmental axes
n_axis <- 10

# Perfect of partial knowledge
if(perf_know==FALSE){part_know<-TRUE} else{part_know<-FALSE}

#Random optimum of species
randomOptSp<-TRUE

#Initialization of the landscape
start_full_landscape<-FALSE
start_one_ind_per_species<-FALSE
start_ten_ind_per_species<-TRUE

#Dependence of the dispersion of seeds to the abundance of species
disp_dep_abund<-TRUE

#Probability of mortality as a function of species performance
mortality_fixed<-TRUE

# Niche width (used only if randomOptSp == FALSE)
niche_width <- 0.25

# Basal mortality
theta <- 0.1

# Strength of unsuitability for mortality
b <- -0.5

# Fecundity (used only if the disp_dep_abund==TRUE)
fecundity <- 0.5

# Number of repetitions
nrep <- 10

# Number of generations
ngen <-10000

############
#Model name#
############

#Perfect knowledge

##Stochastic mortality

###Start with full landscape

if(perf_know==TRUE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_full_mort_stocha"
}

if(perf_know==TRUE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_full_mort_stocha_disp_abund"
}

###Start with one individual per species

if(perf_know==TRUE&&disp_dep_abund==FALSE&&start_one_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_start_1_mort_stocha"
}

if(perf_know==TRUE&&disp_dep_abund==TRUE&&start_one_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_start_1_mort_stocha_disp_abund"
}

###Start with ten individual per species

if(perf_know==TRUE&&disp_dep_abund==FALSE&&start_ten_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_start_10_mort_stocha"
}

if(perf_know==TRUE&&disp_dep_abund==TRUE&&start_ten_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Perf_know_start_10_mort_stocha_disp_abund"
}

##Fixed mortality

###Start with full landscape

####Dispersal is not dependent of species abundance

if(perf_know==TRUE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Perf_know_full_mort_fixed"
}

####Dispersal is dependent of species abundance

if(perf_know==TRUE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Perf_know_full_mort_fixed_disp_abund"
}

###Start with ten individuals per species

####Dispersal is not dependent of species abundance

if(perf_know==TRUE&&disp_dep_abund==FALSE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Perf_know_start_10_mort_fixed"
}

####Dispersal is dependent of species abundance

if(perf_know==TRUE&&disp_dep_abund==TRUE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Perf_know_start_10_mort_fixed_disp_abund"
}


#Partial knowledge

##With IV

###Stochastic mortality

####Start with full landscape

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_IV_full_mort_stocha"
  model_perf <- "Perf_know_full_mort_stocha"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_IV_full_mort_stocha_disp_abund"
  model_perf <- "Perf_know_full_mort_stocha_disp_abund"
}

####Start with one individual per species

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==FALSE&&start_one_ind_per_species==TRUE){
  model <- "Part_know_IV_start_1_mort_stocha"
  model_perf <- "Perf_know_start_1_mort_stocha"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==TRUE&&start_one_ind_per_species==TRUE){
  model <- "Part_know_IV_start_1_mort_stocha_disp_abund"
  model_perf <- "Perf_know_start_1_mort_stocha_disp_abund"
}

###Start with ten individual per species

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==FALSE&&start_ten_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_IV_start_10_mort_stocha"
  model_perf <- "Perf_know_start_10_mort_stocha"
}

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==TRUE&&start_ten_ind_per_species==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_IV_start_10_mort_stocha_disp_abund"
  model_perf <- "Perf_know_start_10_mort_stocha_disp_abund"
}

###Fixed mortality

####Start with full landscape

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_IV_full_mort_fixed"
  model_perf <- "Perf_know_full_mort_fixed"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_IV_full_mort_fixed_disp_abund"
  model_perf <- "Perf_know_full_mort_fixed_disp_abund"
}

####Start with ten individual per species

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==FALSE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_IV_start_10_mort_fixed"
  model_perf <- "Perf_know_start_10_mort_fixed"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==TRUE&&disp_dep_abund==TRUE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_IV_start_10_mort_fixed_disp_abund"
  model_perf <- "Perf_know_start_10_mort_fixed_disp_abund"
}

##No IV

###Stochastic mortality

####Start with full landscape

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_full_mort_stocha"
  model_perf <- "Perf_know_full_mort_stocha"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==FALSE){
  model <- "Part_know_full_mort_stocha_disp_abund"
  model_perf <- "Perf_know_full_mort_stocha_disp_abund"
}

####Start with one individual per species

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==FALSE&&start_one_ind_per_species==TRUE){
  model <- "Part_know_start_1_mort_stocha"
  model_perf <- "Perf_know_start_1_mort_stocha"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==TRUE&&start_one_ind_per_species==TRUE){
  model <- "Part_know_start_1_mort_stocha_disp_abund"
  model_perf <- "Perf_know_start_1_mort_stocha_disp_abund"
}

###Fixed mortality

####Start with full landscape

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==FALSE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_full_mort_fixed"
  model_perf <- "Perf_know_full_mort_fixed"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==TRUE&&start_full_landscape==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_full_mort_fixed_disp_abund"
  model_perf <- "Perf_know_full_mort_fixed_disp_abund"
}

####Start with ten individual per species

#####Dispersal is not dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==FALSE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_start_10_mort_fixed"
  model_perf <- "Perf_know_start_10_mort_fixed"
}

#####Dispersal is dependent of species abundance

if(part_know==TRUE&&IV==FALSE&&disp_dep_abund==TRUE&&start_ten_ind_per_species==TRUE&&mortality_fixed==TRUE){
  model <- "Part_know_start_10_mort_fixed_disp_abund"
  model_perf <- "Perf_know_start_10_mort_fixed_disp_abund"
}

# Add the number of axes and observed axes for partial knowledge only
if(part_know==FALSE){
  model <- paste0(model, "_", n_axis, "_axes")
}
if(part_know==TRUE){
  model <- paste0(model, "_", n_axis, "_axes_", n_observed_axis, "_obs")
  model_perf <- paste0(model_perf, "_", n_axis, "_axes")
}

# Add the seed

model <- paste0(model, "_seed_", seed)
if(part_know==TRUE){
  model_perf <- paste0(model_perf, "_seed_", seed)
}

##############
# Load files #
##############

#Load species parameters and environmental variables
if(part_know==TRUE){
  # Charge files from the model without IV
  # Use model run just before
  load(file=here::here("outputs", model_perf, glue::glue("lm_fit_{n_observed_axis}_obs_axes.RData")))
  load(here::here("outputs", model_perf, "sites.RData"))
  load(here::here("outputs", model_perf, "env.RData"))
  load(file=here::here("outputs", model_perf, "niche_optimum.RData"))
  load(file=here::here("outputs", model_perf, glue::glue("Inferred_species_parameters_{n_observed_axis}_obs_axes.RData")))
  load(file=here::here("outputs", model_perf, "env_entrelac.RData"))
}

#Number of species
if(part_know==TRUE){
  nsp <- length(unique(model.frame(lm_fit)$Species))
} else {nsp <- 20}

#Load IV
if(IV==TRUE){
  load(file=here::here("outputs", model_perf, glue::glue("V_intra_{n_observed_axis}_obs_axes.RData")))
}

#Size of the side of the environmental matrix
if(part_know==TRUE){
  nsite_side <- sqrt(nrow(sites))
} else {nsite_side <- 25}

nsite <- nsite_side^2
