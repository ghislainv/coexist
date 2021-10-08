# Seed for reproducibility
seed <- 1234

# Figure width
fig_width <- 16.6 # in cm

nsite_side=25
n_axis=3
nsp=20
IV=TRUE
#part_know <- FALSE
randomOptSp=TRUE
Disp_dep_abund=TRUE
start_full_landscape=TRUE

#Model name

#Perfect knowledge, 3 axis
if(n_axis==3&&IV==FALSE&&Disp_dep_abund==FALSE&&start_full_landscape==TRUE){
  model <- "Perf_know"
}

if(n_axis==3&&IV==FALSE&&Disp_dep_abund==TRUE&&start_full_landscape==TRUE){
  model <- "Perf_know_disp_abund"
}

if(n_axis==3&&IV==FALSE&&Disp_dep_abund==FALSE&&start_full_landscape==FALSE){
  model <- "Perf_know_start_1_ind_sp"
}

if(n_axis==3&&IV==FALSE&&Disp_dep_abund==TRUE&&start_full_landscape==FALSE){
  model <- "Perf_know_disp_abund_start_1_ind_sp"
}

#Partial knowledge
if(IV==TRUE&&Disp_dep_abund==FALSE&&start_full_landscape==TRUE){
  model <- "Part_know"
}

if(IV==TRUE&&Disp_dep_abund==TRUE&&start_full_landscape==TRUE){
  model <- "Part_know_disp_abund"
}

if(IV==TRUE&&Disp_dep_abund==FALSE&&start_full_landscape==FALSE){
  model <- "Part_know_start_1_ind_sp"
}

if(IV==TRUE&&Disp_dep_abund==TRUE&&start_full_landscape==FALSE){
  model <- "Part_know_disp_abund_start_1_ind_sp"
}

if(IV==TRUE){
  # Charge files from the model without IV
  load(file=here::here("outputs", "m0", "lm_fit.RData"))
  load(file=here::here("outputs", "m0", "V_intra.RData"))
  load(here::here("outputs", "m0", "sites.RData"))
  load(here::here("outputs", "m0", "env.RData"))
}

#Size of the side of the environmental matrix

if(IV==TRUE){
  nsite_side <- sqrt(nrow(sites))
} else {nsite_side <- nsite_side}

nsite <- nsite_side^2

# Niche width (used only if randomOptSp == FALSE)
niche_width <- 0.25

#Number of species (used only if randomOptSp == TRUE)

if(IV==TRUE){
  nsp <- nrow(V_intra)
} else {nsp <- nsp}

# Basal mortality
theta <- 0.1

# Strength of unsuitability for mortality
b <- -0.5

# Fecundity (used only if the fecundity option is activated)
fecundity <- 0.5

# Number of repetitions
nrep <- 10

# Number of generations
ngen <- 1000
