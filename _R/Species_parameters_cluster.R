generate_species_optima <- function(randomOptSp, niche_width, nsp, env, seed_sp){
  if (!randomOptSp){
    
    # Niches per axis
    n_niche <- 1/niche_width
    
    # Number of species
    nsp <- n_niche^n_axes
    
    # Species coordinates on three niche axis (x, y, z)
    base_coord <- seq(0, 1, length.out=n_niche+1)[-(n_niche+1)]+niche_width/2
    sp_x <- rep(rep(base_coord, n_niche), n_niche)
    sp_y <- rep(rep(base_coord, each=n_niche), n_niche)
    sp_z <- rep(base_coord, each=n_niche^2)
    niche_optimum <- as.data.frame(cbind(sp_x, sp_y, sp_z))
  }
  
  # Random optimum for species (unlimited number of species)
  if (randomOptSp) {
    set.seed(seed_sp)
    niche_optimum <- data.frame(matrix(ncol=n_axes, nrow=nsp))
    for(k in 1:n_axes){niche_optimum[,k]<-runif(n=nsp, min = min(env[[k]]), max = max(env[[k]]))}
  }
  
  rownames(niche_optimum)<-1:nrow(niche_optimum)
  save(niche_optimum, file=here::here("outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  
#   # Plot the species niche
#   if(n_axes==3){
#     plot_species_optima(model, fig_width, niche_optimum)
#   }
#   return(niche_optimum)
}