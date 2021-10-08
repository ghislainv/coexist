generate_species_optima <- function(model, randomOptSp, niche_width, nsp, env){
  if (!randomOptSp){
    
    # Niches per axis
    n_niche <- 1/niche_width
    
    # Number of species
    nsp <- n_niche^n_axis
    
    # Species coordinates on the three niche axis (x, y, z)
    base_coord <- seq(0, 1, length.out=n_niche+1)[-(n_niche+1)]+niche_width/2
    sp_x <- rep(rep(base_coord, n_niche), n_niche)
    sp_y <- rep(rep(base_coord, each=n_niche), n_niche)
    sp_z <- rep(base_coord, each=n_niche^2)
    niche_optimum <- as.data.frame(cbind(sp_x, sp_y, sp_z))
  }
  
  # Random optimum for species
  if (randomOptSp) {
    set.seed(seed)
    niche_optimum <- data.frame(sp_x=runif(n=nsp, min = min(env[[1]]), max = max(env[[1]])),
                                sp_y=runif(n=nsp, min = min(env[[2]]), max = max(env[[2]])),
                                sp_z=runif(n=nsp, min = min(env[[3]]), max = max(env[[3]])))
    niche_optimum <- niche_optimum[order(niche_optimum$sp_z, niche_optimum$sp_y, niche_optimum$sp_x), ]
  }
  
  rownames(niche_optimum)<-1:nrow(niche_optimum)
  
  # Plot the species niche
  plot_species_optima(model, fig_width, niche_optimum)
  return(niche_optimum)
}