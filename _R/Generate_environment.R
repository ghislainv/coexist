generate_environment <- function(nsite_side, model, seed){
  # Landscape
  mat <- matrix(0, nrow=nsite_side, ncol=nsite_side)
  r <- raster::raster(mat, crs="+proj=utm +zone=1")
  nsite <- nsite_side^2
  coords <- sp::coordinates(r)
  
  # Neighbourhood matrix
  neighbors.mat <- raster::adjacent(r, cells=c(1:nsite), directions=8,
                            pairs=TRUE, sorted=TRUE)
  # Number of neighbours by site
  n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
  # Adjacent sites
  adj <- neighbors.mat[,2]
  # Generate symmetric adjacency matrix, A
  A <- matrix(0,nsite,nsite)
  index.start <- 1
  for (i in 1:nsite) {
    index.end <- index.start+n.neighbors[i]-1
    A[i,adj[c(index.start:index.end)]] <- 1
    index.start <- index.end+1
  }
  
  # Generate spatial random effects
  Vrho.target <- 1 # Variance of spatial random effects
  d <- 1 # Spatial dependence parameter = 1 for intrinsic CAR
  Q <- diag(n.neighbors)-d*A + diag(.0001,nsite) # Add small constant to make Q non-singular
  covrho <- Vrho.target*solve(Q) # Covariance of rhos
  
  # Environment on each site
  sites <- data.frame(matrix(ncol=n_axes, nrow=nsite))
  colnames(sites) <- c(sprintf("V%d_env", 1:n_axes))
  env <- list()
  for (i in 1:n_axes) {
    seed_env <- seed + i - 1
    rho <- c(rmvn(1, mu=rep(0, nsite), V=covrho, seed=seed_env)) # Spatial Random Effects
    rho <- rho-mean(rho) # Centering rhos on zero
    #rho <- scales::rescale(rho, to=c(0, 1))
    rho <- (rho - min(rho)) / (max(rho) - min(rho))
    sites[,i] <- rho
    env[[i]] <- matrix(rho, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  }
  
  #save the environmental data to avoid simulating them again in the next model
  save(sites, file = here::here("outputs", model, "sites.RData"))
  save(env, file = here::here("outputs", model, "env.RData"))
  
  #Look at the correlations between environmental variables
  Corr_env <- as.data.frame(as.matrix(t(combn(c(1:n_axes), 2))))
  colnames(Corr_env) <- c("Var1", "Var2")
  Corr_env$Corr <- numeric(nrow(Corr_env))
  
  for(k in 1:nrow(Corr_env)){
    Corr_env$Corr[k] <- cor(raster::values(raster::raster(env[[Corr_env$Var1[k]]])), raster::values(raster::raster(env[[Corr_env$Var2[k]]])))
  }
  
  save(Corr_env, file = here::here("outputs", model, "Corr_env.RData"))
  
  # Plot the environment
  plot_environment(model=model, fig_width=fig_width, n_axes=n_axes, env=env, sites=sites)
  
  # Plot the habitat frequency
  plot_hab_freq(n_axes=n_axes, model=model, fig_width=fig_width, env=env)
}