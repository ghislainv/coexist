generate_environment <- function(nsite_side, model){
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
  sites <- data.frame(matrix(ncol=n_axis, nrow=nsite))
  colnames(sites) <- c(sprintf("V%d_env", 1:n_axis))
  env <- list()
  for (i in 1:n_axis) {
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
  print(glue::glue("Variables 1 and 2 have a correlation of {cor(raster::values(raster::raster(env[[1]])), raster::values(raster::raster(env[[2]])))},
                   variables 1 and 3 have a correlation of {cor(raster::values(raster::raster(env[[1]])), raster::values(raster::raster(env[[3]])))}
                   and variables 2 and 3 have a correlation of {cor(raster::values(raster::raster(env[[2]])), raster::values(raster::raster(env[[3]])))}"))
  
  
  # Plot the environment
  plot_environment(model=model, fig_width=fig_width, n_axis=n_axis, env=env, sites=sites)
  
  # Plot the habitat frequency
  plot_hab_freq(n_axis=n_axis, model=model, fig_width=fig_width, env=env)
}