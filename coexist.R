## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# Libraries
library(raster)
library(here)
library(scales)
library(ggplot2)

# Seed for reproducibility
seed <- 1234

# =========================
# Landscape and environment
# =========================

# Landscape
ncell_side <- 25
mat <- matrix(0, nrow=ncell_side, ncol=ncell_side)
r <- raster(mat, crs="+proj=utm +zone=1")
ncell <- ncell(r)
coords <- coordinates(r)

# Neighbourhood matrix
neighbors.mat <- adjacent(r, cells=c(1:ncell), directions=8,
                          pairs=TRUE, sorted=TRUE)
# Number of neighbours by cell
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
# Adjacent cells
adj <- neighbors.mat[,2]
# Generate symmetric adjacency matrix, A
A <- matrix(0,ncell,ncell)
index.start <- 1
for (i in 1:ncell) {
  index.end <- index.start+n.neighbors[i]-1
  A[i,adj[c(index.start:index.end)]] <- 1
  index.start <- index.end+1
}

# Function to draw in a multivariate normal
rmvn <- function(n, mu=0, V=matrix(1), seed=1234) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) {
    stop("Dimension problem!")
  }
  D <- chol(V)
  set.seed(seed)
  t(matrix(rnorm(n*p),ncol=p)%*%D+rep(mu,rep(n,p)))
}

# Generate spatial random effects
Vrho.target <- 1 # Variance of spatial random effects
d <- 1 # Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,ncell) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos
rho <- c(rmvn(1, mu=rep(0,ncell), V=covrho, seed=seed)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero
rho <- scales::rescale(rho, to=c(0, 1))
env <- matrix(rho, nrow=ncell_side, ncol=ncell_side, byrow=TRUE)

# Plot
plot(raster(env), main="Environment")
hist(env, main="Environment", freq=FALSE)

# =========================================
# Species performance given the environment
# =========================================

nsp <- 50
# Remove last species (otherwise last species=first species)
perf <- seq(0, 1, length.out=nsp+1)[-(nsp+1)]

# Habitat frequency for each species
sp_hab_freq <- dnorm(perf, mean=mean(env), sd=sd(env))
plot(sp_hab_freq)

# Matrix of species performance on each cell (distances)
# Cell in rows, Species in columns
dist_E_Sp <- matrix(NA, nrow=ncell, ncol=nsp) 
for (i in 1:ncell) {
  for (j in 1:nsp) {
    dist <- abs(perf[j]-env[i])
    if (dist>0.5) {
      dist <- 1-dist
    }
    dist_E_Sp[i,j] <- dist
  }
}

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 10
# Number of generations
ngen <- 100
# Mortality probability
theta <- 0.2

# Function to identifying the species with the highest performance
high_perf_sp <- function(dist, sp_pres) {
  dist_pres <- dist[sp_pres]
  min_dist <- min(dist_pres)
  sp_pref <- which(dist==min_dist)
  # If two species adapted, selection at random
  if (length(sp_pref)==2) {
    # Random permutation
    sp_pref <- sample(sp_pref)
    sp_pref <- sp_pref[1]
  }
  return(sp_pref)
}

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)

# Species rank at the end of the simulations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)

# Loop on repetitions
for (r in 1:nrep) {

  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Draw species at random in the landscape (one individual per cell)
  sp <- sample(1:nsp, size=ncell, replace=TRUE)
  #hist(sp)
  scene_start <- matrix(sp, nrow=ncell_side, ncol=ncell_side, byrow=TRUE)
  if (r==1) {
    pdf(file=here("outputs", "scene_start.pdf"))
    plot(raster(scene_start), main="Species", col=c(rev(terrain.colors(255))))
    dev.off()
  }
  scene <- scene_start

  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Species richness
  sp_rich[1, r] <- length(unique(scene[scene!=0]))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(scene, levels=1:nsp))
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # ******************
    # Mortality
    # ******************
    
    # Mortality events
    mort_ev <- rbinom(ncell, size=1, prob=theta)
    mortality <- matrix(mort_ev, nrow=ncell_side, ncol=ncell_side, byrow=TRUE)
    
    # Update scene
    scene[mortality==1] <- 0
    # Plot once
    if (r==1 & g==1) {
      pdf(file=here("outputs", "mortality_events.pdf"))
      plot(raster(scene), main="Species", col=c("black", rev(terrain.colors(255))))
      dev.off()
    }

    # *********************
    # Fecundity/Recruitment
    # *********************
    
    # Species present in the community
    sp_present <- sort(unique(scene[scene!=0]))
    nsp_present <- length(sp_present)
    
    # Vacant sites
    scene_rast <- raster(scene)
    sites_vacant <- which(values(scene_rast)==0)
    
    # Performance of species on vacant sites
    dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
    
    # Identify the species with the highest performance
    new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
    
    # Recruitment
    scene_rast[sites_vacant] <- new_ind
    scene <- as.matrix(scene_rast)
    
    # *********************
    # Diversity
    # *********************
    
    sp_rich[g+1, r] <- length(unique(scene[scene!=0]))
    abund[g+1, ] <- table(factor(scene, levels=1:nsp))
    
  } # End ngen
  
  # Species rank
  rank_sp[r, ] <- rank(-abund[ngen, ], ties.method="min")
  
} # End nrep

# =========================
# Diversity analysis
# =========================

sp_rich
rank_sp

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# Mean final rank
sp_mean_rank <- apply(rank_sp, 2, mean)
# Plot
df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  xlab("Species habitat frequency") +
  ylab("Species mean rank (high rank = low abundance)")
ggsave(p, filename=here("outputs", "mean_rank-habitat_freq.pdf"))

# =========================
# End of file
# =========================