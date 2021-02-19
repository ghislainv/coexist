## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# One dimension
# Mortality rate depending on species performance

# Libraries
library(raster)
library(here)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glue)
library(geoR) # variog()

# Create output directories
dir.create(here("outputs"))
dir.create(here("outputs/m3"))

# Seed for reproducibility
seed <- 1234

# Figure width
fig_width <- 16.6 # in cm

# logit/inv_logit functions
logit <- function(x, min=0, max=1) {
  p <- (x-min)/(max-min)
  return(log(p/(1-p)))
}

inv_logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
  return(p * (max-min) + min)
}

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
save(env, file=here("outputs", "m3", "env.rda"))

# Plot
png(file=here("outputs", "m3", "environment.png"))
plot(raster(env), main="Environment", col=topo.colors(255))
dev.off()

# Habitat frequency
png(file=here("outputs", "m3", "habitat_freq.png"))
hist(env, main="Environment", freq=FALSE)
dev.off()

# =========================================
# Species performance given the environment
# =========================================

# Number of species
nsp <- 50
# Remove last species (otherwise last species=first species)
niche_optimum <- seq(0, 1, length.out=nsp+1)[-(nsp+1)]

# Habitat frequency for each species
sp_hab_freq <- dnorm(niche_optimum, mean=mean(env), sd=sd(env))
png(file=here("outputs", "m3", "species_habitat_freq.png"))
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Matrix of species performance on each cell (distances)
# Cell in rows, Species in columns
dist_E_Sp <- matrix(NA, nrow=ncell, ncol=nsp) 
for (i in 1:ncell) {
  env_i <- values(raster(env))[i] 
  for (j in 1:nsp) {
    dist <- abs(niche_optimum[j]-env_i)
    if (dist>0.5) {
      dist <- 1-dist
    }
    dist_E_Sp[i, j] <- dist
  }
}
perf_E_Sp <- 1-dist_E_Sp

# Function to identify the species with the highest performance
high_perf_sp <- function(dist, sp_pres) {
  dist_pres <- dist[sp_pres]
  min_dist <- min(dist_pres)
  sp_high_perf <- sp_pres[which(dist_pres==min_dist)]
  # If more than one species, selection at random
  if (length(sp_high_perf)>1) {
    # Random permutation
    sp_high_perf <- sample(sp_high_perf)
    sp_high_perf <- sp_high_perf[1]
  }
  return(sp_high_perf)
}

# Probability of dying of each species on each cell
# Strength of unsuitability
b <- 0.5
m_dist <- mean(dist_E_Sp)
s_dist <- sd(dist_E_Sp)
scale_dist_E_Sp <- (dist_E_Sp - m_dist)/s_dist
mortality_E_Sp <- inv_logit(logit(0.1) + b * scale_dist_E_Sp)
# Mortality rate distribution
hist(mortality_E_Sp)

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 10
# Number of generations
ngen <- 100

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Species rank at the end of the generations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
env_filt <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen+1, ncol=nrep)

# Loop on repetitions
for (r in 1:nrep) {

  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Draw species at random in the landscape (one individual per cell)
  sp <- sample(1:nsp, size=ncell, replace=TRUE)
  #hist(sp)
  community_start <- matrix(sp, nrow=ncell_side, ncol=ncell_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m3", "community_start.png"))
    plot(raster(community_start), main="Species - Start", zlim=c(0, 50),
         col=c("black", rev(terrain.colors(50))))
    dev.off()
  }
  community <- community_start
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # Environmental filtering
  community_niche_opt <- matrix(niche_optimum[as.vector(community)], ncol=ncell_side)
  env_filt[1, r] <- mean(abs(env-community_niche_opt))
  
  # Mean mortality rate
  theta_cell <- diag(mortality_E_Sp[, as.vector(t(community))])
  theta_comm[1, r] <- mean(theta_cell) 

  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # ******************
    # Mortality
    # ******************
    
    # Mortality rate on each cell
    theta_cell <- diag(mortality_E_Sp[, as.vector(t(community))])
    
    # Mortality events
    mort_ev <- rbinom(ncell, size=1, prob=theta_cell)
    mortality <- matrix(mort_ev, nrow=ncell_side, ncol=ncell_side, byrow=TRUE)
    
    # Update community
    community[mortality==1] <- 0
    # Plot once
    if (r==1 & g==1) {
      png(file=here("outputs", "m3", "mortality_events.png"))
      plot(raster(community), main="Species - with vacant sites", zlim=c(0, 50),
           col=c("black", rev(terrain.colors(50))))
      dev.off()
    }

    # *********************
    # Fecundity/Recruitment
    # *********************
    
    # Species present in the community
    sp_present <- sort(unique(community[community!=0]))
    nsp_present <- length(sp_present)
    
    # Vacant sites
    community_rast <- raster(community)
    sites_vacant <- which(values(community_rast)==0)
    
    # Performance of species on vacant sites
    dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
    
    # Identify the species with the highest performance
    new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
    
    # Recruitment
    community_rast[sites_vacant] <- new_ind
    community <- as.matrix(community_rast)
    
    # *********************
    # Diversity
    # *********************
    
    # Species richness
    sp_rich[g+1, r] <- length(unique(as.vector(community)))
    abund[g+1, ] <- table(factor(as.vector(community), levels=1:nsp))
    
    # Environmental filtering
    community_niche_opt <- matrix(niche_optimum[as.vector(community)], ncol=ncell_side)
    env_filt[g+1, r] <- mean(abs(env-community_niche_opt))
    
    # Mean mortality rate in the community
    theta_cell <- diag(mortality_E_Sp[, as.vector(t(community))])
    theta_comm[g+1, r] <- mean(theta_cell)
    
  } # End ngen
  
  # Plot final community once
  if (r==1) {
    png(file=here("outputs", "m3", "community_end.png"))
    plot(raster(community), main=glue("Species - End (ngen={ngen})"),
         zlim=c(0, 50), col=c("black", rev(terrain.colors(50))))
    dev.off()
  }
  
  # Species rank
  rank_sp[r, ] <- rank(-abund[ngen, ], ties.method="min")
  
} # End nrep

# =========================
# Diversity analysis
# =========================

sp_rich
rank_sp

# ---------------------------------------------
# Plot species richness
# ---------------------------------------------

sp_rich <- data.frame(sp_rich)
sp_rich_long <- sp_rich %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="sp_rich")
p <- ggplot(data=sp_rich_long, aes(x=gen, y=sp_rich, col=rep)) +
  geom_line() + 
  xlab("Generations") + 
  ylab("Species richness")
ggsave(p, filename=here("outputs", "m3", "species_richness_with_time.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

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
ggsave(p, filename=here("outputs", "m3", "mean_rank-habitat_freq.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# ---------------------------------------------
# Environmental filtering
# ---------------------------------------------

env_filt <- data.frame(env_filt)
env_filt_long <- env_filt %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="env_filt")
p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
  geom_line() +
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean env-species perf difference")
ggsave(p, filename=here("outputs", "m3", "environmental_filtering.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ---------------------------------------------
# Theta community
# ---------------------------------------------

theta_comm <- data.frame(theta_comm)
theta_comm_long <- theta_comm %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="theta_comm")
p <- ggplot(data=theta_comm_long, aes(x=gen, y=theta_comm, col=rep)) +
  geom_line() +
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean mortality rate in the community")
ggsave(p, filename=here("outputs", "m3", "mortality_rate_community.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(community)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=values(raster(env)))
# Plot with correlation
png(file=here("outputs", "m3", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2))
plot(vario_sp, main="Species")
plot(vario_env, main="Environment")
plot(vario_env$v, vario_sp$v,
     xlab="Semivariance for environment",
     ylab="Semivariance for species")
m <- lm(vario_sp$v ~ vario_env$v-1)
abline(a=0, b=coef(m), col="red")
dev.off()

## Conclusions

# 1. Because mortality rate is equal for all species, species with the lowest habitat frequency have a higher probability to disappear
# 2. Stable coexistence (cf. species rank from one repetition to the other)
# 3. Species abundance at the end correlated with species habitat frequency
# 4. Species autocorrelation with correspondence between species and environment
# 5. For a given habitat distribution, the number of species at the equilibrium depends on the mortality rate
#    If the mortality rate is higher, a higher number of low abundance species disappear.   

# ===========================
# "Intraspecific variability"
# ===========================

# Data-set
df <- data.frame(perf_E_Sp)
names(df) <- sprintf("Sp_%03d", 1:50)
df_perf <- tibble(df) %>%
  mutate(Env=values(raster(env))) %>%
  mutate(Env2=Env^2) %>%
  pivot_longer(cols=Sp_001:Sp_050, names_to="Species", values_to="Perf")

# Observed niche
df_Sp <- df_perf %>% filter(Species=="Sp_025")
lm_fit <- lm(Perf~Env+Env2, data=df_Sp)
df_Sp_pred <- data.frame(Perf=predict(lm_fit, df_Sp), Env=df_Sp$Env)
ggplot(data=df_Sp, aes(x=Env, y=Perf)) +
  geom_point() +
  geom_line(data=df_Sp_pred, col="red")

# Observed intraspecific variability
lm_fit <- lm(Perf~Species, data=df_perf)
V_intra <- df_perf %>%
  mutate(res=lm_fit$residuals) %>%
  group_by(Species) %>%
  summarise(V=var(res))

# =========================
# End of file
# =========================