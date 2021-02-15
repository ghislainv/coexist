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
library(tidyr)
library(dplyr)
library(glue)

# Create output directories
dir.create("outputs")
dir.create("outputs/m2")

# Seed for reproducibility
seed <- 1234

# Figure width
fig_width <- 16.6 # in cm

# =========================
# Landscape and environment
# =========================

# Landscape
ncell_side <- 25
mat <- matrix(0, nrow=ncell_side, ncol=ncell_side)
r <- raster(mat, crs="+proj=utm +zone=1")
ncell <- ncell(r)
coords <- coordinates(r)

# Environment
load(here("outputs", "m1", "env.rda"))
env_m1 <- data.frame(rasterToPoints(raster(env)))
names(env_m1) <- c("x", "y", "env")
set.seed(seed)
env <- runif(ncell, min=0, max=1)
df <- env_m1[order(env_m1$env),]
df$env <- sort(env)
env <- as.matrix(rasterFromXYZ(df))

# Plot
png(file=here("outputs", "m2", "environment.png"))
plot(raster(env), main="Environment", col=topo.colors(255))
dev.off()

# Habitat frequency
png(file=here("outputs", "m2", "habitat_freq.png"))
hist(env, main="Environment", freq=FALSE)
dev.off()

# =========================================
# Species performance given the environment
# =========================================

# Number of species
nsp <- 50
# Remove last species (otherwise last species=first species)
perf <- seq(0, 1, length.out=nsp+1)[-(nsp+1)]

# Habitat frequency for each species
sp_hab_freq <- dunif(perf, min=0, max=1)
png(file=here("outputs", "m2", "species_habitat_freq.png"))
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Matrix of species performance on each cell (distances)
# Cell in rows, Species in columns
dist_E_Sp <- matrix(NA, nrow=ncell, ncol=nsp) 
for (i in 1:ncell) {
  env_i <- values(raster(env))[i] 
  for (j in 1:nsp) {
    dist <- abs(perf[j]-env_i)
    if (dist>0.5) {
      dist <- 1-dist
    }
    dist_E_Sp[i, j] <- dist
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

# Function to identify the species with the highest performance
high_perf_sp <- function(dist, sp_pres) {
  dist_pres <- dist[sp_pres]
  min_dist <- min(dist_pres)
  sp_high_perf <- sp_pres[which(dist_pres==min_dist)]
  # If two species adapted, selection at random
  if (length(sp_high_perf)==2) {
    # Random permutation
    sp_high_perf <- sample(sp_high_perf)
    sp_high_perf <- sp_high_perf[1]
  }
  return(sp_high_perf)
}

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)

# Species rank at the end of the generations
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
    png(file=here("outputs", "m2", "scene_start.png"))
    plot(raster(scene_start), main="Species - Start", zlim=c(0, 50),
         col=c("black", rev(terrain.colors(50))))
    dev.off()
  }
  scene <- scene_start

  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(scene)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(scene), levels=1:nsp))
  
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
      png(file=here("outputs", "m2", "mortality_events.png"))
      plot(raster(scene), main="Species - with vacant sites", zlim=c(0, 50),
           col=c("black", rev(terrain.colors(50))))
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
    
    sp_rich[g+1, r] <- length(unique(as.vector(scene)))
    abund[g+1, ] <- table(factor(as.vector(scene), levels=1:nsp))
    
  } # End ngen
  
  # Plot final scene once
  if (r==1) {
    png(file=here("outputs", "m2", "scene_end.png"))
    plot(raster(scene), main=glue("Species - End (ngen={ngen})"),
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
ggsave(p, filename=here("outputs", "m2", "species_richness_with_time.png"),
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
ggsave(p, filename=here("outputs", "m2", "mean_rank-habitat_freq.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

library(geoR)

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(scene)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=values(raster(env)))
# Plot with correlation
png(file=here("outputs", "m2", "sp_autocorrelation.png"),
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

# 1. The number of species depends on the environment heterogeneity

# =========================
# End of file
# =========================