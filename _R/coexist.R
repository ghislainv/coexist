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
library(geoR) # variog()
library(plot3D) # scatter3D()
library(Rcpp)
library(RcppArmadillo)
library(viridisLite)
library(purrr)

# Create output directories
dir.create(here("outputs/m0"), recursive=TRUE)

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

# Cpp function to compute distance between Sites and Species
Rcpp::sourceCpp(here("_src", "dist_Site_Sp.cpp"))

# =========================
# Landscape and environment
# =========================

# Landscape
nsite_side <- 25
mat <- matrix(0, nrow=nsite_side, ncol=nsite_side)
r <- raster(mat, crs="+proj=utm +zone=1")
nsite <- ncell(r)
coords <- coordinates(r)

# Environment on each site
Sites <- data.frame(site=1:nsite, V1_env=NA, V2_env=NA, V3_env=NA)

# Neighbourhood matrix
neighbors.mat <- adjacent(r, cells=c(1:nsite), directions=8,
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
Q <- diag(n.neighbors)-d*A + diag(.0001,nsite) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos

# Number of axis for the niche
n_axis <- 3

# Environment on each site
sites <- data.frame(V1_env=rep(NA, nsite), V2_env=NA, V3_env=NA)
env <- list()
for (i in 1:n_axis) {
  seed <- 1234 + i - 1
  rho <- c(rmvn(1, mu=rep(0, nsite), V=covrho, seed=seed)) # Spatial Random Effects
  rho <- rho-mean(rho) # Centering rhos on zero
  #rho <- scales::rescale(rho, to=c(0, 1))
  rho <- (rho - min(rho)) / (max(rho) - min(rho))
  sites[,i] <- rho
  env[[i]] <- matrix(rho, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
}
env_stack <- stack(raster(env[[1]]*255), raster(env[[2]]*255), raster(env[[3]]*255))
crs(env_stack) <- "+proj=utm +zone=1"
#CGT 25/06/2021: save the environmental data to avoid simulating them again in the next model
save(sites, file = here::here("outputs", "m0", "sites.RData"))
save(env, file = here::here("outputs", "m0", "env.RData"))

# Plot
png(file=here("outputs", "m0", "environment.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(raster(env[[1]]), main="Environment var1", col=topo.colors(255), cex.main=1.2)
plot(raster(env[[2]]), main="Environment var2", col=topo.colors(255), cex.main=1.2, legend=FALSE)
plot(raster(env[[3]]), main="Environment var3", col=topo.colors(255), cex.main=1.2, legend=FALSE)
plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
dev.off()

# Habitat frequency
png(file=here("outputs", "m0", "hab_freq_1.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env[[1]], main="", xlab="Environment var1")    
dev.off()

png(file=here("outputs", "m0", "hab_freq_2.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env[[2]], main="", xlab="Environment var2")    
dev.off()

png(file=here("outputs", "m0", "hab_freq_3.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env[[3]], main="", xlab="Environment var3")    
dev.off()


# =========================================
# Species niche
# =========================================

# Niche width
niche_width <- 0.25
# Niches per axis
n_niche <- 1/niche_width
# Number of species
nsp <- n_niche^n_axis
# Species coordinates on the three niche axis (x, y, z)
nsp <- 100
base_coord <- seq(0, 1, length.out=n_niche+1)[-(n_niche+1)]+niche_width/2
sp_x <- rep(rep(base_coord, n_niche), n_niche)
sp_y <- rep(rep(base_coord, each=n_niche), n_niche)
sp_z <- rep(base_coord, each=n_niche^2)
niche_optimum <- as.data.frame(cbind(sp_x, sp_y, sp_z))

# Random optimum for species
randomOptSp <- TRUE
if (randomOptSp) {
  set.seed(seed)
  niche_optimum <- data.frame(sp_x=runif(n=nsp, min = min(env[[1]]), max = max(env[[1]])),
                              sp_y=runif(n=nsp, min = min(env[[2]]), max = max(env[[2]])),
                              sp_z=runif(n=nsp, min = min(env[[3]]), max = max(env[[3]])))
  niche_optimum <- niche_optimum[order(niche_optimum$sp_z, niche_optimum$sp_y, niche_optimum$sp_x), ]
}
#CGT 10/06/2021
rownames(niche_optimum)<-1:nrow(niche_optimum)

# Plot the species niche
png(file=here("outputs", "m0", "species_niche.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mar=c(1,1,2,2))
scatter3D(niche_optimum$sp_x, niche_optimum$sp_y, niche_optimum$sp_z,
          pch=16, 
          colvar=1:nsp, col=viridis(nsp),
          bty = "f", main ="Three-dimensional species optima", phi=0,
          xlim=c(0,1), ylim=c(0,1), zlim=c(0,1))
dev.off()

# Matrix of species performance on each site (distances)
# Sites in rows, Species in columns
dist_E_Sp <- dist_Site_Sp(as.matrix(sites), as.matrix(niche_optimum))
dprim_E_Sp <- (dist_E_Sp-mean(dist_E_Sp))/sd(dist_E_Sp)
perf_E_Sp <- -dprim_E_Sp

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

# Probability of dying of each species on each site
# Strength of unsuitability
b <- -0.5
mortality_E_Sp <- inv_logit(logit(0.1) + b * perf_E_Sp)
# Mortality rate distribution
png(file=here("outputs", "m0", "hist_mortality.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(mortality_E_Sp)
dev.off()
# Mortality probability function
png(file=here::here("outputs", "m0", "function_mort_proba.png"),
    width=fig_width, height = fig_width/1.5, units="cm", res=300)
plot(x=c(perf_E_Sp),
     y=c(mortality_E_Sp),
     xlab="Performance",
     ylab="Mortality probability",
     main="Mortality probability as a function of performance",
     cex.lab=1.5,
     cex.main=1.5)
dev.off()

# Habitat frequency for each species
rank_dist_E <- t(apply(dist_E_Sp, 1, rank, ties.method="min"))
sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
sp_hab_freq <- as.table(sp_hab_freq)
names(sp_hab_freq) <- 1:nsp
save(sp_hab_freq, file = here::here("outputs", "m0", "sp_hab_freq.RData"))
png(file=here("outputs", "m0", "species_habitat_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(cex.lab=1.5)
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Species with no habitat
sp_no_habitat <- as.vector(which(sp_hab_freq==0))
nsp_no_habitat <- length(sp_no_habitat)

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 50
# Number of generations
ngen <- 1000

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Species rank at the end of the generations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
env_filt <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen+1, ncol=nrep)
#CGT 16/06/2021
Shannon <- c(nrep)
#CGT 17/06/2021
#To keep the abundance matrixes in order to infer alpha matrix
Abundances_m0<-list()

# Loop on repetitions
for (r in 1:nrep) {

  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Draw species at random in the landscape (one individual per site)
  sp <- sample(1:nsp, size=nsite, replace=TRUE)
  #hist(sp)
  community_start_m0 <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m0", "community_start.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    par(bty = "n")
    plot(raster(community_start_m0), main="Species - Start", zlim=c(0, nsp),
         col=c("black",  viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
    dev.off()
  }
  community <- community_start_m0
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # Environmental filtering
  dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
  env_filt[1, r] <- mean(dist_site)
  
  # Mean mortality rate
  theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
  theta_comm[1, r] <- mean(theta_site)

  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # ******************
    # Mortality
    # ******************
    
    # Mortality rate on each site
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    
    # Mortality events
    mort_ev <- rbinom(nsite, size=1, prob=theta_site)
    mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    
    # Update community
    community[mortality==1] <- 0
    # Plot once
    if (r==1 & g==1) {
      png(file=here("outputs", "m0", "mortality_events.png"),
          width=fig_width, height=fig_width, units="cm", res=300)
      par(bty = "n")
      plot(raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
           col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
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
    dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
    env_filt[g+1, r] <- mean(dist_site)

    # Mean mortality rate in the community
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    theta_comm[g+1, r] <- mean(theta_site)
    
  } # End ngen
  
  #CGT 28/06/2021 to compare the two models
  community_end_m0 <- community
  
  # Plot final community once
  if (r==1) {
    png(file=here("outputs", "m0", "community_end.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    par(bty = "n")
    plot(raster(community), main=glue("Species - End (ngen={ngen})"),
         zlim=c(0, nsp), col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
    dev.off()
  }
  
  # Species rank
  #CGT 16/06/2021 : ngen --> ngen +1
  rank_sp[r, ] <- rank(-abund[ngen+1, ], ties.method="min")
  
  #CGT 16/06/2021
  df_shannon <- data.frame(Species = 1:nsp,
                           Abundance = abund[ngen+1, ])%>%
    mutate(Proportion = Abundance / sum(Abundance))%>%
    filter(Abundance > 0)%>%
    mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
  
  Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
  
  #CGT 17/06/2021
  #To keep the abundance matrixes in order to infer alpha matrix
  Abundances_m0[[r]] <- abund
  
} # End nrep

save(Abundances_m0, file = here::here("outputs", "m0", "Abundances_m0.RData"))
save(sp_rich, file=here::here("outputs", "m0", "sp_rich_m0.RData"))
save(rank_sp, file=here::here("outputs", "m0", "rank_sp_m0.RData"))

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
  pivot_longer(cols=X1:X50, names_to="rep",
               names_prefix="X", values_to="sp_rich")
p <- ggplot(data=sp_rich_long, aes(x=gen, y=sp_rich, col=rep)) +
  geom_line() +
  scale_colour_viridis_d()+
  xlab("Generations") + 
  ylab("Species richness")+
  theme(legend.position = "none",
        text = element_text(size = 20))
ggsave(p, filename=here("outputs", "m0", "species_richness_with_time.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

#CGT 16/06/2021
sp_rich_final <- sp_rich[ngen+1,]
save(sp_rich_final, file=here::here("outputs", "m0", "Species_richness_m0.RData"))

#CGT 16/06/2021
# ---------------------------------------------
# Shannon index and Shannon equitability index
# ---------------------------------------------
save(Shannon, file = here::here("outputs", "m0", "Shannon_m0.RData"))
Equitability <- Shannon/log(as.numeric(sp_rich[ngen+1,]))
save(Equitability, file = here::here("outputs", "m0", "Equitability_m0.RData"))

#CGT 16/06/2021
# ---------------------------------------------------------------------------------
# pairwise Spearman correlation on the species ranks at the end of each simulation
# ---------------------------------------------------------------------------------

Spearman <- as.dist(round(cor(t(rank_sp), method="spearman"),2))
save(Spearman, file = here::here("outputs", "m0", "Spearman_m0.RData"))

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# Mean final rank
sp_mean_rank <- apply(rank_sp, 2, mean)
# Plot
df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
  geom_point() +
  geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
  xlab("Species suitable habitat frequency") +
  ylab("Species mean rank (higher rank = lower abundance)") +
  theme(axis.title=element_text(size=16))
ggsave(p, filename=here("outputs", "m0", "mean_rank-habitat_freq.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

#CGT 15/06/2021
# To compare with the other model
Df_rank_sp <- data.frame(rank_sp)%>%
  mutate(rep = 1:nrep)%>%
  pivot_longer(cols=X1:X100, values_to="rank_sp")
p <- ggplot(data=Df_rank_sp, aes(x=rep, y=rank_sp)) +
  geom_line(aes(colour=name)) +
  scale_colour_viridis_d()+
  xlab("Repetition") + 
  ylab("Species rank")+
  theme(legend.position = "none",
        text = element_text(size = 20))
ggsave(p, filename=here("outputs", "m0", "species_rank_with_repetitions.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)


# ---------------------------------------------
# Environmental filtering
# ---------------------------------------------

env_filt <- data.frame(env_filt)
env_filt_long <- env_filt %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X50, names_to="rep",
               names_prefix="X", values_to="env_filt")
p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
  geom_line() +
  scale_colour_viridis_d()+
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean env-species perf difference")
ggsave(p, filename=here("outputs", "m0", "environmental_filtering.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Plot
png(file=here("outputs", "m0", "spatial_comp_env_sp.png"), 
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(raster(community_start_m0), main="Species - Start", zlim=c(0, nsp),
     col=c("black", viridis(nsp)), legend=FALSE)
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", viridis(nsp)), legend=FALSE)
plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", viridis(nsp)), legend=FALSE)
dev.off()

# ---------------------------------------------
# Theta community
# ---------------------------------------------

theta_comm <- data.frame(theta_comm)
theta_comm_long <- theta_comm %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X50, names_to="rep",
               names_prefix="X", values_to="theta_comm")
p <- ggplot(data=theta_comm_long, aes(x=gen, y=theta_comm, col=rep)) +
  geom_line() +
  scale_colour_viridis_d()+
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean mortality rate in the community")
ggsave(p, filename=here("outputs", "m0", "mortality_rate_community.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(community)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation
# 3D voxel for each site
x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
# Plot with correlation
png(file=here("outputs", "m0", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(vario_sp, main="Species - End")
plot(vario_env, main="Environment")
plot(vario_env$v, vario_sp$v,
     main = "Regression",
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
names(df) <- c(sprintf("Sp_%03d", 1:(nsp-1)), sprintf("Sp_%d", nsp))
df_perf <- tibble(df) %>%
  mutate(Env=values(raster(env[[1]]))) %>%
  mutate(Env2=Env^2) %>%
  pivot_longer(cols=c(Sp_001:glue("Sp_0{nsp-1}"), glue("Sp_{nsp}")), names_to="Species", values_to="Perf")

# Observed niche
# Select 8 species at random
set.seed(seed)
sp_sel <- sample(unique(df_perf$Species), 9, replace=FALSE)
df_sp_sel <- df_perf %>% filter(Species %in% sp_sel)
p <- ggplot(data=df_sp_sel, aes(x=Env, y=Perf)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE) +
  facet_wrap(vars(Species), nrow=3) +
  xlab("Environment (first axis)") +
  ylab("Performance")
ggsave(p, filename=here("outputs", "m0", "infering_species_niche.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)
  
# Observed intraspecific variability
lm_fit <- lm(Perf~Species+Species*Env+Species*Env2, data=df_perf)
save(lm_fit, file = here::here("outputs", "m0","lm_fit.RData"))

summary(lm_fit)$adj.r.squared
hist(summary(lm_fit)$residuals)
shapiro.test(sample(summary(lm_fit)$residuals, size = 5000))
ks.test(x=summary(lm_fit)$residuals, y='pnorm',alternative='two.sided')

V_intra <- df_perf %>%
  mutate(res=lm_fit$residuals) %>%
  group_by(Species) %>%
  summarise(V=var(res))
save(V_intra, file = here::here("outputs", "m0", "V_intra.RData"))

p <- ggplot(data=V_intra, aes(x=Species, y=V)) +
  geom_col() +
  theme(axis.text.x=element_text(angle=90, size=6),
        text = element_text(size = 20)) +
  ylab("Intraspecific variance")
ggsave(p, filename=here("outputs", "m0", "intraspecific_variance.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# CGT 11/06/2021
#At mean environment, x1 = ~ 0.5 (see generation of the environment)
mean_x1 <- mean(sites[,1])

#Species coefficients
beta_0 <- as.vector(lm_fit$coefficients[1:nsp])
beta_1 <- as.vector(c(lm_fit$coefficients[nsp+1], lm_fit$coefficients[(nsp+3):(2*nsp+1)]))
beta_2 <- as.vector(c(lm_fit$coefficients[nsp+2], lm_fit$coefficients[(2*nsp+2):(3*nsp)]))

x <- seq(-3, 3, 0.01)
n_ind_simul <- length(x)
df_perf_IV <- data.frame(Species = rep(1:nsp, each=n_ind_simul), X = rep(x, nsp), Mean=rep(beta_0+beta_1*mean_x1+beta_2*(mean_x1^2), each=n_ind_simul))
df_perf_IV <- df_perf_IV%>%
  group_by(Species)%>%
  mutate(Perf = dnorm(x, mean=Mean, sd=sqrt(V_intra$V[Species])))%>%
  ungroup()

p <- ggplot(df_perf_IV, aes(x = X, y = Perf, colour = as.factor(Species))) +
  geom_line()+
  scale_colour_viridis_d()+
  theme(legend.position = "none")+
  labs(x = "Performance at mean environment variable X1",
       y = "Density")
ggsave(p, filename=here("outputs", "m0", "Perf_overlap_IV.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Compare species habitat frequency and IV #

Niche_width <- data.frame(Species = c(1:nsp),
                       Optimum = niche_optimum$sp_x,
                       Width = numeric(nsp))
Niche_width <- Niche_width[order(Niche_width$Optimum), ]
rownames(Niche_width) <- c(1:nrow(Niche_width))
# Width of the first niche
Niche_width$Width[1] <- Niche_width$Optimum[1] + 0.5 * (Niche_width$Optimum[2] - Niche_width$Optimum[1])
# Width of the last niche
Niche_width$Width[nsp] <- (1-Niche_width$Optimum[nsp]) + 0.5 * (Niche_width$Optimum[nsp] - Niche_width$Optimum[nsp-1])
#Width of all other niches
for (k in 2:(nsp-1)){
  Niche_width$Width[k] <- (Niche_width$Optimum[k]-Niche_width$Optimum[k-1])/2 + (Niche_width$Optimum[k+1]-Niche_width$Optimum[k])/2
}
Niche_width <- Niche_width[order(Niche_width$Species), ]
rownames(Niche_width) <- c(1:nrow(Niche_width))

#
load(here::here("outputs", "m0", "V_intra.RData"))
load(here::here("outputs", "m0", "sp_hab_freq.RData"))

V_intra_hab_freq <- data.frame(Species = c(1:nsp),
                               IV = V_intra$V,
                               Hab_freq = as.data.frame(sp_hab_freq)$Freq,
                               Pres_m0 = numeric(nsp),
                               Pres_m1 = numeric(nsp),
                               Code_presence = numeric(nsp),
                               Niche_width = Niche_width$Width)

load(here::here("outputs", "m0", "Abundances_m0.RData"))
for (sp in 1:nsp) {
  if (Abundances_m0[[1]][nrow(Abundances_m0[[1]]), sp] == 0) {
    V_intra_hab_freq$Pres_m0[sp] <- 0
  }
  else V_intra_hab_freq$Pres_m0[sp] <- 1
}

load(here::here("outputs", "m1", "Abundances_m1.RData"))
for (sp in 1:nsp) {
  if (Abundances_m1[[1]][nrow(Abundances_m1[[1]]), sp] == 0) {
    V_intra_hab_freq$Pres_m1[sp] <- 0
  }
  else V_intra_hab_freq$Pres_m1[sp] <- 1
}

for (k in 1:nrow(V_intra_hab_freq)){
  if (V_intra_hab_freq$Pres_m0[k]==1 && V_intra_hab_freq$Pres_m1[k]==1){
    V_intra_hab_freq$Code_presence[k] <- 3
  }
  if (V_intra_hab_freq$Pres_m0[k]==0 && V_intra_hab_freq$Pres_m1[k]==0){
    V_intra_hab_freq$Code_presence[k] <- 0
  }
  if (V_intra_hab_freq$Pres_m0[k]==1 && V_intra_hab_freq$Pres_m1[k]==0){
    V_intra_hab_freq$Code_presence[k] <- 1
  }
  if (V_intra_hab_freq$Pres_m0[k]==0 && V_intra_hab_freq$Pres_m1[k]==1){
    V_intra_hab_freq$Code_presence[k] <- 2
  }
}

V_intra_hab_freq$Code_presence <- as.factor(V_intra_hab_freq$Code_presence)

ggplot2::ggplot(data = V_intra_hab_freq, aes(x = Hab_freq, y = IV, colour = Code_presence))+
  geom_point()+
  scale_color_viridis_d("Presence of the species", labels = c("Disappeared in both models", "Disappeared in second model", "Disappeared in first model", "Maintained in both models"))+
  labs(title = "Relationship between the species suitable habitat frequency in m0, \n the IV on axis X1 and the species presence at the end of a simulation",
       x = "Suitable habitat frequency in m0",
       y = "Inferred IV on axis X1")+
  theme(text = element_text(size = 15))

ggplot2::ggplot(data = V_intra_hab_freq, aes(x = Niche_width, y = IV, colour = Code_presence))+
  geom_point()+
  scale_color_viridis_d("Presence of the species", labels = c("Disappeared in both models", "Disappeared in second model", "Maintained in both models"))+
  labs(title = "Relationship between the species suitable habitat frequency in m0, \n the IV on axis X1 and the species presence at the end of a simulation",
       x = "Niche width in m0",
       y = "Inferred IV on axis X1")
  

#Find the time of disppearing of each species
abund_m0 <- Abundances_m0[[1]]
Time_disap_m0 <- c()
for (k in 1:ncol(abund_m0)){
t <- abund_m0[,k]%>%purrr::detect_index(function(x) x == 0)
Time_disap_m0 <- c(Time_disap_m0, t)
}

abund_m1 <- Abundances_m1[[1]]
Time_disap_m1 <- c()
for (k in 1:ncol(abund_m1)){
  t <- abund_m1[,k]%>%purrr::detect_index(function(x) x == 0)
  Time_disap_m1 <- c(Time_disap_m1, t)
}

V_intra_hab_freq$Time_disap_m0 <- Time_disap_m0
V_intra_hab_freq$Time_disap_m1 <- Time_disap_m1

V_intra_hab_freq[which(V_intra_hab_freq$Time_disap_m0==0),]$Time_disap_m0 <- ngen
V_intra_hab_freq[which(V_intra_hab_freq$Time_disap_m1==0),]$Time_disap_m1 <- ngen

ggplot2::ggplot(data=V_intra_hab_freq, ggplot2::aes(x=IV, y=Time_disap_m0))+
  ggplot2::geom_point()+
  ggplot2::labs(x="IV", y="Time of species disappearing", title = "m0, 1 repetition, 1000 generations")

ggplot2::ggplot(data=V_intra_hab_freq, ggplot2::aes(x=Hab_freq, y=Time_disap_m0))+
  ggplot2::geom_point()+
  ggplot2::labs(x="Species habitat frequency", y="Time of species disappearing", title = "m0, 1 repetition, 1000 generations")

ggplot2::ggplot(data=V_intra_hab_freq, ggplot2::aes(x=IV, y=Time_disap_m1))+
  ggplot2::geom_point()+
  ggplot2::labs(x="IV", y="Time of species disappearing", title = "m1, 1 repetition, 1000 generations")

ggplot2::ggplot(data=V_intra_hab_freq, ggplot2::aes(x=Hab_freq, y=Time_disap_m1))+
  ggplot2::geom_point()+
  ggplot2::labs(x="Species habitat frequency", y="Time of species disappearing", title = "m1, 1 repetition, 1000 generations")

#Linear model of time of disppearance as a function of IV, habitat frequency and the interaction of both

lm_fit_m0 <- lm(Time_disap_m0~IV+Hab_freq+IV*Hab_freq, data=V_intra_hab_freq[Time_disap_m0!=0,])
lm_fit_m1 <- lm(Time_disap_m1~IV+Hab_freq+IV*Hab_freq, data=V_intra_hab_freq[Time_disap_m1!=0,])
summary(lm_fit_m0)
summary(lm_fit_m1)


# =========================
# End of file
# =========================