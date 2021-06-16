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

# Create output directories
dir.create(here("outputs/m1"), recursive=TRUE)

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

# Charge files from the other model
load(file=here::here("outputs", "m0", "lm_fit.RData"))
load(file=here::here("outputs", "m0", "V_intra.RData"))


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
Sites <- data.frame(site=1:nsite, V1_env=NA)

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

# Environment on each site
sites <- data.frame(V1_env=rep(NA, nsite))
env <- data.frame()
rho <- c(rmvn(1, mu=rep(0, nsite), V=covrho, seed=seed)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero
rho <- scales::rescale(rho, to=c(0, 1))
sites <- rho
env <- matrix(rho, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
env_RGB <- raster(env*255)
crs(env_RGB) <- "+proj=utm +zone=1"

# Plot
png(file=here("outputs", "m1", "environment.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(1,1))
plot(raster(env), main="Environment var1", col=topo.colors(255))
dev.off()

# Habitat frequency
png(file=here("outputs", "m1", "hab_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env, main="", xlab="Environment var1")    
dev.off()

# =========================================
# Species niche
# =========================================

# # Niche width
# niche_width <- 0.25
# # Niches per axis
# n_niche <- 1/niche_width

# Number of species
 nsp <- nrow(V_intra)

# # Species coordinates on the three niche axis (x, y, z)
# base_coord <- seq(0, 1, length.out=n_niche+1)[-(n_niche+1)]+niche_width/2
# sp_x <- rep(rep(base_coord, n_niche), n_niche)
# niche_optimum <- as.data.frame(sp_x)
# 
# # Random optimum for species
# randomOptSp <- TRUE
# if (randomOptSp) {
#   set.seed(seed)
#   niche_optimum <- data.frame(sp_x=runif(nsp))
#   niche_optimum <- niche_optimum[order(niche_optimum$sp_x)]
# }
# rownames(niche_optimum)<-1:nrow(niche_optimum)
# 
# # Plot the species niche
# png(file=here("outputs", "m1", "species_niche.png"),
#     width=fig_width, height=fig_width, units="cm", res=300)
# par(mar=c(1,1,2,2))
# ggplot(niche_optimum,
#       aes(x=sp_x, colour=),
#       pch=16, 
#       col=viridis(nsp),
#           bty = "f", main ="Species niche", phi=0,
#           xlim=c(0,1))
# dev.off()

# Matrix of species performance on each site
# Sites in rows, Species in columns
beta_0 <- as.vector(lm_fit$coefficients[1:nsp])
beta_1 <- as.vector(c(lm_fit$coefficients[nsp+1], lm_fit$coefficients[(nsp+3):(2*nsp+1)]))
beta_2 <- as.vector(c(lm_fit$coefficients[nsp+2], lm_fit$coefficients[(2*nsp+2):(3*nsp)]))

perf_ind_init <- as.data.frame(matrix(ncol=nsp, nrow=nsite))

for(col in 1:nsp){
  sp <- col
  for (row in 1:nsite) {
    perf_ind_init[row,col] <- beta_0[sp]+beta_1[sp]*sites[row]+beta_2[sp]*(sites[row])^2+rnorm(1, mean=0, sd=sqrt(V_intra$V[sp]))
  }
}

# Identify the species with the highest performance (in each line find the max)
sp_high_perf <- apply(perf_ind_init, 1, which.max)

# Probability of dying of each species on each site
# Strength of unsuitability
b <- -0.5
mortality_E_Sp <- inv_logit(logit(0.1) + b * t(perf_ind_init))
# Mortality rate distribution
png(file=here("outputs", "m1", "hist_mortality.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(mortality_E_Sp)
dev.off()

# Habitat frequency for each species
# rank_dist_E <- t(apply(dist_E_Sp, 1, rank, ties.method="min"))
# sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
# sp_hab_freq <- as.table(sp_hab_freq)
# names(sp_hab_freq) <- 1:nsp
# png(file=here("outputs", "m1", "species_habitat_freq.png"),
#     width=fig_width, height=fig_width, units="cm", res=300)
# plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
# dev.off()
# 
# # Species with no habitat
# sp_no_habitat <- as.vector(which(sp_hab_freq==0))
# nsp_no_habitat <- length(sp_no_habitat)

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 10
# Number of generations
ngen <- 500

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Species rank at the end of the generations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
#env_filt <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen+1, ncol=nrep)

Shannon <- c(nrep)

# Loop on repetitions
for (r in 1:nrep) {
  
  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Draw species at random in the landscape (one individual per site)
  sp <- sample(1:nsp, size=nsite, replace=TRUE)
  #hist(sp)
  community_start <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m1", "community_start.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    plot(raster(community_start), main="Species - Start", zlim=c(0, nsp),
         col=c("black", viridis(nsp)))
    dev.off()
  }
  community <- community_start
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # Environmental filtering
  # dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
  # env_filt[1, r] <- mean(dist_site)
  
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
      png(file=here("outputs", "m1", "mortality_events.png"),
          width=fig_width, height=fig_width, units="cm", res=300)
      plot(raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
           col=c("black", viridis(nsp)))
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
    
    #Environmental conditions on vacant sites
    env_sites_vacant <- sites[sites_vacant]
    
    # Performance of individuals on vacant sites
    #dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
    perf_ind_pres <- as.data.frame(matrix(ncol=nsp_present, nrow=length(sites_vacant)))
    
    for(col in 1:nsp_present){
      sp <- sp_present[col]
      for (row in 1:length(sites_vacant)) {
        perf_ind_pres[row,col] <- beta_0[sp]+beta_1[sp]*env_sites_vacant[row]+beta_2[sp]*(env_sites_vacant[row])^2+rnorm(1, mean=0, sd=sqrt(V_intra$V[sp]))
      } 
    }

    # Identify the species with the highest performance
    
    new_ind <- sp_present[apply(perf_ind_pres, 1, which.max)]
    
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
    # dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
    # env_filt[g+1, r] <- mean(dist_site)
    
    # Mean mortality rate in the community
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    theta_comm[g+1, r] <- mean(theta_site)
    
  } # End ngen
  
  # Plot final community once
  if (r==1) {
    png(file=here("outputs", "m1", "community_end.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    plot(raster(community), main=glue("Species - End (ngen={ngen})"),
         zlim=c(0, nsp), col=c("black", viridis(nsp)))
    dev.off()
  }
  
  # Species rank
  rank_sp[r, ] <- rank(-abund[ngen+1, ], ties.method="min")
  
  df_shannon <- data.frame(Species = 1:64,
                           Abundance = abund[ngen+1, ])%>%
    mutate(Proportion = Abundance / sum(Abundance))%>%
    filter(Abundance > 0)%>%
    mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
  
  Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
  
} # End nrep

#######################################################################################
#Plots#

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
  scale_colour_viridis_d()+
  xlab("Generations") + 
  ylab("Species richness")
ggsave(p, filename=here("outputs", "m1", "species_richness_with_time.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ---------------------------------------------------------------------------------
# pairwise Spearman correlation on the species ranks at the end of each simulation
# ---------------------------------------------------------------------------------

Spearman <- as.dist(round(cor(t(rank_sp), method="spearman"),2))
save(Spearman, file = here::here("outputs", "m1", "Spearman_m1.RData"))

rank_sp <- data.frame(rank_sp)%>%
  mutate(rep = 1:nrep)%>%
  pivot_longer(cols=X1:X64, values_to="rank_sp")
p <- ggplot(data=rank_sp, aes(x=rep, y=rank_sp)) +
  geom_line(aes(colour=name)) + 
  scale_colour_viridis_d()+
  xlab("Repetition") + 
  ylab("Species rank")+
  theme(legend.position = "none")
ggsave(p, filename=here("outputs", "m1", "species_rank_with_repetitions.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)


sp_rich_final <- sp_rich[ngen+1,]
save(sp_rich_final, file=here::here("outputs", "m1", "Species_richness_m1.RData"))


# ---------------------------------------------
# Shannon index and Shannon equitability index
# ---------------------------------------------
save(Shannon, file = here::here("outputs", "m1", "Shannon_m1.RData"))
Equitability <- Shannon/log(as.numeric(sp_rich[ngen+1,]))
save(Equitability, file = here::here("outputs", "m1", "Equitability_m1.RData"))

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# # Mean final rank
# sp_mean_rank <- apply(rank_sp, 2, mean)
# # Plot
# df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
# p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
#   geom_point() +
#   geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
#   xlab("Species habitat frequency") +
#   ylab("Species mean rank (higher rank = lower abundance)") +
#   theme(axis.title=element_text(size=16))
# ggsave(p, filename=here("outputs", "m1", "mean_rank-habitat_freq.png"),
#        width=fig_width, height=fig_width, units="cm", dpi=300)

# ---------------------------------------------
# Environmental filtering
# ---------------------------------------------

# env_filt <- data.frame(env_filt)
# env_filt_long <- env_filt %>%
#   mutate(gen=1:(ngen+1)) %>%
#   pivot_longer(cols=X1:X10, names_to="rep",
#                names_prefix="X", values_to="env_filt")
# p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
#   geom_line() +
#   labs(title="Environmental filtering") +
#   xlab("Generations") + 
#   ylab("Mean env-species perf difference")
# ggsave(p, filename=here("outputs", "m1", "environmental_filtering.png"),
#        width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Plot
png(file=here("outputs", "m1", "spatial_comp_env_sp.png"), 
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2))
plot(raster(community_start), main="Species - Start", zlim=c(0, nsp),
     col=c("black", viridis(nsp)))
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", viridis(nsp)))
plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", viridis(nsp)))
dev.off()

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
  scale_colour_viridis_d()+
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean mortality rate in the community")
ggsave(p, filename=here("outputs", "m1", "mortality_rate_community.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(community)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sites)

# Plot with correlation
png(file=here("outputs", "m1", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2))
plot(vario_sp, main="Species - End")
plot(vario_env, main="Environment")
plot(vario_env$v, vario_sp$v,
     xlab="Semivariance for environment",
     ylab="Semivariance for species")
m <- lm(vario_sp$v ~ vario_env$v-1)
abline(a=0, b=coef(m), col="red")
dev.off()

################################
# Comparison of the two models #
################################

#Species richness

load(here::here("outputs", "m0", "Species_richness_m0.RData"))
sp_rich_final_m0 <- sp_rich_final

load(here::here("outputs", "m1", "Species_richness_m1.RData"))
sp_rich_final_m1 <- sp_rich_final

sp_rich_final <- data.frame(Sp_rich = c(as.numeric(sp_rich_final_m0), as.numeric(sp_rich_final_m1)),
                            Model = as.factor(rep(1:2, each=length(sp_rich_final_m0))))
g <- ggplot(sp_rich_final, aes(x= Model, y=Sp_rich, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the species richness \n at the end of each repetition with both models", y="Species richness")

ggsave(g, filename=here("outputs", "Comparison", "Species_richness.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Shannon diversity and equitability index

load(here::here("outputs", "m0", "Shannon_m0.RData"))
Shannon_m0 <- Shannon

load(here::here("outputs", "m1", "Shannon_m1.RData"))
Shannon_m1 <- Shannon

Shannon <- data.frame(Shannon = c(Shannon_m0, Shannon_m1),
                            Model = as.factor(rep(1:2, each=length(Shannon_m0))))
g <- ggplot(Shannon, aes(x= Model, y=Shannon, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon diversity index \n at the end of each repetition with both models", y="Shannon diversity index")

ggsave(g, filename=here("outputs", "Comparison", "Shannon.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

load(here::here("outputs", "m0", "Equitability_m0.RData"))
Equitability_m0 <- Equitability

load(here::here("outputs", "m1", "Equitability_m1.RData"))
Equitability_m1 <- Equitability

Equitability <- data.frame(Equitability = c(Equitability_m0, Equitability_m1),
                      Model = as.factor(rep(1:2, each=length(Equitability_m0))))

g <- ggplot(Equitability, aes(x= Model, y=Equitability, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon equitability index \n at the end of each repetition with both models", y="Shannon equitability index")

ggsave(g, filename=here("outputs", "Comparison", "Equitability.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Spearman pairwise correlation of species ranks

load(here::here("outputs", "m0", "Spearman_m0.RData"))
Spearman_m0 <- Spearman

load(here::here("outputs", "m1", "Spearman_m1.RData"))
Spearman_m1 <- Spearman

Spearman <- data.frame(Spearman = c(as.vector(Spearman_m0), as.vector(Spearman_m1)),
                           Model = as.factor(rep(1:2, each=length(Spearman_m0))))

g <- ggplot(Spearman, aes(x= Model, y=Spearman, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the pairwise Spearman correlation of species ranks \n at the end of each repetition with both models", y="Shannon equitability index")
ggsave(g, filename=here("outputs", "Comparison", "Spearman_ranks.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

