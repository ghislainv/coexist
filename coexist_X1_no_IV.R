
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
dir.create(here("outputs/m2"), recursive=TRUE)

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
load(here::here("outputs", "m0", "sites.RData"))
load(here::here("outputs", "m0", "env.RData"))

# =========================
# Landscape and environment
# =========================

# Landscape
nsite_side <- sqrt(nrow(sites))
nsite <- nrow(sites)

#Environment (axis 1)
X1 <- as.vector(sites[,1])

# Plot
png(file=here("outputs", "m2", "environment.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(1,1), bty = "n")
plot(raster(env[[1]]), main="Environment var1", col=topo.colors(255), cex.main=2)
dev.off()

# Habitat frequency
png(file=here("outputs", "m2", "hab_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env[[1]], main="", xlab="Environment var1")    
dev.off()

# =========================================
# Species parameters
# =========================================

# Number of species
nsp <- nrow(V_intra)

beta_0 <- as.vector(lm_fit$coefficients[1:nsp])
beta_1 <- as.vector(c(lm_fit$coefficients[nsp+1], lm_fit$coefficients[(nsp+3):(2*nsp+1)]))
beta_2 <- as.vector(c(lm_fit$coefficients[nsp+2], lm_fit$coefficients[(2*nsp+2):(3*nsp)]))

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 50
# Number of generations
ngen <- 1000

# Species richness
sp_rich_m2 <- matrix(NA, nrow=ngen, ncol=nrep)
# Species rank at the end of the generations
rank_sp_m2 <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
env_filt <- matrix(NA, nrow=ngen, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen, ncol=nrep)
#Shannon index of the community at the end of each repetition
Shannon <- c(nrep)
#list of the matrices of the species abundance through the generations for each repetition
Abundances_m2<-list()

# Matrix of mean species performance on each site
# Sites in rows, Species in columns
rank_dist_E <- t(apply(-perf_Sp_mean, 1, rank, ties.method="min"))
sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
sp_hab_freq <- as.table(sp_hab_freq)
names(sp_hab_freq) <- 1:nsp
png(file=here("outputs", "m2", "species_habitat_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Species with no habitat
sp_no_habitat <- as.vector(which(sp_hab_freq==0))
nsp_no_habitat <- length(sp_no_habitat)

# b = strength of unsuitability for mortality
b <- -0.5
# theta = basal mortality
theta <- 0.1

# Environmental filtering
env_filt <- matrix(NA, nrow=ngen, ncol=nrep)

# Loop on repetitions
for (r in 1:nrep) {
  
  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Matrix of species performance on each site
  # Sites in rows, Species in columns
  
  # Mean performance of each species on each site
  beta_0_mat <- matrix(rep(beta_0,each=nsite), ncol=nsp)
  beta_1_mat <- matrix(rep(beta_1,each=nsite), ncol=nsp)
  beta_2_mat <- matrix(rep(beta_2,each=nsite), ncol=nsp)
  X1_mat <- matrix(rep(X1, nsp), ncol=nsp)
  perf_Sp_mean <- beta_0_mat+beta_1_mat*X1_mat+beta_2_mat*X1_mat^2
  
  perf_ind <- perf_Sp_mean
  
  # Draw species at random in the landscape (one individual per site)
  sp <- sample(1:nsp, size=nsite, replace=TRUE)
  #hist(sp)
  community_start_m2 <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m2", "community_start.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    par(bty = "n")
    plot(raster(community_start_m2), main="Species - Start", zlim=c(0, nsp),
         col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
    dev.off()
  }
  community <- community_start_m2
  
  # Environmental filtering
  # distance between site environment and performance
  dist_site <- sqrt((diag(as.matrix(perf_ind)[, as.vector(t(community))])-X1)^2)
  #  env_filt[1, r] <- mean(dist_site, na.rm = TRUE)
  
  # Species richness
  #  sp_rich_m2[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen)
  #  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # ******************
    # Mortality
    # ******************
    
    # Probability of dying of each species on each site
    mortality_ind <- inv_logit(logit(theta) + b * as.matrix(perf_ind))

    # Mortality rate distribution
    # Plot once
    if (r==1 & g==1) {
      png(file=here("outputs", "m2", "hist_mortality.png"),
          width=fig_width, height=fig_width, units="cm", res=300)
      hist(mortality_ind)
      dev.off()
    }
    
    # Mortality rate on each site
    theta_site <- diag(mortality_ind[, as.vector(t(community))])
    
    # Mortality events
    mort_ev <- rbinom(nsite, size=1, prob=theta_site)
    mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    
    # Number of deaths
    n_mort <- sum(mort_ev)
    
    if(n_mort!=0){
      
      # Update community
      community[mortality==1] <- 0
      # Plot once
      if (r==1 & g==1) {
        png(file=here("outputs", "m2", "mortality_events.png"),
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
      nsite_vacant <- length(sites_vacant)
      
      # Identify the present species with the highest performance on vacant sites (maximum in each line)
      sp_high_perf <- sp_present[apply(matrix(perf_ind[sites_vacant, sp_present], ncol=nsp_present), 1, which.max)]
      
      # Recruitment
      community_rast[sites_vacant] <- sp_high_perf
      community <- as.matrix(community_rast)
      
    }  # end condition on n_mort
    
    
    # *********************
    # Diversity
    # *********************
    
    # Species richness
    sp_rich_m2[g, r] <- length(unique(as.vector(community)))
    abund[g, ] <- table(factor(as.vector(community), levels=1:nsp))
    
    # Environmental filtering
    dist_site <- sqrt((diag(as.matrix(perf_ind)[, as.vector(t(community))])-X1)^2)
    env_filt[g, r] <- mean(dist_site, na.rm = TRUE)
    
    # Mean mortality rate in the community
    theta_site <- diag(mortality_ind[, as.vector(t(community))])
    theta_comm[g, r] <- mean(theta_site, na.rm = TRUE)
    
  } # End ngen
  
  #CGT 28/06/2021 to compare the two models
  community_end_m2 <- community
  
  # Plot final community once
  if (r==1) {
    png(file=here("outputs", "m2", "community_end.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    par(bty = "n")
    plot(raster(community), main=glue("Species - End (ngen={ngen})"),
         zlim=c(0, nsp), col=c("black", viridis(nsp)), cex.main = 2, cex.axis = 1.5, legend=FALSE)
    dev.off()
  }
  
  # Species rank
  rank_sp_m2[r, ] <- rank(-abund[ngen, ], ties.method="min")
  
  df_shannon <- data.frame(Species = 1:nsp,
                           Abundance = abund[ngen, ])%>%
    mutate(Proportion = Abundance / sum(Abundance))%>%
    filter(Abundance > 0)%>%
    mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
  
  Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
  
  #To keep the abundance matrices in order to infer alpha matrix
  Abundances_m2[[r]] <- abund
  
} # End nrep

save(Abundances_m2, file = here::here("outputs", "m2", "Abundances_m2.RData"))
save(sp_rich_m2, file = here::here("outputs", "m2", "sp_rich_m2.RData"))
save(rank_sp_m2, file = here::here("outputs", "m2", "rank_sp_m2.RData"))

#######################################################################################
#Plots#

sp_rich_m2
rank_sp_m2

# ---------------------------------------------
# Plot species richness
# ---------------------------------------------

sp_rich_m2 <- data.frame(sp_rich_m2)
sp_rich_long <- sp_rich_m2 %>%
  mutate(gen=1:(ngen)) %>%
  pivot_longer(cols=colnames_long, names_to="rep",
               names_prefix="X", values_to="sp_rich_m2")
p <- ggplot(data=sp_rich_long, aes(x=gen, y=sp_rich_m2, col=rep)) +
  geom_line() +
  scale_colour_viridis_d()+
  xlab("Generations") + 
  ylab("Species richness")+
  theme(legend.position = "none", text = element_text(size = 20))+
  ylim(0, 100)
ggsave(p, filename=here("outputs", "m2", "species_richness_with_time.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# Mean final rank
sp_mean_rank <- apply(rank_sp_m2, 2, mean)

# minus because the max must have rank 1
rank_dist_E <- t(apply(-perf_mean, 1, rank, ties.method="min"))
# find the number of cells on which each species wins.
sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
sp_hab_freq <- as.table(sp_hab_freq)
names(sp_hab_freq) <- 1:nsp
# Plot
df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
  geom_point() +
  geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
  xlab("Species suitable habitat frequency") +
  ylab("Species mean rank (higher rank = lower abundance)") +
  theme(axis.title=element_text(size=16))
ggsave(p, filename=here("outputs", "m2", "mean_rank-habitat_freq.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# ---------------------------------------------------------------------------------
# pairwise Spearman correlation on the species ranks at the end of each simulation
# ---------------------------------------------------------------------------------

Spearman <- as.dist(round(cor(t(rank_sp_m2), method="spearman"),2))
save(Spearman, file = here::here("outputs", "m2", "Spearman_m2.RData"))

Df_rank_sp <- data.frame(rank_sp_m2)%>%
  mutate(rep = 1:nrep)%>%
  pivot_longer(cols=X1:X100, values_to="rank_sp_m2")
p <- ggplot(data=Df_rank_sp, aes(x=rep, y=rank_sp_m2)) +
  geom_line(aes(colour=name)) + 
  scale_colour_viridis_d()+
  xlab("Repetition") + 
  ylab("Species rank")+
  theme(legend.position = "none")
ggsave(p, filename=here("outputs", "m2", "species_rank_with_repetitions.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

sp_rich_final <- sp_rich_m2[ngen,]
save(sp_rich_final, file=here::here("outputs", "m2", "Species_richness_m2.RData"))

# ---------------------------------------------
# Shannon index and Shannon equitability index
# ---------------------------------------------
save(Shannon, file = here::here("outputs", "m2", "Shannon_m2.RData"))
Equitability <- Shannon/log(as.numeric(sp_rich_m2[ngen,]))
save(Equitability, file = here::here("outputs", "m2", "Equitability_m2.RData"))

# ---------------------------------------------
# Environmental filtering
# ---------------------------------------------

colnames_long <- paste0("X", 1:nrep)

env_filt <- data.frame(env_filt)
env_filt_long <- env_filt %>%
  mutate(gen=1:(ngen)) %>%
  pivot_longer(cols=colnames_long, names_to="rep",
               names_prefix="X", values_to="env_filt")
p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
  geom_line() +
  scale_colour_viridis_d() +
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Distance between the performance \n and the environment on each site")
ggsave(p, filename=here("outputs", "m2", "environmental_filtering.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Plot
png(file=here("outputs", "m2", "spatial_comp_env_sp.png"), 
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(raster(community_start_m2), main="Species - Start", zlim=c(0, nsp),
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
  mutate(gen=1:(ngen)) %>%
  pivot_longer(cols=colnames_long, names_to="rep",
               names_prefix="X", values_to="theta_comm")
p <- ggplot(data=theta_comm_long, aes(x=gen, y=theta_comm, col=rep)) +
  geom_line() +
  scale_colour_viridis_d()+
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean mortality rate in the community")
ggsave(p, filename=here("outputs", "m2", "mortality_rate_community.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(community)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation X1
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=X1)

# 3D voxel for each site
niche_width <- 0.25
n_niche <- 1/niche_width
x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
vario_env_all <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)

# Plot with correlation
png(file=here("outputs", "m2", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(vario_sp, main="Species - End")
plot(vario_env, main="Environment (one axis)")
plot(vario_env_all, main="Environment (three axes)")
plot(vario_env$v, vario_sp$v,
     main = "Regression",
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

load(here::here("outputs", "m2", "Species_richness_m2.RData"))
sp_rich_final_m2 <- sp_rich_final

sp_rich_final <- data.frame(sp_rich_m2 = c(as.numeric(sp_rich_final_m0), as.numeric(sp_rich_final_m2)),
                            Model = as.factor(rep(0:1, each=length(sp_rich_final_m0))))
p <- ggplot(sp_rich_final, aes(x= Model, y=sp_rich_m2, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the species richness \n at the end of each repetition \n with both models", y="Species richness")+
  theme(text = element_text(size = 20))

ggsave(p, filename=here("outputs", "Comparison", "Species_richness.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Shannon diversity and equitability index

load(here::here("outputs", "m0", "Shannon_m0.RData"))
Shannon_m0 <- Shannon

load(here::here("outputs", "m2", "Shannon_m2.RData"))
Shannon_m2 <- Shannon

Shannon <- data.frame(Shannon = c(Shannon_m0, Shannon_m2),
                      Model = as.factor(rep(0:1, each=length(Shannon_m0))))
p <- ggplot(Shannon, aes(x= Model, y=Shannon, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon diversity \n index at the end of each repetition \n with both models", y="Shannon diversity index")+
  theme(text = element_text(size = 20))

ggsave(p, filename=here("outputs", "Comparison", "Shannon.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

load(here::here("outputs", "m0", "Equitability_m0.RData"))
Equitability_m0 <- Equitability

load(here::here("outputs", "m2", "Equitability_m2.RData"))
Equitability_m2 <- Equitability

Equitability <- data.frame(Equitability = c(Equitability_m0, Equitability_m2),
                           Model = as.factor(rep(0:1, each=length(Equitability_m0))))

p <- ggplot(Equitability, aes(x= Model, y=Equitability, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon equitability \n index at the end of each repetition \n with both models", y="Shannon equitability index")+
  theme(text = element_text(size = 20))

ggsave(p, filename=here("outputs", "Comparison", "Equitability.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Spearman pairwise correlation of species ranks

load(here::here("outputs", "m0", "Spearman_m0.RData"))
Spearman_m0 <- Spearman

load(here::here("outputs", "m2", "Spearman_m2.RData"))
Spearman_m2 <- Spearman

Spearman <- data.frame(Spearman = c(as.vector(Spearman_m0), as.vector(Spearman_m2)),
                       Model = as.factor(rep(0:1, each=length(Spearman_m0))))

p <- ggplot(Spearman, aes(x= Model, y=Spearman, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the pairwise Spearman \n correlation of species ranks \n at the end of each repetition \n with both models", y="pairwise Spearman correlation")+
  theme(text = element_text(size = 20))
ggsave(p, filename=here("outputs", "Comparison", "Spearman_ranks.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)


png(file=here("outputs", "Comparison", "Community.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2), bty = "n")
plot(raster(community_start_m0), main="m0 - Species at start", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)), legend=FALSE, cex.main=1.5)
plot(raster(community_start_m2), main="m2 - Species at start", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)), legend=FALSE, cex.main=1.5)
plot(raster(community_end_m0), main="m0 - Species at end", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)), legend=FALSE, cex.main=1.5)
plot(raster(community_end_m2), main="m2 - Species at end", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)), legend=FALSE, cex.main=1.5)
dev.off()


