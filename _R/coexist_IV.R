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
png(file=here("outputs", "m1", "environment.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(1,1))
plot(raster(env[[1]]), main="Environment var1", col=topo.colors(255))
dev.off()

# Habitat frequency
png(file=here("outputs", "m1", "hab_freq.png"),
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
nrep <- 10
# Number of generations
ngen <- 500

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Species rank at the end of the generations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
env_filt <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen+1, ncol=nrep)
#Shannon index of the community at the end of each repetition
Shannon <- c(nrep)
#list of the matrices of the species abundance through the generations for each repetition
Abundances_m1<-list()

# Matrix of mean species performance on each site
# Sites in rows, Species in columns
perf_mean <- matrix(nrow = nsite, ncol = nsp)
for (k in 1:nsp) {
  perf_mean[,k] <- beta_0[k]+beta_1[k]*X1+beta_2[k]*(X1^2)
}
rank_dist_E <- t(apply(-perf_mean, 1, rank, ties.method="min"))
sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
sp_hab_freq <- as.table(sp_hab_freq)
names(sp_hab_freq) <- 1:nsp
png(file=here("outputs", "m1", "species_habitat_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Species with no habitat
sp_no_habitat <- as.vector(which(sp_hab_freq==0))
nsp_no_habitat <- length(sp_no_habitat)

# Loop on repetitions
for (r in 1:nrep) {
  
  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Matrix of species performance on each site
  # Sites in rows, Species in columns
  perf_ind_init <- as.data.frame(matrix(ncol=nsp, nrow=nsite))
  
  for(sp in 1:nsp){
    perf_ind_init[,sp] <- beta_0[sp]+beta_1[sp]*X1+beta_2[sp]*(X1^2)+rnorm(1, mean=0, sd=sqrt(V_intra$V[sp]))
  }
  
  # Identify the species with the highest performance (in each line find the max)
  sp_high_perf <- apply(perf_ind_init, 1, which.max)
  
  # Probability of dying of each species on each site
  # b = strength of unsuitability
  b <- -0.5
  mortality_E_Sp <- inv_logit(logit(0.1) + b * as.matrix(perf_ind_init))
  
  # Mortality rate distribution
  png(file=here("outputs", "m1", "hist_mortality_init.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  hist(mortality_E_Sp)
  dev.off()
  
  # Draw species at random in the landscape (one individual per site)
  sp <- sample(1:nsp, size=nsite, replace=TRUE)
  #hist(sp)
  community_start_m1 <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m1", "community_start.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    plot(raster(community_start_m1), main="Species - Start", zlim=c(0, nsp),
         col=c("black", viridis(nsp)))
    dev.off()
  }
  community <- community_start_m1
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # Environmental filtering
  # distance between site environment and performance
  dist_site <- sqrt((diag(as.matrix(perf_ind_init)[, as.vector(t(community))])-X1)^2)
  env_filt[1, r] <- mean(dist_site, na.rm = TRUE)
  
  # Mean mortality rate
  theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
  theta_comm[1, r] <- mean(theta_site, na.rm = TRUE)
  
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
    
    # Performance of individuals of species present in the community
    perf_ind_pres_vacant <- as.data.frame(matrix(ncol=nsp_present, nrow=length(sites_vacant)))
    perf_ind <- as.data.frame(matrix(ncol=nsp, nrow=nsite))
    
    for(sp in sp_present){
      perf_ind[,sp] <- beta_0[sp]+beta_1[sp]*X1+beta_2[sp]*(X1^2)+rnorm(1, mean=0, sd=sqrt(V_intra$V[sp]))
    }
    
    #Species that are absent have no measured performance
    perf_ind[, which(!((1:nsp) %in% sp_present))] <- NA
    
    perf_ind_pres_vacant <- perf_ind[sites_vacant, sp_present]
    rownames(perf_ind_pres_vacant) <- 1:nrow(perf_ind_pres_vacant)

    # Identify the species with the highest performance
    new_ind <- sp_present[apply(perf_ind_pres_vacant, 1, which.max)]
    
    # Recruitment
    community_rast[sites_vacant] <- new_ind
    community <- as.matrix(community_rast)
    
    # Update of mortality each iteration
    mortality_E_Sp <- inv_logit(logit(0.1) + b * as.matrix(perf_ind))
    
    # *********************
    # Diversity
    # *********************
    
    # Species richness
    sp_rich[g+1, r] <- length(unique(as.vector(community)))
    abund[g+1, ] <- table(factor(as.vector(community), levels=1:nsp))
    
    # Environmental filtering
    dist_site <- sqrt((diag(as.matrix(perf_ind)[, as.vector(t(community))])-X1)^2)
    env_filt[g+1, r] <- mean(dist_site, na.rm = TRUE)
    
    # Mean mortality rate in the community
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    theta_comm[g+1, r] <- mean(theta_site, na.rm = TRUE)
    
  } # End ngen
  
  #CGT 28/06/2021 to compare the two models
  community_end_m1 <- community
  
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
  
  df_shannon <- data.frame(Species = 1:nsp,
                           Abundance = abund[ngen+1, ])%>%
    mutate(Proportion = Abundance / sum(Abundance))%>%
    filter(Abundance > 0)%>%
    mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
  
  Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
  
  #To keep the abundance matrices in order to infer alpha matrix
  Abundances_m1[[r]] <- abund
  
} # End nrep

save(Abundances_m1, file = here::here("outputs", "m1", "Abundances_m1.RData"))

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

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# Mean final rank
sp_mean_rank <- apply(rank_sp, 2, mean)

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
ggsave(p, filename=here("outputs", "m1", "mean_rank-habitat_freq.png"),
        width=fig_width, height=fig_width, units="cm", dpi=300)

# ---------------------------------------------------------------------------------
# pairwise Spearman correlation on the species ranks at the end of each simulation
# ---------------------------------------------------------------------------------

Spearman <- as.dist(round(cor(t(rank_sp), method="spearman"),2))
save(Spearman, file = here::here("outputs", "m1", "Spearman_m1.RData"))

Df_rank_sp <- data.frame(rank_sp)%>%
  mutate(rep = 1:nrep)%>%
  pivot_longer(cols=X1:X64, values_to="rank_sp")
p <- ggplot(data=Df_rank_sp, aes(x=rep, y=rank_sp)) +
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
# Environmental filtering
# ---------------------------------------------

 env_filt <- data.frame(env_filt)
 env_filt_long <- env_filt %>%
   mutate(gen=1:(ngen+1)) %>%
   pivot_longer(cols=X1:X10, names_to="rep",
                names_prefix="X", values_to="env_filt")
 p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
   geom_line() +
   scale_colour_viridis_d() +
   labs(title="Environmental filtering") +
   xlab("Generations") + 
   ylab("Distance between the performance \n and the environment on each site")
 ggsave(p, filename=here("outputs", "m1", "environmental_filtering.png"),
        width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Plot
png(file=here("outputs", "m1", "spatial_comp_env_sp.png"), 
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2))
plot(raster(community_start_m1), main="Species - Start", zlim=c(0, nsp),
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
png(file=here("outputs", "m1", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2))
plot(vario_sp, main="Species - End")
plot(vario_env, main="Environment (one axis)")
plot(vario_env_all, main="Environment (three axes)")
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
p <- ggplot(sp_rich_final, aes(x= Model, y=Sp_rich, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the species richness \n at the end of each repetition with both models", y="Species richness")

ggsave(p, filename=here("outputs", "Comparison", "Species_richness.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Shannon diversity and equitability index

load(here::here("outputs", "m0", "Shannon_m0.RData"))
Shannon_m0 <- Shannon

load(here::here("outputs", "m1", "Shannon_m1.RData"))
Shannon_m1 <- Shannon

Shannon <- data.frame(Shannon = c(Shannon_m0, Shannon_m1),
                            Model = as.factor(rep(1:2, each=length(Shannon_m0))))
p <- ggplot(Shannon, aes(x= Model, y=Shannon, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon diversity index \n at the end of each repetition with both models", y="Shannon diversity index")

ggsave(p, filename=here("outputs", "Comparison", "Shannon.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

load(here::here("outputs", "m0", "Equitability_m0.RData"))
Equitability_m0 <- Equitability

load(here::here("outputs", "m1", "Equitability_m1.RData"))
Equitability_m1 <- Equitability

Equitability <- data.frame(Equitability = c(Equitability_m0, Equitability_m1),
                      Model = as.factor(rep(1:2, each=length(Equitability_m0))))

p <- ggplot(Equitability, aes(x= Model, y=Equitability, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the Shannon equitability index \n at the end of each repetition with both models", y="Shannon equitability index")

ggsave(p, filename=here("outputs", "Comparison", "Equitability.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# Spearman pairwise correlation of species ranks

load(here::here("outputs", "m0", "Spearman_m0.RData"))
Spearman_m0 <- Spearman

load(here::here("outputs", "m1", "Spearman_m1.RData"))
Spearman_m1 <- Spearman

Spearman <- data.frame(Spearman = c(as.vector(Spearman_m0), as.vector(Spearman_m1)),
                           Model = as.factor(rep(1:2, each=length(Spearman_m0))))

p <- ggplot(Spearman, aes(x= Model, y=Spearman, colour = Model))+
  geom_boxplot()+
  scale_colour_viridis_d(option="plasma")+
  labs(title = "Boxplot of the pairwise Spearman correlation of species ranks \n at the end of each repetition with both models", y="Shannon equitability index")
ggsave(p, filename=here("outputs", "Comparison", "Spearman_ranks.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)


png(file=here("outputs", "Comparison", "Community.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2))
plot(raster(community_start_m0), main="m0 - Species at start", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)))
plot(raster(community_start_m1), main="m1 - Species at start", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)))
plot(raster(community_end_m0), main="m0 - Species at end", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)))
plot(raster(community_end_m1), main="m1 - Species at end", zlim=c(0, nsp),
     col=c("black",  viridis(nsp)))
dev.off()


perf_ind <- as.data.frame(matrix(ncol=nsp, nrow=nsite))
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
load(here::here("outputs", "m0", "sp_hab_freq.RData"))
sp_hab_freq <- as.data.frame(sp_hab_freq)

#Step 1 : no repetitions
Df_rank_hab <- data.frame(Rank = rank_sp[1,], Freq = sp_hab_freq$Freq)
ggplot(data = Df_rank_hab, aes(x = Freq, y = Rank))+
  geom_point()+
  labs(title="Relationship between the suitable habitat frequency in m0  \n and the species rank at the end of one simulation in m1",
       x = "Suitable habitat frequency in m0",
       y = "Species rank at the end of one m1 simulation")

#Step 2 : repetitions for rank
Df_rank_hab$Mean_rank <- sp_mean_rank
Df_rank_hab$sd_rank <- apply(rank_sp, 2, sd)
ggplot(data = Df_rank_hab, aes(x = Freq, y = Mean_rank))+
  geom_point()+
  geom_point()+
  geom_errorbar(aes(ymin=Mean_rank-sd_rank, ymax=Mean_rank+sd_rank), width=.2,
                position=position_dodge(0.05))+
  labs(title="Relationship between the suitable habitat frequency in m0  \n and the mean and standard deviation of the species rank at the end of ten simulations in m1",
       x = "Suitable habitat frequency in m0",
       y = "Mean rank and standard deviations from 10 m1 simulations")

#Step 3 : repetitions for habitat frequency
Hab_freq <- matrix(ncol=nsp, nrow=10)

for (r in 1:10) {
  for(sp in 1:nsp){
    perf_ind[,sp] <- beta_0[sp]+beta_1[sp]*X1+beta_2[sp]*(X1^2)+rnorm(1, mean=0, sd=sqrt(V_intra$V[sp]))
  }
  rank_dist_E_m1 <- t(apply(-perf_ind, 1, rank, ties.method="min"))
  # find the number of cells on which each species wins in model m0
  sp_hab_freq_m1 <- apply(rank_dist_E_m1, 2, function(x){sum(x==1)})
  sp_hab_freq_m1 <- as.table(sp_hab_freq_m1)      
  names(sp_hab_freq_m1) <- 1:nsp
  sp_hab_freq_m1 <- as.data.frame(sp_hab_freq_m1)
  Hab_freq[r,] <- sp_hab_freq_m1$Freq
}

Df_rank_hab$Mean_freq <- apply(Hab_freq, 2, mean)
Df_rank_hab$sd_freq <- apply(Hab_freq, 2, sd)
ggplot(data = Df_rank_hab, aes(x = Freq, y = Mean_rank))+
  geom_point()+
  geom_point()+
  geom_errorbar(aes(ymin=Mean_rank-sd_rank, ymax=Mean_rank+sd_rank), width=.2,
                position=position_dodge(0.05))+
  geom_errorbarh(aes(xmin=Freq-sd_freq, xmax=Freq+sd_freq))+
  labs(title="Relationship between the suitable habitat frequency in m0 \n with the standard deviation from m1 over 10 repetitions,  \n and the mean and standard deviation of the species rank at the end of ten simulations in m1",
       x = "Suitable habitat frequency in m0 and standard deviation from 10 m1 simulations",
       y = "Mean rank and standard deviations from 10 m1 simulations")

