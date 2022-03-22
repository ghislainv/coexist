## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================
launch_model <- function(){

  source(file = "./call_libraries.R")
  
  source(file = here::here("Basic_parameters.R"))
  
  # Create output directories
  dir.create(here::here("outputs", model), recursive=TRUE)
  
  source(file=here::here("Math_functions.R"))
  
  source(file=here::here("Plot_functions.R"))
  
  source(file=here::here("Generate_environment.R"))
  
  source(file=here::here("Species_parameters.R"))
  
  # Function to identify the species with the highest performance
  source(file=here::here("Functions_perf_mort.R"))
  
  # =========================
  # Landscape and environment
  # =========================
  
  if(perf_know==TRUE){
    generate_environment(nsite_side=nsite_side, model=model)
    load(here::here("outputs", model, "sites.RData"))
    load(here::here("outputs", model, "env.RData"))
  } else {
    X1 <- as.vector(sites[,1])
  }
  # Plot the environment
  plot_environment(model=model, fig_width=fig_width, n_axis=n_axis, env=env)
  # Plot the habitat frequency
  plot_hab_freq(n_axis=n_axis, model=model, fig_width=fig_width, env=env)
  
  # =========================================
  # Species optima
  # =========================================
  
  if(perf_know==TRUE){
    niche_optimum <- NULL
    niche_optimum <- generate_species_optima(model=model, randomOptSp=randomOptSp, niche_width=niche_width, nsp=nsp, env=env)
  } else {
    beta_0 <- as.vector(lm_fit$coefficients[1:nsp])
    beta_1 <- as.vector(c(lm_fit$coefficients[nsp+1], lm_fit$coefficients[(nsp+3):(2*nsp+1)]))
    beta_2 <- as.vector(c(lm_fit$coefficients[nsp+2], lm_fit$coefficients[(2*nsp+2):(3*nsp)]))
    beta_0_mat <- matrix(rep(beta_0,each=nsite), ncol=nsp)
    beta_1_mat <- matrix(rep(beta_1,each=nsite), ncol=nsp)
    beta_2_mat <- matrix(rep(beta_2,each=nsite), ncol=nsp)
    X1_mat <- matrix(rep(X1, nsp), ncol=nsp)
  }
  
  # =========================================
  # Species performance
  # =========================================
  
  # Matrix of species performance on each site (distances)
  # Sites in rows, Species in columns
  if(perf_know==TRUE){
    dist_E_Sp <- dist_Site_Sp(as.matrix(sites), as.matrix(niche_optimum))
    dprim_E_Sp <- (dist_E_Sp-mean(dist_E_Sp))/sd(dist_E_Sp)
    perf_E_Sp <- -dprim_E_Sp
    save(perf_E_Sp, file=here::here("outputs", model, "perf_E_Sp.RData"))
  } else {
    perf_Sp_mean <- beta_0_mat+beta_1_mat*X1_mat+beta_2_mat*X1_mat^2
  }
  
  # =========================================
  # Species mortality probability
  # =========================================
  
  if(perf_know==TRUE){
    if(mortality_fixed==FALSE){
      # Probability of dying of each species on each site
      mortality_E_Sp <- inv_logit(logit(theta) + b * perf_E_Sp)
      # Mortality rate distribution
      plot_histogram_mortality(model=model, fig_width=fig_width, mortality_E_Sp=mortality_E_Sp)
      # Mortality probability function
      plot_function_mort_proba(model=model, fig_width=fig_width, perf_E_Sp=perf_E_Sp, mortality_E_Sp=mortality_E_Sp)
    }
    
    # Habitat frequency for each species
    # Identify the rank of the performance of species
    rank_dist_E <- t(apply(dist_E_Sp, 1, rank, ties.method="min"))
    # Identify the best performing species
    sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
    sp_hab_freq <- as.table(sp_hab_freq)
    names(sp_hab_freq) <- 1:nsp
    save(sp_hab_freq, file = here::here("outputs", model, "sp_hab_freq.RData"))
  }
  
  if(perf_know==FALSE){
    rank_dist_E <- t(apply(-perf_Sp_mean, 1, rank, ties.method="min"))
    sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
    sp_hab_freq <- as.table(sp_hab_freq)
    names(sp_hab_freq) <- 1:nsp
    save(sp_hab_freq, file = here::here("outputs", model=model, "sp_hab_freq.RData"))
    mortality_Sp_mean <- inv_logit(logit(theta) + b * (perf_Sp_mean))
  }
  
  plot_species_habitat_freq(model=model, fig_width=fig_width, sp_hab_freq=sp_hab_freq)
  
  # Species with no habitat
  sp_no_habitat <- as.vector(which(sp_hab_freq==0))
  nsp_no_habitat <- length(sp_no_habitat)
  
  print(glue::glue("There are {nsp_no_habitat} species with no suitable habitat"))
  
  # =========================================
  # Repetitions
  # =========================================
  
  # Species richness
  sp_rich <- matrix(NA, nrow=ngen, ncol=nrep)
  # Species rank at the end of the generations
  rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
  # Environmental filtering
  env_filt <- matrix(NA, nrow=ngen, ncol=nrep)
  # Mean mortality rate in the community
  theta_comm <- matrix(NA, nrow=ngen, ncol=nrep)
  #Distance between site and species optimum
  dist_site <- matrix(NA, nrow=ngen, ncol=nrep)
  #Shannon diversity index
  Shannon <- c(nrep)
  #Abundance matrices
  Abundances <- list()
  
  # Loop on repetitions
  for (r in 1:nrep) {
    
    abund <- matrix(NA, ncol=nsp, nrow=ngen)
    #used if disp_dep_abund==TRUE 
    abund_after_mortality <- NULL
    
    # -----------------------------------------
    # Initial conditions
    # -----------------------------------------
    
    if(start_full_landscape==TRUE){
      # Draw species at random in the landscape (one individual per site)
      sp <- sample(1:nsp, size=nsite, replace=TRUE)
      community <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    }
    if(start_one_ind_per_species==TRUE){
      # One individual per species distributed randomly over the grid
      sites_start <- sample(1:nsite, size=nsp, replace=FALSE)
      sp_start<- sample(1:nsp, size=nsp, replace=FALSE)
      community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      community_rast <- raster(community)
      community_rast[sites_start] <- sp_start
      community <- as.matrix(community_rast)
    }
    if(start_ten_ind_per_species==TRUE){
      # Ten individuals per species distributed randomly over the grid
      sites_start <- sample(1:nsite, size=nsp*10, replace=FALSE)
      sp_start<- rep(1:nsp, each=10)
      community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      community_rast <- raster(community)
      community_rast[sites_start] <- sp_start
      community <- as.matrix(community_rast)
    }
    community_start <- community
    
    # Plot the community at the start of the first repetition
    if (r==1) {
      plot_community_start(model=model, fig_width=fig_width, community=community_start, nsp=nsp)
    }
    
    if(perf_know==FALSE&&IV==TRUE){
      # Adding endogenous IV to compute mortality for the first generation
      epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
      perf_ind <- perf_Sp_mean + epsilon
    }
  
    # -----------------------------------------
    # Dynamics
    # -----------------------------------------
    
    # Simulating generation
    for (g in 1:ngen) {
      
      # Species richness
      sp_rich[g, r] <- length(unique(as.vector(community[community!=0])))
      abund[g, ] <- table(factor(as.vector(community), levels=1:nsp))
      
      # ******************
      # Mortality
      # ******************
      
      if(mortality_fixed==FALSE){
        
        if(perf_know==FALSE&&IV==TRUE){
          
          # Probability of dying of each species on each site
          mortality_ind <- inv_logit(logit(theta) + b * (perf_ind))
          
          # Mean mortality correction
          epsilon_mat <- matrix(epsilon, ncol=nsp)
          theta_var <- inv_logit(logit(theta) + b * epsilon_mat)
          diff <- mean(theta_var)-theta
          mortality_ind <- mortality_ind - diff
          # /!\ can modify mean!!
          mortality_ind[mortality_ind<0] <- 0
          
        } # end condition on IV
        
        # Mortality rate on each site
        # w0 for vacant sites
        theta_site <- rep(NA, nsite)
        w0 <- (as.vector(t(community))==0)
        theta_site[w0] <- 0
        
        if(perf_know==TRUE&&IV==FALSE){
          if(length(c(theta_site[!w0]))>1){
            theta_site[!w0] <- diag(mortality_E_Sp[!w0, as.vector(t(community))[!w0]])
          }
          if(length(c(theta_site[!w0]))==1){
            theta_site[!w0] <- mortality_E_Sp[!w0, as.vector(t(community))[!w0]]
          }
        }
        
        if(perf_know==FALSE&&IV==FALSE){
          if(length(c(theta_site[!w0]))>1){
            theta_site[!w0] <- diag(mortality_Sp_mean[!w0, as.vector(t(community))[!w0]])
          }
          if(length(c(theta_site[!w0]))==1){
            theta_site[!w0] <- mortality_Sp_mean[!w0, as.vector(t(community))[!w0]]
          }
        }
        
        if(perf_know==FALSE&&IV==TRUE){
          if(length(c(theta_site[!w0]))>1){
            theta_site[!w0] <- diag(mortality_ind[!w0, as.vector(t(community))[!w0]])
          }
          if(length(c(theta_site[!w0]))==1){
            theta_site[!w0] <- mortality_ind[!w0, as.vector(t(community))[!w0]]
          }
        }
        
        # Mortality events
        mort_ev <- rbinom(nsite, size=1, prob=theta_site)
        mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      
      }
      
      #/!\ to be tested
      
      if(mortality_fixed==TRUE){
        # No stocasticity: the 10 less performant individuals of the community die each generation
        
        perf_present <- rep(NA, nsite)
        w0 <- (as.vector(t(community))==0)
        perf_present[w0] <- NA
        
        if(perf_know==TRUE&&IV==FALSE){
          #performance of the present individuals
          #empty sites are not taken into account
          perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community))[!w0]])
        }
        
        if(perf_know==FALSE&&IV==FALSE){
          #performance of the present individuals
          perf_present[!w0] <- diag(perf_Sp_mean[!w0, as.vector(t(community))[!w0]])
        }
        
        if(perf_know==FALSE&&IV==TRUE){
          #performance of the present individuals
          perf_present[!w0] <- diag(perf_ind[!w0, as.vector(t(community))[!w0]])
        }
        
        #identify the 10 less performant present individuals and kill them
        mort_ev <- rep(0, length(perf_present))
        mort_ev[which(perf_present<=sort(perf_present)[10])] <- 1 #sort ignores NA
        mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      }
      
      # Number of deaths
      n_mort <- sum(mort_ev)
        
      if(n_mort!=0){
        
        # Update community
        community[mortality==1] <- 0
        
        # Plot once
        if (r==1 & g==1) {
          plot_mortality_events(model=model, fig_width=fig_width, community=community, nsp=nsp)
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
        
        if(perf_know==FALSE&&IV==TRUE){
          # New individual effects for potential new individuals (one per site and species)
          epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
          perf_ind_pot <- perf_Sp_mean + epsilon
          
        }
          
        if(disp_dep_abund==FALSE){
          
          if(perf_know==TRUE&&IV==FALSE){
            
            # Performance of species on vacant sites
            dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
            
            # Identify the species with the highest performance
            new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
            
            # Recruitment
            community_rast[sites_vacant] <- new_ind
            
          } # end condition on no IV
          
          if(perf_know==FALSE&&IV==FALSE){
            # Identify the present species with the highest performance on vacant sites (maximum in each line)
            sp_high_perf <- sp_present[apply(matrix(perf_Sp_mean[sites_vacant, sp_present], ncol=nsp_present), 1, which.max)]
            
            # Recruitment
            community_rast[sites_vacant] <- sp_high_perf
          } # end condition on IV
          
          if(perf_know==FALSE&&IV==TRUE){
            # Update performance matrix
            perf_ind[sites_vacant, sp_present] <- perf_ind_pot[sites_vacant, sp_present]
            
            # Identify the present species with the highest performance on vacant sites (maximum in each line)
            sp_high_perf <- sp_present[apply(matrix(perf_ind_pot[sites_vacant, sp_present], ncol=nsp_present), 1, which.max)]
            
            # Recruitment
            community_rast[sites_vacant] <- sp_high_perf
          } # end condition on IV
          
        } # end condition on abundance dependence
        
        if (disp_dep_abund==TRUE){
          
          abund_after_mortality <- as.data.frame(table(factor(as.vector(community), levels=1:nsp)))$Freq
          nb_seeds_sp <- round(abund_after_mortality*fecundity)
          nb_seeds_tot <- sum(nb_seeds_sp)
          
          sites_of_seeds <- sample(sites_vacant, nb_seeds_tot, replace=TRUE)
          
          # Performance of species on vacant sites
          seeds <- data.frame(Species=rep(1:nsp, times=nb_seeds_sp), Sites=sites_of_seeds)
          
          if (nrow(seeds)>1) {
            if(perf_know==TRUE&&IV==FALSE){
              seeds$Perf <- diag(perf_E_Sp[seeds$Sites, seeds$Species])
            }
            if(perf_know==FALSE&&IV==FALSE){
              seeds$Perf <- diag(perf_Sp_mean[seeds$Sites, seeds$Species])
            }
            if(perf_know==FALSE&&IV==TRUE){
              seeds$Perf <- diag(perf_ind_pot[seeds$Sites, seeds$Species])
            }
          } else {
            if(perf_know==TRUE&&IV==FALSE){
              seeds$Perf <- perf_E_Sp[seeds$Sites, seeds$Species]
            }
            if(perf_know==FALSE&&IV==FALSE){
              seeds$Perf <- perf_Sp_mean[seeds$Sites, seeds$Species]
            }
            if(perf_know==FALSE&&IV==TRUE){
              seeds$Perf <- perf_ind_pot[seeds$Sites, seeds$Species]
            }
          } # end else number of seeds
          
          new_ind <- plyr::ddply(seeds, c("Sites"), f_sp)
          colnames(new_ind) <- c("Sites", "Species")
          
          # Recruitment
          community_rast[new_ind$Sites] <- new_ind$Species
          
          if(perf_know==FALSE&&IV==TRUE){
            # /!\ remove loop
            for(k in 1:nrow(new_ind)){
              
            perf_ind[new_ind$Sites[k], new_ind$Species[k]] <- perf_ind_pot[new_ind$Sites[k], new_ind$Species[k]]
            
            }
            
          }
          
        } # end condition on option on abundance dependence of dispersion 
        
        
      }  # end condition on n_mort (contains all cases)
      
      # update community
      community <- as.matrix(community_rast)
      
      # *********************
      # Diversity
      # *********************
      
      # Environmental filtering
      if(IV==FALSE&&perf_know==TRUE){
        #Some sites have a 0 but then they are not taken into account since there is no distance associated with "species 0".
        dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
      }
      if(perf_know==FALSE){
        
        # /!\ to check
        
        # Environmental filtering for the partial knowledge model
        #distance between the observed optimum of the species and the environmental conditions
        E_seq <- seq(min(X1), max(X1), length.out=100)
        E_seq_mat <- matrix(rep(E_seq, nsp), ncol=nsp)
        
        beta_0_mat_E_seq <- matrix(rep(beta_0,each=length(E_seq)), ncol=nsp)
        beta_1_mat_E_seq <- matrix(rep(beta_1,each=length(E_seq)), ncol=nsp)
        beta_2_mat_E_seq <- matrix(rep(beta_2,each=length(E_seq)), ncol=nsp)
        
        Mat_perf_inferred <- beta_0_mat_E_seq+beta_1_mat_E_seq*E_seq_mat+beta_2_mat_E_seq*E_seq_mat^2
        Mat_perf_inferred <- as.data.frame(cbind(E_seq, Mat_perf_inferred))
        
        colnames(Mat_perf_inferred) <- c("E", paste0("Sp",1:nsp))
        
        Optimum_Sp_inferred <- c()
        for (k in 1:nsp) {
          Optimum_Sp_inferred[k] <- Mat_perf_inferred[which.max(Mat_perf_inferred[,k+1]),]$E
        }
        
        Optimum_Sp_inferred_community <- Optimum_Sp_inferred[community]
        
        dist_site <- sqrt((Optimum_Sp_inferred_community-X1[which(as.vector(community)!=0)])^2)
        
      }
      
      env_filt[g, r] <- mean(dist_site)
  
      # Mean mortality rate in the community
      #/!\ do we want this mortality rate to take empty sites into account?
      if(mortality_fixed==FALSE){
        if (IV==FALSE&&perf_know==TRUE){
          theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
        }
        if (IV==TRUE&&perf_know==FALSE){
          theta_site <- diag(mortality_ind[, as.vector(t(community))])
        }
        theta_comm[g, r] <- mean(theta_site)
      }
      
    } # End ngen
    
    community_end <- community
    
    # Plot final community once
    if (r==1) {
      plot_community_end(model=model, fig_width=fig_width, community=community_end, nsp=nsp)
      save(community_end, file=here::here("outputs", model, glue::glue("community_end_{model}.RData")))
    }
    
    # Species rank
    rank_sp[r, ] <- rank(-abund[ngen, ], ties.method="min")
    
    df_shannon <- data.frame(Species = 1:nsp,
                             Abundance = abund[ngen, ])%>%
      mutate(Proportion = Abundance / sum(Abundance))%>%
      filter(Abundance > 0)%>%
      mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
    
    Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
    
    #To keep the abundance matrixes in order to infer alpha matrix
    Abundances[[r]] <- abund
    
    
  } # End nrep
  
  sp_rich <- data.frame(sp_rich)
  env_filt <- data.frame(env_filt)
  theta_comm <- data.frame(theta_comm)
  
  save(Abundances, file = here::here("outputs", model, glue::glue("Abundances_{model}.RData")))
  save(sp_rich, file=here::here("outputs", model, glue::glue("sp_rich_{model}.RData")))
  save(rank_sp, file=here::here("outputs", model, glue::glue("rank_sp_{model}.RData")))
  
  # =========================
  # Diversity analysis
  # =========================
  
  # ---------------------------------------------
  # Plot species abundance
  # ---------------------------------------------

  plot_abundance_species(Abundances, model, fig_width)
  
  # ---------------------------------------------
  # Plot species richness
  # ---------------------------------------------
  
  plot_species_richness(nrep=nrep, sp_rich=sp_rich, model=model, fig_width=fig_width)
  
  sp_rich_final <- sp_rich[ngen,]
  save(sp_rich_final, file=here::here("outputs", model, glue::glue("Species_richness_{model}.RData")))
  
  # ---------------------------------------------
  # Shannon index and Shannon equitability index
  # ---------------------------------------------
  save(Shannon, file = here::here("outputs", model, glue::glue("Shannon_{model}.RData")))
  Equitability <- Shannon/log(as.numeric(sp_rich[ngen,]))
  save(Equitability, file = here::here("outputs", model, glue::glue("Equitability_{model}.RData")))
  
  # ---------------------------------------------------------------------------------
  # pairwise Spearman correlation on the species ranks at the end of each simulation
  # ---------------------------------------------------------------------------------
  
  Spearman <- as.dist(round(cor(t(rank_sp), method="spearman"),2))
  save(Spearman, file = here::here("outputs", model, glue::glue("Spearman_{model}.RData")))
  
  # ---------------------------------------------
  # Link between final rank and habitat frequency
  # ---------------------------------------------
  
  # Mean final rank
  sp_mean_rank <- apply(rank_sp, 2, mean)
  # Plot
  plot_mean_rank_hab_freq(sp_mean_rank=sp_mean_rank, sp_hab_freq=sp_hab_freq, model=model, fig_width=fig_width)
  
  # ---------------------------------------------
  # Environmental filtering
  # ---------------------------------------------
  
  plot_env_filt(nrep, env_filt, model, fig_width)
  plot_env_species(model, fig_width, community_start, community_end, env)
  
  # ---------------------------------------------
  # Theta community
  # ---------------------------------------------
  
  plot_theta_community(theta_comm, ngen, model, fig_width)
  
  # ----------------------------------
  # Spatial autocorrelation of species
  # ----------------------------------
  
  plot_spatial_autocorr(community_end, sites, niche_width, model, fig_width)
  
} # End function

load(here::here("outputs", model, glue::glue("Abundances_{model}.RData")))
load(here::here("outputs", model, glue::glue("sp_rich_{model}.RData")))
load(here::here("outputs", model, glue::glue("rank_sp_{model}.RData")))
load(here::here("outputs", model, glue::glue("community_end_{model}.RData")))

# ========================================
# Infer observed intraspecific variability
# ========================================

infer_IV <- function(){
  load(here::here("outputs", model, "perf_E_Sp.RData"))
  load(here::here("outputs", model, "env.RData"))
  nsp<-ncol(perf_E_Sp)

  # Data-set
  df <- data.frame(perf_E_Sp)
  names(df) <- c(sprintf("Sp_%03d", 1:(nsp-1)), sprintf("Sp_%d", nsp))
  df_perf <- tibble(df) %>%
    mutate(Env_1=values(raster(env[[1]])), Env_2=values(raster(env[[2]])), Env_3=values(raster(env[[3]]))) %>%
    mutate(Env_1_sq=Env_1^2, Env_2_sq=Env_2^2, Env_3_sq=Env_3^2) %>%
    pivot_longer(cols=c(Sp_001:glue("Sp_0{nsp-1}"), glue("Sp_{nsp}")), names_to="Species", values_to="Perf")
  
  # Observed niche
  plot_random_species_niche(seed, df_perf, model, fig_width)
  
  #Check that the model well fits the data if all environmental variables are included
  lm_all_env <- lm(Perf~Species+Species*Env_1+Species*Env_1_sq+Species*Env_2+Species*Env_2_sq+Species*Env_3+Species*Env_3_sq, data=df_perf)
  
  print(glue::glue("the r-squared with all environmental variables is {summary(lm_all_env)$adj.r.squared}"))
    
  # Observed intraspecific variability
  lm_fit <- lm(Perf~Species+Species*Env_1+Species*Env_1_sq, data=df_perf)
  save(lm_fit, file = here::here("outputs", model,"lm_fit.RData"))
  
  print(glue::glue("the r-squared with the first environmental variable is {summary(lm_fit)$adj.r.squared}"))
  
  #Histogram of the residuals to check normality
  hist(summary(lm_fit)$residuals)
  
  V_intra <- df_perf %>%
    mutate(res=lm_fit$residuals) %>%
    group_by(Species) %>%
    summarise(V=var(res))
  save(V_intra, file = here::here("outputs", model, "V_intra.RData"))
  
  plot_IV(V_intra, model, fig_width)
  
  return(V_intra)
  
} # End function inferring IV


plot_inferred_perf_environment(env, lm_fit, nsp, model, fig_width)

plot_comp_IV_niche_width

plot_comp_IV_hab_freq


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
  scale_color_manual(values=c("red", "blue", "green", "orange"), "Presence of the species", labels = c("Disappeared in both models", "Disappeared in second model", "Maintained in both models"))+
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