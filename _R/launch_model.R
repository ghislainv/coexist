launch_model <- function(){
  
  source(file = here::here("_R", "call_libraries.R"))
  
  # Create output directories
  dir.create(here::here("outputs", model), recursive=TRUE)
  
  source(file=here::here("_R", "Math_functions.R"))
  
  source(file=here::here("_R", "Plot_functions.R"))
  
  source(file=here::here("_R", "Generate_environment.R"))
  
  source(file=here::here("_R", "Species_parameters.R"))
  
  # =========================
  # Landscape and environment
  # =========================
  print("Generating or loading the environment...")
  
  if(perf_know==TRUE){
    generate_environment(nsite_side=nsite_side, model=model)
    load(here::here("outputs", model, "sites.RData"))
    load(here::here("outputs", model, "env.RData"))
    load(here::here("outputs", model, "env_entrelac.RData"))
  } else {
    Obs_env <- sites[,1:n_observed_axes]
  }
  
  # =========================================
  # Species optima
  # =========================================
  print("Generating or loading species parameters...")
  
  if(perf_know==TRUE){
    niche_optimum <- NULL
    niche_optimum <- generate_species_optima(model=model, randomOptSp=randomOptSp, niche_width=niche_width, nsp=nsp, env=env, seed=seed)
  } else {
    
    Inferred_species_parameters_mat <-list()
    
    for(k in 1:ncol(Inferred_species_parameters)){
      Inferred_species_parameters_mat[[k]] <- matrix(rep(Inferred_species_parameters[,k],each=nsite), ncol=nsp)
    }
    
    Obs_env_mat <- list()
    
    if(is.null(ncol(Obs_env))==FALSE){
      
      for(k in 1:ncol(Obs_env)){
        Obs_env_mat[[k]] <- matrix(rep(Obs_env[,k], nsp), ncol=nsp)
        Obs_env_mat[[k+ncol(Obs_env)]] <- matrix(rep((Obs_env[,k])^2, nsp), ncol=nsp)
      }
    }
    
    else{
      Obs_env_mat[[1]] <- matrix(rep(Obs_env, nsp), ncol=nsp)
      Obs_env_mat[[2]] <- matrix(rep(Obs_env^2, nsp), ncol=nsp)
    }

  }
  
  # =========================================
  # Species performance
  # =========================================
  print("Computing species performance on each site...")
  
  # Matrix of species performance on each site (distances)
  # Sites in rows, Species in columns
  if(perf_know==TRUE){
    dist_E_Sp <- dist_Site_Sp(as.matrix(sites), as.matrix(niche_optimum))
    dprim_E_Sp <- (dist_E_Sp-mean(dist_E_Sp))/sd(dist_E_Sp)
    perf_E_Sp <- -dprim_E_Sp
    save(perf_E_Sp, file=here::here("outputs", model, "perf_E_Sp.RData"))
  } else {
    
    perf_Sp_mean <- Inferred_species_parameters_mat[[1]]
    
    for(k in 1:length(Obs_env_mat)){
      perf_Sp_mean <- perf_Sp_mean + Inferred_species_parameters_mat[[k+1]]*Obs_env_mat[[k]]
    }
    save(perf_Sp_mean, file = here::here("outputs", model, "perf_Sp_mean.RData"))
  }
  
  # =========================================
  # Species mortality probability
  # =========================================
  print("Computing species ranks...")
  
  if(perf_know==TRUE){
    if(mortality_stocha==TRUE){
      # Probability of dying of each species on each site
      if(mortality_stocha_basal==TRUE){
        mortality_E_Sp <- inv_logit(logit(theta) + b * perf_E_Sp)
      }
      if(mortality_stocha_basal==FALSE){
        mortality_E_Sp <- inv_logit(b*perf_E_Sp)
      }
      # Mortality rate distribution
      plot_histogram_mortality(model=model, fig_width=fig_width, mortality_E_Sp=mortality_E_Sp)
      # Mortality probability function
      plot_function_mort_proba(model=model, fig_width=fig_width, perf_E_Sp=perf_E_Sp, mortality_E_Sp=mortality_E_Sp)
    }#end condition on stochastic mortality
    
    # Habitat frequency for each species
    # Identify the rank of the performance of species
    rank_dist_E <- t(apply(dist_E_Sp, 1, rank, ties.method="max"))
    # Identify the best performing species
    sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
    sp_hab_freq <- as.table(sp_hab_freq)
    names(sp_hab_freq) <- 1:nsp
    save(sp_hab_freq, file = here::here("outputs", model, "sp_hab_freq.RData"))
  }#end condition on perfect knowledge
  
  if(perf_know==FALSE){
    rank_dist_E <- t(apply(-perf_Sp_mean, 1, rank, ties.method="max"))
    sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
    sp_hab_freq <- as.table(sp_hab_freq)
    names(sp_hab_freq) <- 1:nsp
    save(sp_hab_freq, file = here::here("outputs", model=model, "sp_hab_freq.RData"))
    if(mortality_stocha==TRUE){
      if(mortality_stocha_basal==TRUE){
        mortality_Sp_mean <- inv_logit(logit(theta) + b * (perf_Sp_mean))
      }
      if(mortality_stocha_basal==FALSE){
        mortality_Sp_mean <- inv_logit(b*perf_Sp_mean)
      }
    }#end condition on stochastic mortality
  }#end condition on partial knowledge
  
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
  #Community at the beginning and the end of the simulation
  community_end <- list()
  community_start <- list()
  
  # Loop on repetitions
  #print("Entering the repetition loop")
  
  for (r in 1:nrep) {
    
    print(paste0("Repetition ", r, "/", nrep))
    
    abund <- matrix(NA, ncol=nsp, nrow=ngen)
    #used if nb_seeds_dep_abund==TRUE 
    abund_after_mortality <- NULL
    
    # -----------------------------------------
    # Initial conditions
    # -----------------------------------------
    print("Initialising the landscape...")
    
    if(start_full_landscape==TRUE){
      # Draw species at random in the landscape (one individual per site)
      set.seed(Seeds[r])
      sp <- sample(1:nsp, size=nsite, replace=TRUE)
      community <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    }
    if(start_one_ind_per_species==TRUE){
      # One individual per species distributed randomly over the grid
      set.seed(Seeds[r])
      sites_start <- sample(1:nsite, size=nsp, replace=FALSE)
      set.seed(Seeds[r])
      sp_start<- sample(1:nsp, size=nsp, replace=FALSE)
      community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      community_rast <- raster::raster(community)
      community_rast[sites_start] <- sp_start
      community <- as.matrix(community_rast)
    }
    if(start_ten_ind_per_species==TRUE){
      # Ten individuals per species distributed randomly over the grid
      set.seed(Seeds[r])
      sites_start <- sample(1:nsite, size=nsp*10, replace=FALSE)
      sp_start<- rep(1:nsp, each=10)
      community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      community_rast <- raster::raster(community)
      community_rast[sites_start] <- sp_start
      community <- as.matrix(community_rast)
    }
    community_start[[r]] <- community
    
    # Plot the community at the start of the first repetition
    if (r==1) {
      plot_community_start(model=model, fig_width=fig_width, community=community_start[[r]], nsp=nsp)
    }
    
    if(perf_know==FALSE&&IV==TRUE){
      # Adding endogenous IV to compute mortality for the first generation
      epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
      perf_ind <- perf_Sp_mean + epsilon
    }
    
    # -----------------------------------------
    # Dynamics
    # -----------------------------------------
    
    #print("Entering the loop on generations")
    
    # Simulating generation
    for (g in 1:ngen) {
      
      print(paste0("Generation ", g, "/", ngen, " of repetition ", r, "/", nrep))
      
      # Species richness
      sp_rich[g, r] <- length(unique(as.vector(community[community!=0])))
      abund[g, ] <- table(factor(as.vector(community), levels=1:nsp))
      
      # ******************
      # Mortality
      # ******************
      
      #print("Computing mortality...")
      
      #Total abundance in the community, used for mortality rate as a proportion of abundance
      abund_tot <- length(community[community!=0])
      
      if(mortality_stocha==TRUE){
        
        if(perf_know==FALSE&&IV==TRUE){
          
          if(mortality_stocha_basal==TRUE){
            # Probability of dying of each species on each site
            mortality_ind <- inv_logit(logit(theta) + b * (perf_ind))
            
            # Mean mortality correction
            epsilon_mat <- matrix(epsilon, ncol=nsp)
            theta_var <- inv_logit(logit(theta) + b * epsilon_mat)
            diff <- mean(theta_var)-theta
            mortality_ind <- mortality_ind - diff
            # /!\ can modify mean!
            mortality_ind[mortality_ind<0] <- 0
          }#end condition on basal mortality
          
          if(mortality_stocha_basal==FALSE){
            mortality_ind <- inv_logit(b*perf_ind)
          }
          
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
        if(mortality_stocha_basal==TRUE){
          mort_ev <- rbinom(nsite, size=1, prob=theta_site)
        }
        if(mortality_stocha_basal==FALSE){
          mort_ev <- rmultinom(1, size=round(0.01*abund_tot), prob=theta_site)
        }
        mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
        
      }#end condition on mortality not stochastic
      
      if(mortality_stocha==FALSE){
        # No stochasticity: the 10 less performing individuals of the community OR a fixed proportion of the total abundance die each generation
        
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
        
        mort_ev <- rep(0, length(perf_present))
        #identify the 10 less performant present individuals and kill them
        if(mortality_fixed==TRUE){
          mort_ev[which(perf_present<=sort(perf_present)[10])] <- 1 #sort ignores NA
        }
        #identify the 0.01*total abudance less performant present individuals and kill them
        if(mortality_proportion==TRUE){
          mort_ev[which(perf_present<=sort(perf_present)[round(abund_tot*0.01)])] <- 1 #sort ignores NA
        }
        mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      }
      
      # Number of deaths
      n_mort <- sum(mort_ev)
      
      if(n_mort!=0){
        
        # Plot once
        if (r==1 & g==1) {
          plot_mortality_events(model, fig_width, community, mortality, nsp)
        }
        
        # Update community
        community[mortality==1] <- 0
        
        
        # *********************
        # Fecundity/Recruitment
        # *********************
        
        #print("Computing recruitment...")
        
        # Species present in the community
        sp_present <- sort(unique(community[community!=0]))
        nsp_present <- length(sp_present)
        
        # Vacant sites
        community_rast <- raster::raster(community)
        sites_vacant <- which(raster::values(community_rast)==0)
        nsite_vacant <- length(sites_vacant)
        
        if(perf_know==FALSE&&IV==TRUE){
          # New individual effects for potential new individuals (one per site and species)
          epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
          perf_ind_pot <- perf_Sp_mean + epsilon
        }
        
        if(nb_seeds_dep_abund==FALSE){
          
          if(perf_know==TRUE&&IV==FALSE){
            
            # Performance of species on vacant sites
            dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
            
            # Identify the present species with the highest performance on vacant sites
            new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
            
            # Recruitment
            community_rast[sites_vacant] <- new_ind
            
          } # end condition perfect knowledge and no IV
          
          if(perf_know==FALSE&&IV==FALSE){
            
            # Identify the present species with the highest performance on vacant sites
            #sp_high_perf <- sp_present[apply(matrix(perf_Sp_mean[sites_vacant, sp_present], ncol=nsp_present), 1, which.max)]
            sp_high_perf <- apply(-matrix(perf_Sp_mean[sites_vacant, ], ncol=nsp), 1, high_perf_sp, sp_pres=sp_present)
            
            # Recruitment
            community_rast[sites_vacant] <- sp_high_perf
            
          } # end condition partial knowledge and no IV
          
          if(perf_know==FALSE&&IV==TRUE){
            
            # Identify the present species with the highest performance on vacant sites (maximum in each line)
            #sp_high_perf <- sp_present[apply(matrix(perf_ind_pot[sites_vacant, sp_present], ncol=nsp_present), 1, which.max)]
            sp_high_perf <- apply(-matrix(perf_ind_pot[sites_vacant, ], ncol=nsp), 1, high_perf_sp, sp_pres=sp_present)
            
            # Recruitment
            community_rast[sites_vacant] <- sp_high_perf
            
            # Update performance matrix
            # not OK if we want to change only recruited species:
            #perf_ind[sites_vacant, sp_present] <- perf_ind_pot[sites_vacant, sp_present]
            # OK:
            # /!\ remove loop?
            for(k in 1:nsite_vacant){
              perf_ind[sites_vacant[k], sp_high_perf[k]] <- perf_ind_pot[sites_vacant[k], sp_high_perf[k]]
            }
            
          } # end condition partial knowledge + IV
          
        } # end condition no abundance dependence
        
        if (nb_seeds_dep_abund==TRUE){
          
          abund_after_mortality <- as.data.frame(table(factor(as.vector(community), levels=1:nsp)))$Freq
          nb_seeds_sp <- round(abund_after_mortality*fecundity)
          nb_seeds_tot <- sum(nb_seeds_sp)
          
          # Each seed is dispersed to a random vacant site; several seeds can land on the same site.
          sites_of_seeds <- sample(sites_vacant, nb_seeds_tot, replace=TRUE)
          
          # Performance on vacant sites of species which have seeds
          seeds <- data.frame(Species=rep(1:nsp, times=nb_seeds_sp), Sites=sites_of_seeds)
          
          if (nb_seeds_tot>1) {
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
          } # end else (one seed)
          
          new_ind <- plyr::ddply(seeds, c("Sites"), f_sp)
          colnames(new_ind) <- c("Sites", "Species")
          
          # Recruitment
          community_rast[new_ind$Sites] <- new_ind$Species
          
          # Update performance matrix with recruited seeds
          if(perf_know==FALSE&&IV==TRUE){
            # /!\ remove loop?
            for(k in 1:nrow(new_ind)){
              perf_ind[new_ind$Sites[k], new_ind$Species[k]] <- perf_ind_pot[new_ind$Sites[k], new_ind$Species[k]]
            }
          }
          
        } # end condition on abundance dependence of number of seeds
        
      }  # end condition on n_mort (contains all cases)
      
      # update community
      community <- as.matrix(community_rast)
      
      # *********************
      # Diversity
      # *********************
      
      #print("Computing abundances and ecological indices...")
      
      # Environmental filtering
      if(IV==FALSE&&perf_know==TRUE){
        #Some sites have a 0 but then they are not taken into account since there is no distance associated with "species 0".
        dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
      }
      
      if(perf_know==FALSE){
        
        # Environmental filtering for the partial knowledge model
        #distance between the observed optimum of the species and the environmental conditions
        
        E_seq_mat <- list()
        
        if(n_observed_axes>1){
          
          E_seq <- matrix(nrow=100, ncol=n_observed_axes)
          
          for(k in 1:ncol(Obs_env)){
            E_seq[,k] <- seq(min(Obs_env[,k]), max(Obs_env[,k]), length.out=nrow(E_seq))
            E_seq_mat[[k]] <- matrix(rep(E_seq[,k], nsp), ncol=nsp)
            E_seq_mat[[k+ncol(Obs_env)]] <- E_seq_mat[[k]]^2
          }
          
        }else{
          E_seq <- seq(min(Obs_env), max(Obs_env), length.out=100)
          E_seq_mat[[1]] <- matrix(rep(E_seq, nsp), ncol=nsp)
          E_seq_mat[[2]] <- E_seq_mat[[1]]^2
        }
        
        Inferred_parameters_mat_E_seq <-list()
        
        if(n_observed_axes>1){
          for(k in 1:ncol(Inferred_species_parameters)){
            Inferred_parameters_mat_E_seq[[k]] <- matrix(rep(Inferred_species_parameters[,k],each=nrow(E_seq)), ncol=nsp)
          }
        }else{
          for(k in 1:ncol(Inferred_species_parameters)){
            Inferred_parameters_mat_E_seq[[k]] <- matrix(rep(Inferred_species_parameters[,k],each=length(E_seq)), ncol=nsp)
          }
        }
        
        Mat_perf_inferred <- Inferred_parameters_mat_E_seq[[1]]
        
        for(k in 1:length(E_seq_mat)){
          Mat_perf_inferred <- Mat_perf_inferred + Inferred_parameters_mat_E_seq[[k+1]]*E_seq_mat[[k]]
        }
        
        if(r==1&&g==1){
          plot_inferred_perf_environment(E_seq, Mat_perf_inferred, nsp, model, fig_width)
        }
        
        if(n_observed_axes>1){
          Optimum_Sp_inferred <- matrix(nrow=n_observed_axes, ncol=nsp)
          for(sp in 1:nsp){
            for(axis in 1:n_observed_axes){
              Optimum_Sp_inferred[axis,sp] <- (-Inferred_species_parameters[sp,axis+1])/(2*Inferred_species_parameters[sp,2*axis+1])
            }
          }
          Optimum_Sp_inferred <- t(Optimum_Sp_inferred)
        }else{
          Optimum_Sp_inferred <- c()
          for(sp in 1:nsp){
            Optimum_Sp_inferred[sp] <- (-Inferred_species_parameters[sp,2])/(2*Inferred_species_parameters[sp,3])
          }
        }
        
        dist_site <- dist_Site_Sp(as.matrix(Obs_env), as.matrix(Optimum_Sp_inferred))
        dist_site <- diag(dist_site[, as.vector(t(community))])
        
      }#end condition perf_know==FALSE
      
      env_filt[g, r] <- mean(dist_site)
      
      # Mean mortality rate in the community
      #/!\ do we want this mortality rate to take empty sites into account?
      if(mortality_stocha==TRUE){
        if (IV==FALSE&&perf_know==TRUE){
          theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
        }
        if (IV==TRUE&&perf_know==FALSE){
          theta_site <- diag(mortality_ind[, as.vector(t(community))])
        }
        theta_comm[g, r] <- mean(theta_site)
      }
      
    } # End ngen
    
    # store final community and plot it once
    community_end[[r]] <- community
    
    if (r==1) {
      plot_community_end(model=model, fig_width=fig_width, community=community_end[[r]], nsp=nsp)
    }
    
    # Species rank
    rank_sp[r, ] <- rank(-abund[ngen, ], ties.method="max")
    
    df_shannon <- data.frame(Species = 1:nsp,
                             Abundance = abund[ngen, ])%>%
      mutate(Proportion = Abundance / sum(Abundance))%>%
      filter(Abundance > 0)%>%
      mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
    
    Shannon[r] <- -sum(df_shannon$prop_times_ln_prop)
    
    #To keep the abundance matrixes in order to infer alpha matrix
    Abundances[[r]] <- abund
    
  } # End nrep
  
  save(community_start, file=here::here("outputs", model, "community_start.RData"))
  save(community_end, file=here::here("outputs", model, "community_end.RData"))
  
  sp_rich <- data.frame(sp_rich)
  env_filt <- data.frame(env_filt)
  theta_comm <- data.frame(theta_comm)
  
  save(Abundances, file = here::here("outputs", model, "Abundances.RData"))
  save(sp_rich, file=here::here("outputs", model, "sp_rich.RData"))
  save(rank_sp, file=here::here("outputs", model, "rank_sp.RData"))
  save(env_filt, file = here::here("outputs", model, "env_filt.RData"))
  save(theta_comm, file = here::here("outputs", model, "theta_comm.RData"))
  
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
  save(sp_rich_final, file=here::here("outputs", model, "Species_richness.RData"))
  
  # ---------------------------------------------
  # Shannon index and Shannon equitability index
  # ---------------------------------------------
  save(Shannon, file = here::here("outputs", model, "Shannon.RData"))
  Equitability <- Shannon/log(as.numeric(sp_rich[ngen,]))
  save(Equitability, file = here::here("outputs", model, "Equitability.RData"))
  
  # ---------------------------------------------------------------------------------
  # pairwise Spearman correlation on the species ranks at the end of each simulation
  # ---------------------------------------------------------------------------------
  
  Spearman <- as.dist(round(cor(t(rank_sp), method="spearman"),2))
  save(Spearman, file = here::here("outputs", model, "Spearman.RData"))
  
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
  plot_env_species(model, fig_width, community_start[[1]], community_end[[1]], class_site)
  
  # ---------------------------------------------
  # Theta community
  # ---------------------------------------------
  
  if(mortality_stocha==TRUE){
    plot_theta_community(theta_comm, ngen, nrep, model, fig_width)
  }
  
  # ----------------------------------
  # Spatial autocorrelation of species
  # ----------------------------------
  #plot_spatial_autocorr(nrep, community_end, n_axes, sites, niche_optimum, niche_width, model, fig_width)

  # ------------------------------------------------------
  # Performance of species that *should* win vs. *do* win
  # in the partial knowledge model
  # ------------------------------------------------------
  if(perf_know==FALSE){
    plot_perf_suitable_habitat(perf_Sp_mean, sites, fig_width)
    plot_perf_community_end(community_end, perf_Sp_mean, sites, fig_width)
  }
  
}