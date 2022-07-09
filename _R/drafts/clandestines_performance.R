mortality <- "prop"
fecundity <- "fixed"

Comparison_clandestines <- data.frame(
  Mod=character(),
  Nb_obs=integer(),
  Seed=integer(),
  Seed_r=integer(),
  perf_clandestines=numeric(),
  perf_not_clandestines=numeric(),
  perf_tot=numeric())

for(seed in 1:10){
  for(seed_r in 1:10){
    mod <- "Perf_know"
    load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
    
    #clandestines <- unique(as.vector(t(community_end))[which((winner==as.vector(t(community_end)))==FALSE)])
    #mean_perf_clandestines <- mean(apply(perf_E_Sp[,clandestines], 2, mean))
    #nb_sites_clandestines <- length(which((winner==as.vector(t(community_end)))==FALSE))
    
    
    perf_present <- rep(NA, length(community_end))
    w0 <- (as.vector(t(community_end))==0)
    perf_present[w0] <- NA
    perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
    
    
    Comparison_clandestines_temp <- data.frame(
      Mod=mod,
      Nb_obs=NA,
      Seed=seed,
      Seed_r=seed_r,
      perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
      perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
      perf_tot=round(mean(perf_present, na.rm = TRUE), digits = 2))
    
    Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
    
    mod <- "Part_know_IV"
    for (nb_obs in 1:15){
      
      load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
      
      #clandestines <- unique(as.vector(t(community_end))[which((winner==as.vector(t(community_end)))==FALSE)])
      #mean_perf_clandestines <- mean(apply(perf_E_Sp[,clandestines], 2, mean))
      #nb_sites_clandestines <- length(which((winner==as.vector(t(community_end)))==FALSE))
      
      
      perf_present <- rep(NA, length(community_end))
      w0 <- (as.vector(t(community_end))==0)
      perf_present[w0] <- NA
      perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
      
      Comparison_clandestines_temp <- data.frame(
        Mod=mod,
        Nb_obs=nb_obs,
        Seed=seed,
        Seed_r=seed_r,
        perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
        perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
        perf_tot=round(mean(perf_present, na.rm = TRUE), digits = 2))
      
      Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
      
    }
  }
}