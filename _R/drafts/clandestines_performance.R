
load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_perf_E_Sp.RData")))
winner <- apply(X=perf_E_Sp, MARGIN=1, FUN=which.max)

Comparison_clandestines <- data.frame(
  Mod=character(),
  Nb_obs=integer(),
  Seed=integer(),
  Seed_r=integer(),
  Abund_clandestines=numeric(),
  perf_clandestines=numeric(),
  perf_not_clandestines=numeric(),
  perf_tot=numeric())

for(seed in 1:10){
  for(seed_r in 1:10){
    mod <- "Perf_know"
    load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
    
    
    perf_present <- rep(NA, length(community_end))
    w0 <- (as.vector(t(community_end))==0)
    perf_present[w0] <- NA
    perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
    
    
    Comparison_clandestines_temp <- data.frame(
      Mod=mod,
      Nb_obs=NA,
      Seed=seed,
      Seed_r=seed_r,
      Abund_clandestines=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
      perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)], na.rm=TRUE),
      perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)], na.rm=TRUE),
      perf_tot=mean(perf_present, na.rm = TRUE))
    
    Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
    
    for (mod in c("Part_know", "Part_know_IV")){
      for (nb_obs in 1:15){
        
        load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        
        
        perf_present <- rep(NA, length(community_end))
        w0 <- (as.vector(t(community_end))==0)
        perf_present[w0] <- NA
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
        
        Comparison_clandestines_temp <- data.frame(
          Mod=mod,
          Nb_obs=nb_obs,
          Seed=seed,
          Seed_r=seed_r,
          Abund_clandestines=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
          perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)], na.rm=TRUE),
          perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)], na.rm=TRUE),
          perf_tot=mean(perf_present, na.rm = TRUE))
        
        Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
        
      }
    }
  }
}