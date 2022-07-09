for(seed in Seeds){
  
  for(mort in c("mort_stocha", "mort_fixed")) {

    model <- glue::glue("Perf_know_start_10_{mort}_seeds_dep_abund_10_axes_seed_{seed}")
    
    load(here::here("outputs", model, "perf_E_Sp.RData"))
    load(here::here("outputs", model, "env.RData"))
    # Data-set
    df <- data.frame(perf_E_Sp)
    names(df) <- c(sprintf("Sp %02d", 1:nsp))
    
    df_perf <- data.frame(matrix(nrow=nrow(df), ncol=ncol(df)+2*n_axes))
    df_perf[,1:ncol(df)] <- df
    for(k in 1:n_axes){
      df_perf[,ncol(df)+k] <- raster::values(raster::raster(env[[k]]))
      df_perf[,ncol(df)+(k+n_axes)] <- (df_perf[,ncol(df)+k])^2
    }
    
    colnames(df_perf) <- c(sprintf("Sp %02d", 1:nsp), sprintf("Env_%d", 1:n_axes), sprintf("Env_%d_sq", 1:n_axes))
    
    df_perf <- df_perf %>%
      tidyr::pivot_longer(cols=c("Sp 01":glue::glue("Sp {nsp}")), names_to="Species", values_to="Perf")
    
    # Observed niche
    plot_species_niche(seed, df_perf, model, fig_width)
  }
}