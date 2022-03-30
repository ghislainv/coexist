infer_IV <- function(model, n_observed_axis){
  
  source(file=here::here("_R", "Plot_functions.R"))
  
  load(here::here("outputs", model, "perf_E_Sp.RData"))
  load(here::here("outputs", model, "env.RData"))
  load(here::here("outputs", model, "sites.RData"))
  
  nsp <- ncol(perf_E_Sp)
  n_axis <- length(env)
  Obs_env <- sites[,1:n_observed_axis]
  
  # Data-set
  df <- data.frame(perf_E_Sp)
  names(df) <- c(sprintf("Sp_%02d", 1:nsp))
  
  df_perf <- data.frame(matrix(nrow=nrow(df), ncol=ncol(df)+2*n_axis))
  df_perf[,1:ncol(df)] <- df
  for(k in 1:n_axis){
    df_perf[,ncol(df)+k] <- raster::values(raster::raster(env[[k]]))
    df_perf[,ncol(df)+(k+n_axis)] <- (df_perf[,ncol(df)+k])^2
  }
  
  colnames(df_perf) <- c(sprintf("Sp_%02d", 1:nsp), sprintf("Env_%d", 1:n_axis), sprintf("Env_%d_sq", 1:n_axis))
  
  df_perf <- df_perf %>%
    tidyr::pivot_longer(cols=c(Sp_01:glue::glue("Sp_{nsp}")), names_to="Species", values_to="Perf")
  
  # Observed niche
  plot_species_niche(seed, df_perf, model, fig_width)
  
  #Check that the model well fits the data if all environmental variables are included
  # formula <- as.formula(paste0("Perf~-1+Species+Species:", paste0(colnames(df_perf)[1:(2*n_axis)], collapse= "+Species:")))
  # lm_all_env <- lm(formula, data=df_perf)
  
  # Observed intraspecific variability
  
  formula <- as.formula(paste0("Perf~-1+Species+Species:",
                               paste0(colnames(df_perf)[1:n_observed_axis], collapse= "+Species:"),
                               "+Species:",
                               paste0(colnames(df_perf)[(n_axis+1):(n_axis+n_observed_axis)], collapse= "+Species:")))
  
  lm_fit <- lm(formula, data=df_perf)
  save(lm_fit, file = here::here("outputs", model, glue::glue("lm_fit_{n_observed_axis}_obs_axes.RData")))
  
  #Histogram of the residuals to check normality
  hist(summary(lm_fit)$residuals)
  
  Inferred_species_parameters <- data.frame(matrix(nrow=nsp, ncol=1+2*n_observed_axis))
  Inferred_species_parameters[,1] <- as.vector(lm_fit$coefficients[1:nsp])
  colnames(Inferred_species_parameters)[1]<-"beta_0"
  for(k in 1:(2*n_observed_axis)){
    colnames(Inferred_species_parameters)[k+1] <- glue::glue("beta_{k}")
  }
  
  for(k in 1:n_observed_axis){
    #with complete interactions
    #Inferred_species_parameters[,k+1] <- as.vector(c(lm_fit$coefficients[nsp+k], lm_fit$coefficients[nsp+k]+lm_fit$coefficients[(nsp+2*n_observed_axis+k*(nsp-1)-(nsp-2)):(nsp+2*n_observed_axis+k*(nsp-1))]))
    #Inferred_species_parameters[,k+1+n_observed_axis] <- as.vector(c(lm_fit$coefficients[nsp+n_observed_axis+k], lm_fit$coefficients[nsp+n_observed_axis+k]+lm_fit$coefficients[(nsp+2*n_observed_axis+k*(nsp-1)-(nsp-2)+n_observed_axis*(nsp-1)):(nsp+2*n_observed_axis+k*(nsp-1)+n_observed_axis*(nsp-1))]))
    Inferred_species_parameters[,k+1] <- as.vector(c(lm_fit$coefficients[(k*nsp+1):((k+1)*nsp)]))
    Inferred_species_parameters[,k+1+n_observed_axis] <- as.vector(c(lm_fit$coefficients[(nsp*(n_observed_axis+k)+1):(nsp*(n_observed_axis+k+1))]))
  }
  
  save(Inferred_species_parameters, file=here::here("outputs", model, glue::glue("Inferred_species_parameters_{n_observed_axis}_obs_axes.RData")))
  
  load(file=here::here("outputs", model, glue::glue("Inferred_species_parameters_{n_observed_axis}_obs_axes.RData")))
  load(file=here::here("outputs", model, "niche_optimum.RData"))
  
  plot_optima_real_estim(nsp, n_observed_axis, niche_optimum, Inferred_species_parameters, model, fig_width)
  
  V_intra <- df_perf %>%
    mutate(res=lm_fit$residuals) %>%
    group_by(Species) %>%
    summarise(V=var(res))
  save(V_intra, file = here::here("outputs", model, glue::glue("V_intra_{n_observed_axis}_obs_axes.RData")))
  
  load(here::here("outputs", model, glue::glue("V_intra_{n_observed_axis}_obs_axes.RData")))
  
  plot_IV(V_intra, model, fig_width)
  
  plot_inferred_perf_IV(n_observed_axis, Obs_env=sites[,1:n_observed_axis], nsp, Inferred_species_parameters, V_intra, model, fig_width)
  
  return(V_intra)
  
}