infer_IV <- function(model, n_observed_axis){
  
  source(file=here::here("Plot_functions.R"))
  
  load(here::here("outputs", model, "perf_E_Sp.RData"))
  
  load(here::here("outputs", model, "env.RData"))
  load(here::here("outputs", model, "sites.RData"))
  nsp <- ncol(perf_E_Sp)
  # CGT 01/12/2021
  n_axis <- length(env)
  Obs_env <- sites[,1:n_observed_axis]
  
  # Data-set
  df <- data.frame(perf_E_Sp)
  names(df) <- c(sprintf("Sp_%02d", 1:nsp))
  
  #df_perf <- tibble(df) %>%
  #mutate(Env_1=values(raster(env[[1]])), Env_2=values(raster(env[[2]])), Env_3=values(raster(env[[3]]))) %>%
  #mutate(Env_1_sq=Env_1^2, Env_2_sq=Env_2^2, Env_3_sq=Env_3^2) %>%
  #pivot_longer(cols=c(Sp_001:glue("Sp_0{nsp-1}"), glue("Sp_{nsp}")), names_to="Species", values_to="Perf")
  
  df_perf <- data.frame(matrix(nrow=nrow(df), ncol=ncol(df)+2*n_axis))
  df_perf[,1:ncol(df)] <- df
  for(k in 1:n_axis){
    df_perf[,ncol(df)+k] <- values(raster(env[[k]]))
    df_perf[,ncol(df)+(k+n_axis)] <- (df_perf[,ncol(df)+k])^2
  }
  
  colnames(df_perf) <- c(sprintf("Sp_%02d", 1:nsp), sprintf("Env_%d", 1:n_axis), sprintf("Env_%d_sq", 1:n_axis))
  
  df_perf <- df_perf %>%
    pivot_longer(cols=c(Sp_01:glue("Sp_{nsp}")), names_to="Species", values_to="Perf")
  
  # Observed niche
  plot_species_niche(seed, df_perf, model, fig_width)
  
  #Check that the model well fits the data if all environmental variables are included
  
  #lm_all_env <- lm(Perf~Species+Species*Env_1+Species*Env_1_sq+Species*Env_2+Species*Env_2_sq+Species*Env_3+Species*Env_3_sq, data=df_perf)
  
  formula <- as.formula(paste0("Perf~-1+Species+Species:", paste0(colnames(df_perf)[1:(2*n_axis)], collapse= "+Species:")))
  
  lm_all_env <- lm(formula, data=df_perf)
  
  print(glue::glue("the r-squared with all environmental variables is {summary(lm_all_env)$adj.r.squared}"))
  
  # Observed intraspecific variability
  
  formula <- as.formula(paste0("Perf~-1+Species+Species:",
                               paste0(colnames(df_perf)[1:n_observed_axis], collapse= "+Species:"),
                               "+Species:",
                               paste0(colnames(df_perf)[(n_axis+1):(n_axis+n_observed_axis)], collapse= "+Species:")))
  
  #lm_fit <- lm(Perf~Species+Species*Env_1+Species*Env_1_sq, data=df_perf)
  lm_fit <- lm(formula, data=df_perf)
  save(lm_fit, file = here::here("outputs", model,"lm_fit.RData"))
  
  print(glue::glue("the r-squared with the first environmental variable is {summary(lm_fit)$adj.r.squared}"))
  
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
  
  save(Inferred_species_parameters, file=here::here("outputs", model, "Inferred_species_parameters.RData"))
  
  plot_optima_real_estim(nsp, n_observed_axis, niche_optimum, Inferred_species_parameters, model, fig_width)
  
  V_intra <- df_perf %>%
    mutate(res=lm_fit$residuals) %>%
    group_by(Species) %>%
    summarise(V=var(res))
  save(V_intra, file = here::here("outputs", model, "V_intra.RData"))
  
  plot_IV(V_intra, model, fig_width)
  
  plot_inferred_perf_IV(n_observed_axis, Obs_env, nsp, Inferred_species_parameters, V_intra, model, fig_width)
  
  return(V_intra)
  
}