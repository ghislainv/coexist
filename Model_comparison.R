compare_models<-function(){
  dir.create(here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs")))
  # models_all <- c("Perf_know_full_mort_stocha",
  #                 "Perf_know_full_mort_stocha_disp_abund",
  #                 "Perf_know_start_1_mort_stocha",
  #                 "Perf_know_start_1_mort_stocha_disp_abund",
  #                 "Perf_know_start_10_mort_stocha",
  #                 "Perf_know_start_10_mort_stocha_disp_abund",
  #                 "Perf_know_full_mort_fixed",
  #                 "Perf_know_full_mort_fixed_disp_abund",
  #                 "Perf_know_start_10_mort_fixed",
  #                 "Perf_know_start_10_mort_fixed_disp_abund",
  #                 "Part_know_IV_full_mort_stocha",
  #                 "Part_know_IV_full_mort_stocha_disp_abund",
  #                 "Part_know_IV_start_1_mort_stocha",
  #                 "Part_know_IV_start_1_mort_stocha_disp_abund",
  #                 "Part_know_IV_start_10_mort_stocha",
  #                 "Part_know_IV_start_10_mort_stocha_disp_abund",
  #                 "Part_know_IV_full_mort_fixed",
  #                 "Part_know_IV_full_mort_fixed_disp_abund",
  #                 "Part_know_IV_start_10_mort_fixed",
  #                 "Part_know_IV_start_10_mort_fixed_disp_abund",
  #                 "Part_know_full_mort_stocha",
  #                 "Part_know_full_mort_stocha_disp_abund",
  #                 "Part_know_start_1_mort_stocha",
  #                 "Part_know_start_1_mort_stocha_disp_abund",
  #                 "Part_know_full_mort_fixed",
  #                 "Part_know_full_mort_fixed_disp_abund",
  #                 "Part_know_start_10_mort_fixed",
  #                 "Part_know_start_10_mort_fixed_disp_abund")
  # 
  # models_choices <- c("Perf_know_start_10_mort_fixed",
  #                     "Perf_know_start_10_mort_fixed_disp_abund",
  #                     "Part_know_start_10_mort_fixed",
  #                     "Part_know_start_10_mort_fixed_disp_abund",
  #                     "Part_know_IV_start_10_mort_fixed",
  #                     "Part_know_IV_start_10_mort_fixed_disp_abund")
  # 
  # models_choices <- c("Perf_know_start_10_mort_fixed_disp_abund_10_axes_1_obs",
  #                     "Perf_know_start_10_mort_fixed_disp_abund_10_axes_5_obs",
  #                     "Part_know_start_10_mort_fixed_disp_abund_10_axes_1_obs",
  #                     "Part_know_start_10_mort_fixed_disp_abund_10_axes_5_obs",
  #                     "Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_1_obs",
  #                     "Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_5_obs")
  
  models_choices <- c(glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_{n_observed_axis}_obs_seed_{seed}"),
                      glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_{n_observed_axis}_obs_seed_{seed}"),
                      glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_{n_observed_axis}_obs_seed_{seed}"))
  
  models <- models_choices
  
  # 1: Compare the species diversity at the end of the simulations within a model (Shannon diversity index)
  
  Shannon_all_models <- data.frame(Shannon=numeric(), Model=factor())
  
  for (m in 1:length(models)) {
    model <- models[m]
    load(here::here("outputs", model, glue::glue("Shannon_{model}.RData")))
    Shannon_all_models <- rbind(Shannon_all_models, data.frame(Shannon=Shannon, Model=rep(m, length(Shannon))))
  }
  
  p <- ggplot2::ggplot(data=Shannon_all_models, ggplot2::aes(x=as.factor(Model), y=Shannon))+
    ggplot2::geom_boxplot()+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Model)), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#80002D", "#008071", "#088000"))+
    #ggplot2::scale_colour_manual(values=c("#80002D", "#BF0043", "#008071", "#00BFA9", "#088000", "#0DBF00"))+
    ggplot2::labs(title = "Beeswarmplot of the Shannon diversity index \n at the end of each simulation",
                  x = "Model",
                  y = "Shannon diversity index")+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Shannon_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  
  # 2: Compare the species ranks in the end of the simulations within a model (Spearman pairwise correlation)
  
  Spearman_all_models <- data.frame(Spearman=numeric(), Model=factor())
  
  for (m in 1:length(models)) {
    model <- models[m]
    load(here::here("outputs", model, glue::glue("Spearman_{model}.RData")))
    Spearman_all_models <- rbind(Spearman_all_models, data.frame(Spearman=c(Spearman), Model=rep(m, length(Spearman))))
  }
  
  p <- ggplot2::ggplot(data=Spearman_all_models, ggplot2::aes(x=as.factor(Model), y=Spearman))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Model)), alpha=0.1)+
    #ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Model)), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#80002D", "#008071", "#088000"))+
    #ggplot2::scale_colour_manual(values=c("#80002D", "#BF0043", "#008071", "#00BFA9", "#088000", "#0DBF00"))+
    # ggplot2::labs(title = "Beeswamplot of the Spearman pairwise correlation \n of species ranks at the end of each simulation",
    #               x = "Model",
    #               y = "Spearman pairwise correlation")+
    ggplot2::labs(title = "Jitterplot of the Spearman pairwise correlation \n of species ranks at the end of each simulation",
                  x = "Model",
                  y = "Spearman pairwise correlation")+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Spearman_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  # 3: Compare the composition of the community at the end of the simulations between models
  # Jaccard index (similarity matrix)
  
  # 4:  Percentage similarity of species abudances
  
  jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
  }
  
  percentage_similarity <- function(a, b) {
    A <- sum(a)
    B <- sum(b)
    W <- sum(pmin(a, b))
    return((2*W)/(A+B))
  }
  
  Jaccard_all_models <- data.frame(Jaccard=numeric(), Model=factor())
  Percentage_similarity_all_models <- data.frame(Percentage_similarity=numeric(), Model=factor())
  
  combi_models <- gtools::combinations(n = length(c(1:length(models))), r = 2, v = c(1:length(models)), repeats.allowed = TRUE)
  
  combi_rep <- expand.grid(A=c(1:nrep), B=c(1:nrep))
  
  combi_rep_same_model <- t(combn(c(1:nrep), 2))
  
  for (m in 1:length(models)) {
    model <- models[m]
    load(here::here("outputs", model, glue::glue("Abundances_{model}.RData")))
    abund_end <- matrix(nrow=nrep, ncol=nsp)
    species_end <- list()
    for (r in 1:nrep) {
      abund_end[r,] <- Abundances[[r]][ngen,]
      species_end[[r]] <- unique(which(abund_end[r,]!=0))
    }
    assign(glue::glue("abund_end_{m}"), abund_end)
    assign(glue::glue("species_end_{m}"), species_end)
  }
  
  list_jaccard_matrix <- list()
  mean_jaccard <- c()
  
  list_percentage_similarity_matrix <- list()
  mean_percentage_similarity <- c()
  
  vec_jaccard <- c()
  vec_combi_models_jaccard <- c()
  vec_combi_rep_jaccard <- c()
  
  vec_percentage_similarity <- c()
  vec_combi_models_percentage_similarity <- c()
  vec_combi_rep_percentage_similarity <- c()
  
  for (l in 1:nrow(combi_models)){
    mod1 <- combi_models[l, 1]
    mod2 <- combi_models[l, 2]
    if(mod1==mod2){
      list_jaccard_matrix[[l]] <- matrix(nrow=nrow(combi_rep_same_model), ncol=1)
      list_percentage_similarity_matrix[[l]] <- matrix(nrow=nrow(combi_rep_same_model), ncol=1)
      for (c in 1:nrow(combi_rep_same_model)){
        list_jaccard_matrix[[l]][c,1] <- jaccard(a=get(glue::glue("species_end_{mod1}"))[[combi_rep_same_model[c,1]]], b=get(glue::glue("species_end_{mod2}"))[[combi_rep_same_model[c,2]]])
        list_percentage_similarity_matrix[[l]][c,1] <- percentage_similarity(a=get(glue::glue("abund_end_{mod1}"))[combi_rep_same_model[c,1],], b=get(glue::glue("abund_end_{mod2}"))[combi_rep_same_model[c,2],])
      }
      vec_percentage_similarity <- c(vec_percentage_similarity, c(list_percentage_similarity_matrix[[l]]))
      vec_combi_models_percentage_similarity <- c(vec_combi_models_percentage_similarity, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep_same_model)))
      vec_combi_rep_percentage_similarity <- c(vec_combi_rep_percentage_similarity, paste0(combi_rep_same_model[,1], "-", combi_rep_same_model[,2]))
      vec_jaccard <- c(vec_jaccard, c(list_jaccard_matrix[[l]]))
      vec_combi_models_jaccard <- c(vec_combi_models_jaccard, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep_same_model)))
      vec_combi_rep_jaccard <- c(vec_combi_rep_jaccard, paste0(combi_rep_same_model[,1], "-", combi_rep_same_model[,2]))
    }
    else {
      list_jaccard_matrix[[l]] <- matrix(nrow=nrep, ncol=nrep)
      list_percentage_similarity_matrix[[l]] <- matrix(nrow=nrep, ncol=nrep)
      for (c in 1:nrow(combi_rep)){
        list_jaccard_matrix[[l]][combi_rep[c,1], combi_rep[c,2]] <- jaccard(a=get(glue::glue("species_end_{mod1}"))[[combi_rep[c,1]]], b=get(glue::glue("species_end_{mod2}"))[[combi_rep[c,2]]])
        list_percentage_similarity_matrix[[l]][combi_rep[c,1], combi_rep[c,2]] <- percentage_similarity(a=get(glue::glue("abund_end_{mod1}"))[combi_rep[c,1],], b=get(glue::glue("abund_end_{mod2}"))[combi_rep[c,2],])
      }
      vec_percentage_similarity <- c(vec_percentage_similarity, c(list_percentage_similarity_matrix[[l]]))
      vec_combi_models_percentage_similarity <- c(vec_combi_models_percentage_similarity, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep)))
      vec_combi_rep_percentage_similarity <- c(vec_combi_rep_percentage_similarity, paste0(combi_rep[,1], "-", combi_rep[,2]))
      vec_jaccard <- c(vec_jaccard, c(list_jaccard_matrix[[l]]))
      vec_combi_models_jaccard <- c(vec_combi_models_jaccard, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep)))
      vec_combi_rep_jaccard <- c(vec_combi_rep_jaccard, paste0(combi_rep[,1], "-", combi_rep[,2]))
    }
    mean_jaccard[l] <- mean(list_jaccard_matrix[[l]])
    mean_percentage_similarity[l] <- mean(list_percentage_similarity_matrix[[l]])
  }
  
  combi_models_jaccard <- as.data.frame(combi_models)
  combi_models_jaccard$Jaccard <- mean_jaccard
  colnames(combi_models_jaccard)[1:2] <- c("mod1", "mod2")
  
  combi_models_percentage_similarity <- as.data.frame(combi_models)
  combi_models_percentage_similarity$Percentage_similarity <- mean_percentage_similarity
  colnames(combi_models_percentage_similarity)[1:2] <- c("mod1", "mod2")
  
  p <- ggplot2::ggplot(data = combi_models_jaccard, ggplot2::aes(x=mod1, y=mod2, fill=Jaccard))+ 
    ggplot2::geom_tile()+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::coord_fixed()+
    ggplot2::labs(title="Jaccard index of the composition of the final community \n (mean of all repetitions)",
                  x="Model 1",
                  y = "Model 2",
                  fill="Jaccard index")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Jaccard_tileplot.png"),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  p <- ggplot2::ggplot(data = combi_models_percentage_similarity, ggplot2::aes(x=mod1, y=mod2, fill=Percentage_similarity))+ 
    ggplot2::geom_tile()+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::coord_fixed()+
    ggplot2::labs(title="Percentage similarity index of the species abundances \n of the final community (mean of all repetitions)",
                  x="Model 1",
                  y = "Model 2",
                  fill="Percentage similarity")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Percentage_similarity_tileplot.png"),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  Jaccard_all_models <- data.frame(jaccard=vec_jaccard,
                                   combi_models=vec_combi_models_jaccard,
                                   combi_rep=vec_combi_rep_jaccard)
  
  p <- ggplot2::ggplot(data=Jaccard_all_models, ggplot2::aes(x=as.factor(combi_models), y=jaccard))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(combi_models)), alpha=0.1)+
    #ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(combi_models)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    # ggplot2::labs(title = "Beeswarmplot of the Jaccard index of the composition \n of the final community; \n each point is a comparison between two repetitions ",
    #               x = "Model combination",
    #               y = "Jaccard index of the composition of the community")+
    ggplot2::labs(title = "Jitterplot of the Jaccard index of the composition \n of the final community; \n each point is a comparison between two repetitions ",
                  x = "Model combination",
                  y = "Jaccard index of the composition of the community")+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Jaccard_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Percentage_similarity_all_models <- data.frame(percentage_similarity=vec_percentage_similarity,
                                                 combi_models=vec_combi_models_percentage_similarity,
                                                 combi_rep=vec_combi_rep_percentage_similarity)
  
  p <- ggplot2::ggplot(data=Percentage_similarity_all_models, ggplot2::aes(x=as.factor(combi_models), y=percentage_similarity))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(combi_models)), alpha=0.1)+
    #ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(combi_models)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    # ggplot2::labs(title = "Beeswarmplot of the percentage similarity of species abundance; \n each point is a comparison between two repetitions ",
    #               x = "Model combination",
    #               y = "Percentage similarity of species abundance")+
    ggplot2::labs(title = "Jitterplot of the percentage similarity of species abundance; \n each point is a comparison between two repetitions ",
                  x = "Model combination",
                  y = "Percentage similarity of species abundance")+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}_{n_observed_axis}_obs"), "Percentage_similarity_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
}

Compare_IV_axis_nb <- function(seed){
  dir.create(here::here("outputs", glue::glue("Comparison_seed_{seed}")))
  
  models <- c(glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_1_obs_seed_{seed}"),
              glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_3_obs_seed_{seed}"),
              glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_5_obs_seed_{seed}"),
              glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_7_obs_seed_{seed}"))
  
  nb_obs_axes <- c(1, 3, 5, 7)
  
  IV_all_models <- data.frame(IV=numeric(), Model=factor())
  
  for (m in 1:length(models)) {
    model <- models[m]
    load(here::here("outputs", model, "V_intra.RData"))
    IV_all_models <- rbind(IV_all_models, data.frame(IV=V_intra$V, Model=rep(m, nrow(V_intra))))
  }
  
  load(here::here("outputs", models[1], "sites.RData"))
  sum_var <- 0
  for(k in 1:ncol(sites)){sum_var <- sum_var + var(sites[,k])}
  percentage_inertia <- c()
  percentage_inertia_cum <- c()
  for(k in 1:ncol(sites)){
    percentage_inertia[k]<-(var(sites[k])/sum_var)*100
    if(k == 1){percentage_inertia_cum[k]<-percentage_inertia[k]}
    if(k > 1){percentage_inertia_cum[k]<-percentage_inertia_cum[k-1]+percentage_inertia[k]}
  }
  
  p <- ggplot2::ggplot(data=IV_all_models, ggplot2::aes(x=as.factor(Model), y=IV))+
    ggplot2::geom_boxplot()+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Model)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Beeswarmplot of the intraspecific variability \n inferred with different levels of knowledge",
                  x = "Number of observed axes",
                  y = "IV")+
    ggplot2::scale_x_discrete(labels=nb_obs_axes)+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")+
    ggplot2::geom_text(ggplot2::aes(x=1, y=1), label=paste(round(percentage_inertia_cum[nb_obs_axes])[1], "%", "inertia")) +
    ggplot2::geom_text(ggplot2::aes(x=2, y=1), label=paste(round(percentage_inertia_cum[nb_obs_axes])[2], "%"))+
    ggplot2::geom_text(ggplot2::aes(x=3, y=1), label=paste(round(percentage_inertia_cum[nb_obs_axes])[3], "%"))+
    ggplot2::geom_text(ggplot2::aes(x=4, y=1), label=paste(round(percentage_inertia_cum[nb_obs_axes])[4], "%"))
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_seed_{seed}"), "IV_nb_axes.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
}

Compare_spatial_structure <- function(){
  models <- c(glue::glue("Perf_know_start_10_mort_fixed_disp_abund_10_axes_1_obs_seed_{seed}"),
              glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_1_obs_seed_{seed}"),
              glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_3_obs_seed_{seed}"),
              glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_5_obs_seed_{seed}"),
              glue::glue("Part_know_start_10_mort_fixed_disp_abund_10_axes_7_obs_seed_{seed}"),
              glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_1_obs_seed_{seed}"),
              glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_3_obs_seed_{seed}"),
              glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_5_obs_seed_{seed}"),
              glue::glue("Part_know_IV_start_10_mort_fixed_disp_abund_10_axes_7_obs_seed_{seed}"))
  model_names <- c("Perfect knowledge",
                   "Partial knowledge \n (1 observed axis)",
                   "Partial knowledge \n (3 observed axis)",
                   "Partial knowledge \n (5 observed axis)",
                   "Partial knowledge \n (7 observed axis)",
                   "Partial knowledge + IV \n (1 observed axis)",
                   "Partial knowledge + IV \n (3 observed axis)",
                   "Partial knowledge + IV \n (5 observed axis)",
                   "Partial knowledge + IV \n (7 observed axis)")
  # Models <- c("Perf", rep("Part", 4), rep("Part_IV", 4))
  # Nb_axes <- c(" ", rep(c(1, 3, 5, 7), 2))
  nb_obs_axes <- c(1, 3, 5, 7)
  
  # semivar_all_models <- data.frame(Semivar_sp=numeric(), Semivar_env=numeric(), Model=factor(), Nb_axes=factor())
  # 
  # for (m in 1:length(models)) {
  #   model <- models[m]
  #   load(here::here("outputs", model, glue::glue("semivar_multidim.RData")))
  #   if(semivar_multidim$Sample_size[nrow(semivar_multidim)]<500){
  #     semivar_multidim <- semivar_multidim[1:(nrow(semivar_multidim)-1),]
  #   }
  #   semivar_all_models <- rbind(semivar_all_models, data.frame(Semivar_sp=semivar_multidim$vario_sp, Semivar_env=semivar_multidim$Vario_env, Model=rep(Models[m], nrow(semivar_multidim)), Nb_axes=rep(Nb_axes[m], nrow(semivar_multidim))))
  # }
  
  png(file=here::here("outputs", glue::glue("Comparison_seed_{seed}"), "Semivar_nb_axes.png"),
      width=fig_width, height=fig_width*0.8, units="cm", res=300)
  par(mfrow=c(length(nb_obs_axes),3), bty = "n")
  
  for (m in 1:length(models)) {
    model <- models[m]
    load(here::here("outputs", model, glue::glue("semivar_multidim.RData")))
    # if(semivar_multidim$Sample_size[nrow(semivar_multidim)]<500){
    #   semivar_multidim <- semivar_multidim[1:(nrow(semivar_multidim)-1),]
    # }
    plot(semivar_multidim$Vario_env, semivar_multidim$Vario_sp,
         main=model_names[m],
         xlab="Semivariance for environment",
         ylab="Semivariance for species")
    m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)
    abline(a=as.numeric(coef(m)[1]), b=as.numeric(coef(m)[2]), col="#008071")
  }
  
  # for(i in unique(semivar_all_models$Model)){
  #   for(j in nb_obs_axes){
  #     
  #   }
  # }
  
  dev.off()
  
}
