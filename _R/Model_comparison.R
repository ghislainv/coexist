source(file=here::here("_R", "Math_functions.R"))

Compare_IV_axis_nb <- function(Seeds, nsp, nb_obs_axes){
  
  dir.create(here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}")))
  
  IV_all_models <- data.frame(Sp=(rep(1:nsp, length(Seeds)*length(nb_obs_axes))),
                              IV=numeric(nsp*length(nb_obs_axes)*length(Seeds)),
                              Seed=rep(Seeds, each=nsp*length(nb_obs_axes)),
                              Nb_obs_axes=rep(rep(nb_obs_axes, each=nsp), length(Seeds)))
  Percentage_inertia <- data.frame(Seed=rep(Seeds, each=n_axes),
                                   Nb_obs_axes=rep(1:n_axes, length(Seeds)),
                                   PI=numeric(n_axes*length(Seeds)),
                                   PI_cum=numeric(n_axes*length(Seeds)))
  Level_explanation_axes_nb <- data.frame(Seed=rep(Seeds, each=length(nb_obs_axes)),
                                          Nb_obs_axes=rep(nb_obs_axes, length(Seeds)),
                                          R2=numeric(length(Seeds)*length(nb_obs_axes)))
  
  for(seed in Seeds){
  
    model <- glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}")
    
    for (n_observed_axes in nb_obs_axes) {
      load(here::here("outputs", model, glue::glue("V_intra_{n_observed_axes}_obs_axes.RData")))
      IV_all_models[IV_all_models$Seed==seed&IV_all_models$Nb_obs_axes==n_observed_axes,]$IV <- V_intra$V
    }
    
    # Compute % of inertia of each environmental axis
    load(here::here("outputs", model, "sites.RData"))
    sum_var <- 0
    for(k in 1:ncol(sites)){sum_var <- sum_var + var(sites[,k])}
    percentage_inertia <- c()
    percentage_inertia_cum <- c()
    for(k in 1:ncol(sites)){
      percentage_inertia[k]<-(var(sites[k])/sum_var)*100
      if(k == 1){percentage_inertia_cum[k]<-percentage_inertia[k]}
      if(k > 1){percentage_inertia_cum[k]<-percentage_inertia_cum[k-1]+percentage_inertia[k]}
    }
    Percentage_inertia[Percentage_inertia$Seed==seed,]$PI <- percentage_inertia
    Percentage_inertia[Percentage_inertia$Seed==seed,]$PI_cum <- percentage_inertia_cum
    
    # Retrieve R2 of associated with each number of environmental axes
    for(n_observed_axes in nb_obs_axes){
      load(here::here("outputs", model, glue::glue("lm_fit_{n_observed_axes}_obs_axes.RData")))
      Level_explanation_axes_nb[Level_explanation_axes_nb$Seed==seed,]$R2[n_observed_axes] <- summary(lm_fit)$r.squared
    }
  }
  
  Summary_percentage_inertia <- Percentage_inertia%>%
    dplyr::group_by(Nb_obs_axes)%>%
    dplyr::mutate(Mean_PI_cum=mean(PI_cum), Sd_PI_cum=mean(PI_cum))%>%
    dplyr::slice(1)%>%
    dplyr::ungroup()%>%
    dplyr::select(-Seed, -PI, -PI_cum)
  
  Summary_level_explanation_axes_nb <- Level_explanation_axes_nb%>%
    dplyr::group_by(Nb_obs_axes)%>%
    dplyr::mutate(Mean_explanation=mean(R2), Sd=sd(R2))%>%
    dplyr::slice(1)%>%
    dplyr::ungroup()%>%
    dplyr::select(-Seed, -R2)
  
  p <- ggplot2::ggplot(data=IV_all_models, ggplot2::aes(x=as.factor(Nb_obs_axes), y=IV))+
    ggplot2::geom_ribbon(data=Summary_level_explanation_axes_nb, ggplot2::aes(x=Nb_obs_axes, y=Mean_explanation, ymin=Mean_explanation-Sd, ymax=Mean_explanation+Sd), colour="deeppink3", fill="hotpink3", alpha = 0.3)+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::geom_point(data=Summary_level_explanation_axes_nb, ggplot2::aes(x=Nb_obs_axes, y=Mean_explanation), colour="deeppink3")+
    ggplot2::geom_line(data=Summary_level_explanation_axes_nb, ggplot2::aes(x=Nb_obs_axes, y=Mean_explanation), colour="deeppink3")+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(x = "Number of observed axes",
                  y = "IV")+
    ggplot2::scale_x_discrete(labels=nb_obs_axes)+
    ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")+
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . * 1 / 1 , name = "Proportion of variance explained by the axes"))
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "IV_nb_axes.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
}


compare_models<-function(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes){
  
  dir.create(here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}")))
  
  # 1: Compare species richness and diversity at the end of the simulations within a model (Shannon diversity index)
  Species_richness_all_models <- data.frame(Species_richness=numeric(), Model=factor(), Nb_obs_axes=factor(), Seed=factor())
  Shannon_all_models <- data.frame(Shannon=numeric(), Model=factor(), Nb_obs_axes=factor(), Seed=factor())
  # 2: Compare the species ranks in the end of the simulations within a model (Spearman pairwise correlation)
  Spearman_all_models <- data.frame(Spearman=numeric(), Model=factor(), Nb_obs_axes=factor(), Seed=factor())
  # 3: Compare the composition of the community at the end of the simulations between models:
  #Jaccard index (similarity matrix)
  Jaccard_all_models <- data.frame(Jaccard=numeric(),
                                   Combi_model=factor(),
                                   Combi_rep=factor(),
                                   Nb_obs_axes=factor(),
                                   Seed=factor())
  # 4:  Percentage similarity of species abundances
  Percentage_similarity_all_models <- data.frame(Percentage_similarity=numeric(),
                                                 Combi_model=factor(),
                                                 Combi_rep=factor(),
                                                 Nb_obs_axes=factor(),
                                                 Seed=factor())
  
  for (seed in Seeds) {
    
    for(n_observed_axes in nb_obs_axes){
      
      models <- c(glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}"),
                  glue::glue("Part_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_{n_observed_axes}_obs_seed_{seed}"),
                  glue::glue("Part_know_IV_{start}_{mort}_{nb_seeds}_{n_axes}_axes_{n_observed_axes}_obs_seed_{seed}"))
      
      Mod_type <- c("Perf_know", "Part_know", "Part_know_IV")
      
      combi_models <- gtools::combinations(n = length(c(1:length(models))), r = 2, v = c(1:length(models)), repeats.allowed = TRUE)
      combi_rep <- expand.grid(A=c(1:nrep), B=c(1:nrep))
      combi_rep_same_model <- t(combn(c(1:nrep), 2))
      
      for (m in 1:length(models)) {
        model <- models[m]
        load(here::here("outputs", model, "Species_richness.RData"))
        sp_rich_final <- unlist(unname(sp_rich_final))
        Species_richness_all_models <- rbind(Species_richness_all_models,
                                    data.frame(Species_richness=sp_rich_final,
                                               Model=rep(Mod_type[m], length(sp_rich_final)),
                                               Nb_obs_axes=rep(n_observed_axes, length(sp_rich_final)),
                                               Seed=rep(seed, length(sp_rich_final))
                                    )
        )
        load(here::here("outputs", model, "Shannon.RData"))
        Shannon_all_models <- rbind(Shannon_all_models,
                                    data.frame(Shannon=Shannon,
                                               Model=rep(Mod_type[m], length(Shannon)),
                                               Nb_obs_axes=rep(n_observed_axes, length(Shannon)),
                                               Seed=rep(seed, length(Shannon))
                                               )
                                    )
        load(here::here("outputs", model, "Spearman.RData"))
        Spearman_all_models <- rbind(Spearman_all_models,
                                     data.frame(Spearman=c(Spearman),
                                                Model=rep(Mod_type[m], length(Spearman)),
                                                Nb_obs_axes=rep(n_observed_axes, length(Spearman)),
                                                Seed=rep(seed, length(Spearman))
                                                )
                                     )
        load(here::here("outputs", model, "Abundances.RData"))
        abund_end <- matrix(nrow=nrep, ncol=nsp)
        species_end <- list()
        
        for (r in 1:nrep) {
          abund_end[r,] <- Abundances[[r]][ngen,]
          species_end[[r]] <- unique(which(abund_end[r,]!=0))
        }
        assign(glue::glue("abund_end_{m}"), abund_end)
        assign(glue::glue("species_end_{m}"), species_end)
      } #loop on models
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
          
          vec_jaccard <- c(vec_jaccard, c(list_jaccard_matrix[[l]]))
          vec_combi_models_jaccard <- c(vec_combi_models_jaccard, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep_same_model)))
          vec_combi_rep_jaccard <- c(vec_combi_rep_jaccard, paste0(combi_rep_same_model[,1], "-", combi_rep_same_model[,2]))
          
          vec_percentage_similarity <- c(vec_percentage_similarity, c(list_percentage_similarity_matrix[[l]]))
          vec_combi_models_percentage_similarity <- c(vec_combi_models_percentage_similarity, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep_same_model)))
          vec_combi_rep_percentage_similarity <- c(vec_combi_rep_percentage_similarity, paste0(combi_rep_same_model[,1], "-", combi_rep_same_model[,2]))
          
        } #end if mod1==mod2
        else {
          list_jaccard_matrix[[l]] <- matrix(nrow=nrep, ncol=nrep)
          list_percentage_similarity_matrix[[l]] <- matrix(nrow=nrep, ncol=nrep)
          
          for (c in 1:nrow(combi_rep)){
            list_jaccard_matrix[[l]][combi_rep[c,1], combi_rep[c,2]] <- jaccard(a=get(glue::glue("species_end_{mod1}"))[[combi_rep[c,1]]], b=get(glue::glue("species_end_{mod2}"))[[combi_rep[c,2]]])
            list_percentage_similarity_matrix[[l]][combi_rep[c,1], combi_rep[c,2]] <- percentage_similarity(a=get(glue::glue("abund_end_{mod1}"))[combi_rep[c,1],], b=get(glue::glue("abund_end_{mod2}"))[combi_rep[c,2],])
          }
          
          vec_jaccard <- c(vec_jaccard, c(list_jaccard_matrix[[l]]))
          vec_combi_models_jaccard <- c(vec_combi_models_jaccard, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep)))
          vec_combi_rep_jaccard <- c(vec_combi_rep_jaccard, paste0(combi_rep[,1], "-", combi_rep[,2]))
          
          vec_percentage_similarity <- c(vec_percentage_similarity, c(list_percentage_similarity_matrix[[l]]))
          vec_combi_models_percentage_similarity <- c(vec_combi_models_percentage_similarity, rep(paste0(combi_models[l,1], "-", combi_models[l,2]), nrow(combi_rep)))
          vec_combi_rep_percentage_similarity <- c(vec_combi_rep_percentage_similarity, paste0(combi_rep[,1], "-", combi_rep[,2]))
          
        } #end if mod1!=mod2
        mean_jaccard[l] <- mean(list_jaccard_matrix[[l]])
        mean_percentage_similarity[l] <- mean(list_percentage_similarity_matrix[[l]])
      } #loop nrow combi models
      
      combi_models_jaccard <- as.data.frame(combi_models)
      combi_models_jaccard$Jaccard <- mean_jaccard
      colnames(combi_models_jaccard)[1:2] <- c("mod1", "mod2")
      
      combi_models_percentage_similarity <- as.data.frame(combi_models)
      combi_models_percentage_similarity$Percentage_similarity <- mean_percentage_similarity
      colnames(combi_models_percentage_similarity)[1:2] <- c("mod1", "mod2")
      
      
      Jaccard_all_models <- rbind(Jaccard_all_models,
                                  data.frame(Jaccard=vec_jaccard,
                                             Combi_model=vec_combi_models_jaccard,
                                             Combi_rep=vec_combi_rep_jaccard,
                                             Nb_obs_axes=rep(n_observed_axes, length(vec_jaccard)),
                                             Seed=rep(seed, length(vec_jaccard))
                                  )
      )
      
      Percentage_similarity_all_models <- rbind(Percentage_similarity_all_models,
                                  data.frame(Percentage_similarity=vec_percentage_similarity,
                                             Combi_model=vec_combi_models_percentage_similarity,
                                             Combi_rep=vec_combi_rep_percentage_similarity,
                                             Nb_obs_axes=rep(n_observed_axes, length(vec_percentage_similarity)),
                                             Seed=rep(seed, length(vec_percentage_similarity))
                                  )
      )
      
    } #loop n_observed_axes
    
  } #loop seeds
  
  Species_richness_A <- ggplot2::ggplot(data=Species_richness_all_models[Species_richness_all_models$Model=="Perf_know"&Species_richness_all_models$Nb_obs_axes==1,], ggplot2::aes(x=factor(0), y=Species_richness))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Perfect knowledge",
                  y = "Species richness")+
    ggplot2::theme(text=ggplot2::element_text(size = 16),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())
  
  Species_richness_B <- ggplot2::ggplot(data=Species_richness_all_models[Species_richness_all_models$Model=="Part_know",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Species_richness))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Partial knowledge",
                  x = "Number of axes",
                  y = "Species richness")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")+
    ggplot2::ylim(range(c(Species_richness_all_models[Species_richness_all_models$Model=="Part_know",]$Species_richness, Species_richness_all_models[Species_richness_all_models$Model=="Part_know_IV",]$Species_richness)))
  
  
  Species_richness_C <- ggplot2::ggplot(data=Species_richness_all_models[Species_richness_all_models$Model=="Part_know_IV",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Species_richness))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Species richness")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")+
    ggplot2::ylim(range(c(Species_richness_all_models[Species_richness_all_models$Model=="Part_know",]$Species_richness, Species_richness_all_models[Species_richness_all_models$Model=="Part_know_IV",]$Species_richness)))
  
  Species_richness_all_models_together <- Species_richness_all_models[-which(Species_richness_all_models$Model == "Perf_know" & Species_richness_all_models$Nb_obs_axes != 1),]
  Species_richness_all_models_together[Species_richness_all_models_together$Model=="Perf_know",]$Nb_obs_axes <- "Perfect knowledge"
  
  Species_richness_one_plot <- ggplot2::ggplot(Species_richness_all_models_together, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Species_richness))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Model), group=Model), alpha=0.6, position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Model), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Species richness")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Species_richness_one_plot, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Species_richness_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Species_richness_arranged <- ggpubr::ggarrange(Species_richness_B, Species_richness_C, Species_richness_A, nrow=1, ncol=3)
  
  ggplot2::ggsave(Species_richness_arranged, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Species_richness_boxplot_arranged.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Shannon_A <- ggplot2::ggplot(data=Shannon_all_models[Shannon_all_models$Model=="Perf_know"&Shannon_all_models$Nb_obs_axes==1,], ggplot2::aes(x=factor(0), y=Shannon))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Perfect knowledge",
                  y = "Shannon diversity index")+
    ggplot2::theme(text=ggplot2::element_text(size = 16),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())
  
  Shannon_B <- ggplot2::ggplot(data=Shannon_all_models[Shannon_all_models$Model=="Part_know",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Shannon))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Partial knowledge",
                  x = "Number of axes",
                  y = "Shannon diversity index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")+
    ggplot2::ylim(range(c(Shannon_all_models[Shannon_all_models$Model=="Part_know",]$Shannon, Shannon_all_models[Shannon_all_models$Model=="Part_know_IV",]$Shannon)))
  
  
  Shannon_C <- ggplot2::ggplot(data=Shannon_all_models[Shannon_all_models$Model=="Part_know_IV",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Shannon))+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::labs(title = "Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Shannon diversity index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")+
    ggplot2::ylim(range(c(Shannon_all_models[Shannon_all_models$Model=="Part_know",]$Shannon, Shannon_all_models[Shannon_all_models$Model=="Part_know_IV",]$Shannon)))
  
  Shannon_all_models_together <- Shannon_all_models[-which(Shannon_all_models$Model == "Perf_know" & Shannon_all_models$Nb_obs_axes != 1),]
  Shannon_all_models_together[Shannon_all_models_together$Model=="Perf_know",]$Nb_obs_axes <- "Perfect knowledge"
  
  Shannon_one_plot <- ggplot2::ggplot(Shannon_all_models_together, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Shannon))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Model), group=Model), alpha=0.6, position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Model), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Shannon diversity index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Shannon_one_plot, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Shannon_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Shannon_arranged <- ggpubr::ggarrange(Shannon_B, Shannon_C, Shannon_A, nrow=1, ncol=3)
  
  ggplot2::ggsave(Shannon_arranged, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Shannon_boxplot_arranged.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Richness_Shannon <- ggpubr::ggarrange(Species_richness_one_plot, Shannon_one_plot, nrow=2, ncol=1, labels = c("A", "B"))
  ggplot2::ggsave(Richness_Shannon, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Richness_Shannon.png"),
                  width=fig_width*2, height=fig_width*2, units="cm", dpi=300)
  
  
  Spearman_A <- ggplot2::ggplot(data=Spearman_all_models[Spearman_all_models$Model=="Perf_know"&Spearman_all_models$Nb_obs_axes==1,], ggplot2::aes(x=factor(0), y=Spearman))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.5)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge",
                  y = "Spearman pairwise correlation index")+
    ggplot2::theme(text=ggplot2::element_text(size = 16),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())
  
  Spearman_B <- ggplot2::ggplot(data=Spearman_all_models[Spearman_all_models$Model=="Part_know",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Spearman))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.5)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge",
                  x = "Number of axes",
                  y = "Spearman pairwise correlation index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")+
    ggplot2::ylim(range(c(Spearman_all_models[Spearman_all_models$Model=="Part_know",]$Spearman, Spearman_all_models[Spearman_all_models$Model=="Part_know_IV",]$Spearman)))
  
  
  Spearman_C <- ggplot2::ggplot(data=Spearman_all_models[Spearman_all_models$Model=="Part_know_IV",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Spearman))+
    ggplot2::geom_boxplot()+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.5)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Spearman pairwise correlation index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())+
    ggplot2::ylim(range(c(Spearman_all_models[Spearman_all_models$Model=="Part_know",]$Spearman, Spearman_all_models[Spearman_all_models$Model=="Part_know_IV",]$Spearman)))
  
  Spearman_all_models_together <- Spearman_all_models[-which(Spearman_all_models$Model == "Perf_know" & Spearman_all_models$Nb_obs_axes != 1),]
  Spearman_all_models_together[Spearman_all_models_together$Model=="Perf_know",]$Nb_obs_axes <- "Perfect knowledge"
  
  Spearman_one_plot <- ggplot2::ggplot(Spearman_all_models_together, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Spearman))+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Model))+
    ggplot2::scale_colour_manual(values=c("#B2BEB5","black", "#54626F"))+
    ggnewscale::new_scale("colour")+
    ggbeeswarm::geom_beeswarm(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Model), group=Model), alpha=0.6, dodge.width = 0.755)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Spearman pairwise correlation index")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Spearman_one_plot, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Spearman_boxplot.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Spearman_arranged <- ggpubr::ggarrange(Spearman_B, Spearman_C, Spearman_A)
  
  ggplot2::ggsave(Spearman_arranged, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Spearman_boxplot_arranged.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)

  
  Jaccard_all_models_together <- Jaccard_all_models[-which(Jaccard_all_models$Combi_model == "1-1" & Jaccard_all_models$Nb_obs_axes != 1),]
  Jaccard_all_models_together[Jaccard_all_models_together$Combi_model=="1-1",]$Nb_obs_axes <- "Perfect knowledge"
  
  Jaccard_all_models_together_within <- Jaccard_all_models_together[Jaccard_all_models_together$Combi_model%in%c("1-1", "2-2", "3-3"),]
  
  Jaccard_one_plot_within <- ggplot2::ggplot(Jaccard_all_models_together_within, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Combi_model), group=Combi_model),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_shape_manual(values=c(15, 16, 17))+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Combi_model))+
    ggplot2::scale_colour_manual(values=c("#54626F", "#B2BEB5", "black"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Jaccard similarity index of the composition of the community")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Jaccard_one_plot_within, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Jaccard_boxplot_within.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  set.seed(0)
  Jaccard_A <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                               ggplot2::aes(x=factor(0), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text=ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  Jaccard_B <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="2-2",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::geom_jitter(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Jaccard, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Jaccard), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge",
                  x = "Number of axes",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  Jaccard_C <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="3-3",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::geom_jitter(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Jaccard, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Jaccard), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  Jaccard_D <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-2",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha=0.6)+
    ggplot2::geom_jitter(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Jaccard, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Jaccard), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge - Partial knowledge",
                  x = "Number of axes",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  Jaccard_E <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-3",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Jaccard, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Jaccard), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge - Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  Jaccard_F <- ggplot2::ggplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="2-3",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Jaccard))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Jaccard, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Jaccard_all_models[Jaccard_all_models$Combi_model=="1-1"&Jaccard_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Jaccard), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge - Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Jaccard similarity index")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Jaccard_all_models_together$Jaccard))
  
  ggplot2::ggsave(Jaccard_F, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Jaccard_boxplot_between_Partials.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Jaccard_arranged_all <- ggpubr::ggarrange(Jaccard_A,
                                            Jaccard_B,
                                            Jaccard_C,
                                            Jaccard_D,
                                            Jaccard_E,
                                            Jaccard_F)
  
  ggplot2::ggsave(Jaccard_arranged_all, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Jaccard_boxplot_arranged_all.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Jaccard_arranged_between <- ggpubr::ggarrange(Jaccard_D,
                                                Jaccard_E,
                                                nrow=1,
                                                ncol=2)
  
  ggplot2::ggsave(Jaccard_arranged_between, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Jaccard_boxplot_arranged_between.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Percentage_similarity_all_models_together <- Percentage_similarity_all_models[-which(Percentage_similarity_all_models$Combi_model == "1-1" & Percentage_similarity_all_models$Nb_obs_axes != 1),]
  Percentage_similarity_all_models_together[Percentage_similarity_all_models_together$Combi_model=="1-1",]$Nb_obs_axes <- "Perfect knowledge"
  
  Percentage_similarity_all_models_together_within <- Percentage_similarity_all_models_together[Percentage_similarity_all_models_together$Combi_model%in%c("1-1", "2-2", "3-3"),]
  Percentage_similarity_all_models_together_between <- Percentage_similarity_all_models_together[!(Percentage_similarity_all_models_together$Combi_model%in%c("1-1", "2-2", "3-3")),]
  
  Percentage_similarity_one_plot_within <- ggplot2::ggplot(Percentage_similarity_all_models_together_within, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Combi_model), group=Combi_model),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_shape_manual(values=c(15, 16, 17))+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Combi_model), alpha = 0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Percentage similarity of the final species abundances")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Percentage_similarity_one_plot_within, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Percentage_similarity_boxplot_within.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Percentage_similarity_one_plot_between <- ggplot2::ggplot(Percentage_similarity_all_models_together_between[Percentage_similarity_all_models_together_between$Combi_model!="2-3",], ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Combi_model), group=Combi_model),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed), shape=as.factor(Combi_model)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_shape_manual(values=c(15, 16, 17))+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=as.factor(Combi_model)), alpha = 0.6)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Combi_model)), alpha = 0.6)+
    ggplot2::scale_colour_manual(values=c("#4C5866", "#808080", "black"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Percentage similarity of the final species abundances")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(Percentage_similarity_one_plot_between, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Percentage_similarity_boxplot_between.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Percentage_similarity_A <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                                             ggplot2::aes(x=factor(0), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text=ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  Percentage_similarity_B <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="2-2",],
                                             ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge",
                  x = "Number of axes",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  Percentage_similarity_C <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="3-3",],
                                             ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  Percentage_similarity_D <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-2",],
                                             ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge - \n Partial knowledge",
                  x = "Number of axes",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  Percentage_similarity_E <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-3",],
                                             ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Perfect knowledge - \n Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none",
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank())+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  Percentage_similarity_F <- ggplot2::ggplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="2-3",],
                                             ggplot2::aes(x=as.factor(Nb_obs_axes), y=Percentage_similarity))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(alpha = 0.6)+
    ggplot2::geom_jitter(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                         ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
    ggplot2::geom_boxplot(data=Percentage_similarity_all_models[Percentage_similarity_all_models$Combi_model=="1-1"&Percentage_similarity_all_models$Nb_obs_axes==1,],
                          ggplot2::aes(x="Perfect \n knowledge", y=Percentage_similarity), alpha = 0.6)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title = "Partial knowledge - \n Partial knowledge + IV",
                  x = "Number of axes",
                  y = "Percentage similarity")+
    ggplot2::theme(plot.title = ggplot2::element_text(size=14),
                   text = ggplot2::element_text(size = 14),
                   legend.position = "none")+
    ggplot2::ylim(range(Percentage_similarity_all_models_together$Percentage_similarity))
  
  ggplot2::ggsave(Percentage_similarity_F, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Percentage_similarity_boxplot_between_Partials.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  
  Percentage_similarity_arranged_all <- ggpubr::ggarrange(Percentage_similarity_A,
                                                      Percentage_similarity_B,
                                                      Percentage_similarity_C,
                                                      Percentage_similarity_D,
                                                      Percentage_similarity_E,
                                                      Percentage_similarity_F)
  
  Percentage_similarity_arranged_all <- ggpubr::annotate_figure(Percentage_similarity_arranged_all,
                  bottom = ggpubr::text_grob("Number of axes", face = "bold", size = 14),
                  left = ggpubr::text_grob("Percentage similarity", face = "bold", size = 14, rot=90))
  
  ggplot2::ggsave(Percentage_similarity_arranged_all, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Percentage_similarity_boxplot_arranged_all.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  Percentage_similarity_arranged_between <- ggpubr::ggarrange(Percentage_similarity_D,
                                                      Percentage_similarity_E,
                                                      Percentage_similarity_F,
                                                      nrow=1,
                                                      ncol=3)
  
  Percentage_similarity_arranged_between <- ggpubr::annotate_figure(Percentage_similarity_arranged_between,
                          bottom = ggpubr::text_grob("Number of axes", face = "bold", size = 14),
                          left = ggpubr::text_grob("Percentage similarity", face = "bold", size = 14, rot=90))
  
  ggplot2::ggsave(Percentage_similarity_arranged_between, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Percentage_similarity_boxplot_arranged_between.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  for (seed in Seeds) {
    
    colourCount = nsp
    getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    
    png(file=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), glue::glue("Final_communities_{seed}.png")),
        width=fig_width, height=fig_width, units="cm", res=300)
    
    par(mfrow=c(2,2), bty = "n")
    
    load(here::here("outputs", glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}"), glue::glue("community_end.RData")))
    raster::plot(raster::raster(community_end[[1]]), main="Species - Perfect knowledge", zlim=c(0, nsp),
                 col=c("black", getPalette(colourCount)), legend=FALSE)
    load(here::here("outputs", glue::glue("Part_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_1_obs_seed_{seed}"), glue::glue("community_end.RData")))
    raster::plot(raster::raster(community_end[[1]]), main="Species - Partial knowledge", 
         zlim=c(0, nsp), col=c("black", getPalette(colourCount)), legend=FALSE, yaxt = "n", xaxt = "n")
    load(here::here("outputs", glue::glue("Part_know_IV_{start}_{mort}_{nb_seeds}_{n_axes}_axes_1_obs_seed_{seed}"), glue::glue("community_end.RData")))
    raster::plot(raster::raster(community_end[[1]]), main="Species - Partial knowledge + IV", 
         zlim=c(0, nsp), col=c("black", getPalette(colourCount)), legend=FALSE, yaxt = "n", xaxt = "n")
    load(here::here("outputs", glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}"), "env_entrelac.RData"))
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), main="Environment summary", col=viridisLite::viridis(255^3), yaxt = "n", xaxt = "n")
    
    dev.off()
  }
  
  Correlation_env_sp <- data.frame(Model = c(rep("Perf_know", length(Seeds)*nrep), rep(c("Part_know", "Part_know_IV"), each=length(nb_obs_axes)*length(Seeds)*nrep)),
                                   Seed = c(rep(Seeds, each=nrep), rep(rep(Seeds, each=length(nb_obs_axes)*nrep),2)),
                                   Nb_obs_axes = c(rep("Perfect knowledge", length(Seeds)*nrep), rep(rep(rep(nb_obs_axes, each=nrep), length(Seeds)), 2)),
                                   Rep = c(rep(1:nrep, length(Seeds)), rep(rep(1:nrep, length(Seeds)*length(nb_obs_axes)),2)),
                                   Correlation = numeric((2*length(nb_obs_axes)+1)*length(Seeds)*nrep))
  
  for(k in 1:nrow(Correlation_env_sp)){
    if(Correlation_env_sp$Model[k]=="Perf_know"){
      model <- glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{Correlation_env_sp$Seed[k]}")
    }else{
      model <- glue::glue("{Correlation_env_sp$Model[k]}_{start}_{mort}_{nb_seeds}_{n_axes}_axes_{Correlation_env_sp$Nb_obs_axes[k]}_obs_seed_{Correlation_env_sp$Seed[k]}")
    }
    
    load(here::here("outputs", model, "semivar_multidim.RData"))
    
    semivar_multidim <- semivar_multidim[[Correlation_env_sp$Rep[k]]]
    
    semivar_multidim<-semivar_multidim%>%
      filter(Sample_size>500)
    
    m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)
    Correlation_env_sp$Correlation[k] <- round(sqrt(summary(m)$r.squared), digits = 2)
  }
  
  Summary_correlation_env_sp<-Correlation_env_sp%>%
    dplyr::group_by(Model, Nb_obs_axes)%>%
    dplyr::mutate(Mean_corr=mean(Correlation, na.rm=TRUE), Sd=sd(Correlation, na.rm=TRUE))%>%
    dplyr::slice(1)%>%
    dplyr::ungroup()%>%
    dplyr::select(-Correlation, -Seed, -Rep)
  
  save(Summary_correlation_env_sp, file=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Mean_correlation_env_sp.RData"))
  
  p <- ggplot2::ggplot(Correlation_env_sp, ggplot2::aes(x=as.factor(Nb_obs_axes), y=Correlation))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Model), group=Model), alpha=0.6, position=ggplot2::position_jitterdodge(jitter.width=0.5))+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Model), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = "Number of observed axes",
                  y = "Correlation between environment and species semivariance")+
    ggplot2::theme(text = ggplot2::element_text(size = 16), legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Corr_env_sp.png"),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  #Correlation between environmental variables
  Corr_env_all_configs <- data.frame()
    
  for (seed in Seeds) {
    model <- glue::glue("Perf_know_{start}_{mort}_{nb_seeds}_{n_axes}_axes_seed_{seed}")
    load(here::here("outputs", model, "Corr_env.RData"))
    Corr_env_all_configs <- rbind(Corr_env_all_configs, Corr_env)
  }
  
  Corr_env_all_configs$Seed <- rep(Seeds, each=(n_axes^2-n_axes)/2)
  save(Corr_env_all_configs, file=here::here("outputs", glue::glue("Comparison_{mort}_{nb_seeds}"), "Corr_env_all_configs.RData"))
}