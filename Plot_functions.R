colourCount = nsp
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

plot_environment <- function(model, fig_width, n_axis, env){
  
  png(file=here::here("outputs", model, "environment.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
  
  par(mfrow=c(4,ceiling(n_axis/4)), bty = "n")
  
  for(k in 1:n_axis){
    plot(raster::raster(env[[k]]), main=glue::glue("Environment var{k}"), col=topo.colors(255), cex.main=1.2)
  }
  if(grepl(pattern="Perf_know", model)==TRUE){
    #RGB environment
    env_stack <- stack(raster::raster(env[[1]]*255), raster::raster(env[[2]]*255), raster::raster(env[[3]]*255))
    crs(env_stack) <- "+proj=utm +zone=1"
    class_site <- entrelac(
      raster::values(env_stack@layers[[1]]),
      raster::values(env_stack@layers[[2]]),
      raster::values(env_stack@layers[[3]])
    )
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), col=viridisLite::viridis(255^3))
    #plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
  }
  
  dev.off()
}

plot_hab_freq <- function(n_axis, model, fig_width, env){
  
  for(k in 1:n_axis){
    png(file=here::here("outputs", model, glue::glue("hab_freq_{k}.png")),
        width=fig_width, height=fig_width, units="cm", res=300)
    hist(env[[k]], main="", xlab=glue::glue("Environment var{k}")) 
    dev.off()
  }
}

plot_species_optima <- function(model, fig_width, niche_optimum){
  png(file=here::here("outputs", model, "species_niche.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mar=c(1,1,2,2))
  scatter3D(niche_optimum[,1], niche_optimum[,2], niche_optimum[,3],
            pch=16, 
            colvar=1:nsp, col=getPalette(colourCount),
            bty = "f", main ="Three-dimensional species optima", phi=0,
            xlim=c(0,1), ylim=c(0,1), zlim=c(0,1))
  dev.off()
}

plot_histogram_mortality <- function(model, fig_width, mortality_E_Sp){
  png(file=here::here("outputs", model, "hist_mortality.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  hist(mortality_E_Sp)
  dev.off()
}

plot_function_mort_proba <- function(model, fig_width, perf_E_Sp, mortality_E_Sp){
  png(file=here::here("outputs", model, "function_mort_proba.png"),
      width=fig_width, height=fig_width/1.5, units="cm", res=300)
  plot(x=c(perf_E_Sp),
       y=c(mortality_E_Sp),
       xlab="Performance",
       ylab="Mortality probability",
       main="Mortality probability as a function of performance",
       cex.lab=1.5,
       cex.main=1.5)
  dev.off()
}

plot_species_habitat_freq <- function(model, fig_width, sp_hab_freq){
  png(file=here::here("outputs", model, "species_habitat_freq.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(cex.lab=1.5)
  plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency", col=getPalette(colourCount))
  dev.off()
}

plot_community_start <- function(model, fig_width, community, nsp){
  png(file=here::here("outputs", model, "community_start.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  # plot(raster::raster(community), main="Species - Start", zlim=c(0, nsp),
  #      col=c("black",  viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  plot(raster::raster(community), main="Species - Start", zlim=c(0, nsp),
       col=c("black",  getPalette(colourCount)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_mortality_events <- function(model, fig_width, community, mortality, nsp){
  png(file=here::here("outputs", model, "mortality_events.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  community_mortality <- community
  community_mortality[mortality==1] <- nsp+1
  plot(raster::raster(community_mortality), main="Species (colours) and vacant sites \n (black = old, white = new)", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), "white", legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_community_end <- function(model, fig_width, community, nsp){
png(file=here::here("outputs", model, "community_end.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(bty = "n")
# plot(raster::raster(community), main=glue::glue("Species - End (ngen={ngen})"),
#      zlim=c(0, nsp), col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
plot(raster::raster(community), main=glue::glue("Species - End (ngen={ngen})"), 
     zlim=c(0, nsp), col=c("black", getPalette(colourCount)), legend=TRUE, cex.main=2, cex.axis=1.5)
dev.off()
}

plot_species_richness <- function(nrep, sp_rich, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  sp_rich_long <- sp_rich %>%
    dplyr::mutate(gen=1:(ngen)) %>%
    tidyr::pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="sp_rich")
  #One curve per repetition
  p <- ggplot2::ggplot(data=sp_rich_long, ggplot2::aes(x=gen, y=sp_rich, col=rep)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_viridis_d()+
    ggplot2::xlab("Generations") + 
    ggplot2::ylab("Species richness")+
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 20))+
    ggplot2::ylim(0,nsp)
  ggplot2::ggsave(p, filename=here::here("outputs", model, "species_richness_with_time.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  #Mean and 95% interval
  if(nrep>1){
    sp_rich_mean <- data.frame(Gen = 1:nrow(sp_rich), Mean=rowMeans(sp_rich), Low=apply(sp_rich, 1, quantile, probs=0.025), High=apply(sp_rich, 1, quantile, probs=0.975))
    
    p <- ggplot2::ggplot(data=sp_rich_mean, ggplot2::aes(x=Gen, y=Mean)) +
      ggplot2::geom_line(col="#008071") +
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
    ggplot2::ggsave(p, filename=here::here("outputs", model, "species_richness_with_time_mean.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
  
  # Log_10 plot
  
  if(nrep==1){
    sp_rich_log <- data.frame(Sp_rich=sp_rich, Gen=c(1:ngen))
    p <- ggplot2::ggplot(data=sp_rich_log, ggplot2::aes(x=as.numeric(Gen), y=Sp_rich)) +
      ggplot2::geom_line(col="#008071")+
      scale_x_continuous(trans='log10')+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, "species_richness_log_10.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    }
  
  else{
    p <- ggplot2::ggplot(data=sp_rich_mean, ggplot2::aes(x=as.numeric(Gen), y=Mean)) +
      ggplot2::geom_line(col="#008071")+
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      scale_x_continuous(trans='log10')+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, "species_richness_log_10.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    }
}

plot_mean_rank_hab_freq <- function(sp_mean_rank, sp_hab_freq, model, fig_width){
  mean_rank_hab_freq <- data.frame(cbind(Species=1:nsp, Hab_freq=sp_hab_freq, Rank=sp_mean_rank))
  
  p <- ggplot2::ggplot(data=mean_rank_hab_freq, ggplot2::aes(x=Hab_freq, y=Rank, colour=as.factor(Species))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="#008071", fill="#69b3a2", se=TRUE) +
    ggplot2::xlab("Species suitable habitat frequency") +
    ggplot2::ylab("Species mean rank (higher rank = lower abundance)") +
    ggplot2::theme(axis.title=ggplot2::element_text(size=16))+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "mean_rank-habitat_freq.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_env_filt<- function(nrep, env_filt, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  env_filt_long <- env_filt %>%
    dplyr::mutate(gen=1:(ngen)) %>%
    tidyr::pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="env_filt")
  
  p <- ggplot2::ggplot(data=env_filt_long, ggplot2::aes(x=gen, y=env_filt, col=rep)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title="Environmental filtering") +
    ggplot2::xlab("Generations") + 
    ggplot2::ylab("Mean env-species perf difference")+
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 20))
  ggplot2::ggsave(p, filename=here::here("outputs", model, "environmental_filtering.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  #Mean and 95% interval
  if(nrep>1){
    env_filt_mean <- data.frame(Gen = 1:nrow(env_filt), Mean=rowMeans(env_filt), Low=apply(env_filt, 1, quantile, probs=0.025), High=apply(env_filt, 1, quantile, probs=0.975))
    
    p <- ggplot2::ggplot(data=env_filt_mean, ggplot2::aes(x=Gen, y=Mean)) +
      ggplot2::geom_line(col="#008071") +
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, "environmental_filtering_mean.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}
  
plot_env_species<- function(model, fig_width, community_start, community_end, env){
  png(file=here::here("outputs", model, "spatial_comp_env_sp.png"), 
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  plot(raster::raster(community_start), main="Species - Start", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), legend=FALSE)
  plot(raster::raster(community_end), main="Species - End", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), legend=FALSE)
  if(n_axis==3){
    env_stack <- raster::stack(raster::raster(env[[1]]*255), raster::raster(env[[2]]*255), raster::raster(env[[3]]*255))
    raster::crs(env_stack) <- "+proj=utm +zone=1"
    save(env_stack, file = here::here("outputs", model, "env_stack.RData"))
    raster::plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
    dev.off()
  }
  else{
    pca_env <- prcomp(sites, scale = TRUE)
    pca_env$rotation <- -pca_env$rotation
    pca_env$x <- -pca_env$x
    env_stack <- raster::stack(
      raster::raster(matrix(range_0_255(pca_env$x[,1]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE)),
      raster::raster(matrix(range_0_255(pca_env$x[,2]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE)),
      raster::raster(matrix(range_0_255(pca_env$x[,3]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE))
      )
    raster::crs(env_stack) <- "+proj=utm +zone=1"
    save(env_stack, file = here::here("outputs", model, "env_stack.RData"))
    raster::plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
    dev.off()
    var_explained <- pca_env$sdev^2 / sum(pca_env$sdev^2)
    p <- ggplot2::ggplot(data.frame(Prin_comp=c(1:length(var_explained)),
                                    Prop=cumsum(var_explained)),
                         ggplot2::aes(x=Prin_comp, y=Prop))+
      ggplot2::geom_col()+
      ggplot2::labs(x="Principal component",
                    y="Cumulative proportion of explained variance")
    ggplot2::ggsave(p, filename=here::here("outputs", model, "pca_prop_var.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_theta_community<-function(theta_comm, ngen, nrep, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  theta_comm_long <- theta_comm %>%
    dplyr::mutate(gen=1:(ngen)) %>%
    tidyr::pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="theta_comm")
  p <- ggplot2::ggplot(data=theta_comm_long, ggplot2::aes(x=gen, y=theta_comm, col=rep)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(title="Environmental filtering") +
    ggplot2::xlab("Generations") + 
    ggplot2::ylab("Mean mortality rate in the community")+
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 20))
  ggplot2::ggsave(p, filename=here::here("outputs", model, "mortality_rate_community.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  #Mean and 95% interval
  if(nrep>1){
    theta_comm_mean <- data.frame(Gen = 1:nrow(theta_comm), Mean=rowMeans(theta_comm), Low=apply(theta_comm, 1, quantile, probs=0.025), High=apply(theta_comm, 1, quantile, probs=0.975))
    
    p <- ggplot2::ggplot(data=theta_comm_mean, ggplot2::aes(x=Gen, y=Mean)) +
      ggplot2::geom_line(col="#008071") +
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, "mortality_rate_community_mean.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_spatial_autocorr <- function(community_end, sites, niche_width, env_stack, model, fig_width){
  # Species autocorrelation
  sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
  names(sp_XY) <- c("x", "y", "sp")
  vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
  #plot(vario_sp[["bins.lim"]][-length(vario_sp[["bins.lim"]])], vario_sp[["v"]])
  
  # Environment autocorrelation
  # 3D voxel for each site
  if(n_axis==3){
    if(randomOptSp==FALSE){
      x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
      y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
      z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
      n_niche <- 1/niche_width
      # This is done to avoid having the same class for different combinations (the function is not continuous)
      class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
      vario_env <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
    }
    else{
      x_site <- floor(sites$V1_env*10)
      y_site <- floor(sites$V2_env*10)
      z_site <- floor(sites$V3_env*10)
      class_site <- z_site*10^2+y_site*10+x_site
      vario_env <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
    }
  }
  else{
    class_site <- entrelac(
      raster::values(env_stack@layers[[1]]),
      raster::values(env_stack@layers[[2]]),
      raster::values(env_stack@layers[[3]])
      )
    
    png(file=here::here("outputs", model, "pca_rgb.png"),
        width=fig_width, height=fig_width*0.8, units="cm", res=300)
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), col=viridisLite::viridis(255^3))
    dev.off()
    
    #vario_env <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
    
    semivar_multidim_env <- semivar_mutlidim(sites, n_axis, sp_XY, vario_sp)
  }
  # Plot with correlation
  png(file=here::here("outputs", model, "sp_autocorrelation.png"),
      width=fig_width, height=fig_width*0.8, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  plot(vario_sp, main="Species - End")
  #plot(vario_env, main="Environment")
  plot(vario_sp[["bins.lim"]][-length(vario_sp[["bins.lim"]])], semivar_multidim_env,
       main="Environment",
       xlab="distance",
       ylab="semivariance")
  # plot(vario_env$v, vario_sp$v,
  #      main = "Regression",
  #      xlab="Semivariance for environment",
  #      ylab="Semivariance for species")
  # m <- lm(vario_sp$v ~ vario_env$v)
  # /!\ a=intercept and b=slope
  # abline(a=0, b=coef(m), col="red")
  plot(semivar_multidim_env, vario_sp[["v"]],
       main = "Regression",
       xlab="Semivariance for environment",
       ylab="Semivariance for species")
  m <- lm(vario_sp$v ~ semivar_multidim_env)
  abline(a=as.numeric(coef(m)[1]), b=as.numeric(coef(m)[2]), col="#008071")
  dev.off()
}

plot_species_niche <- function(seed, df_perf, model, fig_width){
  # Select 8 species at random
  #set.seed(seed)
  #sp_sel <- sample(unique(df_perf$Species), 9, replace=FALSE)
  #df_sp_sel <- df_perf %>% filter(Species %in% sp_sel)
  
  df_perf_sp_niche <- df_perf
  load(here::here("outputs", model, "niche_optimum.RData"))
  df_perf_sp_niche$optimum <- rep(niche_optimum$X1, nsite)
  df_perf_sp_niche$dist <- sqrt((df_perf_sp_niche$Env_1-df_perf_sp_niche$optimum)^2)
  df_perf_sp_niche$dist <- df_perf_sp_niche$dist - mean(df_perf_sp_niche$dist) / sd(df_perf_sp_niche$dist)

  p <- ggplot2::ggplot(data=df_perf_sp_niche, ggplot2::aes(x=Env_1, y=Perf)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE, col="#008071") +
    #ggplot2::stat_function(fun=function(Env_1) Env_1-optimum)+
    ggplot2::geom_line(ggplot2::aes(x=Env_1, y=-dist), col="#088000")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=optimum), col="#80002D")+
    ggplot2::facet_wrap(ggplot2::vars(Species), nrow=3) +
    ggplot2::xlab("Environment (first axis)") +
    ggplot2::ylab("Performance")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "infering_species_niche.png"),
         width=fig_width*2, height=fig_width, units="cm", dpi=300)
}

plot_IV <- function(V_intra, model, fig_width){
  p <- ggplot2::ggplot(data=V_intra, ggplot2::aes(x=Species, y=V, fill=as.factor(Species))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, size=6),
                   text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size=14),
                   legend.text =  ggplot2::element_text(size=10),
                   legend.key.size = unit(0.3, 'cm'),
                   plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                   legend.position=c(1.1, 0.4)) +
    ggplot2::ylab("Intraspecific variance")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "intraspecific_variance.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_relationship_IV_inferred_optima <- function(V_intra,
                                                 n_observed_axis,
                                                 Inferred_species_parameters,
                                                 community_end,
                                                 sp_hab_freq,
                                                 model,
                                                 fig_width){

  Freq_sp_model_IV <- as.data.frame(table(unlist(community_end)))
  
  IV_opt_sp_estim <- data.frame(Species = 1:nrow(V_intra),
                                Opt_estim=numeric(nrow(V_intra)),
                                Perf_opt=numeric(nrow(V_intra)),
                                IV=V_intra$V,
                                Frequency=numeric(nrow(V_intra)),
                                Presence=numeric(nrow(V_intra)),
                                Hab_freq=as.data.frame(sp_hab_freq)$Freq,
                                Winner=character(nrow(V_intra)))
  
  for(sp in 1:nrow(V_intra)){
    for(axis in 1:n_observed_axis){
      IV_opt_sp_estim$Opt_estim[sp] <- (-Inferred_species_parameters[sp,axis+1])/(2*Inferred_species_parameters[sp,2*axis+1])
    }
    if(sp%in%unlist(community_end)){
      IV_opt_sp_estim$Frequency[sp] <- Freq_sp_model_IV[Freq_sp_model_IV$Var1==sp,]$Freq
      IV_opt_sp_estim$Presence[sp] <- 1
    }else{
      IV_opt_sp_estim$Frequency[sp]<-0
      IV_opt_sp_estim$Presence[sp] <- 0
    }
    if(IV_opt_sp_estim$Hab_freq[sp]!=0){
      IV_opt_sp_estim$Winner[sp] <- "theoretical"
    }else{IV_opt_sp_estim$Winner[sp] <- "not present"}
  }
  
  if(length(IV_opt_sp_estim[IV_opt_sp_estim$Opt_estim>1,]$Opt_estim)!=0){
    IV_opt_sp_estim[IV_opt_sp_estim$Opt_estim>1,]$Opt_estim <- 1
  }
  if(length(IV_opt_sp_estim[IV_opt_sp_estim$Opt_estim<0,]$Opt_estim)!=0){
    IV_opt_sp_estim[IV_opt_sp_estim$Opt_estim<0,]$Opt_estim <- 0
  }
  
  for(sp in 1:nrow(V_intra)){
    for(axis in 1:n_observed_axis){
      IV_opt_sp_estim$Perf_opt[sp] <- Inferred_species_parameters[sp,axis] + Inferred_species_parameters[sp,axis+1]*IV_opt_sp_estim$Opt_estim[sp] + Inferred_species_parameters[sp,2*axis+1]*(IV_opt_sp_estim$Opt_estim[sp])^2
    }
  }
  
  IV_opt_sp_estim[which(IV_opt_sp_estim$Winner!="theoretical"&IV_opt_sp_estim$Presence==1),]$Winner<-"IV"
  
  IV_opt_sp_estim$Freq_IV_winner <- c(0)
  IV_opt_sp_estim[IV_opt_sp_estim$Winner=="IV",]$Freq_IV_winner <- IV_opt_sp_estim[IV_opt_sp_estim$Winner=="IV",]$Frequency
  
  ggplot2::ggplot(IV_opt_sp_estim, aes(Opt_estim, IV))+
    ggplot2::geom_point(aes(shape=as.factor(Winner), colour=Freq_IV_winner), size=3)+
    ggplot2::scale_colour_viridis_c(option = "plasma")+
    ggplot2::labs(x="Estimated optimum",
                  y="Estimated IV",
                  shape="Species presence",
                  colour="Frequency of presence \n in final community \n (only species maintained thanks to IV)")+
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size=14),
                   legend.text =  ggplot2::element_text(size=14))
  
  ggplot2::ggplot(IV_opt_sp_estim, aes(Opt_estim, IV))+
    ggplot2::geom_point(aes(shape=as.factor(Winner), colour=Perf_opt), size=3)+
    ggplot2::scale_colour_viridis_c(option = "plasma")+
    ggplot2::labs(x="Estimated optimum",
                  y="Estimated IV",
                  shape="Species presence",
                  colour="Performance at optimum")+
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size=16),
                   legend.text =  ggplot2::element_text(size=14))
}

plot_optima_real_estim <- function(nsp, n_observed_axis, niche_optimum, Inferred_species_parameters, model, fig_width){
  
  opt_sp_estim_vs_real <- data.frame(Species = 1:nsp, Real=numeric(nsp), Estim=numeric(nsp))
  
  for(sp in 1:nsp){

    for(axis in 1:n_observed_axis){
      opt_sp_estim_vs_real$Real[sp] <- niche_optimum[sp, axis]
      opt_sp_estim_vs_real$Estim[sp] <- (-Inferred_species_parameters[sp,axis+1])/(2*Inferred_species_parameters[sp,2*axis+1])
    }

  }
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Real, y=Estim, colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::labs(x="Real species optima on axis X1", y="Estimated species optima")

  ggplot2::ggsave(p, filename=here::here("outputs", model, "optima_real_estim.png"),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  opt_sp_estim_vs_real$Horizontal <- rep(0, nrow(opt_sp_estim_vs_real))
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Estim, y=Horizontal,  colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "optima_estim.png"),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  if(length(opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim>1,]$Estim)!=0){
    opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim>1,]$Estim <- 1}
  if(length(opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim<0,]$Estim)!=0){
    opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim<0,]$Estim <- 0
  }
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Estim, y=Horizontal,  colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "optima_estim_0_1.png"),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_inferred_perf_environment <- function(E_seq, Mat_perf_inferred, nsp, model, fig_width){

  if(n_observed_axis>1){
    
    for(k in 1:n_observed_axis){
  
      Mat_perf_inferred_plot <- as.data.frame(cbind(E_seq[,k], Mat_perf_inferred))
      
      colnames(Mat_perf_inferred_plot) <- c("E", paste0("Sp",1:nsp))
      
      Mat_perf_inferred_plot <- Mat_perf_inferred_plot %>%
        tidyr::pivot_longer(cols=2:(nsp+1), names_to="Sp",
                            names_prefix="Sp", values_to="Perf")
      
      Mat_perf_inferred_plot$Sp <- as.numeric(Mat_perf_inferred_plot$Sp)
      
      p <- ggplot2::ggplot(data=Mat_perf_inferred_plot, ggplot2::aes(x=E, y=Perf, col=as.factor(Sp))) +
        ggplot2::geom_line() +
        #ggplot2::scale_colour_viridis_d()+
        ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
        ggplot2::xlab(glue::glue("Environment axis {k}")) + 
        ggplot2::ylab("Species performance")+
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       legend.title = ggplot2::element_text(size=14),
                       legend.text =  ggplot2::element_text(size=10),
                       legend.key.size = unit(0.3, 'cm'),
                       plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                       legend.position=c(1.1, 0.4))
      
      ggplot2::ggsave(p, filename=here::here("outputs", model, glue::glue("inferred_perf_environment_{k}.png")),
                      width=fig_width, height=fig_width/2, units="cm", dpi=300)
      
    }
  } else {
      
      Mat_perf_inferred_plot <- as.data.frame(cbind(E_seq, Mat_perf_inferred))
      
      colnames(Mat_perf_inferred_plot) <- c("E", paste0("Sp",1:nsp))
      
      Mat_perf_inferred_plot <- Mat_perf_inferred_plot %>%
        tidyr::pivot_longer(cols=2:(nsp+1), names_to="Sp",
                            names_prefix="Sp", values_to="Perf")
      
      Mat_perf_inferred_plot$Sp <- as.numeric(Mat_perf_inferred_plot$Sp)
      
      p <- ggplot2::ggplot(data=Mat_perf_inferred_plot, ggplot2::aes(x=E, y=Perf, col=as.factor(Sp))) +
        ggplot2::geom_line() +
        #ggplot2::scale_colour_viridis_d()+
        ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
        ggplot2::xlab("Environment axis 1") + 
        ggplot2::ylab("Species performance")+
        ggplot2::theme(text = ggplot2::element_text(size = 20),
                       legend.title = ggplot2::element_text(size=14),
                       legend.text =  ggplot2::element_text(size=10),
                       legend.key.size = unit(0.3, 'cm'),
                       plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                       legend.position=c(1.1, 0.4))
      ggplot2::ggsave(p, filename=here::here("outputs", model, "inferred_perf_environment.png"),
                      width=fig_width, height=fig_width/2, units="cm", dpi=300)
      }
}

plot_inferred_perf_IV <- function(n_observed_axis, Obs_env, nsp, Inferred_species_parameters, V_intra, model, fig_width){

  x <- seq(-3, 3, 0.01)
  n_ind_simul <- length(x)
  
  if(n_observed_axis>1){
    mean_env <- mean(Obs_env[,k])
  } else{
    mean_env <- mean(c(Obs_env))
    }
  
  for(k in 1:n_observed_axis){
    
    Env_mat <- matrix(nrow=nsp, ncol=2*n_observed_axis+1)
    
    Env_mat[,1] <- rep(1, nsp)
    
    for (l in 1:n_observed_axis){
      Env_mat[,l+1] <- rep(mean_env, nsp)
      Env_mat[,l+1+n_observed_axis] <- rep(mean_env^2, nsp)
    }
    
    perf_mean_env <- rowSums(Env_mat * Inferred_species_parameters)
    
    perf_IV <- data.frame(Species = rep(1:nsp, each=n_ind_simul), X = rep(x, nsp), Mean=rep(perf_mean_env, each=n_ind_simul))
    
    perf_IV <- perf_IV%>%
      dplyr::group_by(Species)%>%
      dplyr::mutate(Density = dnorm(x, mean=Mean, sd=sqrt(V_intra$V[Species])))%>%
      dplyr::ungroup()
    
    p <- ggplot2::ggplot(perf_IV, ggplot2::aes(x = X, y = Density, colour = as.factor(Species))) +
      ggplot2::geom_line()+
      #ggplot2::scale_colour_viridis_d()+
      ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
      ggplot2::labs(x = glue::glue("Performance at mean environment variable {k}"),
                    y = "Density")+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     legend.title = ggplot2::element_text(size=14),
                     legend.text =  ggplot2::element_text(size=10),
                     legend.key.size = unit(0.3, 'cm'),
                     plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                     legend.position=c(1.1, 0.4))
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, glue::glue("Perf_overlap_IV_variable_{k}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
  }
}

plot_abundance_species <- function(Abundances, model, fig_width){
  Tidy_abund <- data.frame(Abund = c(do.call(rbind, Abundances)),
                           Sp = rep(c(1:nsp), each=nrep*ngen),
                           Gen = rep(c(1:ngen), nsp*nrep),
                           Rep = rep(rep(c(1:nrep), each = ngen), nsp))
  
  Tidy_abund$Sp <- as.factor(Tidy_abund$Sp)
  Tidy_abund$Gen <- as.factor(Tidy_abund$Gen)
  Tidy_abund$Rep <- as.factor(Tidy_abund$Rep)
  
  Tidy_abund <- Tidy_abund %>%
    dplyr::group_by(Sp, Gen)%>%
    dplyr::mutate(Mean_sp_gen = mean(Abund), CI_025 = quantile(Abund, probs=0.025), CI_975 = quantile(Abund, probs=0.975))%>%
    dplyr::ungroup()
  
  Abund_plot <- Tidy_abund[Tidy_abund$Rep==1,]
  
  Abund_plot$Gen <- as.numeric(Abund_plot$Gen)
  
  p <- ggplot2::ggplot(data=Abund_plot, ggplot2::aes(x=Gen, y=Mean_sp_gen, colour=as.factor(Sp)))+
    ggplot2::geom_line()+
    #ggplot2::geom_ribbon(ggplot2::aes(y=Mean_sp_gen, ymin=CI_025, ymax=CI_975, fill=Sp), alpha=0.5)+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::labs(x="Generations", y="Species abundance")+
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size=14),
                   legend.text =  ggplot2::element_text(size=10),
                   legend.key.size = unit(0.3, 'cm'),
                   plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                   legend.position=c(1.1, 0.4))
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "Abundances.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_perf_community_end <- function(community_end, perf_Sp_mean, sites, fig_width){
  nrep <- length(community_end)
  nsite <- nrow(sites)
  for(r in 1:nrep){
    community <- community_end[[r]]
    perf_present <- rep(NA, nsite)
    w0 <- (as.vector(t(community))==0)
    perf_present[w0] <- NA
    perf_present[!w0] <- diag(perf_Sp_mean[!w0, as.vector(t(community))[!w0]])
    df_perf_sp_present_x1 <- data.frame(Species = as.vector(t(community)), Perf = perf_present, X1=sites[,1])
    
    p <- ggplot2::ggplot(df_perf_sp_present_x1, aes(x=X1, y=Perf, colour=as.factor(Species)))+
      ggplot2::geom_point()+
      ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount)[sort(unique(df_perf_sp_present_x1$Species))])+
      ggplot2::labs(x="Environment (first axis)", y="Performance of final community")
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, glue::glue("perf_community_end_{r}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_perf_suitable_habitat <- function(perf_Sp_mean, sites, fig_width){
  nsite <- nrow(sites)
  
  df_perf_suitable_x1 <- data.frame(Species = apply(perf_Sp_mean,1,which.max), Perf = apply(perf_Sp_mean,1,max), X1=sites[,1])
    
  p <- ggplot2::ggplot(df_perf_suitable_x1, aes(x=X1, y=Perf, colour=as.factor(Species)))+
    ggplot2::geom_point()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount)[sort(unique(df_perf_suitable_x1$Species))])+
    ggplot2::labs(x="Environment (first axis)", y="Performance of the theoretical winner")
    
    
    ggplot2::ggsave(p, filename=here::here("outputs", model, "perf_suitable_habitat.png"),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
}


plot_perf_opt_clandestine <- function(
                              n_observed_axis,
                              Inferred_species_parameters,
                              community_end,
                              sp_hab_freq,
                              model,
                              fig_width){
  
  Freq_sp_model_part_know <- as.data.frame(table(unlist(community_end)))
  
  opt_estim_winner <- data.frame(Species = 1:nrow(Inferred_species_parameters),
                             Opt_estim=numeric(nrow(Inferred_species_parameters)),
                             Perf_opt=numeric(nrow(Inferred_species_parameters)),
                             Frequency=numeric(nrow(Inferred_species_parameters)),
                             Presence=numeric(nrow(Inferred_species_parameters)),
                             Hab_freq=as.data.frame(sp_hab_freq)$Freq,
                             Winner=character(nrow(Inferred_species_parameters)))
  
  for(sp in 1:nrow(Inferred_species_parameters)){
    for(axis in 1:n_observed_axis){
      opt_estim_winner$Opt_estim[sp] <- (-Inferred_species_parameters[sp,axis+1])/(2*Inferred_species_parameters[sp,2*axis+1])
    }
    if(sp%in%unlist(community_end)){
      opt_estim_winner$Frequency[sp] <- Freq_sp_model_part_know[Freq_sp_model_part_know$Var1==sp,]$Freq
      opt_estim_winner$Presence[sp] <- 1
    }else{
      opt_estim_winner$Frequency[sp]<- 0
      opt_estim_winner$Presence[sp] <- 0
    }
    if(opt_estim_winner$Hab_freq[sp]!=0){
      opt_estim_winner$Winner[sp] <- "theoretical"
    }
  }
  
  if(length(opt_estim_winner[opt_estim_winner$Opt_estim>1,]$Opt_estim)!=0){
    opt_estim_winner[opt_estim_winner$Opt_estim>1,]$Opt_estim <- 1
  }
  if(length(opt_estim_winner[opt_estim_winner$Opt_estim<0,]$Opt_estim)!=0){
    opt_estim_winner[opt_estim_winner$Opt_estim<0,]$Opt_estim <- 0
  }
  
  for(sp in 1:nrow(Inferred_species_parameters)){
    for(axis in 1:n_observed_axis){
      opt_estim_winner$Perf_opt[sp] <- Inferred_species_parameters[sp,axis] + Inferred_species_parameters[sp,axis+1]*opt_estim_winner$Opt_estim[sp] + Inferred_species_parameters[sp,2*axis+1]*(opt_estim_winner$Opt_estim[sp])^2
    }
  }
  
  opt_estim_winner[which(opt_estim_winner$Winner!="theoretical"&opt_estim_winner$Presence==1),]$Winner<-"clandestine"
  opt_estim_winner[which(opt_estim_winner$Winner=="theoretical"&opt_estim_winner$Presence==1),]$Winner<-"theoretical - present"
  if(length(which(opt_estim_winner$Winner=="theoretical"&opt_estim_winner$Presence==0))!=0){
    opt_estim_winner[which(opt_estim_winner$Winner=="theoretical"&opt_estim_winner$Presence==0),]$Winner<-"theoretical - absent"
  }
  opt_estim_winner[which(opt_estim_winner$Winner!="theoretical"&opt_estim_winner$Presence==0),]$Winner<-"disappeared"
  
  ggplot2::ggplot(opt_estim_winner, aes(Opt_estim, Perf_opt))+
    ggplot2::geom_point(aes(colour=as.factor(Winner)), size=3)+
    ggplot2::scale_colour_viridis_d(option = "plasma")+
    ggplot2::labs(x="Estimated optimum",
                  y="Performance at optimum",
                  colour="Presence in final community")+
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.title = ggplot2::element_text(size=16),
                   legend.text =  ggplot2::element_text(size=14))
}
