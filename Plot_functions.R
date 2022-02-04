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
    plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
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

plot_mortality_events <- function(model, fig_width, community, nsp){
  png(file=here::here("outputs", model, "mortality_events.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  plot(raster::raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), legend=FALSE, cex.main=2, cex.axis=1.5)
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
}

plot_species_richness_log_10 <- function(sp_rich, ngen, model, fig_width){
  
  load(here::here("outputs", model, glue::glue("sp_rich_{model}.RData")))
  
  if(ncol(sp_rich)>1){
    sp_rich_gen <- data.frame(sp_rich=sp_rich[,1], gen=c(1:ngen))
  }else{sp_rich_gen <- data.frame(sp_rich=sp_rich, gen=c(1:ngen))}
  
  p <- ggplot2::ggplot(data=sp_rich_gen, ggplot2::aes(x=as.numeric(gen), y=sp_rich)) +
    ggplot2::geom_line() +
    scale_x_continuous(trans='log10')+
    ggplot2::xlab("Generations") + 
    ggplot2::ylab("Species richness")+
    ggplot2::theme(legend.position = "none",
                   text = ggplot2::element_text(size = 20))+
    ggplot2::ylim(0,nsp)
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "species_richness_log_10.png"),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_mean_rank_hab_freq <- function(sp_mean_rank, sp_hab_freq, model, fig_width){
  df <- data.frame(cbind(Species=1:nsp, Hab_freq=sp_hab_freq, Rank=sp_mean_rank))
  
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=Hab_freq, y=Rank, colour=as.factor(Species))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
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
    ggplot2::ylab("Mean env-species perf difference")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "environmental_filtering.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}
  
plot_env_species<- function(model, fig_width, community_start, community_end, env){
  png(file=here::here("outputs", model, "spatial_comp_env_sp.png"), 
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  plot(raster::raster(community_start), main="Species - Start", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), legend=FALSE)
  plot(raster::raster(community_end), main="Species - End", zlim=c(0, nsp),
       col=c("black", getPalette(colourCount)), legend=FALSE)
  env_stack <- stack(raster::raster(env[[1]]*255), raster::raster(env[[2]]*255), raster::raster(env[[3]]*255))
  crs(env_stack) <- "+proj=utm +zone=1"
  plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
  dev.off()
}

plot_theta_community<-function(theta_comm, ngen, model, fig_width){
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
    ggplot2::ylab("Mean mortality rate in the community")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "mortality_rate_community.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_spatial_autocorr <- function(community_end, sites, niche_width, model, fig_width){
  # Species autocorrelation
  sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
  names(sp_XY) <- c("x", "y", "sp")
  vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
  # Environment autocorrelation
  # 3D voxel for each site
  x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
  y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
  z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
  n_niche <- 1/niche_width
  class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
  vario_env <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
  # Plot with correlation
  png(file=here::here("outputs", model, "sp_autocorrelation.png"),
      width=fig_width, height=fig_width*0.8, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  plot(vario_sp, main="Species - End")
  plot(vario_env, main="Environment")
  plot(vario_env$v, vario_sp$v,
       main = "Regression",
       xlab="Semivariance for environment",
       ylab="Semivariance for species")
  m <- lm(vario_sp$v ~ vario_env$v-1)
  abline(a=0, b=coef(m), col="red")
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
    ggplot2::geom_line(aes(x=Env_1, y=-dist), col="#088000")+
    geom_vline(aes(xintercept=optimum), col="#80002D")+
    ggplot2::facet_wrap(vars(Species), nrow=3) +
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
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))
    #ggplot2::xlim(c(0,1))+
    #ggplot2::ylim(c(0,1))

  ggplot2::ggsave(p, filename=here::here("outputs", model, "optima_real_estim.png"),
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