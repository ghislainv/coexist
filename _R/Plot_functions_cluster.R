colourCount = nsp
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

plot_environment <- function( fig_width, n_axes, env, sites){
  
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_environment.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  
  par(mfrow=c(4,ceiling(n_axes/4)), bty = "n")
  
  for(k in 1:n_axes){
    raster::plot(raster::raster(env[[k]]), main=glue::glue("Environment var {k}"), col=topo.colors(255), cex.main=1.2)
  }
  #Summary of environment
  if(n_axes==3){
    env_stack <- raster::stack(
      raster::raster(matrix(range_0_255(env$x[,1]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE)),
      raster::raster(matrix(range_0_255(env$x[,2]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE)),
      raster::raster(matrix(range_0_255(env$x[,3]), nrow=nsite_side, ncol=nsite_side, byrow=TRUE))
    )
    raster::crs(env_stack) <- "+proj=utm +zone=1"
    class_site <- entrelac(
      raster::values(env_stack@layers[[1]]),
      raster::values(env_stack@layers[[2]]),
      raster::values(env_stack@layers[[3]])
    )
    save(class_site, file = paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_env_entrelac.RData")))
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), main="Environment summary", col=viridisLite::viridis(255^3))
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
    class_site <- entrelac(
      raster::values(env_stack@layers[[1]]),
      raster::values(env_stack@layers[[2]]),
      raster::values(env_stack@layers[[3]])
    )
    save(class_site, file = paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_env_entrelac.RData")))
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), main="Environment summary", col=viridisLite::viridis(255^3))
    dev.off()
    
    png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_env_pca.png")),
        width=fig_width, height=fig_width*0.8, units="cm", res=300)
    raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)), main="Environment summary", col=viridisLite::viridis(255^3))
    dev.off()
    
    var_explained <- pca_env$sdev^2 / sum(pca_env$sdev^2)
    p <- ggplot2::ggplot(data.frame(Prin_comp=c(1:length(var_explained)),
                                    Prop=cumsum(var_explained)),
                         ggplot2::aes(x=Prin_comp, y=Prop))+
      ggplot2::geom_col()+
      ggplot2::labs(x="Principal component",
                    y="Cumulative proportion \n of explained variance")
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_pca_prop_var.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_hab_freq <- function(n_axes, fig_width, env){
  
  for(k in 1:n_axes){
    png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_hab_freq_{k}.png")),
        width=fig_width, height=fig_width, units="cm", res=300)
    hist(env[[k]], main="", xlab=glue::glue("Environment var{k}")) 
    dev.off()
  }
}

plot_species_optima <- function( fig_width, niche_optimum){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_niche.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mar=c(1,1,2,2))
  scatter3D(niche_optimum[,1], niche_optimum[,2], niche_optimum[,3],
            pch=16, 
            colvar=1:nsp, col=getPalette(colourCount),
            bty = "f", main ="Three-dimensional species optima", phi=0,
            xlim=c(0,1), ylim=c(0,1), zlim=c(0,1))
  dev.off()
}

plot_histogram_mortality <- function( fig_width, mortality_E_Sp){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_hist_mortality.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  hist(mortality_E_Sp)
  dev.off()
}

plot_function_mort_proba <- function( fig_width, perf_E_Sp, mortality_E_Sp){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_function_mort_proba.png")),
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

plot_species_habitat_freq <- function( fig_width, sp_hab_freq){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_habitat_freq.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(cex.lab=1.5)
  plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency", col=getPalette(colourCount))
  dev.off()
}

plot_community_start <- function( fig_width, community, nsp){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_community_start.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  # raster::plot(raster::raster(community), main="Species - Start", zlim=c(0, nsp),
  #      col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  raster::plot(raster::raster(community), main="Species - Start", zlim=c(0, nsp),
               col=c("black", getPalette(colourCount)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_mortality_events <- function( fig_width, community, mortality, nsp){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mortality_events.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  community_mortality <- community
  community_mortality[mortality==1] <- nsp+1
  raster::plot(raster::raster(community_mortality), main="Species (colours) and vacant sites \n (black = old, white = new)", zlim=c(0, nsp),
               col=c("black", getPalette(colourCount)), "white", legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_community_end <- function( fig_width, community, nsp){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  # raster::plot(raster::raster(community), main=glue::glue("Species - End (ngen={ngen})"),
  #      zlim=c(0, nsp), col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  raster::plot(raster::raster(community), main=glue::glue("Species - End (ngen={ngen})"), 
               zlim=c(0, nsp), col=c("black", getPalette(colourCount)), legend=TRUE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_species_richness <- function(nrep, sp_rich, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  if (nrep > 1){
    sp_rich_long <- sp_rich %>%
      dplyr::mutate(gen=1:(ngen)) %>%
      tidyr::pivot_longer(cols=colnames_long, names_to="rep",
                          names_prefix="X", values_to="sp_rich")
  }else{sp_rich_long<-data.frame(gen=1:ngen, rep=rep(nrep, ngen), sp_rich=sp_rich)}
  
  #One curve per repetition
  if(nrep > 1){
    p <- ggplot2::ggplot(data=sp_rich_long, ggplot2::aes(x=gen, y=sp_rich, col=rep)) +
      ggplot2::geom_line() +
      ggplot2::scale_colour_viridis_d()+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
  }else{
    p <- ggplot2::ggplot(data=sp_rich_long, ggplot2::aes(x=gen, y=sp_rich)) +
      ggplot2::geom_line() +
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
  }
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_richness_with_time.png")),
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
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_richness_with_time_mean.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
  
  # Log_10 plot
  
  if(nrep==1){
    sp_rich_log <- data.frame(Sp_rich=sp_rich$sp_rich, Gen=c(1:ngen))
    p <- ggplot2::ggplot(data=sp_rich_log, ggplot2::aes(x=as.numeric(Gen), y=Sp_rich)) +
      ggplot2::geom_line(col="#008071")+
      scale_x_continuous(trans='log10')+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Species richness")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))+
      ggplot2::ylim(0,nsp)
    
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_richness_log_10.png")),
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
    
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_species_richness_log_10.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
  }
}

plot_mean_rank_hab_freq <- function(sp_mean_rank, sp_hab_freq, fig_width){
  mean_rank_hab_freq <- data.frame(cbind(Species=1:nsp, Hab_freq=sp_hab_freq, Rank=sp_mean_rank))
  
  p <- ggplot2::ggplot(data=mean_rank_hab_freq, ggplot2::aes(x=Hab_freq, y=Rank, colour=as.factor(Species))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="#008071", fill="#69b3a2", se=TRUE) +
    ggplot2::xlab("Species suitable habitat frequency") +
    ggplot2::ylab("Species mean rank (higher rank = lower abundance)") +
    ggplot2::theme(axis.title=ggplot2::element_text(size=16))+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mean_rank-habitat_freq.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_env_filt<- function(nrep, env_filt, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  if(nrep>1){
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
                     text = ggplot2::element_text(size = 16))
  }else{
    
    env_filt_long <- data.frame(gen=1:ngen, rep=rep(nrep, ngen), env_filt=env_filt)
    
    p <- ggplot2::ggplot(data=env_filt_long, ggplot2::aes(x=gen, y=env_filt)) +
      ggplot2::geom_line() +
      ggplot2::labs(title="Environmental filtering") +
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Mean env-species perf difference")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 16))
  }
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_environmental_filtering.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  #Mean and 95% interval
  if(nrep>1){
    env_filt_mean <- data.frame(Gen = 1:nrow(env_filt), Mean=rowMeans(env_filt), Low=apply(env_filt, 1, quantile, probs=0.025), High=apply(env_filt, 1, quantile, probs=0.975))
    
    p <- ggplot2::ggplot(data=env_filt_mean, ggplot2::aes(x=Gen, y=Mean)) +
      ggplot2::geom_line(col="#008071") +
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Environmental filtering")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 20))
    
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_environmental_filtering_mean.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_env_species <- function( fig_width, community_start, community_end, class_site){
  png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_spatial_comp_env_sp.png")), 
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  raster::plot(raster::raster(community_start), main="Species - Start", zlim=c(0, nsp),
               col=c("black", getPalette(colourCount)), legend=FALSE)
  raster::plot(raster::raster(community_end), main="Species - End", zlim=c(0, nsp),
               col=c("black", getPalette(colourCount)), legend=FALSE)
  raster::plot(raster::raster(matrix(class_site, ncol=nsite_side, nrow=nsite_side, byrow=TRUE)),
               main="Environment summary", col=viridisLite::viridis(255^3))
  dev.off()
}

plot_theta_community<-function(theta_comm, ngen, nrep, fig_width){
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
                   text = ggplot2::element_text(size = 18))
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mortality_rate_community.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  #Mean and 95% interval
  if(nrep>1){
    theta_comm_mean <- data.frame(Gen = 1:nrow(theta_comm), Mean=rowMeans(theta_comm), Low=apply(theta_comm, 1, quantile, probs=0.025), High=apply(theta_comm, 1, quantile, probs=0.975))
    
    p <- ggplot2::ggplot(data=theta_comm_mean, ggplot2::aes(x=Gen, y=Mean)) +
      ggplot2::geom_line(col="#008071") +
      ggplot2::geom_ribbon(aes(ymin=Low, ymax=High), col="#008071", fill="#008071", alpha=0.5)+
      ggplot2::xlab("Generations") + 
      ggplot2::ylab("Mean mortality")+
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 18))
    
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mortality_rate_community_mean.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_spatial_autocorr <- function(nrep, community_end, n_axes, sites, niche_optimum, niche_width, fig_width){
  
  semivar_multidim <- list()
  
  for(rep in 1:nrep){
    sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end[[rep]])))
    names(sp_XY) <- c("x", "y", "sp")
    vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
    
    
    if(randomOptSp==FALSE&n_axes==3){
      # 3D voxel for each site
      x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
      y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
      z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
      n_niche <- 1/niche_width
      # This is done to avoid having the same class for different combinations (the function is not continuous)
      class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
      vario_env <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
      plot(vario_sp)
      plot(vario_env)
      plot(vario_env$v, vario_sp$v)
    }else{
      
      semivar_multidim[[rep]] <- compute_semivar_multidim(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end[[rep]])
      semivar_multidim[[rep]]$Vario_sp_geoR <- vario_sp$u
      semivar_multidim[[rep]]$Distance <- vario_sp$bins.lim[-length(vario_sp$bins.lim)]
      semivar_multidim[[rep]]$Sample_size <- vario_sp$n
      
      if(rep==1){
        # Remove last point if the number of pairs is too low
        if(semivar_multidim[[rep]]$Sample_size[nrow(semivar_multidim[[rep]])]<500){
          semivar_multidim_plot <- semivar_multidim[[rep]][1:(nrow(semivar_multidim[[rep]])-1),]
        }
        # Plot with correlation
        png(file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_sp_autocorrelation.png")),
            width=fig_width, height=fig_width*0.8, units="cm", res=300)
        par(mfrow=c(2,2), bty = "n")
        #Species
        plot(semivar_multidim_plot$Distance, semivar_multidim_plot$Vario_sp,
             main="Species - end",
             xlab="distance",
             ylab="semivariance")
        #Environment
        plot(semivar_multidim_plot$Distance, semivar_multidim_plot$Vario_env,
             main="Environment",
             xlab="distance",
             ylab="semivariance")
        #Regression
        plot(semivar_multidim_plot$Vario_env, semivar_multidim_plot$Vario_sp,
             main = "Regression",
             xlab="Semivariance for environment",
             ylab="Semivariance for species")
        m <- lm(semivar_multidim_plot$Vario_sp ~ semivar_multidim_plot$Vario_env)
        abline(a=as.numeric(coef(m)[1]), b=as.numeric(coef(m)[2]), col="#008071")
        text(semivar_multidim_plot$Vario_env[3], 0.95*max(semivar_multidim_plot$Vario_sp), paste("R =", round(sqrt(summary(m)$r.squared), digits = 2)))
        dev.off()
      }
    }
    semivar_multidim[[rep]]$Rep <- rep(rep, nrow(semivar_multidim[[rep]]))
  }
  save(semivar_multidim, file=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_semivar_multidim.RData")))
}

plot_species_niche <- function(seed, df_perf, Inferred_species_parameters, V_intra, sites, fig_width){
  # Select 8 species at random
  #set.seed(seed)
  #sp_sel <- sample(unique(df_perf$Species), 9, replace=FALSE)
  #df_sp_sel <- df_perf %>% filter(Species %in% sp_sel)
  
  df_perf_sp_niche <- df_perf
  load(paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  df_perf_sp_niche$optimum <- rep(niche_optimum$X1, nsite)
  df_perf_sp_niche$dist <- sqrt((df_perf_sp_niche$Env_1-df_perf_sp_niche$optimum)^2)
  df_perf_sp_niche$dist <- df_perf_sp_niche$dist - mean(df_perf_sp_niche$dist) / sd(df_perf_sp_niche$dist)
  
  p <- ggplot2::ggplot(data=df_perf_sp_niche, ggplot2::aes(x=Env_1, y=Perf))+
    ggplot2::geom_point(size=0.5, alpha=0.5) +
    ggplot2::geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE, col="#008071")+
    ggplot2::scale_x_continuous(labels = function(x) ifelse(x == 0, "0", sub("^0+", "", x)))+
    ggplot2::geom_line(ggplot2::aes(x=Env_1, y=-dist), col="#088000")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=optimum), col="#80002D")+
    ggplot2::facet_wrap(ggplot2::vars(Species), nrow=4)+
    ggplot2::xlab("Environment (first axis)") +
    ggplot2::ylab("Performance")+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size=14),
                   strip.text = element_text(size = 16))+
    ggplot2::coord_fixed(ratio = round((max(df_perf_sp_niche$Env_1)-min(df_perf_sp_niche$Env_1))/(max(df_perf_sp_niche$Perf)-min(df_perf_sp_niche$Perf)), digits=2))
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_infering_species_niche.png")),
                  width=fig_width*2, height=fig_width, units="cm", dpi=300)
  
  p <- ggplot2::ggplot(data=df_perf_sp_niche, ggplot2::aes(x=Env_1, y=Perf))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE, col="#008071")+
    ggplot2::scale_x_continuous(labels = function(x) ifelse(x == 0, "0", sub("^0+", "", x)))+
    ggplot2::facet_wrap(ggplot2::vars(Species), nrow=4)+
    ggplot2::xlab("Environment (first axis)")+
    ggplot2::ylab("Performance")+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size=14),
                   strip.text = element_text(size = 16))+
    ggplot2::coord_fixed(ratio = round((max(df_perf_sp_niche$Env_1)-min(df_perf_sp_niche$Env_1))/(max(df_perf_sp_niche$Perf)-min(df_perf_sp_niche$Perf)), digits=2))
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_infering_species_niche_simple.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  
  X1_mat <- list()
  X1_mat[[1]] <- matrix(rep(sites[,1], nrow(Inferred_species_parameters)), ncol=nrow(Inferred_species_parameters))
  X1_mat[[2]] <- X1_mat[[1]]^2
  
  Inferred_parameters_mat_X1 <-list()
  
  for(k in 1:ncol(Inferred_species_parameters)){
    Inferred_parameters_mat_X1[[k]] <- matrix(rep(Inferred_species_parameters[,k],each=length(sites[,1])), ncol=nrow(Inferred_species_parameters))
  }
  
  Mat_perf_inferred <- Inferred_parameters_mat_X1[[1]]
  
  for(k in 1:length(X1_mat)){
    Mat_perf_inferred <- Mat_perf_inferred + Inferred_parameters_mat_X1[[k+1]]*X1_mat[[k]]
  }
  
  Mat_perf_inferred_X1 <- as.data.frame(cbind(sites[,1], Mat_perf_inferred))
  
  colnames(Mat_perf_inferred_X1) <- c("Env_1", sprintf("Sp %02d", 1:nrow(Inferred_species_parameters)))
  
  Mat_perf_inferred_X1 <- Mat_perf_inferred_X1 %>%
    tidyr::pivot_longer(cols=2:(nrow(Inferred_species_parameters)+1), names_to="Species",
                        names_prefix="Species", values_to="Perf")
  
  Mat_perf_inferred_X1 <- Mat_perf_inferred_X1%>%
    dplyr::mutate(sd=rep(sqrt(V_intra$V), nrow(sites)))
  
  p <- ggplot2::ggplot(data=df_perf_sp_niche, ggplot2::aes(x=Env_1, y=Perf))+
    ggplot2::geom_point(size=0.5, alpha=0.5)+
    ggplot2::geom_line(data=Mat_perf_inferred_X1, ggplot2::aes(x=Env_1, y=Perf), col="#008071")+
    ggplot2::geom_ribbon(data=Mat_perf_inferred_X1, aes(ymin=Perf-sd, ymax=Perf+sd), col="#008071", fill="#008071", alpha=0.2, linetype="dotted")+
    ggplot2::scale_x_continuous(labels = function(x) ifelse(x == 0, "0", sub("^0+", "", x)))+
    ggplot2::facet_wrap(ggplot2::vars(Species), nrow=4)+
    ggplot2::xlab("Environment (first axis)")+
    ggplot2::ylab("Performance")+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size=14),
                   strip.text = element_text(size = 16))+
    ggplot2::coord_fixed(ratio = round((max(df_perf_sp_niche$Env_1)-min(df_perf_sp_niche$Env_1))/(max(df_perf_sp_niche$Perf)-min(df_perf_sp_niche$Perf)), digits=2))
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_infering_species_niche_real_params.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
}

plot_IV <- function(V_intra, fig_width){
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
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_IV.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_relationship_IV_inferred_optima <- function(V_intra,
                                                 n_observed_axes,
                                                 Inferred_species_parameters,
                                                 community_end,
                                                 sp_hab_freq,
                                                 
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
    for(axis in 1:n_observed_axes){
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
    for(axis in 1:n_observed_axes){
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

plot_optima_real_estim <- function(nsp, n_observed_axes, niche_optimum, Inferred_species_parameters, fig_width){
  
  opt_sp_estim_vs_real <- data.frame(Species = 1:nsp, Real=numeric(nsp), Estim=numeric(nsp))
  
  for(sp in 1:nsp){
    
    for(axis in 1:n_observed_axes){
      opt_sp_estim_vs_real$Real[sp] <- niche_optimum[sp, axis]
      opt_sp_estim_vs_real$Estim[sp] <- (-Inferred_species_parameters[sp,axis+1])/(2*Inferred_species_parameters[sp,2*axis+1])
    }
    
  }
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Real, y=Estim, colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::labs(x="Real species optima on axis X1", y="Estimated species optima")
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_optima_real_estim.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  opt_sp_estim_vs_real$Horizontal <- rep(0, nrow(opt_sp_estim_vs_real))
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Estim, y=Horizontal, colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_optima_estim.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
  
  if(length(opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim>1,]$Estim)!=0){
    opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim>1,]$Estim <- 1}
  if(length(opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim<0,]$Estim)!=0){
    opt_sp_estim_vs_real[opt_sp_estim_vs_real$Estim<0,]$Estim <- 0
  }
  
  p <- ggplot2::ggplot(data=opt_sp_estim_vs_real, ggplot2::aes(x=Estim, y=Horizontal, colour=factor(Species)))+
    ggplot2::geom_point()+
    #ggplot2::scale_colour_viridis_d()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount))+
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank())
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_optima_estim_0_1.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_inferred_perf_environment <- function(E_seq, Mat_perf_inferred, nsp, fig_width){
  
  if(n_observed_axes>1){
    
    for(k in 1:n_observed_axes){
      
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
      
      ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_inferred_perf_environment_{k}.png")),
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
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_inferred_perf_environment.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}

plot_inferred_perf_IV <- function(n_observed_axes, Obs_env, nsp, Inferred_species_parameters, V_intra, fig_width){
  
  x <- seq(-3, 3, 0.01)
  n_ind_simul <- length(x)
  
  for(k in 1:n_observed_axes){
    
    if(n_observed_axes>1){
      mean_env <- mean(Obs_env[,k])
    } else{
      mean_env <- mean(c(Obs_env))
    }
    
    Env_mat <- matrix(nrow=nsp, ncol=2*n_observed_axes+1)
    
    Env_mat[,1] <- rep(1, nsp)
    
    for (l in 1:n_observed_axes){
      Env_mat[,l+1] <- rep(mean_env, nsp)
      Env_mat[,l+1+n_observed_axes] <- rep(mean_env^2, nsp)
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
      ggplot2::theme(text = ggplot2::element_text(size = 18),
                     legend.title = ggplot2::element_text(size=14),
                     legend.text =  ggplot2::element_text(size=10),
                     legend.key.size = unit(0.3, 'cm'),
                     plot.margin=unit(c(1, 5, 1, 1), 'lines'),
                     legend.position=c(1.1, 0.4))
    
    ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_Perf_overlap_IV_variable_{k}_{n_observed_axes}_obs_axes.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
  }
}

plot_abundance_species <- function(Abundances, fig_width){
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
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_perf_community_end <- function(community_end, perf_Sp_mean, sites, fig_width){
  nrep <- length(community_end)
  nsite <- nrow(sites)
  # Only one plot or one for each repetition
  #for(r in 1:nrep){
  r <- 1
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
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_perf_community_end_{seed_r}.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
  # }
}

plot_perf_suitable_habitat <- function(perf_Sp_mean, sites, fig_width){
  nsite <- nrow(sites)
  
  df_perf_suitable_x1 <- data.frame(Species = apply(perf_Sp_mean,1,which.max), Perf = apply(perf_Sp_mean,1,max), X1=sites[,1])
  
  p <- ggplot2::ggplot(df_perf_suitable_x1, aes(x=X1, y=Perf, colour=as.factor(Species)))+
    ggplot2::geom_point()+
    ggplot2::scale_colour_manual(name="Species", values = getPalette(colourCount)[sort(unique(df_perf_suitable_x1$Species))])+
    ggplot2::labs(x="Environment (first axis)", y="Performance of the theoretical winner")
  
  
  ggplot2::ggsave(p, filename=paste0(directory_writing, "outputs", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_perf_suitable_habitat.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}


plot_perf_opt_clandestine <- function(
  n_observed_axes,
  Inferred_species_parameters,
  community_end,
  sp_hab_freq,
  
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
    for(axis in 1:n_observed_axes){
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
    for(axis in 1:n_observed_axes){
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

###########################################################################

source("~/Code/coexist/_R/call_libraries.R")
source("~/Code/coexist/_R/Math_functions.R")

fig_width <- 35

load(here::here("Array_simulations.RData"))
ngen = 10000
nsp = 20
n_axes = 15
dir.create(here::here("outputs", glue::glue("Comparison")))

# Build dataset

Abundances_all <- data.frame(Mortality = numeric(),
                             Fecundity = numeric(),
                             Seed = numeric(),
                             Seed_r = numeric(),
                             Mod = character(),
                             Nb_obs = numeric(),
                             Species = numeric(),
                             Abundance = numeric())

Species_all <- data.frame(Mortality = numeric(),
                          Fecundity = numeric(),
                          Seed = numeric(),
                          Seed_r = numeric(),
                          Mod = character(),
                          Nb_obs = numeric(),
                          N_sp = numeric(),
                          Shannon = numeric())

for (simu in c(1:nrow(Simulations))){
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
    
    if(mod == "Perf_know"){
      
      load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
      
      Abundances_tmp <- data.frame(Mortality=rep(mortality, 20),
                                   Fecundity=rep(fecundity, 20),
                                   Mod=rep(mod, 20),
                                   Nb_obs=rep(0, 20),
                                   Seed=rep(seed, 20),
                                   Seed_r=rep(seed_r, 20),
                                   Species=c(1:20),
                                   Abundance=Abundances[ngen,])
      Abundances_all <- rbind(Abundances_all, Abundances_tmp)
      
      nsp_final <- length(which(Abundances[ngen,]!=0))
      
      df_shannon <- data.frame(Species = 1:nsp,
                               Abundance = Abundances[ngen,])%>%
        dplyr::mutate(Proportion = Abundance / sum(Abundance))%>%
        dplyr::filter(Abundance > 0)%>%
        dplyr::mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
      
      shannon <- -sum(df_shannon$prop_times_ln_prop)
      
      Species_tmp <- data.frame(Mortality=mortality,
                                Fecundity=fecundity,
                                Mod=mod,
                                Nb_obs=0,
                                Seed=seed,
                                Seed_r=seed_r,
                                N_sp=nsp_final,
                                Shannon=shannon)
                          
      Species_all <- rbind(Species_all, Species_tmp)
      
    }else{
      
      for (nb_obs in c(0:15)){
        
        load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
        Abundances_tmp <- data.frame(Mortality=rep(mortality, 20),
                                     Fecundity=rep(fecundity, 20),
                                     Mod=rep(mod, 20),
                                     Nb_obs=rep(nb_obs, 20),
                                     Seed=rep(seed, 20),
                                     Seed_r=rep(seed_r, 20),
                                     Species=c(1:20),
                                     Abundance=Abundances[ngen,])
        Abundances_all <- rbind(Abundances_all, Abundances_tmp)
        
        nsp_final <- length(which(Abundances[ngen,]!=0))
        
        df_shannon <- data.frame(Species = 1:nsp,
                                 Abundance = Abundances[ngen,])%>%
          dplyr::mutate(Proportion = Abundance / sum(Abundance))%>%
          dplyr::filter(Abundance > 0)%>%
          dplyr::mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
        
        shannon <- -sum(df_shannon$prop_times_ln_prop)
        
        Species_tmp <- data.frame(Mortality=mortality,
                                  Fecundity=fecundity,
                                  Mod=mod,
                                  Nb_obs=nb_obs,
                                  Seed=seed,
                                  Seed_r=seed_r,
                                  N_sp=nsp_final,
                                  Shannon=shannon)
        
        Species_all <- rbind(Species_all, Species_tmp)
      }
    }
  }
}

save(Abundances_all, file = here::here("outputs", "Comparison", "Abundances_all.RData"))
save(Species_all, file = here::here("outputs", "Comparison", "Species_all.RData"))

Percentage_similarity <- data.frame(
  Mortality = numeric(),
  Fecundity = numeric(),
  Seed = numeric(),
  Seed_r = numeric(),
  Mod_comp = character(),
  Nb_obs = numeric(),
  PS = numeric())

for (simu in c(1:nrow(Simulations))){
  
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
  Abundances_perf <- Abundances[ngen,]
  
  Percentage_similarity <- rbind(Percentage_similarity,
                                 data.frame(
                                   Mortality = mortality,
                                   Fecundity = fecundity,
                                   Seed = seed,
                                   Seed_r = seed_r,
                                   Mod_comp = "Perf_know",
                                   Nb_obs = NA,
                                   PS = 1
                                 ))
  
  for(nb_obs in 0:15){
    
    load(here::here("outputs_cluster", glue::glue("Part_know_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
    Abundances_part <- Abundances[ngen,]
    
    load(here::here("outputs_cluster", glue::glue("Part_know_IV_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
    Abundances_part_IV <- Abundances[ngen,]
    
    Percentage_similarity <- rbind(Percentage_similarity,
                                   data.frame(
                                     Mortality = rep(mortality,2),
                                     Fecundity = rep(fecundity,2),
                                     Seed = rep(seed,2),
                                     Seed_r = rep(seed_r,2),
                                     Mod_comp = c("Part_know", "Part_know_IV"),
                                     Nb_obs = rep(nb_obs,2),
                                     PS = c(percentage_similarity(Abundances_perf, Abundances_part), percentage_similarity(Abundances_perf, Abundances_part_IV))
                                   ))
  }
}

save(Percentage_similarity, file = here::here("outputs", "Comparison", "Percentage_similarity.RData"))


#### Performance on sites where the species is present in the final community ####

Perf_on_sites <- data.frame(
  Mortality = numeric(),
  Fecundity = numeric(),
  Seed = numeric(),
  Seed_r = numeric(),
  Mod = character(),
  Nb_obs = numeric(),
  Perf = numeric()
)

for (simu in c(1:nrow(Simulations))){
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  #Load the species performance once per configuration
  if(seed_r==1){
    load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_perf_E_Sp.RData")))
  }
  
  for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
    
    if(mod == "Perf_know"){
      
      load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
      
      perf_present <- rep(NA, length(community_end))
      w0 <- (as.vector(t(community_end))==0)
      perf_present[w0] <- NA
      perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
      
      Perf_on_sites_temp <- data.frame(
        Mortality = mortality,
        Fecundity = fecundity,
        Seed = seed,
        Seed_r = seed_r,
        Mod = mod,
        Nb_obs = NA,
        Perf = round(mean(perf_present, na.rm = TRUE), digits = 2)
      )
      
      Perf_on_sites <- rbind(Perf_on_sites, Perf_on_sites_temp)
      
    }else{
      
      for (nb_obs in c(0:15)){
        
        load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        
        perf_present <- rep(NA, length(community_end))
        w0 <- (as.vector(t(community_end))==0)
        perf_present[w0] <- NA
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
        
        Perf_on_sites_temp <- data.frame(
          Mortality = mortality,
          Fecundity = fecundity,
          Seed = seed,
          Seed_r = seed_r,
          Mod = mod,
          Nb_obs = nb_obs,
          Perf = round(mean(perf_present, na.rm = TRUE), digits = 2)
        )
        
        Perf_on_sites <- rbind(Perf_on_sites, Perf_on_sites_temp)
        
      }#for n_obs
    }#else Part_know
  }#for mod
}#for simu

save(Perf_on_sites, file = here::here("outputs", "Comparison", "Perf_on_sites.RData"))

# Retrieve R2 of statistical model #

R2_df <- data.frame(Mortality=character(),
                    Fecundity=character(),
                    Seed=integer(),
                    Nb_obs=integer(),
                    R2=numeric())

IV_all <- data.frame(Mortality=character(),
                     Fecundity=character(),
                     Seed=integer(),
                     Nb_obs=integer(),
                     IV=numeric())

for (simu in c(1:nrow(Simulations[Simulations[,4]==1,]))){
  
  print(paste("Simu", simu))
  
  mortality <- Simulations[Simulations[,4]==1,][simu,1]
  fecundity <- Simulations[Simulations[,4]==1,][simu,2]
  seed <- Simulations[Simulations[,4]==1,][simu,3]
  
  load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_perf_E_Sp.RData")))
  load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_env.RData")))
  load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_sites.RData")))
  load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_niche_optimum.RData")))
  
  nsp <- ncol(perf_E_Sp)
  n_axes <- length(env)
  
  for(n_observed_axes in c(0:n_axes)){
    
    if(n_observed_axes>0){
      Obs_env <- sites[,1:n_observed_axes]
    }
    
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
    
    
    if(n_observed_axes>0){
      
      formula <- as.formula(paste0("Perf~-1+Species+Species:",
                                   paste0(colnames(df_perf)[1:n_observed_axes], collapse= "+Species:"),
                                   "+Species:",
                                   paste0(colnames(df_perf)[(n_axes+1):(n_axes+n_observed_axes)], collapse= "+Species:")))
    }
    if(n_observed_axes==0){
      formula <- as.formula(paste0("Perf~-1+Species"))
    }
    
    lm_fit <- lm(formula, data=df_perf)
    
    R2_df_temp <- data.frame(Mortality=mortality,
                             Fecundity=fecundity,
                             Seed=seed,
                             Nb_obs=n_observed_axes,
                             R2=summary(lm_fit)$r.squared)
    R2_df <- rbind(R2_df, R2_df_temp)
    
    load(here::here("outputs_cluster", glue::glue("V_intra_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_1.RData")))
    IV_all_temp <- data.frame(Mortality=rep(mortality, nsp),
                             Fecundity=rep(fecundity, nsp),
                             Seed=rep(seed, nsp),
                             Nb_obs=rep(n_observed_axes, nsp),
                             IV=V_intra$V)
    IV_all <- rbind(IV_all, IV_all_temp)
    
  }
}

save(R2_df, file = here::here("outputs", "Comparison", "R2_df.RData"))
save(IV_all, file = here::here("outputs", "Comparison", "IV_all.RData"))

## PLOTS ##

load(here::here("outputs", "Comparison", "Species_all.RData"))
load(here::here("outputs", "Comparison", "Percentage_similarity.RData"))
load(here::here("outputs", "Comparison", "Perf_on_sites.RData"))
load(here::here("outputs", "Comparison", "R2_df.RData"))

Summary_level_explanation_axes_nb <- R2_df%>%
  dplyr::group_by(Mortality, Fecundity, Nb_obs)%>%
  dplyr::mutate(Mean_explanation=mean(R2), Sd=sd(R2))%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  dplyr::select(-Seed, -R2)

#One plot per Mortality * Fecundity option
for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")) {
    
    # Level of IV and of explained variation #
    
    data_figure <- IV_all[which(IV_all$Mortality==mortality&IV_all$Fecundity==fecundity),]
    
    data_figure_2 <- Summary_level_explanation_axes_nb[which(Summary_level_explanation_axes_nb$Mortality==mortality&Summary_level_explanation_axes_nb$Fecundity==fecundity),]
    
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=IV))+
      ggplot2::geom_ribbon(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation, ymin=Mean_explanation-Sd, ymax=Mean_explanation+Sd, group=1), colour="deeppink3", fill="hotpink3", alpha=0.3)+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
      ggplot2::geom_boxplot(alpha=0.6)+
      ggplot2::geom_point(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation), colour="deeppink3")+
      ggplot2::geom_line(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation, group=1), colour="deeppink3")+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = "Number of observed axes",
                    y = "Observed uIV")+
      ggplot2::scale_x_discrete(labels=c(0:n_axes))+
      ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")+
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . * 1 / 1 , name = "Proportion of variance explained by the axes"))
    
    ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("IV_nb_axes_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Species richness and Shannon diversity index #
    
    data_figure <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta_SR = N_sp - dplyr::lag(N_sp),
                    Delta_Shannon = Shannon - dplyr::lag(Shannon),
                    ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta_SR)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta_SR, Delta_Shannon)
    
    p5 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta_SR))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Species richness")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p5, filename=here::here("outputs", "Comparison", glue::glue("Delta_species_richness_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    p6 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta_Shannon))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Shannon index")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p6, filename=here::here("outputs", "Comparison", glue::glue("Delta_Shannon_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    #without deltas
    data_figure <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                                 levels = c(as.character(c(0:n_axes)), "PK"))
    
    data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
    data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    p1 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=N_sp))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = "Species richness with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p1, filename=here::here("outputs", "Comparison", glue::glue("Species_richness_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    p2 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Shannon))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = "Shannon index with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p2, filename=here::here("outputs", "Comparison", glue::glue("Shannon_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Percentage similarity #
    
    data_figure <- Percentage_similarity[which(Percentage_similarity$Mortality==mortality&Percentage_similarity$Fecundity==fecundity),]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta= PS - lag(PS), ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::select(Nb_obs, Seed, Seed_r, ID_delta, Delta)
    
    p7 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Similarity with perfect knowledge")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p7, filename=here::here("outputs", "Comparison", glue::glue("Delta_percentage_similarity_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    #without deltas
    data_figure <- Percentage_similarity[which(Percentage_similarity$Mortality==mortality&Percentage_similarity$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    
    # data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    # data_figure$Nb_obs <- factor(data_figure$Nb_obs,
    #                              levels = c(as.character(c(0:n_axes)), "PK"))
    
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- 16
    data_figure$Nb_obs <- as.integer(data_figure$Nb_obs)
    data_1 <- data_figure[data_figure$Mod!="Perf_know",]
    data_2 <- data_figure[data_figure$Mod=="Perf_know",]
    
    p3 <- ggplot2::ggplot(data=data_1, ggplot2::aes(x=Nb_obs, y=PS))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::geom_boxplot(alpha=0.6, colour="black", aes(group=Nb_obs))+
      geom_point(data=data_2, ggplot2::aes(x=Nb_obs, y=PS), colour="darkred")+
      ggplot2::scale_x_continuous(breaks=c(0:(n_axes+1)), labels=c(0:n_axes, "PK"))+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = "Similarity between PK and IP with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p3, filename=here::here("outputs", "Comparison", glue::glue("Percentage_similarity_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Theoretical performance (from the perfect knowledge model) of the final species community #
    
    data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality==mortality&Perf_on_sites$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta = Perf - dplyr::lag(Perf) , ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta)
    
    p8 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Theoretical performance")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p8, filename=here::here("outputs", "Comparison", glue::glue("Delta_performance_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality==mortality&Perf_on_sites$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                                 levels = c(as.character(c(0:n_axes)), "PK"))
    data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
    data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    #without deltas
    p4 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Perf))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = "Theoretical performance with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p4, filename=here::here("outputs", "Comparison", glue::glue("Performance_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    
    data_level_explanation <- Summary_level_explanation_axes_nb[Summary_level_explanation_axes_nb$Mortality==mortality&Summary_level_explanation_axes_nb$Fecundity==fecundity,]
    
    fun <- function(x) ifelse(x == 0, "0", sub("^0+", "", x))
    
    labels <- fun(round(data_level_explanation$Mean_explanation, digits=2))
    
    plot_level_explanation <- ggplot(data_level_explanation, aes(Nb_obs, Mean_explanation))+
      geom_blank()+
      theme_classic()+
      scale_x_continuous(name="Level of explained variance", breaks=data_level_explanation$Nb_obs, labels=labels)+
      theme(text = ggplot2::element_text(size = 16),
            axis.text = ggplot2::element_text(size=13),
            axis.line.x = ggplot2::element_line(arrow = arrow(length=unit(0.15, "inches")), colour="deeppink3"),
            axis.ticks.x=ggplot2::element_line(colour="deeppink3"),
            axis.text.x=ggplot2::element_text(colour="deeppink3"),
            axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.major.y=element_blank(),
            plot.margin = margin(l = 50,
                                 r = 0,
                                 t = 0,
                                 b = 0.5))
    
    arrange_SR_Shannon <- ggpubr::ggarrange(p1, p2, p5, p6,
                                   nrow=2, ncol=2, align = "v")
    arrange_PS_Perf <- ggpubr::ggarrange(p3, p4, p7, p8,
                                            nrow=2, ncol=2, align = "v")
    arrange_arrow <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation,
                                       nrow=1, ncol=2, align="v")
    arrange_results_1 <- ggpubr::ggarrange(arrange_SR_Shannon, arrange_arrow, nrow=2, ncol=1, heights=c(15, 1))
    arrange_results_2 <- ggpubr::ggarrange(arrange_PS_Perf, arrange_arrow, nrow=2, ncol=1, heights=c(15, 1))
    ggplot2::ggsave(arrange_results_1, filename=here::here("outputs", "Comparison", glue::glue("Results_1_{mortality}_{fecundity}.png")),
                    width=40, height=30, units="cm", dpi=300)
    ggplot2::ggsave(arrange_results_2, filename=here::here("outputs", "Comparison", glue::glue("Results_2_{mortality}_{fecundity}.png")),
                    width=40, height=30, units="cm", dpi=300)
    
    arrange_a <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                      nrow=2, ncol=4, align = "v")
    arrange_b <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation, plot_level_explanation, plot_level_explanation,
                           nrow=1, ncol=4, align="v")
    final_plot <- ggpubr::ggarrange(arrange_a, arrange_b, nrow=2, ncol=1, heights=c(15, 1))
    
    ggplot2::ggsave(final_plot, filename=here::here("outputs", "Comparison", glue::glue("Results_all_{mortality}_{fecundity}.png")),
                    width=50, height=50/2, units="cm", dpi=300)
    
    }
}

# # Compute multidimensional semivariance
# 
# Correlation_env_sp <- data.frame(
#   Mortality = numeric(),
#   Fecundity = numeric(),
#   Seed = numeric(),
#   Seed_r = numeric(),
#   Mod = character(),
#   Nb_obs = numeric(),
#   Correlation = numeric()
#   )
# 
# for (simu in c(1:nrow(Simulations))){
#   print(paste("Simu", simu))
#   
#   mortality <- Simulations[simu,1]
#   fecundity <- Simulations[simu,2]
#   seed <- Simulations[simu,3]
#   seed_r <- Simulations[simu,4]
#   
#   #Load the environment and species optima once per configuration (seed)
#   if(seed_r==1){
#     load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_sites.RData")))
#     load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
#   }
#   
#   for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
#     
#     if(mod == "Perf_know"){
#       
#       load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
#       sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
#       names(sp_XY) <- c("x", "y", "sp")
#       vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
#       
#       semivar_multidim <- compute_semivar_multidim(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end)
#       semivar_multidim$Vario_sp_geoR <- vario_sp$u
#       semivar_multidim$Distance <- vario_sp$bins.lim[-length(vario_sp$bins.lim)]
#       semivar_multidim$Sample_size <- vario_sp$n
#       
#       semivar_multidim<-semivar_multidim%>%
#         filter(Sample_size>500)
#       
#       m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)
#       
#       Correlation_env_sp_temp <- data.frame(
#         Mortality = mortality,
#         Fecundity = fecundity,
#         Seed = seed,
#         Seed_r = seed_r,
#         Mod = mod,
#         Nb_obs = NA,
#         Correlation = round(sqrt(summary(m)$r.squared), digits = 2)
#       )
#       
#       Correlation_env_sp <- rbind(Correlation_env_sp, Correlation_env_sp_temp)
#       
#     }else{
#       
#       for (nb_obs in c(0:15)){
#         
#         load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
#         sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
#         names(sp_XY) <- c("x", "y", "sp")
#         vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
#         
#         semivar_multidim <- compute_semivar_multidim(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end)
#         semivar_multidim$Vario_sp_geoR <- vario_sp$u
#         semivar_multidim$Distance <- vario_sp$bins.lim[-length(vario_sp$bins.lim)]
#         semivar_multidim$Sample_size <- vario_sp$n
#         
#         semivar_multidim<-semivar_multidim%>%
#           filter(Sample_size>500)
#         
#         m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)
#         
#         Correlation_env_sp_temp <- data.frame(
#           Mortality = mortality,
#           Fecundity = fecundity,
#           Seed = seed,
#           Seed_r = seed_r,
#           Mod = mod,
#           Nb_obs = nb_obs,
#           Correlation = round(sqrt(summary(m)$r.squared), digits = 2)
#         )
#         
#         Correlation_env_sp <- rbind(Correlation_env_sp, Correlation_env_sp_temp)
#         
#       }#for n_obs
#     }#else Part_know
#   }#for mod
# }#for simu
# 
# save(Correlation_env_sp, file=here::here("outputs", "Comparison", "Correlation_env_sp.RData"))
# 
# load(file=here::here("outputs", "Comparison", "Correlation_env_sp.RData"))
# 
# #One plot per Mortality * Fecundity option
# for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
#   for (fecundity in c("abund", "fixed")){
#     
#     data_figure <- Correlation_env_sp[which(Correlation_env_sp$Mortality==mortality&Correlation_env_sp$Fecundity==fecundity),]
# 
#     Summary_correlation_env_sp<-data_figure%>%
#       dplyr::group_by(Mod, Nb_obs)%>%
#       dplyr::mutate(Mean_corr=mean(Correlation, na.rm=TRUE), Sd=sd(Correlation, na.rm=TRUE))%>%
#       dplyr::slice(1)%>%
#       dplyr::ungroup()%>%
#       dplyr::select(-Correlation, -Seed, -Seed_r)
#     
#     save(Summary_correlation_env_sp, file=here::here("outputs", "Comparison", glue::glue("Mean_correlation_env_sp_{mortality}_{fecundity}.RData")))
#     
#     data_figure <- data_figure[data_figure$Mod!="Perf_know",]
#     
#     #Compute delta between with and without
#     data_figure <- data_figure%>%
#       dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
#       dplyr::mutate(Delta = Correlation - dplyr::lag(Correlation) , ID_delta=Nb_obs)%>%
#       dplyr::filter(is.na(Delta)==FALSE)%>%
#       dplyr::ungroup()%>%
#       dplyr::select(Seed, Seed_r, ID_delta, Delta)
#     
#     p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
#       ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
#       ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
#       ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
#       ggplot2::scale_colour_viridis_d()+
#       ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
#                     y = expression(paste(Delta, " environment-species correlation")))+
#       ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
#                      axis.text = ggplot2::element_text(size=7),
#                      legend.position = "none")
#     
#     ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_corr_env_sp_{mortality}_{fecundity}.png")),
#                     width=fig_width, height=fig_width/2, units="cm", dpi=300)
#   }
# }


## Understanding why delta performance is higher than 0 with the proportional mortality ##

Comparison_clandestines <- data.frame(
  Mortality=character(),
  Fecundity=character(),
  Mod=character(),
  Nb_obs=integer(),
  Seed=integer(),
  Seed_r=integer(),
  Abund_clandestines=numeric(),
  perf_clandestines=numeric(),
  perf_not_clandestines=numeric(),
  perf_tot=numeric())

for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){
    
    for(seed in 1:10){
      
      #same perf_E_Sp for all repetitions of the same configuration
      load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_perf_E_Sp.RData")))
      winner <- apply(X=perf_E_Sp, MARGIN=1, FUN=which.max)
      
      for(seed_r in 1:10){
        mod <- "Perf_know"
        load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        
        perf_present <- rep(NA, length(community_end))
        w0 <- (as.vector(t(community_end))==0)
        perf_present[w0] <- NA
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
        
        Comparison_clandestines_temp <- data.frame(
          Mortality=mortality,
          Fecundity=fecundity,
          Mod=mod,
          Nb_obs=NA,
          Seed=seed,
          Seed_r=seed_r,
          Abund_clandestines=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
          perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
          perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
          perf_tot=mean(perf_present, na.rm = TRUE))
        
        Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
        
        for (mod in c("Part_know", "Part_know_IV")){
          for (nb_obs in 1:n_axes){
            
            load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
            
            
            perf_present <- rep(NA, length(community_end))
            w0 <- (as.vector(t(community_end))==0)
            perf_present[w0] <- NA
            perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
            
            Comparison_clandestines_temp <- data.frame(
              Mortality=mortality,
              Fecundity=fecundity,
              Mod=mod,
              Nb_obs=nb_obs,
              Seed=seed,
              Seed_r=seed_r,
              Abund_clandestines=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
              perf_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
              perf_not_clandestines=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
              perf_tot=mean(perf_present, na.rm = TRUE))
            
            Comparison_clandestines <- rbind(Comparison_clandestines, Comparison_clandestines_temp)
          }
        }
      }
    }
  }
}

save(Comparison_clandestines, file = here::here("outputs", "Comparison", "Comparison_clandestines.RData"))

load(here::here("outputs", "Comparison", "Comparison_clandestines.RData"))


for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){
    
    data_figure <- Comparison_clandestines[which(Comparison_clandestines$Mortality==mortality&Comparison_clandestines$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta = Abund_clandestines - dplyr::lag(Abund_clandestines) , ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta)
    
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Sub-optimal species abundance")))+
      ggplot2::ggtitle(paste0("Mortality ", mortality, ", Fecundity ", fecundity))+
      ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                     axis.text = ggplot2::element_text(size=7),
                     legend.position = "none")
    
    ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abundance_clandestines_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}


mortality <- c("prop", "stocha")
fecundity <- "abund"

data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality%in%mortality&Perf_on_sites$Fecundity==fecundity),]

data_figure <- data_figure[data_figure$Mod!="Perf_know",]

#Compute delta between with and without
data_figure <- data_figure%>%
  dplyr::group_by(Mortality, Seed, Seed_r, Nb_obs)%>%
  dplyr::mutate(Delta = Perf - dplyr::lag(Perf) , ID_delta=Nb_obs)%>%
  dplyr::filter(is.na(Delta)==FALSE)%>%
  dplyr::ungroup()%>%
  dplyr::select(Mortality, Seed, Seed_r, ID_delta, Delta)

p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
  ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                    shape=as.factor(Mortality),
                                    group=Mortality),
                       alpha=0.6,
                       position=ggplot2::position_jitterdodge(jitter.width=0.5),
                       show.legend=F)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
  ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
  ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                y = expression(paste(Delta, " Theoretical performance")))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_performance_comparison_mortality.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)


mortality <- c("prop", "stocha")
fecundity <- "abund"

for (mod in c("Part_know", "Part_know_IV")){
  
  data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality%in%mortality&Perf_on_sites$Fecundity==fecundity),]
  
  data_figure <- data_figure[data_figure$Mod!="Perf_know",]
  
  data_figure <- data_figure[data_figure$Mod==mod,]
  
  p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Perf))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                      shape=as.factor(Mortality),
                                      group=Mortality),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5),
                         show.legend=F)+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                  y = "Theoretical performance")+
    ggplot2::ggtitle(paste0(mod))+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                   axis.text = ggplot2::element_text(size=7))
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Performance_comparison_mortality_{mod}.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}


mortality <- c("prop", "stocha")
fecundity <- "abund"

data_figure <- Comparison_clandestines[which(Comparison_clandestines$Mortality%in%mortality&Comparison_clandestines$Fecundity==fecundity),]

data_figure <- data_figure[data_figure$Mod!="Perf_know",]

#Compute delta between with and without
data_figure <- data_figure%>%
  dplyr::group_by(Mortality, Seed, Seed_r, Nb_obs)%>%
  dplyr::mutate(Delta = Abund_clandestines - dplyr::lag(Abund_clandestines) , ID_delta=Nb_obs)%>%
  dplyr::filter(is.na(Delta)==FALSE)%>%
  dplyr::ungroup()%>%
  dplyr::select(Mortality, Seed, Seed_r, ID_delta, Delta)

p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
  ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                    shape=as.factor(Mortality),
                                    group=Mortality),
                       alpha=0.6,
                       position=ggplot2::position_jitterdodge(jitter.width=0.5),
                       show.legend=F)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
  ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
  ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                y = expression(paste(Delta, " Sub-optimal species abundance")))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abundance_clandestines_comparison_mortality.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)


mortality <- c("prop", "stocha")
fecundity <- "abund"

for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
  
  data_figure <- Comparison_clandestines[which(Comparison_clandestines$Mortality%in%mortality&Comparison_clandestines$Fecundity==fecundity),]
  
  data_figure <- data_figure[data_figure$Mod==mod,]
  
  if(mod=="Perf_know"){
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Abund_clandestines))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                        shape=as.factor(Mortality),
                                        group=Mortality),
                           alpha=0.6,
                           position=ggplot2::position_jitterdodge(jitter.width=0.5),
                           show.legend=F)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
      ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
      ggplot2::labs(x="", y = "Sub-optimal species abundance")+
      ggplot2::ggtitle(paste0(mod))+
      ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                     axis.text = ggplot2::element_text(size=7),
                     axis.text.x=element_blank(), #remove x axis labels
                     axis.ticks.x=element_blank(),)
  }else{
  
  p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Abund_clandestines))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                      shape=as.factor(Mortality),
                                      group=Mortality),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5),
                         show.legend=F)+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = expression(paste("Number of observed axes ( ~ ", frac(sIV,uIV), " )")),
                  y = "Sub-optimal species abundance")+
    ggplot2::ggtitle(paste0(mod))+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                   axis.text = ggplot2::element_text(size=7))
  }
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Abund_clandestines_comparison_mortality_{mod}.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}


fecundity <- "abund"
for (mortality in c("prop", "stocha")){

  Perf_vs_clandestine_abundance <- Comparison_clandestines[which(Comparison_clandestines$Mortality==mortality&Comparison_clandestines$Fecundity==fecundity),]
  Perf_vs_clandestine_abundance <- Perf_vs_clandestine_abundance[Perf_vs_clandestine_abundance$Mod!="Perf_know",]
  
  Perf_vs_clandestine_abundance <- Perf_vs_clandestine_abundance%>%
    dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
    dplyr::mutate(Delta_perf = perf_tot - dplyr::lag(perf_tot) , ID_delta=Nb_obs)%>%
    dplyr::mutate(Delta_abund_clandestine = Abund_clandestines - dplyr::lag(Abund_clandestines) , ID_delta=Nb_obs)%>%
    dplyr::filter(is.na(Delta_perf)==FALSE, is.na(Delta_abund_clandestine)==FALSE)%>%
    dplyr::ungroup()%>%
    dplyr::select(Seed, Seed_r, ID_delta, Delta_perf, Delta_abund_clandestine)
  
  p <- ggplot2::ggplot(data=Perf_vs_clandestine_abundance, ggplot2::aes(x=Delta_abund_clandestine, y=Delta_perf))+
    ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
    ggplot2::geom_point(ggplot2::aes(colour=as.factor(ID_delta)), alpha=0.6)+
    ggplot2::scale_colour_viridis_d("Number of observed axes", option="inferno")+
    ggplot2::labs(x = expression(paste(Delta, "Sub-optimal species abundance")),
                  y = expression(paste(Delta, " Theoretical performance")))+
    ggplot2::ggtitle(paste0("Mortality ", mortality))+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                   axis.text = ggplot2::element_text(size=7))
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abund_clandestines_vs_delta_perf_{mortality}.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}


mortality <- "prop"
fecundity <- "abund"

Perf_vs_clandestine_abundance <- Comparison_clandestines[which(Comparison_clandestines$Mortality==mortality&Comparison_clandestines$Fecundity==fecundity),]
Perf_vs_clandestine_abundance <- Perf_vs_clandestine_abundance[Perf_vs_clandestine_abundance$Mod!="Perf_know",]

ggplot2::ggplot(data=Perf_vs_clandestine_abundance, ggplot2::aes(x=Abund_clandestines, y=perf_tot))+
  ggplot2::geom_point(ggplot2::aes(colour=as.factor(Nb_obs)), alpha=0.6)+
  ggplot2::scale_colour_viridis_d("Number of observed axes", option="inferno")+
  ggplot2::labs(x = "Abundance of sub-optimal species",
                y = "Theoretical performance")+
  ggplot2::ggtitle(paste0("Mortality ", mortality))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

ggplot2::ggplot(data=Perf_vs_clandestine_abundance, ggplot2::aes(x=Abund_clandestines, y=perf_tot))+
  ggplot2::geom_point(ggplot2::aes(colour=as.factor(Mod)), alpha=0.6)+
  ggplot2::scale_colour_manual("Model", values=c("darkblue","gold"))+
  ggplot2::labs(x = "Abundance of sub-optimal species",
                y = "Theoretical performance")+
  ggplot2::ggtitle(paste0("Mortality ", mortality))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

