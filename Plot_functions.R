plot_environment <- function(model, fig_width, n_axis, env){
  
  png(file=here::here("outputs", model, "environment.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
  
  par(mfrow=c(2,2), bty = "n")
  
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
  scatter3D(niche_optimum$sp_x, niche_optimum$sp_y, niche_optimum$sp_z,
            pch=16, 
            colvar=1:nsp, col=viridis(nsp),
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
  plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
  dev.off()
}

plot_community_start <- function(model, fig_width, community, nsp){
  png(file=here::here("outputs", model, "community_start.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  plot(raster::raster(community), main="Species - Start", zlim=c(0, nsp),
       col=c("black",  viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_mortality_events <- function(model, fig_width, community, nsp){
  png(file=here::here("outputs", model, "mortality_events.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  plot(raster::raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
       col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_community_end <- function(model, fig_width, community, nsp){
png(file=here::here("outputs", model, "community_end.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(bty = "n")
plot(raster::raster(community), main=glue::glue("Species - End (ngen={ngen})"),
     zlim=c(0, nsp), col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
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

plot_mean_rank_hab_freq <- function(sp_mean_rank, sp_hab_freq, model, fig_width){
  df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=sp_hab_freq, y=sp_mean_rank)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
    ggplot2::xlab("Species suitable habitat frequency") +
    ggplot2::ylab("Species mean rank (higher rank = lower abundance)") +
    ggplot2::theme(axis.title=ggplot2::element_text(size=16))
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
       col=c("black", viridis(nsp)), legend=FALSE)
  plot(raster::raster(community_end), main="Species - End", zlim=c(0, nsp),
       col=c("black", viridis(nsp)), legend=FALSE)
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

plot_random_species_niche <- function(seed, df_perf, model, fig_width){ 
  # Select 8 species at random
  set.seed(seed)
  sp_sel <- sample(unique(df_perf$Species), 9, replace=FALSE)
  df_sp_sel <- df_perf %>% filter(Species %in% sp_sel)
  p <- ggplot2::ggplot(data=df_sp_sel, ggplot2::aes(x=Env_1, y=Perf)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE) +
    ggplot2::facet_wrap(vars(Species), nrow=3) +
    ggplot2::xlab("Environment (first axis)") +
    ggplot2::ylab("Performance")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "infering_species_niche.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_IV <- function(V_intra, model, fig_width){
  p <- ggplot2::ggplot(data=V_intra, ggplot2::aes(x=Species, y=V)) +
    ggplot2::geom_col() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, size=6),
          text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab("Intraspecific variance")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "intraspecific_variance.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}


plot_inferred_perf_environment <- function(env, lm_fit, nsp, model, fig_width){

  X1 <- env[[1]]
  
  beta_0 <- as.vector(lm_fit$coefficients[1:nsp])
  beta_1 <- as.vector(c(lm_fit$coefficients[nsp+1], lm_fit$coefficients[(nsp+3):(2*nsp+1)]))
  beta_2 <- as.vector(c(lm_fit$coefficients[nsp+2], lm_fit$coefficients[(2*nsp+2):(3*nsp)]))
  
  E_seq <- seq(min(X1), max(X1), length.out=100)
  E_seq_mat <- matrix(rep(E_seq, nsp), ncol=nsp)
  
  beta_0_mat_E_seq <- matrix(rep(beta_0,each=length(E_seq)), ncol=nsp)
  beta_1_mat_E_seq <- matrix(rep(beta_1,each=length(E_seq)), ncol=nsp)
  beta_2_mat_E_seq <- matrix(rep(beta_2,each=length(E_seq)), ncol=nsp)
  
  Mat_perf_inferred <- beta_0_mat_E_seq+beta_1_mat_E_seq*E_seq_mat+beta_2_mat_E_seq*E_seq_mat^2
  Mat_perf_inferred <- as.data.frame(cbind(E_seq, Mat_perf_inferred))
  
  colnames(Mat_perf_inferred) <- c("E", paste0("Sp",1:nsp))
  
  Mat_perf_inferred_long <- Mat_perf_inferred %>%
    tidyr::pivot_longer(cols=2:(nsp+1), names_to="Sp",
                 names_prefix="Sp", values_to="Perf")
  
  p <- ggplot2::ggplot(data=Mat_perf_inferred_long, ggplot2::aes(x=E, y=Perf, col=Sp)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_viridis_d()+
    ggplot2::xlab("Environment") + 
    ggplot2::ylab("Species performance")+
    ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = 20))
  ggplot2::ggsave(p, filename=here::here("outputs", model, "inferred_perf_environment.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  x <- seq(-3, 3, 0.01)
  mean_x1 <- mean(X1)
  n_ind_simul <- length(x)
  df_perf_IV <- data.frame(Species = rep(1:nsp, each=n_ind_simul), X = rep(x, nsp), Mean=rep(beta_0+beta_1*mean_x1+beta_2*(mean_x1^2), each=n_ind_simul))
  df_perf_IV <- df_perf_IV%>%
    dplyr::group_by(Species)%>%
    dplyr::mutate(Density = dnorm(x, mean=Mean, sd=sqrt(V_intra$V[Species])))%>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(df_perf_IV, ggplot2::aes(x = X, y = Density, colour = as.factor(Species))) +
    ggplot2::geom_line()+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::theme(legend.position = "none")+
    ggplot2::labs(x = "Performance at mean environment variable X1",
         y = "Density")
  ggplot2::ggsave(p, filename=here::here("outputs", model, "Perf_overlap_IV.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
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
  
  p <- ggplot2::ggplot(data=Abund_plot, ggplot2::aes(x=Gen, y=Mean_sp_gen, colour=Sp))+
    ggplot2::geom_line()+
    #ggplot2::geom_ribbon(ggplot2::aes(y=Mean_sp_gen, ymin=CI_025, ymax=CI_975, fill=Sp), alpha=0.5)+
    ggplot2::scale_colour_viridis_d()+
    ggplot2::labs(x="Generations", y="Species abundance")+
    ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(p, filename=here::here("outputs", model, "Abundances.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}