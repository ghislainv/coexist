plot_environment <- function(model, fig_width, n_axis, env){
  
  png(file=here("outputs", model, "environment.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
  
  par(mfrow=c(2,2), bty = "n")
  
  for(k in 1:n_axis){
    plot(raster(env[[k]]), main=glue::glue("Environment var{k}"), col=topo.colors(255), cex.main=1.2)
  }
  if(grepl(pattern="Perf_know", model)==TRUE){
    #RGB environment
    env_stack <- stack(raster(env[[1]]*255), raster(env[[2]]*255), raster(env[[3]]*255))
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
  png(file=here("outputs", model, "species_niche.png"),
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
  png(file=here("outputs", model, "hist_mortality.png"),
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
  png(file=here("outputs", model, "species_habitat_freq.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(cex.lab=1.5)
  plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
  dev.off()
}

plot_community_start <- function(model, fig_width, community, nsp){
  png(file=here("outputs", model, "community_start.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  plot(raster(community), main="Species - Start", zlim=c(0, nsp),
       col=c("black",  viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_mortality_events <- function(model, fig_width, community, nsp){
  png(file=here("outputs", model, "mortality_events.png"),
      width=fig_width, height=fig_width, units="cm", res=300)
  par(bty = "n")
  plot(raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
       col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
  dev.off()
}

plot_community_end <- function(model, fig_width, community, nsp){
png(file=here("outputs", model, "community_end.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(bty = "n")
plot(raster(community), main=glue("Species - End (ngen={ngen})"),
     zlim=c(0, nsp), col=c("black", viridis(nsp)), legend=FALSE, cex.main=2, cex.axis=1.5)
dev.off()
}

plot_species_richness <- function(nrep, sp_rich, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  sp_rich_long <- sp_rich %>%
    mutate(gen=1:(ngen)) %>%
    pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="sp_rich")
  p <- ggplot(data=sp_rich_long, aes(x=gen, y=sp_rich, col=rep)) +
    geom_line() +
    scale_colour_viridis_d()+
    xlab("Generations") + 
    ylab("Species richness")+
    theme(legend.position = "none",
          text = element_text(size = 20))+
    ylim(0,nsp)
  ggsave(p, filename=here("outputs", model, "species_richness_with_time.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_mean_rank_hab_freq <- function(sp_mean_rank, sp_hab_freq, model, fig_width){
  df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
  p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
    geom_point() +
    geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
    xlab("Species suitable habitat frequency") +
    ylab("Species mean rank (higher rank = lower abundance)") +
    theme(axis.title=element_text(size=16))
  ggsave(p, filename=here("outputs", model, "mean_rank-habitat_freq.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_env_filt<- function(nrep, env_filt, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  env_filt_long <- env_filt %>%
    mutate(gen=1:(ngen)) %>%
    pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="env_filt")
  
  p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
    geom_line() +
    scale_colour_viridis_d()+
    labs(title="Environmental filtering") +
    xlab("Generations") + 
    ylab("Mean env-species perf difference")
  ggsave(p, filename=here("outputs", model, "environmental_filtering.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}
  
plot_env_species<- function(model, fig_width, community_start, community_end, env){
  png(file=here("outputs", model, "spatial_comp_env_sp.png"), 
      width=fig_width, height=fig_width, units="cm", res=300)
  par(mfrow=c(2,2), bty = "n")
  plot(raster(community_start), main="Species - Start", zlim=c(0, nsp),
       col=c("black", viridis(nsp)), legend=FALSE)
  plot(raster(community_end), main="Species - End", zlim=c(0, nsp),
       col=c("black", viridis(nsp)), legend=FALSE)
  env_stack <- stack(raster(env[[1]]*255), raster(env[[2]]*255), raster(env[[3]]*255))
  crs(env_stack) <- "+proj=utm +zone=1"
  plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
  dev.off()
}

plot_theta_community<-function(theta_comm, ngen, model, fig_width){
  colnames_long <- paste0("X", 1:nrep)
  
  theta_comm_long <- theta_comm %>%
    mutate(gen=1:(ngen)) %>%
    pivot_longer(cols=colnames_long, names_to="rep",
                 names_prefix="X", values_to="theta_comm")
  p <- ggplot(data=theta_comm_long, aes(x=gen, y=theta_comm, col=rep)) +
    geom_line() +
    scale_colour_viridis_d()+
    labs(title="Environmental filtering") +
    xlab("Generations") + 
    ylab("Mean mortality rate in the community")
  ggsave(p, filename=here("outputs", model, "mortality_rate_community.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

plot_spatial_autocorr <- function(community_end, sites, niche_width, model, fig_width){
  # Species autocorrelation
  sp_XY <- data.frame(rasterToPoints(raster(community_end)))
  names(sp_XY) <- c("x", "y", "sp")
  vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
  # Environment autocorrelation
  # 3D voxel for each site
  x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
  y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
  z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
  n_niche <- 1/niche_width
  class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
  vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
  # Plot with correlation
  png(file=here("outputs", model, "sp_autocorrelation.png"),
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
  p <- ggplot(data=df_sp_sel, aes(x=Env_1, y=Perf)) +
    geom_point() +
    geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE) +
    facet_wrap(vars(Species), nrow=3) +
    xlab("Environment (first axis)") +
    ylab("Performance")
  ggsave(p, filename=here("outputs", model, "infering_species_niche.png"),
         width=fig_width, height=fig_width, units="cm", dpi=300)
}

plot_IV <- function(V_intra, model, fig_width){
  p <- ggplot(data=V_intra, aes(x=Species, y=V)) +
    geom_col() +
    theme(axis.text.x=element_text(angle=90, size=6),
          text = element_text(size = 20)) +
    ylab("Intraspecific variance")
  ggsave(p, filename=here("outputs", model, "intraspecific_variance.png"),
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
    pivot_longer(cols=2:(nsp+1), names_to="Sp",
                 names_prefix="Sp", values_to="Perf")
  
  p <- ggplot(data=Mat_perf_inferred_long, aes(x=E, y=Perf, col=Sp)) +
    geom_line() +
    scale_colour_viridis_d()+
    xlab("Environment") + 
    ylab("Species performance")+
    theme(legend.position = "none", text = element_text(size = 20))
  ggsave(p, filename=here("outputs", model, "inferred_perf_environment.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
  
  x <- seq(-3, 3, 0.01)
  mean_x1 <- mean(X1)
  n_ind_simul <- length(x)
  df_perf_IV <- data.frame(Species = rep(1:nsp, each=n_ind_simul), X = rep(x, nsp), Mean=rep(beta_0+beta_1*mean_x1+beta_2*(mean_x1^2), each=n_ind_simul))
  df_perf_IV <- df_perf_IV%>%
    group_by(Species)%>%
    mutate(Density = dnorm(x, mean=Mean, sd=sqrt(V_intra$V[Species])))%>%
    ungroup()
  
  p <- ggplot(df_perf_IV, aes(x = X, y = Density, colour = as.factor(Species))) +
    geom_line()+
    scale_colour_viridis_d()+
    theme(legend.position = "none")+
    labs(x = "Performance at mean environment variable X1",
         y = "Density")
  ggsave(p, filename=here("outputs", model, "Perf_overlap_IV.png"),
         width=fig_width, height=fig_width/2, units="cm", dpi=300)
}