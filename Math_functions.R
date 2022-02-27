# logit/inv_logit functions
logit <- function(x, min=0, max=1) {
  p <- (x-min)/(max-min)
  return(log(p/(1-p)))
}

inv_logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
  return(p * (max-min) + min)
}

# Cpp function to compute distance between Sites and Species
Rcpp::sourceCpp(here::here("_src", "dist_Site_Sp.cpp"))

# Function to draw in a multivariate normal
rmvn <- function(n, mu=0, V=matrix(1), seed=1234) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) {
    stop("Dimension problem!")
  }
  D <- chol(V)
  set.seed(seed)
  t(matrix(rnorm(n*p),ncol=p)%*%D+rep(mu,rep(n,p)))
}

# Function to rescale between 0 and 255 for RGB plots
range_0_255 <- function(x, newMax=255, newMin=0){
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin }

# Function to combine the three RGB values in a synthetic value
entrelac <- function(r, g, b){
  acc <- 1
  total <- 0
  for (i in 0:7){
    total <- total + (bitwAnd(bitwShiftR(b, i), 1))*acc
    acc <- acc * 2
    total <- total + (bitwAnd(bitwShiftR(g, i), 1))*acc
    acc <- acc * 2
    total <- total + (bitwAnd(bitwShiftR(r, i),1))*acc
    acc <- acc * 2
  }
  return (total)
}

# Function to compute a multidimensional semivariance
semivar_mutlidim <- function(sites, n_axis, sp_XY, vario_sp){
  #Compute the sum of squared differences 
  Dist_env <- as.vector(dist(sites[, 1]))^2
  for (k in 2:n_axis){
    Dist_env <- Dist_env + as.vector(dist(sites[, k]))^2
  }
  #Compute the geographical distance between all cells
  Dist_cells <- as.vector(dist(sp_XY[,-3]))
  
  Dist_env_cells <- data.frame(Env=Dist_env, Cells=Dist_cells)
  
  # Using the distance bins made by geoR, compute the environmental semivariance of each bin
  Vario_env <- c()
  
  for(k in 1:(length(vario_sp[["bins.lim"]])-1)){
    Dist_bin <- Dist_env_cells%>%
      dplyr::filter(Cells >= vario_sp[["bins.lim"]][k] & Cells < vario_sp[["bins.lim"]][k+1])
    Vario_env[k] <- sum(Dist_bin$Env)/(2*nrow(Dist_bin))
  }
  
  return(Vario_env)
}