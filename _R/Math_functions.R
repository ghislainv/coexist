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

# Function to identify the species with the highest performance
high_perf_sp <- function(dist, sp_pres) {
  dist_pres <- dist[sp_pres]
  min_dist <- min(dist_pres)
  sp_high_perf <- sp_pres[which(dist_pres==min_dist)]
  # If more than one species, selection at random
  if (length(sp_high_perf)>1) {
    # Random permutation
    sp_high_perf <- sample(sp_high_perf)
    sp_high_perf <- sp_high_perf[1]
  }
  return(sp_high_perf)
}

# Other function to identify the species with the highest performance
# For use in plyr::ddply
f_sp <- function(x) {
  w <- which(x[,3]==max(x[,3]))
  return(x[sample(w, 1),1])
}

# Function to compute a multidimensional semivariance
compute_semivar_multidim <- function(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end){
  
  ## COMPUTE DISTANCES
  
  # Compute the sum of squared differences of environmental conditions
  dist_env <- as.vector(dist(sites[, 1]))^2
  for (k in 2:n_axes){
    dist_env <- dist_env + as.vector(dist(sites[, k]))^2
  }
  
  # Compute the sum of squared differences of species optima
  dist_sp <- as.vector(dist(niche_optimum[,1]))^2
  for (k in 2:n_axes){
    dist_sp <- dist_sp + as.vector(dist(niche_optimum[, k]))^2
  }
  
  # Compute the geographical distance between all cells
  dist_spatial <- as.vector(dist(sp_XY[,-3]))
  
  ## SELECT SPECIES COUPLES
  
  # Associate the species couple to each species distance
  Dist_sp <- data.frame(Sp=t(combn(1:nsp, 2)), Dist_sp=dist_sp)
  # Add diagonal without computing it
  Dist_sp <- rbind(Dist_sp, data.frame(Sp.1=1:nsp, Sp.2=1:nsp, Dist_sp=rep(0, nsp)))
  
  # Associate each geographical distance to the cell number of the pair of cells
  Dist_neighbours_cells <- data.frame(Neighbour=t(combn(1:nrow(sites), 2)), Dist_spatial=dist_spatial)
  
  # List of the species and the cells on which they are present
  List_present_species <- data.frame(Cell=1:nrow(sites), Sp=as.vector(raster::raster(community_end)))
  
  # Associate the species that are present and their cells
  Dist_neighbours_cells <- Dist_neighbours_cells%>%
    dplyr::mutate(Cell=Neighbour.1)%>%
    dplyr::inner_join(List_present_species, by="Cell")%>%
    dplyr::rename(Sp1=Sp)%>%
    dplyr::mutate(Cell=Neighbour.2)%>%
    dplyr::inner_join(List_present_species, by="Cell")%>%
    dplyr::rename(Sp2=Sp)%>%
    dplyr::select(- Cell)%>%
    dplyr::mutate(Sp1_temp = dplyr::case_when(Sp1>Sp2~Sp2, Sp1<=Sp2~Sp1),
                  Sp2_temp = dplyr::case_when(Sp1>Sp2~Sp1, Sp1<=Sp2~Sp2),
                  Neighbour.1_temp = dplyr::case_when(Sp1>Sp2~Neighbour.2, Sp1<=Sp2~Neighbour.1),
                  Neighbour.2_temp = dplyr::case_when(Sp1>Sp2~Neighbour.1, Sp1<=Sp2~Neighbour.2))%>%
    dplyr::select(-Sp1, -Sp2, -Neighbour.1, -Neighbour.2)%>%
    dplyr::rename(Neighbour.1=Neighbour.1_temp, Neighbour.2=Neighbour.2_temp, Sp1=Sp1_temp, Sp2=Sp2_temp)%>%
    dplyr::mutate(Diff_sp = dplyr::case_when(Sp1==Sp2~0, Sp1!=Sp2~1))%>%
    dplyr::mutate(Sp.1=Sp1, Sp.2=Sp2)%>%
    dplyr::inner_join(Dist_sp, by=c("Sp.1", "Sp.2"))%>%
    dplyr::select(-Sp.1, -Sp.2)
  
  # Associate each geographical distance to an environmental distance
  Dist_env_cells <- data.frame(Env=dist_env, Dist_spatial=dist_spatial)
  
  # Using the distance bins made by geoR, compute the environmental semivariance of each bin
  Vario_env <- c()
  Vario_sp <- c()
  Vario_sp_0_1 <- c()
  
  for(k in 1:(length(vario_sp$bins.lim)-1)){
    Dist_bin <- Dist_env_cells%>%
      dplyr::filter(Dist_spatial >= vario_sp$bins.lim[k] & Dist_spatial < vario_sp$bins.lim[k+1])
    Vario_env[k] <- sum(Dist_bin$Env)/(2*nrow(Dist_bin))
    
    Dist_bin <- Dist_neighbours_cells%>%
      dplyr::filter(Dist_spatial >= vario_sp$bins.lim[k] & Dist_spatial < vario_sp$bins.lim[k+1])
    Vario_sp[k] <- sum(Dist_bin$Dist_sp)/(2*nrow(Dist_bin))
    
    Dist_bin <- Dist_neighbours_cells%>%
      dplyr::filter(Dist_spatial >= vario_sp[["bins.lim"]][k] & Dist_spatial < vario_sp$bins.lim[k+1])
    Vario_sp_0_1[k] <- sum(Dist_bin$Diff_sp)/(2*nrow(Dist_bin))
  }
  
  return(data.frame(Vario_sp=Vario_sp, Vario_sp_0_1=Vario_sp_0_1, Vario_env=Vario_env))
}

#Compute Jaccard similarity index
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#Compute percentage similarity
percentage_similarity <- function(a, b) {
  A <- sum(a)
  B <- sum(b)
  W <- sum(pmin(a, b))
  return((2*W)/(A+B))
}