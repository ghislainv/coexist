#A classical Lotka-Volterra competition model is :
# ${\frac{dx_{i}}{dt}}=r_{i}x_{i}\left(1-{\frac{\sum_{{j=1}}^{N}\alpha_{{ij}}x_{j}}{K_{i}}}\right)$

#Here we want to estimate the alpha matrix i.e. to solve the Ordinary Differential Equation system.

#Build abundance matrix and compute abundance delta
load(here::here("outputs", "m0", "Abundances_m0.RData"))
nsp <- ncol(Abundances_m0[[1]])
r = 1
for (r in 1:nrep){
  Abund_matrix_m0_t0 <- as.data.frame(Abundances_m0[[r]])
  colnames(Abund_matrix_m0_t0) <- sprintf("Sp_%03d", 1:nsp)
  Abund_matrix_m0_t1 <- rbind(Abund_matrix_m0_t0[2:nrow(Abund_matrix_m0_t0),], NA)
  rownames(Abund_matrix_m0_t1) <- c(1:nrow(Abund_matrix_m0_t1))
  Abund_diff_matrix_m0 <- Abund_matrix_m0_t1 - Abund_matrix_m0_t0
  colnames(Abund_diff_matrix_m0) <- sprintf("Diff_Sp_%03d", 1:nsp)
}

########################################################################################
########################################## lm ##########################################
########################################################################################

data_sp1 <- data.frame(DiffAbund_sp = Abund_diff_matrix_m0[1:500,1],
                       Abund_sp = Abund_matrix_m0_t0[1:500, 1])
data_sp1 <- cbind(data_sp1, Abund_matrix_m0_t0[1:500,])
data_sp1[,3:ncol(data_sp1)] <- data_sp1[,3:ncol(data_sp1)] * data_sp1$Abund_sp
colnames(data_sp1)[3:ncol(data_sp1)] <- sprintf("Abund_sp%03d", 1:nsp)

Lotka_Volterra_lm <- lm(DiffAbund_sp ~ 0 + Abund_sp + Abund_sp001 + Abund_sp002 + Abund_sp003 
                     + Abund_sp004 + Abund_sp005 + Abund_sp006 + Abund_sp007 + Abund_sp008 
                     + Abund_sp009 + Abund_sp010 + Abund_sp011 + Abund_sp012 + Abund_sp013
                     + Abund_sp014 + Abund_sp015 + Abund_sp016 + Abund_sp017 + Abund_sp018
                     + Abund_sp019 + Abund_sp020 + Abund_sp021 + Abund_sp022 + Abund_sp023
                     + Abund_sp024 + Abund_sp025 + Abund_sp026 + Abund_sp027 + Abund_sp028
                     + Abund_sp029 + Abund_sp030 + Abund_sp031 + Abund_sp032 + Abund_sp033
                     + Abund_sp034 + Abund_sp035 + Abund_sp036 + Abund_sp037 + Abund_sp038
                     + Abund_sp039 + Abund_sp040 + Abund_sp041 + Abund_sp042 + Abund_sp043
                     + Abund_sp044 + Abund_sp045 + Abund_sp046 + Abund_sp047 + Abund_sp048
                     + Abund_sp049 + Abund_sp050 + Abund_sp051 + Abund_sp052 + Abund_sp053
                     + Abund_sp054 + Abund_sp055 + Abund_sp056 + Abund_sp057 + Abund_sp058
                     + Abund_sp059 + Abund_sp060 + Abund_sp061 + Abund_sp062 + Abund_sp063
                     + Abund_sp064, data=data_sp1)

########################################################################################
######################################### brms #########################################
########################################################################################

Lotka_Volterra_brms <- brms::brm(
  formula = DiffAbund_sp ~ 0 + Abund_sp + Abund_sp001 + Abund_sp002 + Abund_sp003 
                     + Abund_sp004 + Abund_sp005 + Abund_sp006 + Abund_sp007 + Abund_sp008 
                     + Abund_sp009 + Abund_sp010 + Abund_sp011 + Abund_sp012 + Abund_sp013
                     + Abund_sp014 + Abund_sp015 + Abund_sp016 + Abund_sp017 + Abund_sp018
                     + Abund_sp019 + Abund_sp020 + Abund_sp021 + Abund_sp022 + Abund_sp023
                     + Abund_sp024 + Abund_sp025 + Abund_sp026 + Abund_sp027 + Abund_sp028
                     + Abund_sp029 + Abund_sp030 + Abund_sp031 + Abund_sp032 + Abund_sp033
                     + Abund_sp034 + Abund_sp035 + Abund_sp036 + Abund_sp037 + Abund_sp038
                     + Abund_sp039 + Abund_sp040 + Abund_sp041 + Abund_sp042 + Abund_sp043
                     + Abund_sp044 + Abund_sp045 + Abund_sp046 + Abund_sp047 + Abund_sp048
                     + Abund_sp049 + Abund_sp050 + Abund_sp051 + Abund_sp052 + Abund_sp053
                     + Abund_sp054 + Abund_sp055 + Abund_sp056 + Abund_sp057 + Abund_sp058
                     + Abund_sp059 + Abund_sp060 + Abund_sp061 + Abund_sp062 + Abund_sp063
                     + Abund_sp064,
  data=data_sp1,
  prior = brms::prior(normal(0, 10^6), class = "b"),
  iter = 10000 ,
  warmup = 5000,
  thin = 5,
  chains = 2,
  cores=2,
  control = list(max_treedepth = 11) 
)

summary(Lotka_Volterra_brms)

brms::mcmc_plot(Lotka_Volterra_brms, type="trace")
brms::mcmc_plot(Lotka_Volterra_brms, type="dens")

# #Tests to see if their is a problem with the NA
# data_sp1 <- data.frame(DiffAbund_sp = sample(-2:0, nrow(data_sp1), replace = T),
# Abund_sp = sample(0:15, ncol(data_sp1), replace = T))
# data_sp1 <- cbind(data_sp1, Abund_matrix_m0_t0[1:500,])
# for (k in 3:ncol(data_sp1)) {
#   data_sp1[,k] <- sample(0:200, nrow(data_sp1), replace = T)
# }
# --> the problem is really in the data
# Probably a colinearity problem

########################################################################################
######################################### GJAM #########################################
########################################################################################

library(gjam)
library(devtools)
d <- "https://github.com/jimclarkatduke/gjam/blob/master/gjamTimeFunctions.R?raw=True"
source_url(d)


#Simulated data to see what the data should look like#


S     <- 6        # no. species
nsite <- 1       # no. time series
ntime <- 100      # mean no. of time steps in a series
obsEffort <- 1    # full census
termB <- FALSE #No environmental effect
termR <- TRUE
termA <- TRUE
predPrey = NULL
zeroAlpha = NULL
PLOT = TRUE
seed <- 999
set.seed( seed )

#Tried to make the model run without termR
# But the problem persists later
#gam <- runif(S, .01, .1)
#rhoTrue <- gam

tmp <- gjamSimTime(S, Q = 0, nsite, ntime, termB, termR, termA, 
                      obsEffort = obsEffort, PLOT = T)
xdata  <- tmp$xdata
ydata  <- tmp$ydata
edata  <- tmp$edata
groups <- tmp$groups
times  <- tmp$times
trueValues <- tmp$trueValues
formula    <- tmp$formula

Data_before_fill <- tmp

# Fill missing times (here only the 0 observation)

timeCol   <- 'times'
groupCol  <- 'groups'
groupVars <- c( 'groups' )

tmp <- gjamFillMissingTimes(xdata, ydata, edata, groupCol, timeCol,
                            FILLMEANS = T, groupVars = groupVars,
                            typeNames = 'DA', missingEffort = .1)    
xdata  <- tmp$xdata
ydata  <- tmp$ydata
edata  <- tmp$edata
tlist  <- tmp$timeList
snames <- colnames(ydata)

effort <- list(columns = 1:S, values = edata)

#Model

rhoPrior  <- list(lo = list(intercept = -.3), 
                  hi = list(intercept = .3) ) 

priorList <- list( formulaRho = as.formula(~ 1), rhoPrior = rhoPrior)
tmp <- gjamTimePrior( xdata, ydata, edata, priorList)

timeList <- mergeList(tlist, tmp)

modelList <- list(typeNames = 'DA', ng = 4000, burnin = 1000,  
                  timeList = timeList, effort = effort) 

outputAR <- gjam(formula, xdata=xdata, ydata=ydata, modelList=modelList)

plotPars  <- list(PLOTALLY=T, trueValues = trueValues, 
                  SAVEPLOTS = T, outFolder = outFolder)
gjamPlot(outputAR, plotPars)


# Using my Data #
load(here::here("outputs", "m0", "Abundances_m0.RData"))
r = 1 #choose which repetition of the dynamics we will use
ydata <- as.data.frame(Abundances_m0[[r]])

S = ncol(ydata)
  
colnames(ydata) <- sprintf("sp%d", 1:ncol(ydata))
xdata <- data.frame(groups = rep(1, nrow(ydata)),
                   times = c(1:nrow(ydata)),
                   intercept = rep(1, nrow(ydata)))
#groups : number of environments, here only one
#times : number of iterations, here 501

#environmental data
edata <-as.data.frame(matrix(1, nrow=nrow(ydata), ncol=ncol(ydata)))
colnames(edata) <- colnames(ydata)

timeCol   <- 'times'
groupCol  <- 'groups'
groupVars <- c( 'groups' )
tmp <- gjamFillMissingTimes(xdata, ydata, edata, groupCol, timeCol,
                            FILLMEANS = T, groupVars = groupVars,
                            typeNames = 'DA', missingEffort = .1)    
xdata  <- tmp$xdata
ydata  <- tmp$ydata
edata  <- tmp$edata
tlist  <- tmp$timeList
snames <- colnames(ydata)

effort <- list(columns = 1:S, values = edata)

rhoPrior  <- list(lo = list(intercept = -.3), 
                  hi = list(intercept = .3) ) 

priorList <- list( formulaRho = as.formula(~ 1), rhoPrior = rhoPrior)

tmp <- gjamTimePrior( xdata, ydata, edata, priorList)

timeList <- mergeList(tlist, tmp)

modelList <- list(typeNames = 'DA', ng = 4000, burnin = 1000,  
                  timeList = timeList, effort = effort)

outputAR <- gjam(formula, xdata=xdata, ydata=ydata, modelList=modelList)

outFolder <- 'gjamOutputAR'
plotPars  <- list(PLOTALLY=T, trueValues = trueValues, 
                  SAVEPLOTS = T, outFolder = outFolder)
gjamPlot(outputAR, plotPars)
save(outputAR, file = paste( outFolder, '/output.rdata', sep=''))

#########################################################################"

S <- ncol(ydata)

effMat <- tmp$edata

alpha <- matrix(NA,S,S)

wstar <- apply(ydata/effMat,2,max,na.rm=T)  # carrying capacity based on observed
wstar <- wstar*10                           # assume it is higher

astar <- (1 - 1.2)/wstar                             # reasonable lambda

z <- sqrt( crossprod( as.matrix(ydata) ) )
alpha[z > 1000] <- 2*matrix(astar,S,S,byrow=T)[z > 1000]
diag(alpha) <- astar

loAlpha <- alpha
hiAlpha <- loAlpha*0
hiAlpha[3,2] <- 1          # a predator and prey
loAlpha[2,3] <- -1

alphaPrior <- list(lo = loAlpha, hi = hiAlpha)

timeList <- list(times = 'times', rowInserts = rowInserts,
                 alphaPrior = alphaPrior, betaPrior = betaPrior, 
                 lambdaPrior = lambdaPrior )
effort <- list(columns = 1:S, values = effMat)

rl <- list(N = 8, r = 5)
modelList <- list(typeNames = 'DA',ng=2000, burnin=500, reductList = rl, 
                  timeList=timeList, effort = effort)

output <- gjam(formula, xdata=xdata, ydata=ydata, modelList=modelList)

#####################################################################
############################### Stan ################################
#####################################################################

library(rstan)

load(here::here("outputs", "m0", "Abundances_m0.RData"))
abundances_increments <- c(as.matrix(na.omit(Abund_diff_matrix_m0)))
#abundances_increments <- as.matrix(na.omit(Abund_diff_matrix_m0))
#One must remove last abundances measures which are not associated to an abundance increment
abundance_matrix <- as.matrix(Abundances_m0[[1]])[1:(nrow(Abund_diff_matrix_m0)-1),]
abundance_vector <- c(abundance_matrix)
abundance_matrix_rep <- do.call(rbind, replicate(nsp, abundance_matrix, simplify=FALSE))
nsp <- ncol(abundance_matrix)
#species <- (1:nsp)
species <- rep(1:nsp, each = nrow(Abund_diff_matrix_m0)-1)
nobs <- length(abundance_vector)

data <- list(N = nobs,
             K = nsp,
             J = nsp,
             L = 1,
             jj = species,
             x = abundance_matrix_rep,
             u = as.matrix(rep(1, nsp)),
             y = abundances_increments)

fit <- stan(file = here::here('Model.stan'), data = data)

data_chol <- list(N = nobs,
             K = nsp,
             J = nsp,
             L = 1,
             jj = species,
             x = abundance_matrix,
             u = t(as.matrix(rep(1, nsp))),
             y = abundances_increments)

fit <- stan(file = here::here('Model_chol.stan'), data = data_chol)

data_Jeanne <- list(N = nobs,
             K = nsp,
             J = nsp,
             jj = species,
             x = abundance_matrix_rep,
             y = abundances_increments)

options(mc.cores = 2)
fit <- stan(file = here::here('Stan_Jeanne.stan'), data = data_Jeanne, chains=2)

save(fit, file=here::here("MCMC_Stan_Jeanne.RData"))


