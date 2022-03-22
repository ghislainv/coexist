
perf_E_Sp
V_intra
mean(perf_E_Sp)

epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
epsilon_mat <- matrix(epsilon, ncol=nsp)
apply(epsilon_mat, 2, var) - V_intra$V
perf_E_Sp_bis <- perf_E_Sp +  epsilon_mat

b <- -0.5
mortality_E_Sp_bis <- inv_logit(logit(0.1) + b * perf_E_Sp_bis)
mean(mortality_E_Sp_bis)
#The mean is higher due to the back-transformation of the inverse logit

#Correction

Mort_diff <- mean(mortality_E_Sp_bis) - mean(mortality_E_Sp)

mortality_E_Sp_bis_corrected <- mortality_E_Sp_bis - Mort_diff

mean(mortality_E_Sp_bis_corrected)
