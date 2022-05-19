#### Create an array where each line correspond to one launch/simulation set-up ####

Simulations <- cbind(
  c(rep("fixed", 200), rep("prop", 200), rep("stocha", 200), rep("stocha_basal", 200)),
  rep(c(rep("abund", 100), rep("fixed", 100)), 4),
  rep(rep(1:10, each=10), 8),
  rep(1:10, 80)
)

save(Simulations, file=here::here("Array_simulations.RData"))
