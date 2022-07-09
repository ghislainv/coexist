seed <- 640748
perf_know <- TRUE
IV <- FALSE
source(file = here::here("_R", "Basic_parameters.R"))

load(here::here("outputs", model, "env.RData"))

for(axis in 1:n_axes){
  
  png(file=here::here("outputs", model, glue::glue("environment_{axis}.png")),
      width=fig_width, height=fig_width, units="cm", res=300)
  
  raster::plot(raster::raster(env[[axis]]),
                 col=topo.colors(255),
                 legend=FALSE,
                 axes=FALSE,
                 box=FALSE)
    
  dev.off()
}


