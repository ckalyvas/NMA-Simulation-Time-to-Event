#*******************************************************************************
#*
#*
#*                   Run Simulations for All Bayesian Models   
#*                   <Fractional polynomial - Informative>                                 
#*
#*
#*******************************************************************************


## Remove all objects from the environment ----
rm(list = ls())


## Install necessary libraries ----
list.of.packages <- c("dplyr", "survival", "coda", "R2jags")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages) 


## Load necessary functions ----
# Declare working directory
source("/hpc/leinehome/spinelil/Survival_NMA_simulation/Simulation FP_function_v2.R")
source("/hpc/leinehome/spinelil/Survival_NMA_simulation/Functions_collection.R")


## Define scenarios ----
scenar <- "infcens" 


## Run the simulation n.sims times and combine the results ----
# Number of simulations ...
n.sims <- 200

# ... and run the simulations.
results <- lapply(1:n.sims, function(x) {message(paste(x, "out of", n.sims, "simulations and scenario", scenar));
    run_simulation_FP(lambda = 0.04, 
                      gamma = 0.5, 
                      time0 = NULL, 
                      comorb.prob = 0.6,
                      model = scenar,  
                      P1 = 0,
                      P2 = 0,
                      n.chains = 3, 
                      n.iter = 10000, 
                      n.burnin = 1000, 
                      n.thin = 10)})

# Combine the data-frame 'res' into a single data frame (Each scenario is repeated n.sims times)
combined_results <- cbind(n.sim = 1:n.sims, do.call(rbind, lapply(results, function(x) x$res)))

# Isolate the 'd_draws_ipd'
d_draws_ipd_res <- lapply(results, function(x) x$d_draws_ipd)

# Isolate the 'd_draws_aggd'
d_draws_aggd_res <- lapply(results, function(x) x$d_draws_aggd)


## Save results in Output as .RData
# 'res'
saveRDS(combined_results, file = "/hpc/leinehome/spinelil/Survival_NMA_simulation/Output_FP/Simulation_FP_results_infcens.rds") 

# 'd_draws_ipd'
saveRDS(d_draws_ipd_res, file = "/hpc/leinehome/spinelil/Survival_NMA_simulation/Output_FP/Simulation_FP_draws_ipd_infcens.rds") 

# 'd_draws_aggd'
saveRDS(d_draws_aggd_res, file = "/hpc/leinehome/spinelil/Survival_NMA_simulation/Output_FP/Simulation_FP_draws_aggd_infcens.rds") 
