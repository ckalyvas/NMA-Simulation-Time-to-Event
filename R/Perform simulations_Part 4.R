#*******************************************************************************
#*
#*
#*                   Run Simulations for All Bayesian Models  
#*                         <Cox PH & Delayed>                                           
#*
#*
#*******************************************************************************


## Remove all objects from the environment ----
rm(list = ls())


## Install necessary libraries ----
list.of.packages <- c("survival", "coda", "R2jags")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages) 


## Load necessary functions ----
source("/hpc/leinehome/spinelil/Survival_NMA_simulation/Simulation_function.R")
source("/hpc/leinehome/spinelil/Survival_NMA_simulation/Functions_collection.R")


## Define scenarios ----
scenar <- "ECOG" 


## Run the simulation n.sims times and combine the results ----
# Number of simulations ...
n.sims <- 300

# ... and run the simulations.
results <- lapply(1:length(scenar), function(y) {closeAllConnections();
  
  lapply(1:n.sims, function(x) {message(paste(x, "out of", n.sims, "simulations and scenario", y, "out of", length(scenar)));
    run_simulation(lambda = 0.04, 
                   gamma = 0.5, 
                   time0 = NULL, 
                   comorb.prob = NULL,
                   model = scenar[y],    
                   n.chains = 3, 
                   n.iter = 10000, 
                   n.burnin = 1000, 
                   n.thin = 10)})})

# Combine the results into a single data frame (Each scenario is repeated n.sims times)
combined_results <- do.call(rbind, unlist(results, recursive = FALSE))

# Save results in Output as .RData
saveRDS(combined_results, file = "/hpc/leinehome/spinelil/Survival_NMA_simulation/Output_Cox/Simulation_Cox_results_Part 4.rds") 
