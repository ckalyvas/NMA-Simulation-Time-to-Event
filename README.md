# Bayesian Network Meta-Analysis Simulation Framework for Time-to-Event Outcomes

This repository contains R functions for simulating individual patient time-to-event data and fitting Bayesian network meta-analysis models under several survival-data scenarios. It supports both proportional hazards and non-proportional hazards settings, and compares one-stage individual participant data (IPD) models with two-stage aggregate-data (AggD) approaches.

The code is designed for methodological research on survival outcomes in network meta-analysis, particularly to investigate the impact of censoring mechanisms, delayed treatment effects, covariate imbalance, and fractional polynomial hazard models.

## Overview

The repository includes three main components.

`run_simulation()` implements a Bayesian network meta-analysis simulation under proportional hazards-type models using a one-stage Cox model with a piecewise baseline hazard and a two-stage Cox-based aggregate-data model.

`run_simulation_FP()` extends the framework to fractional polynomial survival models, allowing time-varying treatment effects and event-weighted estimands in non-proportional hazards settings.

`Functions_collection` contains the helper functions used to generate simulated data, transform IPD into interval-level aggregate data, define JAGS model code, and run the one-stage and two-stage analyses.

## Supported data-generating scenarios

The simulation framework currently supports four scenarios through the `model` argument.

`"uninfcens"` represents proportional hazards with non-informative censoring.

`"ECOG"` represents proportional hazards with an ECOG-like prognostic covariate imbalance across treatment arms.

`"infcens"` represents informative censoring driven by comorbidity and treatment-related censoring differences.

`"delay"` represents delayed treatment effects, where the experimental treatment has no effect before a specified time point and a treatment effect after that point.

The `model` argument can also be passed as numeric codes `1` to `4`.

## Network structure

The simulated treatment network contains three treatments: `A`, `B`, and `C`.

The pairwise evidence structure is:

- 3 studies comparing A vs B
- 2 studies comparing A vs C
- 1 study comparing C vs B

This creates a triangular treatment network with heterogeneous study-level treatment effects generated around comparison-specific mean log hazard ratios.

## Main functions

### `run_simulation()`

This function simulates study-level IPD and fits two Bayesian network meta-analysis models under a Cox proportional hazards framework.

It returns a data frame with the following outputs for each contrast:

- `Study`: contrast label
- `Model`: data-generating scenario
- `True_LHR`: true log hazard ratio used as reference
- `d_IPD`: posterior median log hazard ratio from the one-stage IPD model
- `d_IPD_lCrI`, `d_IPD_uCrI`: lower and upper 95% credible interval bounds for the one-stage model
- `d_IPD_se`: posterior standard deviation from the one-stage model
- `d_AggD`: posterior median log hazard ratio from the two-stage aggregate-data model
- `d_AggD_lCrI`, `d_AggD_uCrI`: lower and upper 95% credible interval bounds for the two-stage model
- `d_AggD_se`: posterior standard deviation from the two-stage model

### `run_simulation_FP()`

This function simulates study-level IPD and fits one-stage and two-stage Bayesian fractional polynomial network meta-analysis models.

In addition to standard treatment-effect summaries, it also computes event-weighted average treatment effects over time for each contrast. These are especially useful in delayed-effect and non-proportional hazards settings where a single constant hazard ratio is not the natural estimand.

The returned data frame includes:

- standard posterior summaries for IPD and AggD models
- event-weighted posterior summaries for IPD and AggD models
- event-weighted “true” estimands based on the data-generating process
- bias in event-weighted treatment effects
- structural step-function truth for delayed-effect scenarios
- the delayed-effect time point used in the simulation

Posterior draws of the treatment-effect parameters are also attached as attributes:

- `attr(res, "d_draws_ipd")`
- `attr(res, "d_draws_aggd")`

These can be used for custom post-processing, plotting time-varying effects, or computing alternative estimands.

## Core helper functions

### Data generation

`data_generation()` simulates study-level IPD for a two-arm trial. It handles treatment allocation, recruitment timing, administrative censoring, arm-specific opening delays, and scenario-specific event generation.

`model_specific_events()` generates event times, and, where needed, censoring times, for the specified scenario.

### One-stage and two-stage proportional hazards models

`one_stage_model()` returns the JAGS model string for the one-stage IPD Cox model with a piecewise-constant baseline hazard.

`two_stage_model()` returns the JAGS model string for the two-stage normal-normal network meta-analysis on study-specific log hazard ratios.

### Fractional polynomial models

`one_stage_FP_model()` returns the JAGS model string for the one-stage fractional polynomial survival model.

`two_stage_FP_model()` returns the JAGS model string for the two-stage aggregate-data fractional polynomial survival model.

`fp_nma()` is a wrapper that fits the two-stage fractional polynomial model and extracts posterior summaries and parameter draws.

### Aggregate-data transformation utilities

`km_and_nrisk_from_ipd()` computes Kaplan–Meier curves and numbers at risk from IPD.

`guyot_ipd_arm()` reconstructs pseudo-IPD for one arm from Kaplan–Meier and risk-table information using a Guyot-style approach.

`create_interval_counts_by_arm()` converts reconstructed arm-level pseudo-IPD into interval-level counts.

`transform_study_dynamic()` applies the reconstruction and binning pipeline study by study to create the aggregate-data input needed for the fractional polynomial NMA.

## Model assumptions and estimands

For proportional hazards settings, the main estimand is the log hazard ratio for each treatment comparison.

For delayed-effect and fractional polynomial settings, the repository supports both:

- coefficient-based treatment effects from the model parameters
- event-weighted average treatment effects over time

The event-weighted summaries are useful when treatment effects vary with time and a single hazard ratio is not constant.

In the delayed-effect scenario, the repository also computes a “least-false” Cox estimand by fitting a stratified Cox model to the pooled simulated IPD for each contrast.

## Dependencies

The code relies on standard R survival and Bayesian modelling tools. At a minimum, you will need packages such as:

- `survival`
- `R2jags`
- `coda`
- `dplyr`

Depending on how the repository is structured, you may also need a working JAGS installation available on your system.

## Example usage

### Proportional hazards simulation

```r
res_ph <- run_simulation(
  lambda      = 0.01,
  gamma       = 1.2,
  model       = "uninfcens",
  n.chains    = 2,
  n.iter      = 5000,
  n.burnin    = 1000,
  n.thin      = 1
)

print(res_ph)
```

### Delayed-effect simulation

```r
res_delay <- run_simulation(
  lambda      = 0.01,
  gamma       = 1.2,
  time0       = 180,
  model       = "delay",
  n.chains    = 2,
  n.iter      = 5000,
  n.burnin    = 1000,
  n.thin      = 1
)

print(res_delay)
```

### Fractional polynomial simulation

```r
res_fp <- run_simulation_FP(
  lambda       = 0.01,
  gamma        = 1.2,
  time0        = 180,
  comorb.prob  = 0.4,
  model        = "delay",
  P1           = 0,
  P2           = 0,
  n.chains     = 2,
  n.iter       = 5000,
  n.burnin     = 1000,
  n.thin       = 1
)

print(res_fp)
```

## Repository structure

A suggested structure for the repository is:

```text
.
├── R/
│   ├── Simulation_function.R
│   ├── Simulation_FP_function.R
│   └── Functions_collection.R
├── README.md
└── ...
```

If your current structure differs, update this section to reflect the exact file layout in the repository.

## Notes for reproducibility

The simulation code uses `set.seed(20250819)` in the main scripts. For fully reproducible results, make sure that:

- the same random seed is used
- the same JAGS version is installed
- the same R package versions are used
- the same MCMC settings are passed to the simulation functions

## Current outputs and interpretation

The contrast labels used in the output are:

- `NC`: treatment B versus A
- `OC`: treatment C versus A
- `NO`: treatment B versus C

All treatment effects are reported on the log hazard ratio scale.

For delayed-effect and fractional polynomial settings, it is important to distinguish between the coefficient-based model output and the event-weighted average effect. The latter is often the more interpretable summary when hazards are non-proportional.

## Potential extensions

This codebase can be extended in several directions, including:

- additional network geometries
- alternative censoring mechanisms
- other time-varying treatment-effect models
- alternative true estimands such as RMST-based contrasts
- simulation summaries across repeated replications, including bias, MSE, and coverage

## Contact

Chrysostomos Kalyvas, Department of Biostatistics and Medical Informatics, Medical Faculty, University of Ljubljana, Slovenia.
E-mail:  chrysostomos.kal@gmail.com
