#*******************************************************************************
#*
#*
#*             Simulation Function for All Scenarios & Bayesian Models    
#*
#* Author: Chrysostomos Kalyvas
#* Date: October 2024
#*******************************************************************************


## Function to apply the different survival models on simulated time-to-event data ----
#' @param start.year Positive integer for the recruitment year.
#' @param stop.year Positive integer for the year where the study ends.
#' @param n.subjects Positive integer for the number of participants.
#' @param n.studies Positive integer for the number of studies.
#' @param lambda Positive scalar for the scale parameter of the Weibull distribution.
#' @param gamma Positive scalar for the shape parameter of the Weibull distribution.
#' @param time0 Positive scalar for the time of delayed treament effects.
#' @param comorb.prob Numeric in the interval from 0 to 1 for the probability of cormobidity.
#' @param model Character in {"uninfcens","ECOG","infcens","delay"} (or codes 1..4).
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#' 
#' @return A data-frame with the pairwise comparison studied (\code{Study}), the 
#' number of studies per comparison (\code{Number_studies}), the size of the studies 
#' (\code{Study_size}), the assumption around the mechanism for the event time or 
#' censoring (\code{Model}), the true hazard ratio for each pairwise comparison 
#' (\code{True_HR}), the hazard ratio estimated from the one-stage model (\code{d_IPD}),
#' the lower bound of the 95%\ credible interval of the hazard ratio based on the
#' one-stage model (\code{d_IPD_lCrI}), the upper bound of the 95%\ credible interval 
#' of the hazard ratio based on the one-stage model (\code{d_IPD_uCrI}), the hazard 
#' ratio estimated from the two-stage model (\code{d_AggD_IPD}), the lower bound of 
#' the 95%\ credible interval of the hazard ratio based on the two-stage model 
#' (\code{d_AggD_lCrI}), the upper bound of the 95%\ credible interval of the hazard 
#' ratio based on the two-stage model (\code{d_AggD_uCrI}).
#' 
#' @export
# Required packages
set.seed(20250819)

run_simulation <- function(lambda,
                           gamma,
                           time0       = NULL,
                           comorb.prob = NULL,
                           model,                  # required
                           n.chains    = n.chains,
                           n.iter      = n.iter,
                           n.burnin    = n.burnin,
                           n.thin      = n.thin) {
  # ---- normalize 'model' to a single canonical value ----
  if (missing(model) || length(model) == 0L) {
    stop("'model' is missing.", call. = FALSE)
  }
  m_raw <- model[[1L]]                      # take first element if a vector
  if (is.numeric(m_raw)) m_raw <- as.character(as.integer(m_raw))
  key <- tolower(as.character(m_raw))
  map <- c("1"="uninfcens","2"="ECOG","3"="infcens","4"="delay",
           "uninfcens"="uninfcens","ecog"="ECOG","infcens"="infcens","delay"="delay")
  if (!key %in% names(map)) {
    stop("Set 'model' to one of: 'uninfcens','ECOG','infcens','delay' (or codes 1..4).", call. = FALSE)
  }
  model <- map[[key]]                       # canonical label
  
  # ---- basic argument checks ----
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0) {
    stop("'lambda' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || gamma <= 0) {
    stop("'gamma' must be a single positive number.", call. = FALSE)
  }
  if (identical(model, "delay") && (is.null(time0) || !is.numeric(time0) || time0 <= 0)) {
    stop("For model='delay', provide a positive 'time0' (e.g., 180).", call. = FALSE)
  }
  if (identical(model, "infcens") && (is.null(comorb.prob) || comorb.prob < 0 || comorb.prob > 1)) {
    stop("For model='infcens', 'comorb.prob' must be in [0,1].", call. = FALSE)
  }
  
  # ---- pairwise comparisons and base log(HR) ----
  compar <- data.frame(
    ctrl     = c("A", "A", "C"),
    exper    = c("B", "C", "B"),
    true_LHR = c(log(0.6), log(0.8), log(0.6 / 0.8))
  )
  
  # Repeat the rows: 3 A-B, 2 A-C, 1 B-C
  compar_des0 <- compar[rep(seq_len(nrow(compar)), times = c(3, 2, 1)), ]
  
  # Set target tau² (between-study variance)
  tau2 <- 0.04
  tau  <- sqrt(tau2)
  
  # Add random noise to all study-level log HRs
  compar_des <- data.frame(
    ctrl     = compar_des0$ctrl,
    exper    = compar_des0$exper,
    true_LHR = compar_des0$true_LHR + rnorm(nrow(compar_des0), mean = 0, sd = tau)
  )
  
  # Mean "true" logHR per comparison
  true_means <- with(compar_des, tapply(true_LHR, paste0(exper, ctrl), mean))
  
  Study <- c("OC","NC","NO")
  true_map <- c(
    NC = unname(true_means["BA"]),  # B vs A
    OC = unname(true_means["CA"]),  # C vs A
    NO = unname(true_means["BC"])   # B vs C
  )
  True_LHR <- unname(true_map[Study])
  
  n_st  <- nrow(compar_des)
  sizes <- round(runif(n_st, 180, 634), 0)
  
  # ---- study timelines ----
  starts        <- sample(2012:2016, n_st, replace = TRUE)
  dur_years     <- sample(c(1, 1.5, 2, 3, 4), n_st, replace = TRUE, prob = c(.25, .25, .25, .20, .05))
  stops         <- floor(starts + dur_years)
  accrual_years <- sample(c(0.5, 1, 2, 3), n_st, replace = TRUE, prob = c(.3, .4, .2, .1))
  arm_lag_days  <- sample(c(0, 60, 120, 180), n_st, replace = TRUE, prob = c(.5, .2, .2, .1))
  
  # ---- per-study IPD ----
  create_study <- if (model == "uninfcens") {
    lapply(seq_len(nrow(compar_des)), function(x)
      data_generation(start.year   = starts[x],
                      stop.year    = stops[x],
                      n.subjects   = sizes[x],
                      trt.names    = as.character(compar_des[x, 1:2]),
                      deltaB       = compar_des[x, 3],
                      lambda       = lambda,
                      gamma        = gamma,
                      time0        = NULL,
                      comorb.prob  = NULL,
                      model        = model,
                      accrual.years = accrual_years[x],
                      arm.lag.days  = arm_lag_days[x]))
  } else if (model == "ECOG") {
    lapply(seq_len(nrow(compar_des)), function(x)
      data_generation(start.year   = starts[x],
                      stop.year    = stops[x],
                      n.subjects   = sizes[x],
                      trt.names    = as.character(compar_des[x, 1:2]),
                      deltaB       = compar_des[x, 3],
                      lambda       = lambda,
                      gamma        = gamma,
                      time0        = NULL,
                      comorb.prob  = NULL,
                      model        = model,
                      accrual.years = accrual_years[x],
                      arm.lag.days  = arm_lag_days[x]))
  } else if (model == "infcens") {
    lapply(seq_len(nrow(compar_des)), function(x)
      data_generation(start.year   = starts[x],
                      stop.year    = stops[x],
                      n.subjects   = sizes[x],
                      trt.names    = as.character(compar_des[x, 1:2]),
                      deltaB       = compar_des[x, 3],
                      lambda       = lambda,
                      gamma        = gamma,
                      time0        = NULL,
                      comorb.prob  = comorb.prob,
                      delta.comorb = log(3),
                      comorb.lambda = 0.06,
                      comorb.gamma  = gamma,
                      comorb.beta   = log(3),
                      model        = model,
                      accrual.years = accrual_years[x],
                      arm.lag.days  = arm_lag_days[x]))
  } else { # model == "delay"
    lapply(seq_len(nrow(compar_des)), function(x)
      data_generation(start.year   = starts[x],
                      stop.year    = stops[x],
                      n.subjects   = sizes[x],
                      trt.names    = as.character(compar_des[x, 1:2]),
                      deltaB       = compar_des[x, 3],
                      lambda       = lambda,
                      gamma        = gamma,
                      time0        = time0,           # e.g., 180
                      comorb.prob  = NULL,
                      model        = model,
                      accrual.years = accrual_years[x],
                      arm.lag.days  = arm_lag_days[x]))
  }
  
  # ---- least-false Cox "true" logHR map (delay scenario only) ----
  if (model == "delay") {
    bind_comp <- function(lbl) {
      map <- list(
        BA = c(ctrl = "A", exper = "B"),
        CA = c(ctrl = "A", exper = "C"),
        BC = c(ctrl = "C", exper = "B")
      )
      pair <- map[[lbl]]
      idx  <- which(compar_des$ctrl == pair["ctrl"] & compar_des$exper == pair["exper"])
      if (length(idx) == 0) return(NULL)
      do.call(rbind, lapply(idx, function(i) {
        dd <- create_study[[i]]
        dd$study_id <- i
        dd
      }))
    }
    get_lhr <- function(lbl) {
      d <- bind_comp(lbl)
      if (is.null(d) || sum(d$censor == 1) == 0) return(NA_real_)
      cf <- coef(summary(coxph(Surv(event.time, censor) ~ treat + ecog + strata(study_id), data = d)))
      unname(cf["treat", "coef"])
    }
    true_lhr_map <- c(NC = get_lhr("BA"), OC = get_lhr("CA"), NO = get_lhr("BC"))
    True_LHR <- unname(true_lhr_map[Study])
  }
  
  # ---------------- ONE-STAGE (JAGS Cox w/ piecewise baseline) ----------------
  time0_m <- if (!is.null(time0)) time0/30.4375 else NA_real_
  
  # tmax from the created studies (proper extraction)
  tmax <- max(unlist(lapply(create_study, function(d) d[["event.time"]])))
  a_vec <- sort(unique(c(seq(0, ceiling(tmax) + 1e-8, by = 2), time0_m)))  # finer and aligned
  m_val <- length(a_vec) - 1L
  
  data_comp0 <- lapply(seq_along(create_study), function(x) {
    list(
      n     = nrow(create_study[[x]]),
      time  = create_study[[x]]$event.time,
      treat = create_study[[x]]$treat,        # 0/1
      delta = create_study[[x]]$censor,       # 1=event, 0=censored
      ecog  = create_study[[x]]$ecog,
      index = rep(x, nrow(create_study[[x]]))
    )
  })
  
  # Properly extract columns from lists
  time_all  <- unlist(lapply(data_comp0, function(z) z[["time"]]))
  treat_all <- unlist(lapply(data_comp0, function(z) z[["treat"]]))
  delta_all <- unlist(lapply(data_comp0, function(z) z[["delta"]]))
  ecog_all  <- unlist(lapply(data_comp0, function(z) z[["ecog"]]))
  index_all <- unlist(lapply(data_comp0, function(z) z[["index"]]))
  zeros_all <- rep(0, length(time_all))
  
  # helper: map a time to its left interval index within a partition 'alpha'
  get_int_obs <- function(npts, alpha, time) {
    idx <- findInterval(time, vec = alpha, rightmost.closed = TRUE, all.inside = TRUE)
    pmax(1, pmin(idx, length(alpha) - 1))
  }
  intobs_all <- get_int_obs(npts = length(time_all), alpha = a_vec, time = time_all)
  
  combined_data <- list(
    n            = length(time_all),
    m            = m_val,
    time         = time_all,
    a            = a_vec,
    treat        = treat_all,
    delta        = delta_all,
    int.obs      = intobs_all,
    zeros        = zeros_all,
    index        = index_all,
    ecog         = ecog_all,
    treat.index  = match(compar_des$exper, LETTERS),
    compar.index = match(compar_des$ctrl, LETTERS),
    num_studies  = nrow(compar_des)
  )
  
  run_model <- function(data) {
    jm <- jags.model(file = textConnection(one_stage_model(n.studies = nrow(compar_des))),
                     data = data, n.chains = n.chains, n.adapt = min(1000, n.iter))
    update(jm, n.burnin)
    coda.samples(model = jm,
                 variable.names = c("beta", "lambda", "d", "logHR"),
                 n.iter = n.iter)
  }
  
  samples <- run_model(combined_data)
  samples_matrix <- as.matrix(samples)
  
  LHR_OC_IPD  <- median(samples_matrix[, "d[3]"])      # C vs A
  LHR_NC_IPD  <- median(samples_matrix[, "d[2]"])      # B vs A
  LHR_NO_IPD  <- median(samples_matrix[, "logHR"])     # B vs C
  
  se_OC_IPD   <- sd(samples_matrix[, "d[3]"])
  se_NC_IPD   <- sd(samples_matrix[, "d[2]"])
  se_NO_IPD   <- sd(samples_matrix[, "logHR"])
  
  CRI_OC_IPD  <- quantile(samples_matrix[, "d[3]"],  c(.025, .975))
  CRI_NC_IPD  <- quantile(samples_matrix[, "d[2]"],  c(.025, .975))
  CRI_NO_IPD  <- quantile(samples_matrix[, "logHR"], c(.025, .975))
  
  # ---------------- TWO-STAGE (studywise Cox, then RE NMA in JAGS) ------------
  cox_model <- lapply(create_study, function(x)
    coxph(Surv(event.time, censor) ~ treat, data = x))
  
  l_LHR    <- sapply(cox_model, function(fit) coef(summary(fit))["treat", "coef"])
  l_LHR_se <- sapply(cox_model, function(fit) coef(summary(fit))["treat", "se(coef)"])
  
  data_list <- list(
    ns            = nrow(compar_des),
    treat.index   = match(compar_des$exper, LETTERS),
    compar.index  = match(compar_des$ctrl, LETTERS),
    y             = as.numeric(l_LHR),
    se            = as.numeric(l_LHR_se)
  )
  
  inits <- function() list(d = c(NA, rnorm(2)))
  
  jags_model <-
    jags(data = data_list, inits = inits,
         parameters.to.save = c("d", "diff.AggD"),
         model.file = textConnection(two_stage_model(n.studies = nrow(compar_des))),
         n.chains = n.chains, n.iter = n.iter,
         n.burnin = n.burnin, n.thin = n.thin, DIC = TRUE)
  
  aggd_samples <- as.matrix(as.mcmc(jags_model))
  
  LHR_OC_AGGD <- median(aggd_samples[, "d[3]"])
  LHR_NC_AGGD <- median(aggd_samples[, "d[2]"])
  LHR_NO_AGGD <- median(aggd_samples[, "diff.AggD"])
  
  se_OC_AGGD  <- sd(aggd_samples[, "d[3]"])
  se_NC_AGGD  <- sd(aggd_samples[, "d[2]"])
  se_NO_AGGD  <- sd(aggd_samples[, "diff.AggD"])
  
  CRI_OC_AGGD <- quantile(aggd_samples[, "d[3]"],      c(.025, .975))
  CRI_NC_AGGD <- quantile(aggd_samples[, "d[2]"],      c(.025, .975))
  CRI_NO_AGGD <- quantile(aggd_samples[, "diff.AggD"], c(.025, .975))
  
  data.frame(
    Study       = c("OC","NC","NO"),
    Model       = rep(model, 3),
    True_LHR    = True_LHR,
    d_IPD       = c(LHR_OC_IPD,  LHR_NC_IPD,  LHR_NO_IPD),
    d_IPD_lCrI  = c(CRI_OC_IPD[1],CRI_NC_IPD[1],CRI_NO_IPD[1]),
    d_IPD_uCrI  = c(CRI_OC_IPD[2],CRI_NC_IPD[2],CRI_NO_IPD[2]),
    d_IPD_se    = c(se_OC_IPD,    se_NC_IPD,    se_NO_IPD),
    d_AggD      = c(LHR_OC_AGGD,  LHR_NC_AGGD,  LHR_NO_AGGD),
    d_AggD_lCrI = c(CRI_OC_AGGD[1],CRI_NC_AGGD[1],CRI_NO_AGGD[1]),
    d_AggD_uCrI = c(CRI_OC_AGGD[2],CRI_NC_AGGD[2],CRI_NO_AGGD[2]),
    d_AggD_se   = c(se_OC_AGGD,    se_NC_AGGD,    se_NO_AGGD)
  )
}