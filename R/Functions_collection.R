#*******************************************************************************
#*
#*
#*                      Compiling all necessary functions                                                                               
#*
#* Author: Chrysostomos Kalyvas
#* Date: October 2024
#*******************************************************************************



## A: Function to simulate survival data with informative censoring ----
#' @param start.year Positive integer for the recruitment year.
#' @param stop.year Positive integer for the year where the study ends.
#' @param n.subjects Positive integer for the number of participants in the study.
#' @param trt.names A vector of two alphabet letters referring to the compared 
#'   treatments. The first element is the control arm.
#' @param deltaB Numeric for the hazard ratio in the logarithmic scale between 
#'   the experimental and control arms.
#' @param lambda Positive scalar for the scale parameter of the Weibull distribution.
#' @param gamma Positive scalar for the shape parameter of the Weibull distribution.
#' @param comorb.gamma Positive scalar for the shape parameter of the Weibull distribution 
#'   for the comorbidity Cox model.
#' @param comorb.lambda Positive scalar for the scale parameter of the Weibull distribution 
#'   for the comorbidity Cox model.
#' @param time0 Positive scalar for the time of delayed treatment effects.
#' @param comorb.prob Numeric in the interval from 0 to 1 for the probability of cormobidity.
#' @param comorb.beta Effect of comorbidity in the model.
#' @param delta.comorb Numeric for the hazard ratio in the logarithmic scale between 
#'   participants with comorbidity and participants without comorbidity.
#' @param model Character with possible values \code{"uninfcens"}, \code{"infcens"}, and
#'   \code{"delay"} to indicate proportional hazards with non-informative censoring,
#'   informative censoring, and delayed treatment effect with non-informative censoring, respectively.
#' 
#' @return A data-frame with the treatment names as alphabet letter (\code{treat.name}) and 
#' corresponding number (\code{treat}), the recruitment date (\code{start.date}), the observed
#' time to event (\code{event.time}), and censoring status (\code{censor}).
#' 
#' @export
data_generation <- function(start.year, 
                            stop.year, 
                            n.subjects,
                            trt.names, 
                            deltaB, 
                            lambda, 
                            gamma, 
                            time0, 
                            comorb.prob,
                            comorb.beta    = log(3),
                            comorb.gamma   = 0.5,
                            comorb.lambda  = 0.06,
                            delta.comorb   = log(3),
                            model,
                            accrual.years  = 1,
                            arm.lag.days   = 0) { 
  
  trt.names <- if (length(trt.names) != 2L) {
    stop("'trt.names' must have exactly two elements.", call. = FALSE)
  } else trt.names
  
  if (!is.numeric(lambda) || lambda <= 0) stop("'lambda' must be > 0.", call. = FALSE)
  if (!is.numeric(gamma)  || gamma  <= 0) stop("'gamma' must be > 0.",  call. = FALSE)
  
  if (missing(model) || length(model) == 0L) stop("'model' is missing.", call. = FALSE)
  m_raw <- model[[1L]]
  if (is.numeric(m_raw)) m_raw <- as.character(as.integer(m_raw))
  key <- tolower(as.character(m_raw))
  map <- c("1"="uninfcens","2"="ECOG","3"="infcens","4"="delay",
           "uninfcens"="uninfcens","ecog"="ECOG","infcens"="infcens","delay"="delay")
  if (!key %in% names(map)) stop("model ??? {'uninfcens','ECOG','infcens','delay'} (or codes 1..4).", call. = FALSE)
  model <- map[[key]]
  
  if (identical(model, "delay")) {
    if (is.null(time0) || !is.numeric(time0) || time0 <= 0) stop("For model='delay', provide positive 'time0' (e.g., 180).", call. = FALSE)
  } else {
    time0 <- NULL
  }
  if (identical(model, "infcens")) {
    if (is.null(comorb.prob) || comorb.prob < 0 || comorb.prob > 1) {
      stop("For model='infcens', 'comorb.prob' must be in [0,1].", call. = FALSE)
    }
  } else {
    comorb.prob <- NULL
  }
  
  assign     <- runif(n.subjects)
  treat.name <- factor(ifelse(assign < 0.5, trt.names[1], trt.names[2]), levels = trt.names)
  treat      <- as.integer(treat.name == trt.names[2])
  ecog       <- rep(0L, n.subjects)
  
  start0        <- as.Date(sprintf("%d-01-01", start.year))
  entry_offset  <- round(runif(n.subjects, 0, 365 * accrual.years))
  lag_vec       <- ifelse(treat == 1, arm.lag.days, 0)
  start.date    <- start0 + entry_offset + lag_vec
  
  CUTOFFDT <- as.Date(sprintf("%d-01-01", stop.year)) - 1
  T.ind    <- pmax(0, as.numeric(CUTOFFDT - start.date, units = "days"))
  
  if (model == "uninfcens") {
    gen.evnt.time <- model_specific_events(
      treat = treat, deltaB = deltaB, lambda = lambda, gamma = gamma,
      time0 = NULL, comorb.prob = NULL, comorb.beta = NULL, comorb.gamma = NULL,
      comorb.lambda = NULL, ecog = NULL, beta.ecog = log(2.5), delta.comorb = NULL,
      model = model, n.subjects = n.subjects, treat.name = treat.name, trt.names = trt.names
    )
    ev_days <- gen.evnt.time$uncensored.events
    C       <- runif(n.subjects) * T.ind
    obs_days_raw <- pmin(ev_days, C, T.ind)
    status <- as.integer(obs_days_raw == ev_days)
    
  } else if (model == "ECOG") {
    ecog <- rbinom(n.subjects, 1, ifelse(treat == 1, 0.7, 0.3))
    gen.evnt.time <- model_specific_events(
      treat = treat, deltaB = deltaB, lambda = lambda, gamma = gamma,
      time0 = NULL, comorb.prob = NULL, comorb.beta = NULL, comorb.gamma = NULL,
      comorb.lambda = NULL, ecog = ecog, beta.ecog = log(2.5), delta.comorb = NULL,
      model = model, n.subjects = n.subjects, treat.name = treat.name, trt.names = trt.names
    )
    ev_days <- gen.evnt.time$uncensored.events
    C       <- runif(n.subjects) * T.ind
    obs_days_raw <- pmin(ev_days, C, T.ind)
    status <- as.integer(obs_days_raw == ev_days)
    
  } else if (model == "infcens") {
    gen.time <- model_specific_events(
      treat = treat, deltaB = deltaB, lambda = lambda, gamma = gamma,
      time0 = NULL,
      comorb.prob  = comorb.prob,
      comorb.beta  = comorb.beta,
      comorb.gamma = comorb.gamma,
      comorb.lambda = comorb.lambda,
      delta.comorb = delta.comorb,
      model = model,
      n.subjects = n.subjects,
      treat.name = treat.name,
      trt.names  = trt.names,
      comorb.p.ctrl    = 0.40,
      comorb.p.exper   = 0.65,
      alpha.cens.treat = log(1.5),
      delta.int        = log(1.3)
    )
    ev_days <- gen.time$uncensored.events
    ce_days <- gen.time$censored.events
    obs_days_raw <- pmin(ev_days, ce_days, T.ind)
    status <- as.integer(ev_days <= pmin(ce_days, T.ind))
    
  } else { # delay
    gen.time <- model_specific_events(
      treat = treat, deltaB = deltaB, lambda = lambda, gamma = gamma,
      time0 = time0, comorb.prob = NULL, comorb.beta = NULL, comorb.gamma = NULL,
      comorb.lambda = NULL, ecog = NULL, beta.ecog = log(2.5), delta.comorb = NULL,
      model = model, n.subjects = n.subjects, treat.name = treat.name, trt.names = trt.names
    )
    ev_days  <- gen.time$uncensored.events
    C_admin  <- runif(n.subjects) * T.ind
    obs_days_raw <- pmin(ev_days, C_admin, T.ind)
    status <- as.integer(obs_days_raw == ev_days)
  }
  
  # --- Apply the "min 1 day" guard safely ---
  obs_days <- obs_days_raw
  
  # If truly no administrative follow-up, keep 1 day for numerics but force censoring
  zero_fu <- (T.ind <= 0)
  if (any(zero_fu)) {
    obs_days[zero_fu] <- 1
    status[zero_fu]   <- 0
  }
  
  # For those with positive follow-up, only bump tiny times
  tiny <- (T.ind > 0) & (obs_days < 1)
  if (any(tiny)) obs_days[tiny] <- 1
  
  event.time <- obs_days / 30.4375  # months
  
  data.frame(treat.name, treat, start.date, event.time, censor = status, ecog)
}



## B: Generate treatment-specific event (and censoring) times ----
# Purpose: produce arm-specific vector(s) of event times; optionally censoring
# times when model == "infcens" (informative censoring).
# Notes:
# - Uses inverse-CDF for Weibull with linear predictor in hazard scale.
# - For model == "delay": experimental arm has HR=1 before time0, HR=exp(deltaB)
# after time0 using conditional inversion (no memorylessness).
#' @param treat Binary variable for the assigned arm of each participant.
#' @param deltaB Numeric for the hazard ration in the logarithmic scale between the experimental and control arms.
#' @param lambda Positive scalar for the scale parameter of the Weibull distribution.
#' @param gamma Positive scalar for the shape parameter of the Weibull distribution.
#' @param time0 Positive scalar for the time of delayed treatment effects.
#' @param comorb.prob Numeric in the interval from 0 to 1 for the probability of cormobidity.
#' @param comorb.beta Coefficient of comorbidity in the Cox model
#' @param comorb.gamma Positive scalar for the shape parameter of the Weibull distribution in the comorbidity survival function.
#' @param comorb.lambda Positive scalar for the scale parameter of the Weibull distribution in the comorbidity survival function.
#' @param delta.comorb Effect of comorbidity in comorbidity (Yes) vs comorbidity (No).
#' @param model Character with possible values \code{"uninfcens"}, \code{"infcens"}, and
#'   \code{"delay"} to indicate proportional hazards with non-informative censoring,
#'   informative censoring, and delayed treatment effect with non-informative censoring, respectively.
#'
#' @return list with elements:
#' uncensored.events: vector of event times for all subjects (in DAYS)
#' censored.events: vector of censor times (only for model == "infcens")
#' 
#' @export
model_specific_events <- function(treat,
                                  deltaB,
                                  lambda,
                                  gamma,
                                  time0            = NULL,
                                  comorb.prob      = NULL,
                                  comorb.beta      = NULL,
                                  comorb.gamma     = NULL,
                                  comorb.lambda    = NULL,
                                  ecog             = NULL,
                                  beta.ecog        = log(2.5),
                                  delta.comorb     = NULL,
                                  model,
                                  n.subjects,
                                  treat.name,
                                  trt.names,
                                  # NEW knobs (defaults keep old behavior)
                                  comorb.p.ctrl    = NULL,
                                  comorb.p.exper   = NULL,
                                  alpha.cens.treat = 0,          # log-HR on censoring for treat arm
                                  delta.int        = 0) {        # trt×comorb interaction on event
  
  # --- Validate basics ---------------------------------------------------
  if (length(unique(treat)) > 2) stop("'treat' must be binary.", call. = FALSE)
  if (!is.numeric(lambda) || lambda <= 0) stop("'lambda' must be > 0.", call. = FALSE)
  if (!is.numeric(gamma)  || gamma  <= 0) stop("'gamma' must be > 0.",  call. = FALSE)
  
  # Canonicalize 'model'
  if (missing(model) || length(model) == 0L) stop("'model' is missing.", call. = FALSE)
  m_raw <- model[[1L]]; if (is.numeric(m_raw)) m_raw <- as.character(as.integer(m_raw))
  key <- tolower(as.character(m_raw))
  map <- c("1"="uninfcens","2"="ECOG","3"="infcens","4"="delay",
           "uninfcens"="uninfcens","ecog"="ECOG","infcens"="infcens","delay"="delay")
  if (!key %in% names(map)) stop("model ??? {'uninfcens','ECOG','infcens','delay'} (or 1..4).", call. = FALSE)
  model <- map[[key]]
  
  if (identical(model, "delay")) {
    if (is.null(time0) || time0 <= 0) stop("'time0' must be positive for model='delay'.", call. = FALSE)
  }
  
  # --- Draw base uniforms ------------------------------------------------
  u       <- runif(n.subjects)
  u.ctrl  <- u[treat.name == trt.names[1]]
  u.exper <- u[treat.name == trt.names[2]]
  AVAL    <- rep(NA_real_, n.subjects)
  CNSR    <- rep(NA_real_, n.subjects)
  
  if (model == "uninfcens") {
    # PH; no covariates in generator
    AVAL[treat.name == trt.names[1]] <- (-log(u.ctrl)  / lambda)^(1/gamma)
    AVAL[treat.name == trt.names[2]] <- (-log(u.exper) / (lambda * exp(deltaB)))^(1/gamma)
    return(list(uncensored.events = AVAL, censored.events = NULL))
    
  } else if (model == "ECOG") {
    if (is.null(ecog)) stop("'ecog' must be provided for model='ECOG'.", call. = FALSE)
    AVAL[treat.name == trt.names[1]] <- (-log(u.ctrl)  / (lambda * exp(beta.ecog * ecog[treat.name == trt.names[1]])))^(1/gamma)
    AVAL[treat.name == trt.names[2]] <- (-log(u.exper) / (lambda * exp(deltaB + beta.ecog * ecog[treat.name == trt.names[2]])))^(1/gamma)
    return(list(uncensored.events = AVAL, censored.events = NULL))
    
  } else if (model == "infcens") {
    # Preconditions
    if (is.null(comorb.beta) || is.null(comorb.gamma) || is.null(comorb.lambda) || is.null(delta.comorb)) {
      stop("Provide comorb.beta, comorb.gamma, comorb.lambda, delta.comorb for 'infcens'.", call. = FALSE)
    }
    
    # Allow arm-imbalanced comorbidity prevalence; fall back to single comorb.prob if provided
    if (!is.null(comorb.p.ctrl) && !is.null(comorb.p.exper)) {
      p_vec <- ifelse(treat == 0, comorb.p.ctrl, comorb.p.exper)
    } else if (!is.null(comorb.prob)) {
      p_vec <- rep(comorb.prob, n.subjects)
    } else {
      stop("Specify comorb.prob or (comorb.p.ctrl & comorb.p.exper) for 'infcens'.", call. = FALSE)
    }
    comorb <- rbinom(n.subjects, 1, p_vec)
    
    # EVENT hazard: baseline * exp( deltaB*treat + delta.comorb*comorb + delta.int*(treat*comorb) )
    lin_ctrl  <- delta.comorb * comorb[treat == 0]
    lin_exper <- deltaB + delta.comorb * comorb[treat == 1] + delta.int * comorb[treat == 1]
    hz_ctrl   <- lambda * exp(lin_ctrl)
    hz_exper  <- lambda * exp(lin_exper)
    
    u <- runif(n.subjects)
    u.c <- u[treat == 0]; u.e <- u[treat == 1]
    AVAL[treat == 0] <- (-log(u.c) / hz_ctrl)^(1 / gamma)
    AVAL[treat == 1] <- (-log(u.e) / hz_exper)^(1 / gamma)
    
    # CENSORING hazard: comorbidity + (optional) treatment effect on censoring
    # baseline_cens * exp( comorb.beta*comorb + alpha.cens.treat*treat )
    u.com  <- runif(n.subjects)
    hz_cen <- comorb.lambda * exp(comorb.beta * comorb + alpha.cens.treat * treat)
    CNSR   <- (-log(u.com) / hz_cen)^(1 / comorb.gamma)
    
    list(uncensored.events = AVAL, censored.events = CNSR)
    
  } else { # delay
    # Control: PH (no trt effect)
    AVAL[treat.name == trt.names[1]] <- ((-log(u.ctrl)) / lambda)^(1/gamma)
    
    # Experimental: HR=1 up to time0; HR=exp(deltaB) after time0
    idx  <- which(treat.name == trt.names[2])
    S_t0 <- exp(-lambda * (time0^gamma))                 # survival up to time0 with HR=1
    fail_before <- u[idx] > S_t0                         # T < time0
    # before time0
    if (any(fail_before)) {
      AVAL[idx[fail_before]] <- ((-log(u[idx[fail_before]])) / lambda)^(1/gamma)
    }
    # after time0: conditional inversion beyond t0
    if (any(!fail_before)) {
      theta   <- exp(deltaB)
      u_resid <- u[idx[!fail_before]] / S_t0             # Uniform(0,1]
      t_pow   <- (time0^gamma) + (-log(u_resid)) / (theta * lambda)
      AVAL[idx[!fail_before]] <- t_pow^(1/gamma)
    }
    return(list(uncensored.events = AVAL, censored.events = NULL))
  }
}



## C: nt.obs function (Alvares, page 7) ----
# Helper: map a time to its left interval index within a partition 'alpha'.
get_int_obs <- function(npts, alpha, time) {
  idx <- findInterval(time, vec = alpha, rightmost.closed = TRUE, all.inside = TRUE)
  # findInterval returns 0 for times < alpha[1]; shift to 1..(length(alpha)-1)
  pmax(1, pmin(idx, length(alpha) - 1))
}


## D: One-stage (fixed-effect or random-effects) NMA with PH Cox regression ----
# JAGS model string using zeros trick to encode the individual log-likelihood for
# a piecewise-constant baseline hazard Cox model. Between-study random effects on
# treatment effects (beta[l]) with reference-coded d's.
one_stage_model <- function(n.studies) {
  
  stringcode <- "
  model {

    # --- Likelihood (zeros trick with clean -logLik) ---
    for (i in 1:n) {
      # Build baseline cumulative hazard H0 up to time[i] via piecewise-constant lambda
      for (k in 1:int.obs[i]) {
        cond[i,k] <- step(time[i] - a[k + 1])
        # length of the k-th piece that contributes
        seglen[i,k] <- cond[i,k] * (a[k + 1] - a[k]) + (1 - cond[i,k]) * (time[i] - a[k])
        # contribution of piece k using the baseline hazard for this study index[i]
        HH[i,k] <- seglen[i,k] * lambda[k, index[i]]
      }
      H[i] <- sum(HH[i, 1:int.obs[i]])   # baseline cumulative hazard up to time[i]

      # Linear predictor
      linpred[i] <- beta[index[i]] * treat[i] + beta_ecog * ecog[i]

      # Negative log-likelihood contribution for subject i:
      #  nll_i = H0_i * exp(eta_i) - delta_i * ( log(lambda_interval,study) + eta_i )
      loglam[i] <- log(lambda[int.obs[i], index[i]])
      nll[i] <- H[i] * exp(linpred[i]) - delta[i] * (loglam[i] + linpred[i])

      # zeros trick
      phi[i] <- nll[i] + 1.0E-6  
      zeros[i] ~ dpois(phi[i])
    }

    # --- Priors ---
    # Random-effects on treatment effects across studies
    tau ~ dunif(0, 5)
    prec <- pow(tau, -2)
    beta_ecog ~ dnorm(0, 0.001)

    for (l in 1:num_studies) {
      beta[l] ~ dnorm(md[l], prec)
      md[l] <- d[treat.index[l]] - d[compar.index[l]]
    }

    # Piecewise baseline hazards per study
    for (l in 1:num_studies) {
      for (k in 1:m) {
        lambda[k, l] ~ dgamma(0.01, 0.01)
      }
    }

    # Treatment effects (reference-coded)
    d[1] <- 0
    for (k in 2:3) {
      d[k] ~ dnorm(0, 1.0E-6)
    }
    logHR <- d[2] - d[3] #in our simulation we define NO as B vs C

    # Constant for zeros trick (data node so it can be passed, default 0)
    #cst <- 1000
  }
  "
  
}
  
## E: Two-stage (fixed-effect or random-effects) NMA of log hazard ratios ----
# Standard normal-normal meta-analysis on study-specific logHR with RE variance tau.
two_stage_model <- function(n.studies) {
  
  stringcode <- "model{ 
                   for (i in 1:ns) { 
                     variance[i] <- pow(se[i], 2) 
                     prec[i] <- 1 / variance[i] 
                     y[i] ~ dnorm(delta[i], prec[i])\n"
  
  stringcode <- paste(stringcode, "delta[i] ~ dnorm(md[i], tau) 
                                   md[i] <- d[treat.index[i]] - d[compar.index[i]]}\n")

  stringcode <- paste(stringcode, "d[1] <- 0 # treatment effect is zero for reference treatment
                                   for (k in 2:3) { 
                                     d[k] ~ dnorm(0, .0001) 
                                   }
                                   diff.AggD <- d[2] - d[3]\n") #in our simulation we define NO as B vs C
  
  stringcode <- paste(stringcode, "sd ~ dunif(0, 5) 
                                   tau <- pow(sd, -2)}\n")

}

## F: Create time intervals per arm (aggregate fractional polynomial) ----
# Purpose: convert (pseudo-)IPD into interval-wise counts by arm: events r, censored m,
# at-risk M at the start of interval. Used as input for FP AggD models.
# --- SAFE binning of pseudo-IPD into intervals (expects YEARS) ----------
create_interval_counts_by_arm <- function(data, interval_length) {
  # data must have: event.time (YEARS), censor (1=event, 0=censor), treat.name (factor)
  data$event.time <- as.numeric(data$event.time)
  
  tmax  <- max(data$event.time, na.rm = TRUE)
  K     <- max(1L, ceiling(tmax / interval_length))               # at least one interval
  edges <- seq(0, K * interval_length, by = interval_length)      # length K+1
  
  out <- vector("list", length=0)
  arms <- levels(factor(data$treat.name))
  
  for (arm in arms) {
    arm_data <- data[data$treat.name == arm, , drop = FALSE]
    
    # clamp times to stay inside last interval for cut()
    x  <- pmin(arm_data$event.time, tail(edges,1) - .Machine$double.eps)
    kk <- as.integer(cut(x, breaks = edges, right = FALSE, include.lowest = TRUE))
    
    # events/censors per bin
    r_by_k <- tapply(arm_data$censor == 1, kk, sum)
    m_by_k <- tapply(arm_data$censor == 0, kk, sum)
    r_vec  <- integer(K); m_vec <- integer(K)
    if (length(r_by_k)) r_vec[as.integer(names(r_by_k))] <- as.integer(r_by_k)
    if (length(m_by_k)) m_vec[as.integer(names(m_by_k))] <- as.integer(m_by_k)
    
    # risk set at left edges: subjects with time > left_edge
    left_edges <- edges[-length(edges)]
    M_vec <- vapply(left_edges, function(L) sum(arm_data$event.time > L), integer(1))
    
    out[[length(out)+1L]] <- data.frame(
      start_time = left_edges,
      end_time   = left_edges + interval_length,
      m = m_vec, r = r_vec, M = M_vec,
      ARM = arm, stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}


# KM + numbers-at-risk per arm from IPD
# Produces stepwise KM (t,S) and numbers at risk at landmarks every 'nrisk_every'.
km_and_nrisk_from_ipd <- function(ipd_arm, nrisk_every) {
  fit <- survival::survfit(survival::Surv(event.time, censor) ~ 1, data = ipd_arm)
  
  sm  <- summary(fit)
  steps <- data.frame(t = c(0, sm$time), S = c(1, sm$surv))
  steps <- steps[order(steps$t), ]
  steps$S <- cummin(steps$S)
  
  lm_times <- seq(0, max(steps$t), by = nrisk_every)
  sm_lm <- summary(fit, times = lm_times, extend = TRUE)
  nrisk <- data.frame(t = sm_lm$time, n_risk = sm_lm$n.risk)
  
  list(steps = steps, nrisk = nrisk)
}

# Guyot et al. pseudo-IPD reconstruction for ONE arm
# Given KM steps and risk table at landmarks, allocate events/censors per interval,
# then expand into pseudo-records (time, status).
guyot_ipd_arm <- function(steps, nrisk, arm_label = "A") {
  # steps: data.frame t,S with t0=0, S0=1; nrisk: data.frame t, n_risk incl. t=0
  stopifnot(steps$t[1] == 0, abs(steps$S[1] - 1) < 1e-9)
  
  t <- steps$t
  S <- steps$S
  
  # If there are no events, summary(survfit) may give only t = 0
  if (length(t) < 2) {
    # fabricate a tiny second time so downstream code has one interval
    t <- c(0, 1e-8)
    S <- c(1, S[1])  # no drop
  }
  
  H  <- -log(pmax(S, 1e-12))
  dH <- diff(H)
  dt <- diff(t)
  
  # Build landmark breaks and ensure there are at least two
  lm_raw <- sort(unique(c(0, nrisk$t, tail(t, 1))))
  if (length(lm_raw) < 2 || !is.finite(tail(lm_raw, 1))) {
    lm <- c(0, max(t, na.rm = TRUE) + 1e-8)
  } else {
    lm <- lm_raw
  }
  
  # Cut the left endpoints of each KM step into blocks
  t_left <- t[-length(t)]  # interior left edges
  if (length(t_left) == 0L) {
    # should not happen after the guard, but be safe
    t_left <- 0
  }
  block_id <- cut(t_left,
                  breaks = lm, include.lowest = TRUE, right = FALSE, labels = FALSE)
  
  n_int <- length(dH)
  M <- r <- m <- integer(n_int)
  
  for (b in sort(unique(block_id))) {
    idx <- which(block_id == b)
    t0  <- min(t[idx])
    
    # start risk from nearest landmark (at t0)
    N0 <- nrisk$n_risk[match(t0, nrisk$t)]
    if (is.na(N0)) {
      prev <- max(lm[lm <= t0])
      N0   <- nrisk$n_risk[match(prev, nrisk$t)]
    }
    if (is.na(N0)) N0 <- 0L  # last-ditch guard
    
    # first pass: events from hazard increment, no censors yet
    n_i <- N0
    for (j in seq_along(idx)) {
      i <- idx[j]
      M[i] <- n_i
      di   <- round(n_i * (1 - exp(-max(0, dH[i]))))
      di   <- max(0L, min(di, n_i))
      r[i] <- di
      m[i] <- 0L
      n_i  <- n_i - di
    }
    
    # match next landmark exactly by allocating censoring if possible
    t1 <- lm[match(t0, lm) + 1]
    if (!is.na(t1) && t1 %in% nrisk$t) {
      N1   <- nrisk$n_risk[match(t1, nrisk$t)]
      need <- max(0L, n_i - N1)
      if (need > 0L && length(idx) > 0L) {
        w    <- dt[idx]
        w    <- if (sum(w) > 0) w / sum(w) else rep(1 / length(idx), length(idx))
        cand <- floor(w * need)
        rem  <- need - sum(cand)
        
        if (rem > 0) cand[order(-w)][seq_len(rem)] <- cand[order(-w)][seq_len(rem)] + 1L
        
        leftover <- 0L
        for (j in seq_along(idx)) {
          i <- idx[j]
          avail <- M[i] - r[i]
          take  <- min(avail, cand[j] + leftover)
          m[i]  <- take
          leftover <- cand[j] + leftover - take
        }
        if (leftover > 0) {
          for (j in seq_along(idx)) {
            i <- idx[j]
            avail <- M[i] - r[i] - m[i]
            add   <- min(avail, leftover)
            m[i]  <- m[i] + add
            leftover <- leftover - add
            if (leftover == 0) break
          }
        }
      }
    }
  }
  
  # Expand to pseudo-IPD
  recs <- list()
  for (i in seq_len(n_int)) {
    if (dt[i] <= 0) next
    if (r[i] > 0) {
      recs[[length(recs) + 1L]] <- data.frame(
        time   = rep(t[i + 1], r[i]),
        status = 1L,
        arm    = arm_label
      )
    }
    if (m[i] > 0) {
      recs[[length(recs) + 1L]] <- data.frame(
        time   = runif(m[i], min = t[i] + 1e-8, max = t[i + 1] - 1e-8),
        status = 0L,
        arm    = arm_label
      )
    }
  }
  if (length(recs) == 0) {
    return(data.frame(time = numeric(0), status = integer(0), arm = character(0)))
  }
  do.call(rbind, recs)
}


## I: Create dataset for two-stage FP model ----
# Build an AggD-style dataset from IPD by: (1) computing KM per arm, (2) reconstructing
# pseudo-IPD via Guyot, (3) binning into fixed-length intervals with counts (r, m, M).
# --- Build AggD rows from IPD via KM -> Guyot -> binning (YEARS) --------
transform_study_dynamic <- function(ipd,
                                    study_index,
                                    interval_length,
                                    compar_des,
                                    nrisk_every = 3 * interval_length) {
  # 1) KM + risk table per arm, then Guyot pseudo-IPD (times in MONTHS here)
  arms <- levels(factor(ipd$treat.name))
  
  pseudo_list <- lapply(arms, function(a) {
    arm_ipd <- ipd[ipd$treat.name == a, , drop = FALSE]
    kn <- km_and_nrisk_from_ipd(arm_ipd, nrisk_every = nrisk_every)
    guyot_ipd_arm(kn$steps, kn$nrisk, arm_label = a)  # time (months), status, arm
  })
  pseudo <- do.call(rbind, pseudo_list)
  
  # 2) Convert to YEARS before binning
  pseudo_df <- data.frame(
    treat.name = factor(pseudo$arm, levels = arms),
    event.time = pseudo$time / 12,   # MONTHS -> YEARS
    censor     = pseudo$status
  )
  
  # 3) Bin into fixed-length intervals (YEARS)
  interval_counts_by_arm <- create_interval_counts_by_arm(
    data = pseudo_df,
    interval_length = interval_length
  )
  
  # 4) Map arms to A=1/B=2/C=3 using the study's ctrl/exper in compar_des
  arm_mapping <- setNames(c(1, 2),
                          c(compar_des[study_index, "ctrl"],
                            compar_des[study_index, "exper"]))
  a_code <- sapply(as.character(interval_counts_by_arm$ARM),
                   function(x) arm_mapping[[x]])
  
  data.frame(
    time = interval_counts_by_arm$start_time,                # years
    m    = interval_counts_by_arm$m,
    r    = interval_counts_by_arm$r,
    M    = interval_counts_by_arm$M,
    a    = a_code,
    dt   = interval_counts_by_arm$end_time - interval_counts_by_arm$start_time,
    s    = rep(study_index, nrow(interval_counts_by_arm)),
    z    = interval_counts_by_arm$M - interval_counts_by_arm$m
  )
}

## J: One-stage (fixed-effect or random-effects) NMA with fractional polynomials ----
# IPD-based FP model: trapezoidal integration over a fine grid, with ECOG effect
# included and study/arm-specific FP coefficients linked by treatment effects.
one_stage_FP_model <- function(n.studies) {
  
  stringcode <- "
  model {

    for (i in 1:n) {
      for (k in 1:steps) {

        tL[i,k] <- (k - 1) * dt
        tR[i,k] <- k * dt

        geL[i,k] <- step(time[i] - tL[i,k])
        geR[i,k] <- step(time[i] - tR[i,k])

        seg[i,k] <- geL[i,k] * ( (1 - geR[i,k]) * (time[i] - tL[i,k]) + geR[i,k] * (tR[i,k] - tL[i,k]) )
        tU[i,k]  <- tL[i,k] + seg[i,k]

        # FP basis
        t1L[i,k] <- equals(P1,0)*log(tL[i,k] + 1.0E-12) + (1 - equals(P1,0)) * pow(tL[i,k], P1)
        t1U[i,k] <- equals(P1,0)*log(tU[i,k] + 1.0E-12) + (1 - equals(P1,0)) * pow(tU[i,k], P1)

        t2L[i,k] <- (1 - equals(P2,P1)) * ( equals(P2,0)*log(tL[i,k] + 1.0E-12) + (1 - equals(P2,0)) * pow(tL[i,k], P2) ) +
                    equals(P2,P1) * ( equals(P2,0) * pow(log(tL[i,k] + 1.0E-12), 2) +
                    (1 - equals(P2,0)) * pow(tL[i,k], P2) * log(tL[i,k] + 1.0E-12) )

        t2U[i,k] <- (1 - equals(P2,P1)) * ( equals(P2,0)*log(tU[i,k] + 1.0E-12) + (1 - equals(P2,0)) * pow(tU[i,k], P2) ) +
                    equals(P2,P1) * ( equals(P2,0) * pow(log(tU[i,k] + 1.0E-12), 2) +
                    (1 - equals(P2,0)) * pow(tU[i,k], P2) * log(tU[i,k] + 1.0E-12) )

        loghL[i,k] <- Beta[index[i], arm[i], 1] +
                      Beta[index[i], arm[i], 2]*t1L[i,k] +
                      Beta[index[i], arm[i], 3]*t2L[i,k] +
                      beta_ecog * ecog[i]

        loghU[i,k] <- Beta[index[i], arm[i], 1] +
                      Beta[index[i], arm[i], 2]*t1U[i,k] +
                      Beta[index[i], arm[i], 3]*t2U[i,k] +
                      beta_ecog * ecog[i]

        hL[i,k] <- exp(loghL[i,k])
        hU[i,k] <- exp(loghU[i,k])

        cum_h_component[i,k] <- (hL[i,k] + hU[i,k]) * seg[i,k] / 2
      }

      cum_h[i] <- sum(cum_h_component[i,1:steps])

      t1i[i] <- equals(P1,0)*log(time[i] + 1.0E-12) + (1 - equals(P1,0)) * pow(time[i], P1)
      t2i[i] <- (1 - equals(P2,P1)) * ( equals(P2,0)*log(time[i] + 1.0E-12) + (1 - equals(P2,0)) * pow(time[i], P2) ) +
               equals(P2,P1) * ( equals(P2,0) * pow(log(time[i] + 1.0E-12), 2) +
               (1 - equals(P2,0)) * pow(time[i], P2) * log(time[i] + 1.0E-12) )

      log_h[i] <- Beta[index[i], arm[i], 1] +
                  Beta[index[i], arm[i], 2]*t1i[i] +
                  Beta[index[i], arm[i], 3]*t2i[i] +
                  beta_ecog * ecog[i]

      h[i] <- exp(log_h[i])

      nll[i] <- cum_h[i] - event[i] * log(h[i] + 1.0E-12)
      phi[i] <- nll[i] + 1.0E-6
      zeros[i] ~ dpois(phi[i])
    }

    # Study-specific baselines
    for (i in 1:num_studies) {
      delta[i, 1] <- 0

      for (k in 1:3) {
        mu[i, k] ~ dnorm(0, 0.0001)
      }

      for (j in 1:2) {
        Beta[i, j, 1] <- mu[i, 1] + delta[i, j]
        Beta[i, j, 2] <- mu[i, 2] + d[t[i, j], 2] - d[t[i, 1], 2]
        Beta[i, j, 3] <- mu[i, 3] + d[t[i, j], 3] - d[t[i, 1], 3]
      }
    }

    # Random effects only on PH/level component (as you had)
    tau  ~ dunif(0, 5)
    prec <- pow(tau, -2)

    for (i in 1:num_studies) {
      delta[i, 2] ~ dnorm(md[i, 2], prec)
      md[i, 2] <- d[t[i, 2], 1] - d[t[i, 1], 1]
    }

    beta_ecog ~ dnorm(0, 0.001)

    # ----- Shrinkage priors for time-varying treatment effects -----
    # Global SD for the log(t) treatment interaction (k=2)
    sd_time1 ~ dnorm(0, 1) T(0,)          # half-normal(0,1)
    prec_time1 <- pow(sd_time1, -2)

    # Stronger shrinkage for the quadratic term (k=3)
    sd_time2 ~ dnorm(0, 4) T(0,)          # half-normal(0,0.5): because prec=4 => sd=0.5
    prec_time2 <- pow(sd_time2, -2)

    # Treatment effects (reference-coded)
    for (k in 1:3) {
      d[1,k] <- 0
    }

    # Level (PH component): keep weakly-informative as before
    for (tr in 2:3) {
      d[tr,1] ~ dnorm(0, 0.0001)
    }

    # Time-varying components: shrink toward 0 unless data demand otherwise
    for (tr in 2:3) {
      d[tr,2] ~ dnorm(0, prec_time1)   # log(t) interaction
      d[tr,3] ~ dnorm(0, prec_time2)   # (log(t))^2 interaction (more shrinkage)
    }

    # Derived contrast (B vs C)
    for (k in 1:3) {
      logHR[k] <- d[2,k] - d[3,k]
    }
  }
  "
  stringcode
}


## K: Two-stage (fixed-effect or random-effects) NMA with fractional polynomials ----
# AggD FP model: poisson/binomial-like likelihood on interval counts with FP hazards.
# Beta[i,j,k] tie arm-specific FP coefficients to a common treatment effect.
two_stage_FP_model <- function(n.studies) {
  stringcode <- "
  model {
    # Likelihood on AggD
    for (i in 1:n) {
      timen1[i] <- equals(P1, 0)*log(time[i] + 1.0E-12) + (1 - equals(P1, 0))*pow(time[i], P1)
      timen2[i] <- (1 - equals(P2,P1)) * ( equals(P2,0)*log(time[i] + 1.0E-12) + (1 - equals(P2,0)) * pow(time[i], P2) ) +
                   equals(P2,P1) * ( equals(P2,0) * pow(log(time[i] + 1.0E-12), 2) +
                   (1 - equals(P2,0)) * pow(time[i], P2) * log(time[i] + 1.0E-12) )
      event[i] ~ dbin(p[i], M[i])
      p[i] <- max(0.0001, min(0.9999, 1 - exp(-h[i]*interval[i])))
      log(h[i]) <- Beta[index[i], arm[i], 1] +
                   Beta[index[i], arm[i], 2]*timen1[i] +
                   Beta[index[i], arm[i], 3]*timen2[i]
    }

    # Study-specific FP baselines
    for (i in 1:num_studies) {
      for (k in 1:3) {
        mu[i, k] ~ dnorm(0, 0.0001)
      }
    }

    # RE on the 'level' coefficient difference (delta), declare tau/prec first
    tau  ~ dunif(0, 5)
    prec <- pow(tau, -2)

    for (i in 1:num_studies) {
      delta[i, 1] <- 0
      for (j in 1:2) {
        Beta[i, j, 1] <- mu[i, 1] + delta[i, j]
        Beta[i, j, 2] <- mu[i, 2] + d[t[i, j], 2] - d[t[i, 1], 2]
        Beta[i, j, 3] <- mu[i, 3] + d[t[i, j], 3] - d[t[i, 1], 3]
      }
      delta[i, 2] ~ dnorm(md[i, 2], prec)
      md[i, 2] <- d[t[i, 2], 1] - d[t[i, 1], 1]
    }

    # Derived B vs C contrast for each FP term
    for (k in 1:3) {
      logHR[k] <- d[2, k] - d[3, k]
    }

    # Treatment effects (reference-coded)
    for (k in 1:3) {
      d[1, k] <- 0
      for (tr in 2:3) {
        d[tr, k] ~ dnorm(0, 0.0001)
      }
    }
  }
  "
  stringcode
}


## L: FRACTIONAL POLYNOMIAL NMA ----
# Wrapper to run the two-stage FP JAGS model. Assembles the required data list
# from the AggD table and extracts posterior summaries for logHR and d.
fp_nma <- function(MTCData, 
                   MTCAddData, 
                   n.studies,
                   P1, 
                   P2, 
                   n.chains, 
                   n.iter, 
                   n.burnin, 
                   n.thin) {
  
  r     <- c(MTCData[, "r"])
  z     <- c(MTCData[, "z"])
  a     <- c(MTCData[, "a"])
  s     <- c(MTCData[, "s"])
  time  <- c(MTCData[, "time"])
  dt    <- c(MTCData[, "dt"])
  N     <- length(dt)
  ns    <- length(unique(s))
  Mvec  <- c(MTCData[, "M"])
  
  t <- data.matrix(MTCAddData[, c("t1", "t2")])
  treatmeant.col <- factor(paste0(a, s), levels = unique(paste0(a, s)))
  treat <- rep(matrix(t(t), nrow = 1, ncol = 2 * length(t[, 1])), table(treatmeant.col))
  
  Data.multi <- list("event"       = r, 
                     "M"           = Mvec, 
                     "arm"         = a, 
                     "index"       = s, 
                     "time"        = time, 
                     "interval"    = dt, 
                     "n"           = N, 
                     "num_studies" = ns, 
                     "P1"          = P1, 
                     "P2"          = P2, 
                     "t"           = t)
  
  jagsfit <- jags(data               = Data.multi, 
                  parameters.to.save = c("logHR", "d"),
                  model.file         = textConnection(two_stage_FP_model(n.studies = n.studies)),
                  n.chains           = n.chains, 
                  n.iter             = n.iter, 
                  n.burnin           = n.burnin, 
                  n.thin             = n.thin, 
                  DIC                = FALSE, 
                  digits             = 4)
  
  # Summaries (as before)
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))
  logHR      <- t(getResults %>% dplyr::select(starts_with("logHR[")))
  d          <- t(getResults %>% dplyr::select(starts_with("d[")))
  
  # NEW: posterior draws of d  (iterations × 3 treatments × 3 FP coeffs)
  if (!is.null(jagsfit$BUGSoutput$sims$list)) {
    d_draws <- jagsfit$BUGSoutput$sims$list$d
  } else if (!is.null(jagsfit$BUGSoutput$sims.list)) {
    d_draws <- jagsfit$BUGSoutput$sims.list$d
  } else {
    smat <- jagsfit$BUGSoutput$sims.matrix
    arr  <- array(NA_real_, c(nrow(smat), 3, 3))
    for (tr in 1:3) for (k in 1:3) {
      arr[, tr, k] <- smat[, sprintf("d[%d,%d]", tr, k)]
    }
    d_draws <- arr
  }
  
  # Ensure orientation is (iterations, 3, 3)
  if (length(dim(d_draws)) == 3L && dim(d_draws)[1] == 3 && dim(d_draws)[2] == 3) {
    d_draws <- aperm(d_draws, c(3,1,2))
  }
  stopifnot(is.array(d_draws), dim(d_draws)[2] == 3, dim(d_draws)[3] == 3)
  
  list(logHR = logHR, d = d, d_draws = d_draws)
}
