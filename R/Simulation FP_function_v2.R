#*******************************************************************************
#*
#*
#*             Simulation Function for All Scenarios & Bayesian Models   
#*                            <Fractional polynomial>                            
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
#' @param model Character with possible values \code{"uninfcens"}, \code{"infcens"}, and
#'   \code{"delay"} to indicate proportional hazards with non-informative censoring,
#'   informative censoring, and delayed treatment effect with non-informative censoring, respectively.
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
set.seed(20250819)
run_simulation_FP <- function (lambda, 
                               gamma, 
                               time0, 
                               comorb.prob, 
                               model, 
                               P1,
                               P2,
                               n.chains, 
                               n.iter, 
                               n.burnin, 
                               n.thin) {
  
  
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
  if (!is.numeric(lambda) || length(lambda)!=1L || lambda <= 0) {
    stop("'lambda' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma)!=1L || gamma <= 0) {
    stop("'gamma' must be a single positive number.", call. = FALSE)
  }
  if (identical(model, "delay") && (is.null(time0) || !is.numeric(time0) || time0 <= 0)) {
    stop("For model='delay', provide a positive 'time0' (e.g., 180).", call. = FALSE)
  }
  if (identical(model, "infcens") && (is.null(comorb.prob) || comorb.prob < 0 || comorb.prob > 1)) {
    stop("For model='infcens', 'comorb.prob' must be in [0,1].", call. = FALSE)
  }
  
  # Define the pairwise comparisons and their base log(HR)
  compar <- data.frame(
    ctrl = c("A", "A", "C"),
    exper = c("B", "C", "B"),
    true_LHR = c(log(0.6), log(0.8), log(0.6 / 0.8))
  )
  
  # Repeat the rows: 3 A-B, 2 A-C, 1 B-C
  compar_des0 <- compar[rep(seq_len(nrow(compar)), times = c(3, 2, 1)), ]
  
  # Set target tau? (between-study variance)
  tau2 <- 0.04
  tau <- sqrt(tau2)
  
 
  
  # Add random noise to all study-level log HRs
  compar_des <- data.frame(
    ctrl = compar_des0$ctrl,
    exper = compar_des0$exper,
    true_LHR = compar_des0$true_LHR + rnorm(nrow(compar_des0), mean = 0, sd = tau)
  )
  
  # Mean "true" logHR per comparison
  true_means <- with(compar_des, tapply(true_LHR, paste0(exper, ctrl), mean))
  
  Study <- c("OC","NC","NO")
  
  true_map <- c(
    NC = unname(true_means["BA"]),  # B vs A
    OC = unname(true_means["CA"]),  # C vs A
    NO = unname(true_means["BC"])   # B vs C  (new vs old)
  )
  
  True_LHR <- unname(true_map[Study])
  
  n_st <- nrow(compar_des)
  sizes <- round(runif(n_st, 180, 634), 0)
  
  # Staggered START YEAR per study (e.g., 2012-2016)
  starts <- sample(2012:2016, n_st, replace = TRUE)
  
  # Study duration in YEARS (short early cutoffs + some long studies)
  dur_years <- sample(c(1, 1.5, 2, 3, 4), n_st, replace = TRUE, prob = c(.25, .25, .25, .20, .05))
  stops <- floor(starts + dur_years)  # administrative stop year
  
  # Accrual window length (years): short = concentrated pre/post-delay, long = broad mix
  accrual_years <- sample(c(0.5, 1, 2, 3), n_st, replace = TRUE, prob = c(.3, .4, .2, .1))
  
  # Experimental arm opens later in some studies (days)
  arm_lag_days <- sample(c(0, 60, 120, 180), n_st, replace = TRUE, prob = c(.5, .2, .2, .1))
  
  
  # Generate per-study IPD
  create_study <-
    if (model == "uninfcens") {
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
  
  if (model == "delay") {
    # helper to bind studies for a given comparison label like "BA","CA","BC"
    bind_comp <- function(lbl) {
      # map labels to ctrl/exper
      map <- list(
        BA = c(ctrl = "A", exper = "B"),  # NC: B vs A
        CA = c(ctrl = "A", exper = "C"),  # OC: C vs A
        BC = c(ctrl = "C", exper = "B")   # NO: B vs C  <-- CHANGED from CB
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
    
    # per-contrast stratified Cox (study strata), returns logHR or NA
    get_lhr <- function(lbl) {
      d <- bind_comp(lbl)
      if (is.null(d) || sum(d$censor == 1) == 0) return(NA_real_)
      cf <- coef(summary(coxph(Surv(event.time, censor) ~ treat + ecog + strata(study_id), data = d)))
      unname(cf["treat", "coef"])
    }
    
    # Named map; reorder later to match your output row order
    true_lhr_map <- c(
      NC = get_lhr("BA"),  # B vs A
      OC = get_lhr("CA"),  # C vs A
      NO = get_lhr("BC")   # B vs C  <-- CHANGED
    )
    True_LHR <- unname(true_lhr_map[Study])
  }
  
  # FP basis at time t that matches your JAGS code
  fp_terms <- function(t, P1, P2, eps = 1e-12){
    t <- pmax(t, eps)
    t1 <- if (P1 == 0) log(t) else t^P1
    if (P2 == P1) {
      t2 <- (if (P1 == 0) log(t) else t^P1) * log(t)  # repeated power
    } else {
      t2 <- if (P2 == 0) log(t) else t^P2
    }
    cbind(t1 = t1, t2 = t2)
  }
  
  # Inputs:
  #   combined_data: your list used for JAGS (contains time, event, index, arm, t, dt, steps)
  #   study_ids: integer vector of studies contributing the contrast
  #   tr_codes: numeric vector length 2 with the LETTERS codes of (a,b), e.g. c(1,3) for A,C
  # Output: a list with w (length=steps), dt, steps
  surv_weights_ipd <- function(combined_data, study_ids, tr_codes) {
    # Map each IPD row to global treatment code (1=A,2=B,3=C) via (study, arm)
    tr_per_row <- combined_data$t[cbind(combined_data$index, combined_data$arm)]
    
    keep <- combined_data$index %in% study_ids & tr_per_row %in% tr_codes
    df <- data.frame(
      time  = combined_data$time[keep],
      event = combined_data$event[keep],
      tr    = tr_per_row[keep]
    )
    
    K  <- combined_data$steps
    dt <- combined_data$dt
    t_mid <- ((1:K) - 0.5) * dt
    
    # Safe KM extractor for one arm
    get_S <- function(d) {
      if (!nrow(d)) return(rep(1, K))
      # small jitter to avoid exact zeros
      d$time <- pmax(d$time, 1e-10)
      sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
      # summary(..., times=..., extend=TRUE) returns length(t_mid)
      out <- suppressWarnings(summary(sf, times = t_mid, extend = TRUE)$surv)
      if (length(out) != length(t_mid)) {
        # belt-and-suspenders fallback
        ssum <- summary(sf)
        out <- approx(ssum$time, ssum$surv, xout = t_mid, rule = 2)$y
      }
      as.numeric(out)
    }
    
    Sa <- get_S(subset(df, tr == tr_codes[1]))
    Sb <- get_S(subset(df, tr == tr_codes[2]))
    
    Sa_step <- c(1, Sa); Sb_step <- c(1, Sb)
    w_raw <- 0.5*(pmax(-diff(Sa_step),0) + pmax(-diff(Sb_step),0))
    w <- if (sum(w_raw)>0) w_raw / sum(w_raw) else rep(0, length(w_raw))
    
    list(w = w, dt = dt, steps = K)
  }
  
  # Build event-weighted average Delta from posterior draws of d
  avg_delta_from_draws <- function(d_draws, P1, P2, dt, steps, w, contrast){
    t_mid <- ((1:steps) - 0.5) * dt
    TT <- fp_terms(t_mid, P1, P2); t1 <- TT[,1]; t2 <- TT[,2]
    w <- if (sum(w) > 0) w / sum(w) else rep(1/steps, steps)
    
    a <- contrast[1]; b <- contrast[2]
    d1 <- d_draws[, b, 1] - d_draws[, a, 1]
    d2 <- d_draws[, b, 2] - d_draws[, a, 2]
    d3 <- d_draws[, b, 3] - d_draws[, a, 3]
    
    # draws x steps
    Delta <- outer(d1, rep(1, steps)) + outer(d2, t1) + outer(d3, t2)
    Delta_bar <- as.numeric(Delta %*% w)
    Delta_bar <- Delta_bar[is.finite(Delta_bar)]  # be safe
    
    med  <- stats::median(Delta_bar, na.rm = TRUE)
    lCri <- as.numeric(stats::quantile(Delta_bar, 0.025, names = FALSE, na.rm = TRUE))
    uCri <- as.numeric(stats::quantile(Delta_bar, 0.975, names = FALSE, na.rm = TRUE))
    c(med = med, lCrI = lCri, uCrI = uCri)
  }
  
  #'****************************************************************************
  #'                   ENTER ONE-STAGE NMA FP SPECIFICATION               
  #'****************************************************************************
  
  # Prepare data for each study
  data_comp <-
    lapply(1:length(create_study), 
           function(x) data.frame(time = create_study[[x]]$event.time,            # Random survival times
                                  treat = create_study[[x]]$treat.name,                   # Treatment assigned
                                  treat.arm = ifelse(create_study[[x]]$treat == 1, 2, 1),                    # Treatment arm (control: 0, exper: 1)
                                  event = create_study[[x]]$censor,                       # Censoring indicator
                                  index = rep(x, dim(create_study[[x]])[1]),              # Study ID
                                  ecog = create_study[[x]]$ecog))
  
  time_years  <- unlist(lapply(data_comp, function(x) x$time / 12))   # months -> YEARS
  dt_years    <- 2/12
  steps_years <- max(2L, ceiling(max(time_years, na.rm = TRUE) / dt_years))  # at least 2 bins

  combined_data <- list(
    n           = sum(unlist(lapply(create_study, nrow))),
    time        = time_years,                                   # <-- years
    arm         = unlist(lapply(data_comp, `[[`, "treat.arm")),
    event       = unlist(lapply(data_comp, `[[`, "event")),
    index       = unlist(lapply(data_comp, `[[`, "index")),
    t           = apply(compar_des[, -3], 2, function(x) match(x, LETTERS)),
    num_studies = nrow(compar_des),
    P1          = P1,
    P2          = P2,
    dt          = dt_years,                                     # <-- years
    steps       = steps_years,                                  # <-- computed with years
    ecog        = unlist(lapply(data_comp, `[[`, "ecog")),
    zeros       = rep(0, sum(unlist(lapply(create_study, nrow))))
  )
  
  # Run JAGS model (ONE-STAGE MODEL)
  run_model <- 
    jags(data = combined_data,
         parameters.to.save = c("d", "logHR"),
         model.file = textConnection(one_stage_FP_model(n.studies = dim(compar_des)[1])),
         n.chains = n.chains,
         n.iter = n.iter,
         n.burnin = n.burnin,
         n.thin = n.thin,
         DIC = TRUE)
  
  # Convert JAGS output to MCMC object
  samples <- as.mcmc(run_model)
  
  # Extract the MCMC samples from the jags_model
  samples_matrix <- as.matrix(samples)
  
  # --- draws of d[*,*] to build Delta(t) ---
  pick_ipd <- paste0("d[", rep(1:3, each = 3), ",", rep(1:3, times = 3), "]")
  dmat_ipd <- samples_matrix[, pick_ipd, drop = FALSE]
  d_draws_ipd <- array(NA_real_, dim = c(nrow(dmat_ipd), 3, 3))
  col_map_ipd <- matrix(pick_ipd, nrow = 3, byrow = TRUE)
  for (tr in 1:3) for (k in 1:3) d_draws_ipd[, tr, k] <- dmat_ipd[, col_map_ipd[tr, k]]
  
  # study id sets per contrast (as you defined compar_des)
  idx_BA <- which(compar_des$ctrl == "A" & compar_des$exper == "B")  # NC: B vs A
  idx_CA <- which(compar_des$ctrl == "A" & compar_des$exper == "C")  # OC: C vs A
  idx_BC <- which(compar_des$ctrl == "C" & compar_des$exper == "B")  # NO: B vs C
  
  # event weights per interval from the IPD
  w_BA_ipd <- surv_weights_ipd(combined_data, study_ids = idx_BA, tr_codes = c(1,2)) # A vs B
  w_CA_ipd <- surv_weights_ipd(combined_data, study_ids = idx_CA, tr_codes = c(1,3)) # A vs C
  w_BC_ipd <- surv_weights_ipd(combined_data, study_ids = idx_BC, tr_codes = c(3,2)) # C vs B
  
  
  # ---- Build an event-weighted TRUTH that matches the FP estimand ----
  # helper: per-(study,interval) event table from IPD
  ipd_df <- data.frame(
    time  = combined_data$time,
    index = combined_data$index,
    event = combined_data$event
  )
  # Map each IPD row to its grid bin k (1..steps), using the same units as dt and time
  k_from_time <- function(t, dt, K) {
    k <- floor((pmax(t, 0) - 1e-12) / dt) + 1L   # left-closed bins; exact multiples stay in the left bin
    k <- pmax(1L, pmin(K, k))                    # clamp to [1, K]
    k
  }
  
  # ipd_df must already have: time, index, event
  k_idx <- k_from_time(ipd_df$time, combined_data$dt, combined_data$steps)
  ipd_df$k <- k_idx
  
  # true ??_i(t) for a study i at time t_mid
  # convert time0 (days) -> months to match t_mid units
  time0_y <- time0 / 365.25
  
  true_delta_at_time <- function(delta_i, t_mid_y, model, time0_y){
    if (identical(model, "delay")) ifelse(t_mid_y < time0_y, 0, delta_i) else delta_i
  }
  
  true_evAvg_ipd <- function(study_ids){
    df <- subset(ipd_df, index %in% study_ids & event == 1)
    if (nrow(df) == 0) return(NA_real_)
    
    tab <- tapply(df$event, list(df$index, df$k), sum)
    tab <- as.matrix(tab); tab[is.na(tab)] <- 0
    if (ncol(tab) < combined_data$steps) {
      tab <- cbind(tab, matrix(0, nrow(tab), combined_data$steps - ncol(tab)))
      colnames(tab) <- 1:combined_data$steps
    }
    
    t_mid_y <- ((1:combined_data$steps) - 0.5) * combined_data$dt   # years
    deltas  <- compar_des$true_LHR[study_ids]
    Delta_true <- outer(deltas, t_mid_y, function(d, t) true_delta_at_time(d, t, model, time0_y))
    
    num <- sum(tab * Delta_true); den <- sum(tab)
    if (den == 0) NA_real_ else num / den
  }
  
  True_evAvg_OC <- true_evAvg_ipd(idx_CA)  # C vs A
  True_evAvg_NC <- true_evAvg_ipd(idx_BA)  # B vs A
  True_evAvg_NO <- true_evAvg_ipd(idx_BC)  # B vs C
  
  # event-weighted averages (order below will be OC, NC, NO to match your return)
  ipd_OC_ev <- avg_delta_from_draws(d_draws_ipd, P1, P2,
                                    w_CA_ipd$dt, w_CA_ipd$steps, w_CA_ipd$w, c(1,3))
  ipd_NC_ev <- avg_delta_from_draws(d_draws_ipd, P1, P2,
                                    w_BA_ipd$dt, w_BA_ipd$steps, w_BA_ipd$w, c(1,2))
  ipd_NO_ev <- avg_delta_from_draws(d_draws_ipd, P1, P2,
                                    w_BC_ipd$dt, w_BC_ipd$steps, w_BC_ipd$w, c(3,2))
  
  # Extract posterior median of ds (log hazard ratios)
  LHR_d21 <- median(samples_matrix[, "d[2,1]"])   # 2 vs 1
  LHR_d31 <- median(samples_matrix[, "d[3,1]"])   # 3 vs 1
  LHR_d32      <- median(samples_matrix[, "logHR[1]"])
  
  # Extract posterior standard error of ds 
  seLHR_d21 <- sd(samples_matrix[, "d[2,1]"])   
  seLHR_d31 <- sd(samples_matrix[, "d[3,1]"])   
  seLHR_d32    <- sd(samples_matrix[, "logHR[1]"])
  
  # Calculate the 95% credible interval (2.5th and 97.5th percentiles)
  LHR_d21_cri <- quantile(samples_matrix[, "d[2,1]"], probs = c(0.025, 0.975))
  LHR_d31_cri <- quantile(samples_matrix[, "d[3,1]"], probs = c(0.025, 0.975))
  LHR_d32_cri  <- quantile(samples_matrix[, "logHR[1]"], probs = c(0.025, 0.975))
  
  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(run_model$BUGSoutput$summary))
  
  # Extract R-hat
  R_hat <- range(t(get_results)[, "Rhat"])

  
  #'****************************************************************************
  #'                   ENTER AGGREAGATE NMA FP MODEL SPECIFICATION                 
  #'****************************************************************************
  
  # Process all studies dynamically
  data_test_list <- lapply(seq_along(create_study), function(i) {
    transform_study_dynamic(create_study[[i]],
                            study_index = i,
                            interval_length = 2/12,  # YEARS
                            compar_des = compar_des)
  })
  
  # Combine all transformed datasets
  data_test <- do.call(rbind, data_test_list)
  
  # Define the matrix with the comparisons per trial
  AddData <- apply(compar_des[, -3], 2, function(x) match(x, LETTERS))
  
  # Rename the columns
  colnames(AddData) <- c("t1", "t2")
  
  # Inputs:
  #   MTCData: the AggD table (must contain columns: s (study id), arm (1/2), time, r (events), M (at risk))
  #   MTCAddData: AddData matrix (n_studies x 2) mapping (study, arm)-> global tr code 1=A,2=B,3=C
  #   interval_length: dt
  #   tr_codes: vector length 2 with the global codes for arms a,b (e.g., c(1,3))
  #   keep_studies: integer vector of studies contributing the contrast
  surv_weights_aggd <- function(MTCData, MTCAddData, interval_length, tr_codes, keep_studies) {
    df <- subset(MTCData, s %in% keep_studies)
    if (!nrow(df)) return(list(w = numeric(0), dt = interval_length, steps = 0))
    
    # map each row (study,arm) -> global treatment code (1=A,2=B,3=C)
    df$tr <- MTCAddData[cbind(df$s, df$arm)]
    
    # index grid cells (assuming df$time are left edges: 0,2,4,...)
    idx <- floor(df$time / interval_length) + 1L
    K   <- max(1L, max(idx, na.rm = TRUE))   # never 0
    
    S_arm <- function(target_tr) {
      sel  <- df$tr == target_tr
      at_k <- tapply(df$r[sel], idx[sel], sum)
      nk   <- tapply(df$M[sel], idx[sel], sum)
      
      rvec <- rep(0, K); nvec <- rep(0, K)
      if (length(at_k)) rvec[as.integer(names(at_k))] <- as.numeric(at_k)
      if (length(nk))   nvec[as.integer(names(nk))]   <- as.numeric(nk)
      
      # KM by cells (not used for weighting now, but handy for checks)
      S <- numeric(K); S[1] <- 1
      for (k in 1:K) {
        hk <- if (nvec[k] > 0) rvec[k] / nvec[k] else 0
        S[k] <- if (k == 1) (1 - hk) else S[k - 1] * (1 - hk)
      }
      list(S = S, rvec = rvec, nvec = nvec)
    }
    
    Sa <- S_arm(tr_codes[1])
    Sb <- S_arm(tr_codes[2])
    
    # ---- event-weighted pooling (NOT S(t) weighting) ----
    e_pool <- 0.5 * (Sa$rvec + Sb$rvec)     # events per bin, pooled across the two arms
    w_raw  <- e_pool
    w      <- if (sum(w_raw) > 0) as.numeric(w_raw / sum(w_raw)) else rep(0, K)
    
    list(w = w, dt = interval_length, steps = K)
  }
  
  # A=1, B=2, C=3 (from your AddData coding)
  w_CA_aggd <- surv_weights_aggd(data_test, AddData, interval_length = 2/12,
                                 tr_codes = c(1,3), keep_studies = idx_CA)
  w_BA_aggd <- surv_weights_aggd(data_test, AddData, interval_length = 2/12,
                                 tr_codes = c(1,2), keep_studies = idx_BA)
  w_BC_aggd <- surv_weights_aggd(data_test, AddData, interval_length = 2/12,
                                 tr_codes = c(3,2), keep_studies = idx_BC)
  
  aggd_out <- fp_nma(MTCData = data_test, 
                     MTCAddData = AddData, 
                     n.studies = dim(compar_des)[1],
                     P1 = P1, P2 = P2, 
                     n.chains = n.chains, n.iter = n.iter, 
                     n.burnin = n.burnin, n.thin = n.thin)
  
  # 1) Sanity checks + extract draws
  stopifnot("d_draws" %in% names(aggd_out))
  d_draws_aggd <- aggd_out$d_draws
  stopifnot(is.array(d_draws_aggd), dim(d_draws_aggd)[2] == 3, dim(d_draws_aggd)[3] == 3)
  
  # 2) Compute event-weighted averages for each contrast using those draws
  aggd_OC_ev <- avg_delta_from_draws(d_draws_aggd, P1, P2,
                                     w_CA_aggd$dt, w_CA_aggd$steps, w_CA_aggd$w,
                                     c(1,3))  # C vs A
  aggd_NC_ev <- avg_delta_from_draws(d_draws_aggd, P1, P2,
                                     w_BA_aggd$dt, w_BA_aggd$steps, w_BA_aggd$w,
                                     c(1,2))  # B vs A
  aggd_NO_ev <- avg_delta_from_draws(d_draws_aggd, P1, P2,
                                     w_BC_aggd$dt, w_BC_aggd$steps, w_BC_aggd$w,
                                     c(3,2))  # B vs C
  
  # summaries (unchanged, just from aggd_out now)
  LHR_d21_aggd <- aggd_out$d["d[2,1]", "50%"]
  LHR_d31_aggd <- aggd_out$d["d[3,1]", "50%"]
  LHR_d32_aggd <- aggd_out$logHR["logHR[1]", "50%"]
  
  seLHR_d21_aggd <- aggd_out$d["d[2,1]", "sd"]
  seLHR_d31_aggd <- aggd_out$d["d[3,1]", "sd"]
  seLHR_d32_aggd <- aggd_out$logHR["logHR[1]", "sd"]
  
  LHR_aggd_d21_cri <- c(aggd_out$d["d[2,1]", "2.5%"], aggd_out$d["d[2,1]", "97.5%"])
  LHR_aggd_d31_cri <- c(aggd_out$d["d[3,1]", "2.5%"], aggd_out$d["d[3,1]", "97.5%"])
  LHR_aggd_d32_cri <- c(aggd_out$logHR["logHR[1]", "2.5%"], aggd_out$logHR["logHR[1]", "97.5%"])
  
  
  
  # IPD event-weighted average (median and CrI)
  d_IPD_evAvg       = c(as.numeric(ipd_OC_ev["med"]),  as.numeric(ipd_NC_ev["med"]),  as.numeric(ipd_NO_ev["med"]))
  d_IPD_evAvg_lCrI  = c(as.numeric(ipd_OC_ev["lCrI"]), as.numeric(ipd_NC_ev["lCrI"]), as.numeric(ipd_NO_ev["lCrI"]))
  d_IPD_evAvg_uCrI  = c(as.numeric(ipd_OC_ev["uCrI"]), as.numeric(ipd_NC_ev["uCrI"]), as.numeric(ipd_NO_ev["uCrI"]))
  
  # AggD event-weighted average (median and CrI)
  d_AggD_evAvg      = c(as.numeric(aggd_OC_ev["med"]),  as.numeric(aggd_NC_ev["med"]),  as.numeric(aggd_NO_ev["med"]))
  d_AggD_evAvg_lCrI = c(as.numeric(aggd_OC_ev["lCrI"]), as.numeric(aggd_NC_ev["lCrI"]), as.numeric(aggd_NO_ev["lCrI"]))
  d_AggD_evAvg_uCrI = c(as.numeric(aggd_OC_ev["uCrI"]), as.numeric(aggd_NC_ev["uCrI"]), as.numeric(aggd_NO_ev["uCrI"]))

  ## --- structural (data-generating) step truth per contrast (mean ?? across studies) ---
  true_step_OC <- mean(compar_des$true_LHR[idx_CA])  # C vs A
  true_step_NC <- mean(compar_des$true_LHR[idx_BA])  # B vs A
  true_step_NO <- mean(compar_des$true_LHR[idx_BC])  # B vs C
  
  ## -------- assemble result (and ATTACH DRAWS SAFELY) --------
  res <- data.frame(
    Study       = c("OC","NC","NO"),
    Model       = rep(model, 3),
    True_LHR    = True_LHR,
    
    d_IPD       = c(LHR_d31, LHR_d21, LHR_d32),
    d_IPD_lCrI  = c(LHR_d31_cri[1], LHR_d21_cri[1], LHR_d32_cri[1]),
    d_IPD_uCrI  = c(LHR_d31_cri[2], LHR_d21_cri[2], LHR_d32_cri[2]),
    d_IPD_se    = c(seLHR_d31, seLHR_d21, seLHR_d32),
    
    d_AggD      = c(LHR_d31_aggd, LHR_d21_aggd, LHR_d32_aggd),
    d_AggD_lCrI = c(LHR_aggd_d31_cri[1], LHR_aggd_d21_cri[1], LHR_aggd_d32_cri[1]),
    d_AggD_uCrI = c(LHR_aggd_d31_cri[2], LHR_aggd_d21_cri[2], LHR_aggd_d32_cri[2]),
    d_AggD_se   = c(seLHR_d31_aggd, seLHR_d21_aggd, seLHR_d32_aggd),
    
    d_IPD_evAvg      = c(as.numeric(ipd_OC_ev["med"]),  as.numeric(ipd_NC_ev["med"]),  as.numeric(ipd_NO_ev["med"])),
    d_IPD_evAvg_lCrI = c(as.numeric(ipd_OC_ev["lCrI"]), as.numeric(ipd_NC_ev["lCrI"]), as.numeric(ipd_NO_ev["lCrI"])),
    d_IPD_evAvg_uCrI = c(as.numeric(ipd_OC_ev["uCrI"]), as.numeric(ipd_NC_ev["uCrI"]), as.numeric(ipd_NO_ev["uCrI"])),
    
    d_AggD_evAvg      = c(as.numeric(aggd_OC_ev["med"]),  as.numeric(aggd_NC_ev["med"]),  as.numeric(aggd_NO_ev["med"])),
    d_AggD_evAvg_lCrI = c(as.numeric(aggd_OC_ev["lCrI"]), as.numeric(aggd_NC_ev["lCrI"]), as.numeric(aggd_NO_ev["lCrI"])),
    d_AggD_evAvg_uCrI = c(as.numeric(aggd_OC_ev["uCrI"]), as.numeric(aggd_NC_ev["uCrI"]), as.numeric(aggd_NO_ev["uCrI"])),
    
    True_evAvg    = c(True_evAvg_OC, True_evAvg_NC, True_evAvg_NO),
    bias_IPD_ev   = c(ipd_OC_ev["med"], ipd_NC_ev["med"], ipd_NO_ev["med"]) - c(True_evAvg_OC, True_evAvg_NC, True_evAvg_NO),
    bias_AggD_ev  = c(aggd_OC_ev["med"], aggd_NC_ev["med"], aggd_NO_ev["med"]) - c(True_evAvg_OC, True_evAvg_NC, True_evAvg_NO),
    
    ## NEW: structural step truth (for plotting HR(t) step) and the time0 you used
    True_step_LHR = c(true_step_OC, true_step_NC, true_step_NO),
    time0_days    = if (identical(model,"delay")) rep(time0, 3) else rep(NA_real_, 3)
  )
  
  #attr(res, "d_draws_ipd")  <- d_draws_ipd
  #attr(res, "d_draws_aggd") <- d_draws_aggd
  return(list(res = res,
              d_draws_ipd = d_draws_ipd,
              d_draws_aggd = d_draws_aggd))
}
