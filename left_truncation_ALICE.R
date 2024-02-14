
source("fun_simulate_truncated_exposure_data.R")

require("dplyr")
require("tidyr") # expand_grid
require("survival")
require("MixtureRegLTIC")

# Function by Hadley Wickham to silence cat() texts.
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

Fit_AFT <- function(
    n = 500,
    exposure_dist = "exp_growth",#"unif", # 
    true_dist = "Weibull", #"lognormal"#  ## 
    par_1 = 2.453, par_2 = 6.258,
    assumed_dist = "Weibull", #"lognormal" # 
    window_width = 60,
    p_exact = 0.1,
    g = 2,
    truncation_in_analysis = TRUE
){
  
 
  dat <- Simulate_data_windows_truncated(
    n = n,
    window_width = window_width,
    p.exact = p_exact,
    distribution_T = true_dist,
    par_1 = par_1, par_2 = par_2,
    exposure_dist = exposure_dist,
    g = g
  )
  
  percentiles <- c(0.5, 0.9, 0.95, 0.975, 0.99)
  
  shortest <- dat$R0 - dat$L1;
  longest <- dat$R0 - dat$L0
  
  dat_m <- data.frame(shortest = shortest, longest = longest, 
                      time.trunc = ifelse(dat$L0 == dat$L1, 0, dat$L1 - dat$Li), # Checked: same effect as NULL.
                      
                      #  For interval censored data, the status indicator is 3.
                      # Other: 0 = right censored, 1 = event at time, 
                      # 2 = left censored, 3 = interval censored.                    
                      status = ifelse(dat$L0 == dat$L1, 1, 3))
  
  
  dat_m <- dat_m[sample(1:nrow(dat_m), nrow(dat_m), replace = FALSE),]
  
  
  if(assumed_dist == "Weibull"){
    shape_for_dist <- 1 # for Weibull
  } else if (assumed_dist == "lognormal"){
    shape_for_dist <- 0 # for lognormal
  }
  
  if(truncation_in_analysis == TRUE){
    var.entry <- "time.trunc"} else {
      var.entry <- NULL
    }
  
  # Fit the AFT location-scale model.
  
  fit <- try(
    
    quiet(
      MixtureLogitAFT(
        
        # formula: same as Surv function. Status indicator 3 for interval-censored.
        formula = Surv(shortest, longest, status) ~ 1,
        shape = shape_for_dist,
        # to exclude cure part:
        eventprobreg = NULL, 
        # var.entry: specifies each study subject's left-truncated time at entry in the follow-up study.
        var.entry = var.entry, # Left truncation time.
        # data: a data.frame with the variables named in the formula and var.entry section (and in 
        data = dat_m
        
      )
    )
    
  )
  
  if( class(fit) == "try-error" ){ cat("AFT model gave error.")
    invalid_run <- TRUE
  } else {
    if( fit$convergence != 0 ){ cat("AFT model did not converge.")
      invalid_run <- TRUE
    } else {
      invalid_run <- FALSE
      
      
      # Collect parameter estimates generalized gamma distribution.
      gamma <- fit$par$gamma; alpha <- fit$par$alpha
      
      # __________ Weibull distribution ____________
      
      if(assumed_dist == "Weibull"){
        # Parameter estimates
        
        shape = exp(-alpha) # shape = exp(-alpha)
        scale = exp(gamma) # scale = exp(gamma)
        
        ests_pars <- c(shape, scale)
        
        # Quantiles
        
        ests_quantiles <- qweibull(p = percentiles, 
                                   
                                   shape = exp(-alpha), 
                                   scale = exp(gamma) 
                                   
        )
      } else if(assumed_dist == "lognormal"){
        
        # __________ Lognormal distribution __________
        
        # Parameter estimates
        
        ( meanlog = gamma ) # mu = gamma
        ( sdlog = exp(alpha) ) # sd = exp(alpha)
        
        ests_pars <- c(meanlog, sdlog)
        
        # Quantiles
        
        ests_quantiles <- qlnorm(p = percentiles, 
                                 
                                 meanlog = gamma, # mu = gamma
                                 sdlog = exp(alpha) # sd = exp(alpha)
                                 
        )
      }
      # ____________________________________________
      
      if (true_dist == "Weibull"){
        true_p <- qweibull(p = percentiles, shape = par_1, scale = par_2)
      } else if(true_dist == "lognormal"){
        true_p <- qlnorm(p = percentiles, meanlog = par_1, sdlog = par_2)
      } else if(true_dist == "gamma"){
        true_p <- qgamma(p = percentiles, shape = par_1, rate = par_2)
      }
      
      deviation <- ests_quantiles - true_p
      
    } } # Skip until here if the model did not converge or gave error.
  
  if(invalid_run == TRUE){ out <- list(true_dist = true_dist, 
                                       assumed_dist = assumed_dist,
                                       p_exact = p_exact,
                                       par_1 = NA,
                                       par_2 = NA,
                                       quantiles = rep(NA, length(percentiles)), 
                                       ests_q50 = NA, 
                                       ests_q95 = NA, 
                                       ests_q99 = NA,
                                       error_q50 = NA, 
                                       error_q95 = NA, 
                                       error_q99 = NA,
                                       error_par_1 = NA,
                                       error_par_2 = NA,
                                       error_quantiles = rep(NA, length(percentiles)))
  
  } else { out <- list(true_dist = true_dist, 
                       assumed_dist = assumed_dist,
                       p_exact = p_exact,
                       par_1 = ests_pars[1],
                       par_2 = ests_pars[2],
                       quantiles = ests_quantiles, 
                       ests_q50 = ests_quantiles[1], 
                       ests_q95 = ests_quantiles[2], 
                       ests_q99 = ests_quantiles[3],
                       error_q50 = deviation[1], 
                       error_q95 = deviation[2], 
                       error_q99 = deviation[3],
                       error_par_1 = ests_pars[1] - par_1,
                       error_par_2 = ests_pars[2] - par_2,
                       error_quantiles = deviation)
  }
  
  out
  
}

# _______________________________________________________

Run_scenario <- function(n_runs = 10,
                         n = 500,
                         exposure_dist = "unif", # "exp_growth"
                         true_dist = "Weibull", #"lognormal"
                         par_1 = 2.453, par_2 = 6.258,
                         assumed_dist = "Weibull", #"lognormal"
                         window_width = 5,
                         p_exact = 0.1,
                         g = 2,
                         truncation_in_analysis = TRUE
){
  
  array_out <- replicate(n_runs, Fit_AFT(
    n = n,
    exposure_dist = exposure_dist,
    true_dist = true_dist, 
    par_1 = par_1, par_2 = par_2,
    assumed_dist = assumed_dist,
    window_width = window_width,
    p_exact = p_exact,
    g = g,
    truncation_in_analysis = truncation_in_analysis
  ))
  
  out <- list(n,
              exposure_dist,
              true_dist, 
              par_1, par_2,
              assumed_dist, 
              window_width,
              p_exact = p_exact,
              mean( unlist(array_out["error_q50", ]), na.rm = T), # Bias q50
              mean( unlist(array_out["error_q95", ]), na.rm = T), # Bias q95
              mean( unlist(array_out["error_q99", ]), na.rm = T),  # Bias q99
              mean( unlist(array_out["error_par_1", ]), na.rm = T), # Bias par_1
              mean( unlist(array_out["error_par_2", ]), na.rm = T), # Bias par_1# Bias par_2
              mean( is.na(array_out[10, ]) ) * 100 # Percentage erronous/ not converged iterations
  )
  
  out
  
}



# ------------------- # RUN SIMULATIONS # ------------------- #

####################################################
n_runs <- 10
n <- 5000
exposure_dist <- c("exp_growth","unif", "exp_household")
true_dist <- c("Weibull" #"Weibull"#"Weibull" # "lognormal""
)

par_1 <- NA; par_2 <- NA
assumed_dist <- c("Weibull"#"lognormal"#"Weibull"#, "lognormal"
)
window_width <- c(seq(1,30,2))#, 70,85,100)
p_exact <- 0.1
g <- 0.14 # growth factor exponential growth
truncation_in_analysis <- TRUE#FALSE
####################################################

# Combine all options and select rows.
scenarios <- tidyr::expand_grid(n, n_runs, exposure_dist,
                                true_dist, par_1, par_2,
                                assumed_dist,
                                window_width,
                                p_exact, g,   truncation_in_analysis)

colnames(scenarios) <- c("n", "n_runs", "exposure_dist", "true_dist", "par_1", "par_2", 
                         "assumed_dist", "window_width", "p_exact", "g", "truncation_in_analysis")

scenarios <- scenarios %>% filter(true_dist == assumed_dist)

#### Parameters #######################################################
scenarios[which(scenarios$true_dist == "Weibull"), ]$par_1 <- 2.453
scenarios[which(scenarios$true_dist == "Weibull"), ]$par_2 <- 6.258
#scenarios[which(scenarios$true_dist == "lognormal"), ]$par_1 <- 1.621
#scenarios[which(scenarios$true_dist == "lognormal"), ]$par_2 <- 0.418
#######################################################################

#scenarios <- scenarios[1:5, ]

n_scenarios <- nrow(scenarios)

column_names <- c( 
  "n",
  "exposure_dist",
  "true_dist", 
  "par_1", 
  "par_2",
  "assumed_dist", 
  "window_width",
  "p_exact",
  "Bias_q50",
  "Bias_q95",
  "Bias_q99",
  "Bias_par_1",
  "Bias_par_2",
  "Perc_error"
)

dat_results <- data.frame(matrix(nrow = n_scenarios, ncol = length(column_names)))
colnames(dat_results) <- column_names

for(run in 1:n_scenarios){
  
  cat("Starting run ", run, " of ", n_scenarios, ".\n", sep = "")
  
  dat_results[run, ] <- 
    Run_scenario(n_runs = scenarios$n_runs[run],
                 n = scenarios$n[run],
                 exposure_dist = scenarios$exposure_dist[run],
                 true_dist = scenarios$true_dist[run],
                 par_1 = scenarios$par_1[run], scenarios$par_2[run],
                 assumed_dist = scenarios$assumed_dist[run],
                 window_width = scenarios$window_width[run],
                 p_exact = scenarios$p_exact[run],
                 g = scenarios$g[run],
                 truncation_in_analysis = scenarios$truncation_in_analysis[run]
    )
  
  saveRDS(dat_results[run, ], file = paste("/results/Left_trunc_Run", run, ".RDS", sep = ""))
  
}

saveRDS(dat_results, file = paste("/results/Left_trunc_All.RDS", sep = ""))

