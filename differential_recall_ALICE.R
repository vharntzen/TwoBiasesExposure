
source("fun_estimate_single_interval_censored_ALICE.R") # Expects L0, L1, R0


require('survival')
require('tidyr')
require('magrittr')
require('dplyr')

n = 50000

# Function
Simulate_scenario <- function(distanceCutpoints = 1, n = 10000, n_runs = 2,
                              p.exact = 0.1,
                              true_dist = "Weibull", par_1 = 2.453, par_2 = 6.258,
                              assumed_dist = "Weibull",
                              cut_off_subset = 5, method = "B",
                              perc_observed = 100, diff_recall_rate = 1/4){
  
  out <- array(NA, dim = c(5, 3, n_runs))
  out_subset <- array(NA, dim = c(5, 3, n_runs))
  
  rate <- diff_recall_rate
  p_exact <- p.exact
  
  for(i in 1:n_runs){
    
    
    if(true_dist == "Weibull"){
      mydata  <- data.frame(time = rweibull(n, shape = par_1, scale = par_2)) #%>% round(0)
    } else if(true_dist == "gamma"){
      mydata <- data.frame(time = rgamma(n, shape = par_1, scale = par_2)) #%>% round(0)
    } else if(true_dist == "lognormal"){
      mydata <- data.frame(time = rlnorm(n, meanlog = par_1, sdlog = par_2)) #%>% round(0)
    } else if(true_dist == "exponential"){
      mydata <- data.frame(time = rexp(n, rate = par_1))} #%>% round(0)
    
    if(method %in% c("A", "B")){
      obstimes <- matrix(rep(seq(1,20,distanceCutpoints), n), nrow=20) # Vera
    } else if(method == "R"){
      obstimes <- matrix(rexp(n = 10*n, rate = rate), nrow=10) # Ronald
    }
    
    obstimes <- rbind(matrix(rep(0,n), nrow=1), obstimes)
    obstimes <- apply(obstimes, 2, sort)
    
    # Select the 'sandwiching' points.
    
    for (obs in 1:n){
      
      if(method == "R"){  
        # Method R(onald): P = 1  
        c.vec <- obstimes[, obs][!is.na(obstimes[, obs])]
        
      } else if(method == "A"){
        # Method A:
        O.vec <- rbinom(length(obstimes[,obs]), 1, prob = 1 - pexp(mydata$time[obs], rate = rate) ) 
        c.vec <- obstimes[, obs][!is.na(obstimes[, obs])]
        c.vec <- c.vec[O.vec == 1]
        
      } else if(method == "B"){
        # Method B: P(observation depends on time of observation)    
        O.vec <- rbinom(length(obstimes[,obs]), 1, prob = 1 - pexp(obstimes[, obs], rate = rate) ) 
        c.vec <- obstimes[, obs][!is.na(obstimes[, obs])]
        c.vec <- c.vec[O.vec == 1]
      }
      
      if(sum(mydata$time[obs] <= c.vec) %in% c(0,NA)){
        mydata$upper[obs] <- NA
      } else {
        mydata$upper[obs] <- min(c.vec[ mydata$time[obs] <= c.vec ], na.rm = TRUE) 
      }
      
      if(sum(mydata$time[obs] >= c.vec) %in% c(0,NA)){
        mydata$lower[obs] <- NA
      } else {
        mydata$lower[obs] <- max(c.vec[ mydata$time[obs] >= c.vec ], na.rm = TRUE) 
      }
    }
    
    if(p_exact > 0){
      exact_obs <- sample(1:n, size = round(p_exact*n), replace = FALSE)
      mydata$upper[exact_obs] <- mydata$time[exact_obs]
      mydata$lower[exact_obs] <- mydata$time[exact_obs]
    }
    
    dat <- mydata
    
      dat_subset <- dat %>%
      drop_na(upper) %>%
      subset( (upper-lower) < cut_off_subset )
    
    # NB ADD THE COUNT OF NAs
    
    
      temp.mat <- try(Estimate_single_mle_survival(dat, distribution = assumed_dist
    ) %>% as.matrix())
    if(class(temp.mat)[1] == "try-error"){
      out[,,i] <- matrix(NA, nrow = 5, ncol = 3)
    }else{
      out[,,i] <- temp.mat
    }
    
    temp.mat_subset <- try(Estimate_single_mle_survival(dat_subset, distribution = assumed_dist
                                                       
    ) %>% as.matrix())
    if(class(temp.mat_subset)[1] == "try-error"){
      out_subset[,,i] <- matrix(NA, nrow = 5, ncol = 3)
    }else{
      out_subset[,,i] <- temp.mat_subset
    }
    
  } # End of loop over runs
  
  if(true_dist == "Weibull"){
    real_p <- qweibull(percentiles, shape = par_1, scale = par_2)
  } else if(true_dist == "exponential"){
    real_p <- qexp(percentiles, #shape = 2.453, scale = 6.258
                   rate = par_1)
  } else if(true_dist == "lognormal"){
    real_p <- qlnorm(percentiles, #shape = 2.453, scale = 6.258
                     meanlog = par_1, sdlog = sdlog)
  }
  
  
  out.summary <- data.frame(nrow = 5, ncol = 11)
  out.summary_subset <- data.frame(nrow = 5, ncol = 11)
  
  for(p in 1:5){
    
    out.summary[p,1] <- distanceCutpoints
    out.summary[p,2] <- mean(out[p, 1, 1:n_runs] - real_p[p], na.rm = TRUE)
    out.summary[p,c(3,4)] <- quantile(out[p, 1, 1:n_runs] - real_p[p], c(0.25,0.75), na.rm = TRUE)
    #out.summary_subset[p,5] <- mean(real_p[p] >= out_subset[p, 2, 1:n_runs] & real_p[p] <= out_subset[p, 3, 1:n_runs], na.rm = TRUE)
    out.summary[p,5] <- perc_observed
    out.summary[p,6] <- method
    out.summary[p,7] <- diff_recall_rate
    out.summary[p,8] <- true_dist
    out.summary[p,9] <- assumed_dist
    out.summary[p,10] <- NA
    out.summary[p,11] <- p.exact
  }
  
  colnames(out.summary) <- c("dC", "Bias", "Deviations_p25", "Deviations_p75", "Percentage_cutpoints_observed", "Diff_recall_method", "Diff_recall_rate", "True_dist", "Assumed_dist", "Exposure_dist", "P_exact")
  rownames(out.summary) <- percentiles
  
  out.summary$Subset <- "all"
  out.summary$percentile <- percentiles
  
  for(p in 1:5){
    
    out.summary_subset[p,1] <- distanceCutpoints
    out.summary_subset[p,2] <- mean(out_subset[p, 1, 1:n_runs] - real_p[p], na.rm = TRUE)
    out.summary_subset[p,c(3,4)] <- quantile(out_subset[p, 1, 1:n_runs] - real_p[p], c(0.25,0.75), na.rm = TRUE)
     out.summary_subset[p,5] <- perc_observed
    out.summary_subset[p,6] <- method
    out.summary_subset[p,7] <- diff_recall_rate
    out.summary_subset[p,8] <- true_dist
    out.summary_subset[p,9] <- assumed_dist
    out.summary_subset[p,10] <- NA
    out.summary_subset[p,11] <- p.exact
  }
  
  colnames(out.summary_subset) <- c("dC", "Bias", "Deviations_p25", "Deviations_p75", "Percentage_cutpoints_observed", "Diff_recall_method", "Diff_recall_rate", "True_dist", "Assumed_dist", "Exposure_dist", "P_exact")
  rownames(out.summary_subset) <- percentiles
  
  
  out.summary_subset$Subset <- "subset"
  out.summary_subset$percentile <- percentiles
  
  
  
  out_both <- rbind(out.summary, out.summary_subset)
  
  out_both
  
}

# # exponential
# true_dist_preset <- "exponential"
# par_1_preset = 1/4; par_2_preset = NA
# assumed_dist_preset = "exponential"

# Weibull
true_dist_preset <- "Weibull"
par_1_preset = 2.453; par_2_preset = 6.258
assumed_dist_preset = "Weibull"

# Scenarios

#### Define the parameters ####
shape = 2.453; scale = 6.258
#shape = 1.5; scale = 1.5

rates <- seq(0.1, 2, 0.1)
p_exact <- c(0, 0.1)
method <- c("A", "B", "R")
n_runs <- 1000
dC <- 1
n <- 10000
##############################

scenarios <- expand.grid(rates, p_exact, method, n_runs, dC, n)
colnames(scenarios) <- c("rate", "p_exact", "method", "n_runs", "dC", "n")

# Run the simulation

saveRDS(scenarios, file = paste("/results/Scenarios.RDS", sep = ""))

for(run in 1:nrow(scenarios)){
  
  diff_recall_rate <- scenarios[run, "rate"]
  p_exact <- scenarios[run, "p_exact"]
  method <- scenarios[run, "method"]
  n_runs <- scenarios[run, "n_runs"]
  dC <- scenarios[run, "dC"]
  n <- scenarios[run, "n"]
  
  out <- Simulate_scenario(n = n, 
                           n_runs = n_runs, 
                           distanceCutpoints = dC,
                           true_dist = true_dist_preset, 
                           p.exact = p_exact,
                           par_1 = par_1_preset, par_2 = par_2_preset,
                           diff_recall_rate = diff_recall_rate,
                           assumed_dist = assumed_dist_preset,
                           method = method)
  
  saveRDS(out, file = paste("/results/Run", run, ".RDS", sep = ""))
  
  cat("Now starting run ", run + 1, " of ", nrow(scenarios), ".\n")
  
}

