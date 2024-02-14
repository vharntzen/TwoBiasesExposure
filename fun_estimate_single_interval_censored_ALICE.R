
# Estimate_single_mle_survival

percentiles <- c(0.5, 0.9, 0.95, 0.975, 0.99)
row_names_percentiles <- c("q50", "q90", "q95", "q97.5", "q99")


#percentiles <- c(0.025, 0.050, 0.500, 0.950, 0.975)

# using nonparametric approach in survival package.
Estimate_single_mle_survival <- function(dat, distribution = "lognormal", par_estimates = FALSE){  
  
  # The second approach is to think of each observation as a time interval with 
  #(-infinity, t2) for left censored, (t1, infinity) for right censored, 
  #(t,t) for exact and (t1, t2) 
  
  
  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1
  # alternatively: lower and upper interval
  
  if("lower" %in% colnames(dat)){
    shortest <- dat$lower
    longest <- dat$upper
  } else{
  shortest <- dat$R0 - dat$L1
  longest <- dat$R0 - dat$L0
  }
  
  require('survival')
  require('tidyr')
  require('magrittr')
  require('dplyr')
  
  tab <- as.data.frame(matrix(nrow = length(percentiles), ncol = 3))
  
  # MLE
  # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
 shortest[is.na(shortest)] <- 0.000001
  shortest[shortest == 0] <- shortest[shortest == 0] + 0.000001 # To make the function work.
longest[is.na(longest)] <- Inf
  
  #shortest <- pmax(shortest,min(shortest)/2)
  
  if(distribution == "Weibull"){
    fit <- survreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = "weibull")
    summary(fit)
    
  } else if(distribution == "lognormal"){ # gamma not available.
    fit <- survreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = "lognormal")
  }
  
  if(par_estimates == TRUE){
    
    if(distribution == "lognormal"){
      params <- c(fit$coefficients, fit$scale)
      names(params) <- c("meanlog", "sdlog")  
      
    } else if(distribution == "Weibull"){
      params <- c(1/fit$scale, exp(fit$coefficients))
      names(params) <- c("shape", "scale")
    }
    
    params} else{
      
      for(perc in 1:length(percentiles)){ # Based on https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf page 85.
        q2 <- predict(fit, type = 'uquantile',
                      p = percentiles[perc], se.fit = T)
        ci2 <- cbind(q2$fit, q2$fit - 1.96*q2$se.fit, q2$fit + 1.96*q2$se.fit)
        tab[perc,] <- exp(ci2[1,]) # convert from log scale to original y
      }
      
      rownames(tab) <- row_names_percentiles
      colnames(tab) <- c("est", "lower_CI", "upper_CI")
      
      tab
      #  list(tab, params)
    }
  
}
