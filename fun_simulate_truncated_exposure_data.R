
# temp <- Simulate_data_windows(n = 100, 
#                       type = "single", 
#                       window_lengths = c(1, 5), 
#                       window_props = c(0.05, 0.95),
#                       distribution_T = "Weibull", 
#                       par_1 = 1, par_2 = 7, exposure_dist = "unif")

# _____________________________________________________________________________

# Function to generate single- or doubly interval-censored data
# using window length.

Simulate_data_windows_truncated <- function(
                                  n = 10, 
                                  window_width = 3, 
                                  p.exact = 0.3,
                                  distribution_T = "Weibull", 
                                  par_1 = 1, par_2 = 7, exposure_dist = "unif",
                                  g = 0.14
                                  ){

  # Inputs

  # Output

  # Packages
  require(dplyr)
  require(tidyr)
# require(poweRlaw)
  require(rmutil)  # Burr
#  source("functions_Generate_from_truncated_distribution.R")
  
  mat <- matrix(nrow = n, ncol = 8)
  colnames(mat) <- c("Li", "Ri", "Ti", "L0", "L1", "R0", "R1", "dC")
  
  n_exact <- round(p.exact * n)
  n_windows <- n - n_exact
  
  if(n_exact > 0){
    for (i in 1:n_exact){
    
      # Generate left window
      L0 <- 0
      
      L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
      
      mat[i, 4:5] <- c(L0, L1)
      
      # Generate L, R and T
      
      if(exposure_dist == "exp_growth"){
        a <- 1  # a starting value; g growth factor
        x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))
        Li <- (1/g)*log(x/a) # inverse of cumulative incidence
        
      } else if(exposure_dist == "unif"){
        Li <- runif(1, min = L0, max = L1)  
        
      } else if(exposure_dist == "exp_household"){
        p <- 0.2
        x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
        Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
      }

        if(distribution_T == "Weibull"){
          Ti <- rweibull(1, shape = par_1, scale = par_2)
        } else if(distribution_T == "gamma"){
          Ti <- rgamma(1, shape = par_1, scale = par_2)
        } else if(distribution_T == "lognormal"){
          Ti <- rlnorm(1, meanlog = par_1, sdlog = par_2)
        } else if(distribution_T == "heavytail"){
          Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
        }
        
        Ri <- Li + Ti
      
      mat[i, c(1:3)]  <- c(Li, Ri, Ti)
      mat[i, c(4, 5)] <- Li
      
    } # End of loop per individual.
  } # End of loop over exact observations.
  
  for (i in (n_exact + 1):n){ # For every interval censored observation
    
    # Generate left window
    L0 <- 0
  
    L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
    
    mat[i, 4:5] <- c(L0, L1)
    
    # Generate L, R and T
    
    if(exposure_dist == "exp_growth"){
    a <- 1  # a starting value; g growth factor
    x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))  # draw from uniform, 
    Li <- (1/g)*log(x/a) # inverse of cumulative incidence
   
    } else if(exposure_dist == "unif"){
    Li <- runif(1, min = L0, max = L1)  
    
    } else if(exposure_dist == "exp_household"){
    p <- 0.2
    x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
    Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
    }

      if(distribution_T == "Weibull"){
      Ti <- rweibull(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "gamma"){
      Ti <- rgamma(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "lognormal"){
      Ti <- rlnorm(n = 1, meanlog = par_1, sdlog = par_2) #%>% round(0)
    } else if(distribution_T == "heavytail"){
      Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
    }
    
    Ri <- Li + Ti
    
    mat[i, c(1:3)]  <- c(Li, Ri, Ti)
    
} # End of loop per individual.

# Generate/specify R0 and R1    

  mat[, "R0"] <-  mat[, "R1"] <- mat[, "Ri"] # Ri = R0 = R1

  mat[,8] <- window_width

# Exclude the observations for which symptoms occur before end of exposure.
indices_truncation <- which(mat[, "Ri"] <= mat[, "L1"])
mat <- mat[-indices_truncation, ]

    
  return(mat %>% as.data.frame)
  
  # End of function   
  
}

# ------------------------------------------------------------------------------ 

# 
# temp <- Simulate_data_windows(n = 100,
#                       type = "single",
#                       window_lengths = c(1, 5),
#                       window_props = c(0.05, 0.95),
#                       distribution_T = "Weibull",
#                       par_1 = 1, par_2 = 7, exposure_dist = "unif")

# _____________________________________________________________________________

# Function to generate single- or doubly interval-censored data
# using window length.

Simulate_data_windows_p_for_forgotten_start <- function(
    n = 10, 
    window_width = 3, 
    p.exact = 0,
    p.forgotten_start = 0.1,
    distribution_T = "Weibull", 
    par_1 = 1, par_2 = 7, exposure_dist = "unif",
    g = 0.14
){
  
  # Inputs
  
  # Output
  
  # Packages
  require(dplyr)
  require(tidyr)
  # require(poweRlaw)
  require(rmutil)  # Burr
  source("functions_Generate_from_truncated_distribution.R")
  
  mat <- matrix(nrow = n, ncol = 8)
  colnames(mat) <- c("Li", "Ri", "Ti", "L0", "L1", "R0", "R1", "dC")
  
  n_exact <- round(p.exact * n)
  n_windows <- n - n_exact
  
  if(n_exact > 0){
    for (i in 1:n_exact){
      
      # Generate left window
      L0 <- 0
      
      L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
      
      mat[i, 4:5] <- c(L0, L1)
      
      # Generate L, R and T
      
      if(exposure_dist == "exp_growth"){
        a <- 1  # a starting value; g growth factor
        x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))
        Li <- (1/g)*log(x/a) # inverse of cumulative incidence
        
      } else if(exposure_dist == "unif"){
        Li <- runif(1, min = L0, max = L1)  
        
      } else if(exposure_dist == "exp_household"){
        p <- 0.2
        x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
        Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
      }
      
      if(distribution_T == "Weibull"){
        Ti <- rweibull(1, shape = par_1, scale = par_2)
      } else if(distribution_T == "gamma"){
        Ti <- rgamma(1, shape = par_1, scale = par_2)
      } else if(distribution_T == "lognormal"){
        Ti <- rlnorm(1, meanlog = par_1, sdlog = par_2)
      } else if(distribution_T == "heavytail"){
        Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
      }
      
      Ri <- Li + Ti
      
      mat[i, c(1:3)]  <- c(Li, Ri, Ti)
      mat[i, c(4, 5)] <- Li
      
    } # End of loop per individual.
  } # End of loop over exact observations.
  
  for (i in (n_exact + 1):n){ # For every observation
    
    # Generate left window
    L0 <- 0
    
    L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
    
    mat[i, 5] <- L1
    
    # Set the exposure start to NA with probability p.forgotten_start
    mat[i, 4] <- ifelse(rbinom(n = 1, size = 1, prob = p.forgotten_start) == 1, NA, L0)
    
    # Generate L, R and T
    
    if(exposure_dist == "exp_growth"){
      a <- 1  # a starting value; g growth factor
      x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))  # draw from uniform, 
      Li <- (1/g)*log(x/a) # inverse of cumulative incidence
      
    } else if(exposure_dist == "unif"){
      Li <- runif(1, min = L0, max = L1)  
      
    } else if(exposure_dist == "exp_household"){
      p <- 0.2
      x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
      Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
    }
    
    if(distribution_T == "Weibull"){
      Ti <- rweibull(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "gamma"){
      Ti <- rgamma(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "lognormal"){
      Ti <- rlnorm(n = 1, meanlog = par_1, sdlog = par_2) #%>% round(0)
    } else if(distribution_T == "heavytail"){
      Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
    }
    
    Ri <- Li + Ti
    
    mat[i, c(1:3)]  <- c(Li, Ri, Ti)
    
  } # End of loop per individual.
  
  # Generate/specify R0 and R1    
  
  mat <- mat %>% as.data.frame()
  
  mat[, 6:7] <-  mat[, 2] # Ri = R0 = R1
  
  # If exposure ends after symptom onset, set to symptom onset.
  mat[ , "L1"] <- ifelse(mat[ ,"L1"] > mat[ , "Ri"], mat[ , "Ri"], mat[ , "L1"])
  
  mat[,8] <- window_width
  
  return(mat %>% as.data.frame)
  
  # End of function   
  
}

# ------------------------------------------------------------------------------ 


Simulate_data_windows_p_start_end_forgotten <- function(
    n = 10, 
    window_width = 3, 
    p.exact = 0,
    p.forgotten = 0.1,
    distribution_T = "Weibull", 
    par_1 = 1, par_2 = 7, exposure_dist = "unif",
    g = 0.14
){
  
  # Inputs
  
  # Output
  
  # Packages
  require(dplyr)
  require(tidyr)
  # require(poweRlaw)
  require(rmutil)  # Burr
  source("functions_Generate_from_truncated_distribution.R")
  
  mat <- matrix(nrow = n, ncol = 8)
  colnames(mat) <- c("Li", "Ri", "Ti", "L0", "L1", "R0", "R1", "dC")
  
  n_exact <- round(p.exact * n)
  n_windows <- n - n_exact
  
  if(n_exact > 0){
    for (i in 1:n_exact){
      
      # Generate left window
      L0 <- 0
      
      L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
      
      mat[i, 4] <-  L0
      mat[i, 5] <-  L1
      
      # Generate L, R and T
      
      if(exposure_dist == "exp_growth"){
        a <- 1  # a starting value; g growth factor
        x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))
        Li <- (1/g)*log(x/a) # inverse of cumulative incidence
        
      } else if(exposure_dist == "unif"){
        Li <- runif(1, min = L0, max = L1)  
        
      } else if(exposure_dist == "exp_household"){
        p <- 0.2
        x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
        Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
      }
      
      if(distribution_T == "Weibull"){
        Ti <- rweibull(1, shape = par_1, scale = par_2)
      } else if(distribution_T == "gamma"){
        Ti <- rgamma(1, shape = par_1, scale = par_2)
      } else if(distribution_T == "lognormal"){
        Ti <- rlnorm(1, meanlog = par_1, sdlog = par_2)
      } else if(distribution_T == "heavytail"){
        Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
      }
      
      Ri <- Li + Ti
      
      mat[i, c(1:3)]  <- c(Li, Ri, Ti)
      mat[i, c(4, 5)] <- Li
      
    } # End of loop per individual.
  } # End of loop over exact observations.
  
  for (i in (n_exact + 1):n){ # For every observation
    
    # Generate left window
    L0 <- 0
    
    L1 <- L0 + window_width  # interval L0 to L1 is between length 1 and 5.
    
    # Generate L, R and T
    
    if(exposure_dist == "exp_growth"){
      a <- 1  # a starting value; g growth factor
      x <- runif(1, min = exp(g*(L0)), max = exp(g*(L1)))  # draw from uniform, 
      Li <- (1/g)*log(x/a) # inverse of cumulative incidence
      
    } else if(exposure_dist == "unif"){
      Li <- runif(1, min = L0, max = L1)  
      
    } else if(exposure_dist == "exp_household"){
      p <- 0.2
      x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
      Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
    }
    
    if(distribution_T == "Weibull"){
      Ti <- rweibull(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "gamma"){
      Ti <- rgamma(n = 1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "lognormal"){
      Ti <- rlnorm(n = 1, meanlog = par_1, sdlog = par_2) #%>% round(0)
    } else if(distribution_T == "heavytail"){
      Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
    }
    
    Ri <- Li + Ti
    
    mat[i, c(1:3)]  <- c(Li, Ri, Ti)
    
    
    if(rbinom(n = 1, size = 1, prob = p.forgotten) == 1){
      
      if(rbinom(n = 1, size = 1, prob = 0.5) == 1){
        # Set the exposure start to NA with probability p.forgotten/2
        mat[i, 4] <- NA
        mat[i, 5] <- L1
        
      } else {
        # Set the exposure end to NA with probability p.forgotten/2
        # NB: if window end is missing, but in reality ends before symptom
        # onset, the observation should be omitted. REMOVE
        
        #if(L1 < Ri){ mat[i, 4] <- NA
     # } else {  
          mat[i, 4] <- L0# }
        
        mat[i, 5] <- NA
      }
      
    }
    
    
  } # End of loop per individual.
  
  # Generate/specify R0 and R1    
  
  mat <- mat %>% as.data.frame()
  
  mat[, 6:7] <-  mat[, 2] # Ri = R0 = R1
  
  # Remove individuals with exposure start and end NA
  indices_to_remove <- which(is.na(mat[ , "L0"]) & is.na(mat[ , "L1"]))
  mat <- mat[-indices_to_remove,]
  
  # If exposure ends after symptom onset, set to symptom onset.
  mat[ , "L1"] <- ifelse(mat[ ,"L1"] > mat[ , "Ri"], mat[ , "Ri"], mat[ , "L1"])
  
  # If exposure end is unknown, set to symptom onset (only when end was not before symptom onset,
  # which is taken care of automatic now)
  mat[ , "L1"] <- ifelse(is.na(mat[ ,"L1"]), mat[ , "Ri"], mat[ , "L1"])
  
  mat[,8] <- window_width
  
  return(mat %>% as.data.frame)
  
  # End of function   
  
}

