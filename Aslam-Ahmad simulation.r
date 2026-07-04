# Load required packages
rm(list=ls())
library(survival)
library(MASS)
# Function to calculate Mean Squared Error (MSE)
output <- function(Me, Tb) {
  Tb_row <- matrix(Tb, nrow = 1, ncol = length(Tb))
  SM <- sweep(Me, 2, Tb_row, "-")^2
  col_means_sq_error <- colMeans(SM, na.rm = TRUE)
  MSE <- mean(col_means_sq_error, na.rm = TRUE)
  MSE}
Cox_AA <- function(n, p, rho, tau, CP) {
  
m <- 1000# the number of replications
  
  # True value of the parameters
  Tbeta1 <- c(0.5, 1, -0.6, 0.8)
  Tbeta2 <- rep(0, times = (p - 4))
  Tbeta <- c(Tbeta1, Tbeta2)
  
  # Deviation values from the true null hypothesis (H*beta = h)
  delta <- seq(0, 4, by = 0.5) # the deviation from true values
  
  # Restriction matrix H (for H*beta = h)
  H <- matrix(0, nrow = (p - 4), ncol = p)
  H[, 5:p] <- diag(1, nrow = (p - 4), ncol = (p - 4))
  
  # Matrix to store Relative Total MSE results
  RTMSE <- matrix(nrow = length(delta), ncol = 8)
  names_out <- c("MPL", "AA", "RAA", "AALS", "AAS", "AAPS", "AAPT", "AASPT")
  colnames(RTMSE) <- names_out
  
  
  for (l in 1:length(delta)) {
    set.seed(303) 
    h <- rep(0, times = (p - 4)) 
    h[1] <- delta[l]    # setting diparture value 
    
    # Matrices to store estimates for each replication
    MPLE <- matrix(ncol = p, nrow = m)
    AAE <- matrix(ncol = p, nrow = m)
    RAAE <- matrix(ncol = p, nrow = m)   # restricted
    AALSE <- matrix(ncol = p, nrow = m) # linear shrinkage
    AASE <- matrix(ncol = p, nrow = m)   # stein
    AAPSE <- matrix(ncol = p, nrow = m)  # positive stein
    AAPTE <- matrix(ncol = p, nrow = m)  # pretest
    AASPTE <- matrix(ncol = p, nrow = m) # shrinkage pretest
    
    # REPLICATION STARTS
    for (u in 1:m) {
      # Generating the observations of covariates
      Z <- matrix(rnorm(n * (p + 1)), nrow = n, ncol = (p + 1))
      
      X <- matrix(0, nrow = n, ncol = p)
      for (j in 1:p) {
        X[, j] <- sqrt(1 - rho^2) * Z[, j] + rho * Z[, p + 1]
      }
      
      # Generate survival times and censoring indicators
      linear_predictor <- X %*% Tbeta
      rates <- tau * exp(linear_predictor)
      # Handle potential overflow if exp(linear_predictor) is huge
      rates[rates > 1e10] <- 1e10 # Cap rates to avoid Inf/-Inf times
      time <- rexp(n, rate = rates)
      
      # Simple random censoring (user's original method)
      status <- sample(c(0, 1), n, replace = TRUE, prob = c(CP, 1 - CP))
      
      
      # data frame for the generated data
      data <- data.frame(time = time,status = status, X)
      cox_fit <- coxph(Surv(time, status) ~ X, data = data)
      MPLE[u,]<- coef(cox_fit)
      
      # Deriving information matrix (precision matrix)
      InDTD <- vcov(cox_fit)
      DTD <- solve(InDTD) # Information matrix 
      
      H_InDTD_Ht <- H %*% InDTD %*% t(H)
      Inv_H_InDTD_Ht <- solve(H_InDTD_Ht)
      
      # Deviation vector for restriction
      DTV <- H %*% MPLE[u,] - h
      
      # Restricted estimator (RMPLE)
      RMPLE <- MPLE[u,] - InDTD %*% t(H) %*% Inv_H_InDTD_Ht %*% DTV
      
      # Test statistic Tn (Wald test)
      Tn <- crossprod(DTV, Inv_H_InDTD_Ht) %*% DTV 
      Tn <- as.numeric(Tn) 
      
      # Degrees of freedom for the test
      df_h <- length(h)   
      
      # Multiplier for Stein-type estimators
      Multi <- 1 - ((df_h - 2) / Tn)
      
      # Critical value for pre-test
      qq <- qchisq(0.05, df_h, lower.tail = FALSE)
      
      # Pre-test indicator (1 if H0 not rejected, 0 if rejected)
      TTest <- ifelse(Tn <= qq, 1, 0)
      
      # Eigenvalues and eigenvectors of DTD (Information Matrix)
      eigen_DTD <- eigen(DTD, symmetric = TRUE)
      lam <- eigen_DTD$values
      # Ensure eigenvalues are positive (numerical stability)
      lam[lam < 1e-8] <- 1e-8 
      W <- eigen_DTD$vectors
      I <- diag(p)
      
      lambda <- eigen(DTD)$values
      alpha_mp <- coef(cox_fit)
      
      # --- Parameter selection algorithm (Section 4) ---
      # Step 1: Compute d_max
      d_max <- max((lambda * (alpha_mp^2 - 1)) / (1 + lambda * alpha_mp^2))
      
      # Step 2: Compute k_med
      k_j <- ((lambda + d_max) - lambda * (1 - d_max) * alpha_mp^2) / 
        ((1 + d_max) * (lambda + 1) * alpha_mp^2)
      k_med <- median(k_j)
      
      # Step 3: Compute d_opt (fixed vectorization error)
      numerator <- sum(
        (lambda + k_med * (lambda + 1)) * (lambda^2 - k_med * lambda^2 - k_med * lambda) * alpha_mp^2 -
          lambda * (lambda^2 + k_med * lambda^2 - k_med * lambda)
        / ((lambda + 1)^2))
      
      denominator <- sum(
        (lambda^2 + k_med * lambda^2 - k_med * lambda) +
          (lambda - k_med * (lambda + 1)) * (lambda^2 - k_med * lambda^2 - k_med * lambda) * alpha_mp^2
        / ((lambda + 1)^2))
      
      d_opt <- numerator / denominator
      
      # Correction for d_opt (scalar condition)
      if (d_opt < 0 || d_opt > 1) {
        d_opt <- d_max
      }
      
      # --- Compute AAE (Equation 6) ---
      term1 <- solve(DTD + diag(p))
      term2 <- DTD + d_opt * diag(p)
      term3 <- solve(DTD + k_med * (1 + d_opt) * diag(p))
      M <- term1 %*% term2 %*% term3 %*% DTD
      
      # Aslam Ahmad estimator (AA)
      AAE[u,] <- M %*% MPLE[u,]
      
      # Restricted Aslam Ahmad estimator (RAA)
      RAAE[u,] <- M %*% RMPLE   
      
      diff_estim <- AAE[u,] - RAAE[u,]
      
      # Linear shrinkage estimators (AALS)
      AALSE[u, ] <- AAE[u,] - 0.5 * diff_estim 
      
      # Stein estimator / Shrinkage (AAS)
      AASE[u, ] <- AAE[u,] - Multi * diff_estim
      
      # Positive Stein estimator (AAPS)
      AAPSE[u,] <- AAE[u,] - max(0, Multi) * diff_estim
      
      # Pretest estimator (AAPT)
      AAPTE[u,] <- AAE[u,] - TTest * diff_estim
      
      # Shrinkage pretest estimators (AASPT)
      AASPTE[u,] <- AAE[u,] - (0.5 * TTest) * diff_estim
      
    } # End of inner loop (replications u)
    
    # Calculate MSE for each estimator type, ignoring NAs from failed runs
    MSE_MPLE <- output(MPLE, Tbeta)
    MSE_AAE <- output(AAE, Tbeta)
    MSE_RAAE <- output(RAAE, Tbeta)
    MSE_AALSE <- output(AALSE, Tbeta)
    MSE_AASE <- output(AASE, Tbeta)
    MSE_AAPSE <- output(AAPSE, Tbeta)
    MSE_AAPTE <- output(AAPTE, Tbeta)
    MSE_AASPTE <- output(AASPTE, Tbeta)
    
    # Store MSEs
    TMSE <- c(MSE_MPLE, MSE_AAE, MSE_RAAE, MSE_AALSE,
              MSE_AASE, MSE_AAPSE, MSE_AAPTE, MSE_AASPTE)
    
    # Calculate Relative MSE (relative to AA estimator)
    RTMSE[l, ] <- MSE_AAE / TMSE
    
  } # End of outer loop (delta l)
  
  round(RTMSE, 4)
}

####example

#result1 <- Cox_AA(n = 50, p = 18, rho = 0.50, tau = 1, CP = 0.6)
#print(result1)

#result2 <- Cox_AA(n = 200, p = 18, rho = 0.70 , tau = 1, CP = 0.6)
#print(result2)

#result3 <- Cox_AA(n = 200, p = 18, rho = 0.95, tau = 1, CP = 0.6)
#print(result3)
