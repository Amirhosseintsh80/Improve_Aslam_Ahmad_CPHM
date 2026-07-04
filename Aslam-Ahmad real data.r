# ============================================================
# Install and load required packages
# ============================================================
# Uncomment the following lines if packages are not installed:
# install.packages("penalized")
# install.packages("car")
# install.packages("glmnet")
# install.packages("survival")
# install.packages("corrplot")
# install.packages("MASS")

library(penalized)
library(car)
library(glmnet)
library(survival)
library(corrplot)
library(MASS) # Required for the ginv function

# ============================================================
# 1. Data Preparation
# ============================================================
# Load nki70 dataset
data(nki70)

# Create the specified dataframe
data1 <- data.frame(
  time = nki70$time,
  event = nki70$event,
  IGFBP5.1 = nki70$IGFBP5.1,
  DIAPH3.1 = nki70$DIAPH3.1,
  DIAPH3.2 = nki70$DIAPH3.2,
  IGFBP5 = nki70$IGFBP5,
  PECI.1 = nki70$PECI.1,
  PECI = nki70$PECI,
  MELK = nki70$MELK,
  RAB6B = nki70$RAB6B,
  ZNF533 = nki70$ZNF533,
  LGP2 = nki70$LGP2,
  MS4A7 = nki70$MS4A7,
  Age = nki70$Age,
  QSCN6L1 = nki70$QSCN6L1,
  DIAPH3 = nki70$DIAPH3,
  NUSAP1 = nki70$NUSAP1
)

# ============================================================
# 2. Helper Functions
# ============================================================
mse_per_variable <- function(Me, Tb) {
  Tb_reordered <- Tb[colnames(Me)]
  SM <- sweep(Me, 2, Tb_reordered, "-")^2
  colMeans(SM, na.rm = TRUE)
}

mean_per_variable <- function(Me) {
  colMeans(Me, na.rm = TRUE)
}

# Changed from Variance to Standard Deviation (SD)
sd_per_variable <- function(Me) {
  apply(Me, 2, sd, na.rm = TRUE)
}

# ============================================================
# 3. Main Function: Cox_Liu_real_full
# ============================================================
Cox_Liu_real_full <- function(data, null_vars, m = 80, sample_size = 80, seed = 123) {
  set.seed(seed)
  
  time_col <- "time"
  status_col <- "event"
  pred_names <- setdiff(names(data), c(time_col, status_col))
  p <- length(pred_names)
  
  # Step 1: Fit the full Cox model to obtain true coefficients
  full_fit <- coxph(as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ .")),
                    data = data[, c(time_col, status_col, pred_names)])
  tb <- coef(full_fit)
  tb_full <- rep(0, p)
  names(tb_full) <- pred_names
  tb_full[names(tb)] <- tb
  
  # Step 2: Set the null variables to zero
  true_betas_zeroed <- tb_full
  true_betas_zeroed[null_vars] <- 0
  
  # Step 3: Build the restriction matrix H (denoted as C in the paper)
  null_idx <- match(null_vars, pred_names)
  H <- matrix(0, nrow = length(null_idx), ncol = p)
  for (i in seq_along(null_idx)) H[i, null_idx[i]] <- 1
  h <- rep(0, nrow(H))
  
  # Step 4: Initialize storage matrices with Paper's Estimator Names
  estimators <- list(
    MPLE = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AAE  = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    RAAE = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AALS = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AAS  = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AAPS = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AAP  = matrix(NA, m, p, dimnames = list(NULL, pred_names)),
    AASP = matrix(NA, m, p, dimnames = list(NULL, pred_names))
  )
  
  # ============================================================
  # Monte Carlo Simulation
  # ============================================================
  for (u in 1:m) {
    idx <- sample(seq_len(nrow(data)), sample_size, replace = FALSE)
    dsub <- data[idx, , drop = FALSE]
    
    fit <- try(coxph(Surv(time, event) ~ ., data = dsub), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    b <- coef(fit)
    b_full <- rep(NA, p)
    names(b_full) <- pred_names
    b_full[names(b)] <- b
    if (any(is.na(b_full))) next
    estimators$MPLE[u, ] <- b_full
    
    vc <- try(vcov(fit), silent = TRUE)
    if (inherits(vc, "try-error")) next
    DTD <- try(solve(vc), silent = TRUE)
    if (inherits(DTD, "try-error")) DTD <- ginv(vc)
    
    lambda <- eigen(DTD)$values
    alpha_mp <- coef(fit)
    
    # Parameter selection algorithm for Aslam-Ahmad Estimator
    d_max <- max((lambda * (alpha_mp^2 - 1)) / (1 + lambda * alpha_mp^2))
    k_j <- ((lambda + d_max) - lambda * (1 - d_max) * alpha_mp^2) / ((1 + d_max) * (lambda + 1) * alpha_mp^2)
    k_med <- median(k_j)
    
    numerator <- sum(
      (lambda + k_med * (lambda + 1)) * (lambda^2 - k_med * lambda^2 - k_med * lambda) * alpha_mp^2 -
        lambda * (lambda^2 + k_med * lambda^2 - k_med * lambda)
      / sum((lambda + 1)))
    
    denominator <- sum(
      (lambda^2 + k_med * lambda^2 - k_med * lambda) +
        (lambda - k_med * (lambda + 1)) * (lambda^2 - k_med * lambda^2 - k_med * lambda) * alpha_mp^2
      / ((lambda + 1)))
    
    d_opt <- numerator / denominator
    if (d_opt < 0 || d_opt > 1) d_opt <- d_max
    
    term1 <- solve(DTD + diag(p))
    term2 <- DTD + d_opt * diag(p)
    term3 <- solve(DTD + k_med * (1 + d_opt) * diag(p))
    M <- term1 %*% term2 %*% term3 %*% DTD
    
    # 1. Unrestricted Aslam-Ahmad Estimator (AAE)
    AAE <- as.numeric(M %*% b_full)
    estimators$AAE[u, ] <- AAE
    
    # Restricted Partial Likelihood Estimator calculations
    HvcHt <- H %*% vc %*% t(H)
    Inv_HvcHt <- try(solve(HvcHt), silent = TRUE)
    if (inherits(Inv_HvcHt, "try-error")) Inv_HvcHt <- ginv(HvcHt)
    
    DTV <- H %*% b_full - h
    RMPLE <- b_full - vc %*% t(H) %*% Inv_HvcHt %*% DTV
    
    # 2. Restricted Aslam-Ahmad Estimator (RAAE)
    RAAE <- as.numeric(M %*% RMPLE)
    estimators$RAAE[u, ] <- RAAE
    
    # Differences and Test Statistic
    diff_AA_RAA <- AAE - RAAE
    Tn <- as.numeric(t(DTV) %*% Inv_HvcHt %*% DTV)
    df_h <- length(h) # Degrees of freedom (r)
    qq <- qchisq(0.95, df_h) # Alpha = 0.05
    TTest <- ifelse(Tn <= qq, 1, 0)
    
    # Shrinkage factor calculation
    shrinkage_factor <- (df_h - 2) / (Tn + 1e-8)
    
    # 3. Aslam Ahmed Linear Shrinkage (AALS) estimator (Eq 22, lambda = 0.5)
    estimators$AALS[u, ] <- AAE - 0.5 * diff_AA_RAA
    
    # 4. Aslam Ahmed Stein (AAS) estimator (Eq 23)
    estimators$AAS[u, ] <- AAE - shrinkage_factor * diff_AA_RAA
    
    # 5. Aslam-Ahmed Positive Stein (AAPS) estimator (Eq 24)
    positive_part <- max(0, 1 - shrinkage_factor)
    estimators$AAPS[u, ] <- RAAE + positive_part * diff_AA_RAA
    
    # 6. Aslam Ahmed Pretest (AAP) estimator (Eq 25)
    estimators$AAP[u, ] <- AAE - TTest * diff_AA_RAA
    
    # 7. Aslam Ahmed Shrinkage Pretest (AASP) estimator (Eq 26, lambda = 0.5)
    estimators$AASP[u, ] <- AAE - (0.5 * TTest) * diff_AA_RAA
  }
  
  # Calculate summary metrics for ALL variables (Active + Null)
  mse_results_full <- sapply(estimators, mse_per_variable, Tb = true_betas_zeroed)
  mean_results_full <- sapply(estimators, mean_per_variable)
  sd_results_full <- sapply(estimators, sd_per_variable)
  
  # -------------------------------------------------------------
  # IMPORTANT: Calculate Total MSE across ALL parameters 
  # This reveals the theoretical advantage of restricted estimators
  # -------------------------------------------------------------
  Total_MSE <- colSums(mse_results_full, na.rm = TRUE)
  
  # ============================================================
  # Filter out variables assumed to be zero for specific tables
  # ============================================================
  active_vars <- setdiff(pred_names, null_vars)
  
  mse_results_active <- mse_results_full[active_vars, ]
  mean_results_active <- mean_results_full[active_vars, ]
  sd_results_active <- sd_results_full[active_vars, ]
  
  list(
    Total_MSE = Total_MSE,
    MSE_Active = round(mse_results_active, 6),
    Mean_Active = round(mean_results_active, 6),
    SD_Active = round(sd_results_active, 6),
    True_Betas_Active = true_betas_zeroed[active_vars]
  )
}

# ============================================================
# 4. Execute Function and Display Results
# ============================================================
# Variables whose coefficients were zeroed by Lasso
null_vars <- c("IGFBP5.1", "DIAPH3.2", "PECI", "MELK", "DIAPH3")

res_full <- Cox_Liu_real_full(
  data = data1,
  null_vars = null_vars,
  m = 80,
  sample_size = 80,
  seed = 123
)

# Print Total MSE (Sorted to show which estimator is globally best)
cat("\n===== Total Mean Squared Error (Risk across ALL variables) =====\n")
print(sort(res_full$Total_MSE))

# Print active filtered results
cat("\n===== Mean Squared Error (MSE) - Active Variables Only =====\n")
print(res_full$MSE_Active)

cat("\n===== Mean of Estimated Coefficients - Active Variables Only =====\n")
print(res_full$Mean_Active)

cat("\n===== Standard Deviation (SD) of Estimated Coefficients - Active Variables Only =====\n")
print(res_full$SD_Active)
