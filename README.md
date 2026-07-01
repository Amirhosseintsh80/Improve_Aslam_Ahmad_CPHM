# Improved Aslam–Ahmad Estimator for the Cox Proportional Hazards Model

## Overview

This repository contains the R implementation of the **Improved Aslam–Ahmad estimator** for the Cox proportional hazards model, together with the simulation studies and real-data analysis presented in the accompanying paper.

The proposed estimator extends the classical Improved Aslam–Ahmad shrinkage methodology to survival analysis by incorporating **linear equality restrictions** into the Cox proportional hazards model. The estimator combines unrestricted and restricted maximum partial likelihood estimators through a data-driven shrinkage mechanism, leading to improved estimation accuracy when reliable prior information is available.

The simulation study evaluates the performance of the proposed estimator under various sample sizes, correlation structures, sparsity levels, and censoring rates. Estimation performance is assessed using Relative Average Absolute Error (RAAE), Average Absolute Loss (AALS), Average Absolute Shrinkage (AAS), and other performance measures introduced in the paper.

---

# Methodology

Consider the Cox proportional hazards model

\[
h(t|x)=h_0(t)\exp(x^\top\beta),
\]

where

- \(h_0(t)\) is the unspecified baseline hazard function,
- \(x\) is the vector of explanatory variables,
- \(\beta\) denotes the regression coefficient vector.

Suppose that prior information is available in the form of linear equality restrictions

\[
H\beta=h,
\]

where \(H\) is a known restriction matrix and \(h\) is a known constant vector.

The proposed estimation procedure consists of the following steps:

1. Estimate the regression coefficients using the **Maximum Partial Likelihood Estimator (MPLE)**.

2. Estimate the **Restricted Maximum Partial Likelihood Estimator (RMPLE)** under the specified linear equality restrictions.

3. Compute the shrinkage parameter according to the Improved Aslam–Ahmad methodology.

4. Combine the unrestricted and restricted estimators to obtain the **Improved Aslam–Ahmad estimator**.

The proposed estimator exploits available prior information while reducing estimation variance, resulting in improved finite-sample performance compared with conventional Cox regression estimators.

---

# Simulation Study

The simulation study can be reproduced by modifying the following parameters.

## Sample Size

```r
n
```

Typical values

```text
50
100
200
```

---

## Number of Covariates

```r
p
```

Typical values

```text
6
12
18
```

---

## Sparsity Level

```r
cp
```

Typical values

```text
0.30
0.60
```

---

## Correlation Between Covariates

```r
rho
```

Typical values

```text
0.50
0.70
0.95
```

---

## Censoring Level

```r
tau
```

Different values generate different censoring rates.

---

## Number of Monte Carlo Replications

```r
N.sim
```

---

## True Regression Coefficients

```r
beta
```

---

## Random Seed

```r
set.seed(...)
```

---

Simulation performance is evaluated using

- Relative Average Absolute Error (RAAE)
- Average Absolute Loss (AALS)
- Average Absolute Shrinkage (AAS)

and the additional performance measures presented in the paper.

---

# Real Data Analysis

The proposed estimator should only be applied after verifying that the Cox proportional hazards model is appropriate for the observed data.

## Step 1. Fit the Standard Cox Model

Fit the ordinary Cox proportional hazards model using all available predictors.

---

## Step 2. Assess Model Adequacy

Evaluate the overall goodness-of-fit using appropriate diagnostic procedures such as

- Cox–Snell residuals,
- Martingale residuals,
- Deviance residuals,
- or other suitable goodness-of-fit methods.

Proceed only if the Cox model adequately fits the observed data.

---

## Step 3. Check the Functional Form of Continuous Covariates

Assess the linearity assumption using methods such as

- Martingale residual plots,
- Restricted cubic splines,
- Fractional polynomial models.

Transform nonlinear covariates when necessary.

---

## Step 4. Variable Selection

Fit a **LASSO Cox regression model**.

Predictors with coefficients estimated as exactly zero are regarded as inactive variables.

---

## Step 5. Construct Linear Equality Restrictions

Use the inactive variables identified by LASSO to construct the restriction matrix.

Regression coefficients estimated as zero are imposed as linear equality constraints in the restricted Cox model.

---

## Step 6. Estimate the Improved Aslam–Ahmad Estimator

Using the constructed restriction matrix,

- compute the Restricted Maximum Partial Likelihood Estimator (RMPLE),
- estimate the shrinkage parameter,
- obtain the Improved Aslam–Ahmad estimator.

---

## Step 7. Model Evaluation

Report

- Regression coefficient estimates,
- Hazard ratios,
- Confidence intervals,
- Prediction performance,
- Model comparison with competing estimators.

---

# Limitations

The proposed estimator is subject to the following limitations.

- The proportional hazards assumption must be satisfied.
- The validity of the estimator depends on correctly specified linear equality restrictions.
- Incorrect restrictions may introduce estimation bias.
- The proposed methodology is designed only for linear equality constraints.
- Performance may deteriorate under extremely high censoring rates or very small sample sizes.
- The method assumes independent survival observations.
- Clustered survival data, recurrent events, and time-dependent covariates are beyond the scope of the current implementation.

---

# Recommended Workflow

```text
Raw Survival Data
        │
        ▼
Fit Standard Cox Model
        │
        ▼
Goodness-of-Fit Assessment
        │
        ▼
Check Functional Form of Continuous Covariates
        │
        ▼
LASSO Cox Variable Selection
        │
        ▼
Identify Zero Coefficients
        │
        ▼
Construct Linear Equality Restrictions
        │
        ▼
Restricted Maximum Partial Likelihood Estimation (RMPLE)
        │
        ▼
Improved Aslam–Ahmad Estimation
        │
        ▼
Estimate Hazard Ratios and Evaluate Model Performance
```


## Citation

If you use this code in your research, please cite the accompanying paper.
