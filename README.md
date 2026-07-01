# Improved Aslam-Ahmad estimator in Cox regression model

## Introduction

This repository contains the implementation of the simulation studies and real-data analyses presented in the accompanying paper on constrained adaptive shrinkage estimation for the Cox model.

The proposed methodology combines variable selection with constrained estimation to improve estimation accuracy, particularly when prior information is available. The simulation study evaluates the performance of the proposed estimator under different sample sizes, correlation structures, sparsity levels, and censoring rates. Performance is compared using the Relative Average Absolute Error (RAAE), Average Absolute Loss (AALS), Average Absolute Shrinkage (AAS), and other evaluation criteria described in the paper.

The repository also provides a framework for applying the proposed method to real survival datasets after verifying that the assumptions of the Cox model are satisfied.

---

# Reproducing the Simulation Study

The simulation scripts are fully parameterized. To reproduce the results reported in the paper or conduct additional experiments, the following parameters can be modified.

## 1. Sample Size

Adjust the number of observations

```r
n
```

Typical values used in the paper include

```text
50
100
200
```

---

## 2. Number of Covariates

Modify

```r
p
```

Example

```text
p = 6
p = 12
p = 18
```

---

## 3. Sparsity Level

Specify the proportion of truly non-zero regression coefficients

```r
cp
```

Examples

```text
cp = 0.30
cp = 0.60
```

---

## 4. Correlation Between Covariates

Change the correlation coefficient

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

## 5. Censoring Level

Modify the censoring parameter

```r
tau
```
Different values generate different censoring percentages used throughout the simulation study.
---

## 6. Number of Monte Carlo Replications

Increase or decrease

```r
N.sim
```

for more accurate Monte Carlo estimates.

---

## 7. Regression Coefficients

The true regression coefficients are specified through

```r
beta
```

Users may define any sparse or dense coefficient vector.

---

## 8. Random Seed

For reproducibility, change

```r
set.seed(...)
```

to generate independent simulation runs.

---

# Applying the Method to Real Data

Before fitting the proposed constrained Cox model, the standard assumptions of the Cox proportional hazards model should be verified.

## Step 1. Fit the Initial Cox Model

Fit the ordinary Cox proportional hazards model using all available predictors.

---

## Step 2. Assess Overall Model Adequacy

Before variable selection, evaluate whether the Cox model provides an adequate fit to the data.

Possible goodness-of-fit procedures include

- Cox–Snell residual analysis
- Deviance residual plots
- Martingale residual diagnostics
- Other appropriate goodness-of-fit tests

Only proceed if the Cox model adequately represents the data.

---

## Step 3. Check the Functional Form of Continuous Covariates

Verify that continuous predictors satisfy the required functional form.

Possible approaches include

- Martingale residual plots
- Restricted cubic splines
- Fractional polynomial models
- Component-plus-residual diagnostics

If strong nonlinearity is detected, appropriate transformations or spline terms should be incorporated before model selection.

---

## Step 4. Variable Selection

After confirming that the Cox model assumptions are reasonable, perform variable selection.

A recommended approach is to use the LASSO Cox regression model

```text
LASSO Cox Regression
```

The selected variables correspond to coefficients that remain non-zero after penalization.

Variables whose estimated coefficients shrink exactly to zero are considered inactive predictors.

---

## Step 5. Construct the Constraint Matrix

Based on the LASSO results,

- coefficients selected by LASSO remain unrestricted;
- coefficients shrunk to zero are treated as equality constraints.

These zero coefficients define the constraint matrix used in the proposed constrained estimation procedure.

---

## Step 6. Fit the Proposed Constrained Estimator

Run the constrained estimation algorithm using

- the survival outcome,
- the selected predictors,
- and the constraint matrix constructed from the LASSO results.

The resulting estimator incorporates both data-driven variable selection and prior equality constraints.

---

## Step 7. Evaluate the Final Model

Report

- estimated regression coefficients,
- hazard ratios,
- confidence intervals,
- prediction performance,
- and any additional performance measures used in the paper.

Model comparison may also be performed against

- the standard Cox model,
- Ridge Cox regression,
- LASSO Cox regression,
- Adaptive LASSO,
- or other competing estimators.

---

# Recommended Workflow

The complete workflow is summarized below.

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
Check Linearity of Continuous Covariates
        │
        ▼
LASSO Variable Selection
        │
        ▼
Identify Zero Coefficients
        │
        ▼
Construct Equality Constraints
        │
        ▼
Fit Proposed Constrained Estimator
        │
        ▼
Estimate Hazard Ratios and Evaluate Performance
```

This workflow ensures that the proposed constrained estimation procedure is applied only after verifying the validity of the Cox proportional hazards model assumptions and obtaining a reliable set of candidate variables.
