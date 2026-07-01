# Improved Aslam–Ahmad Estimator in the Cox Model

## Introduction

This repository contains the implementation of the simulation studies and real-data analyses presented in the accompanying paper on the Improved Aslam–Ahmad estimator for the Cox model.

The proposed methodology extends the Improved Aslam–Ahmad shrinkage estimator to the Cox model by incorporating equality constraints derived from a preliminary variable selection procedure. The approach combines variable selection with constrained estimation to improve estimation accuracy, particularly in the presence of multicollinearity, sparse regression coefficients, and prior structural information.

The simulation study evaluates the performance of the proposed estimator under different sample sizes, correlation structures, sparsity levels, and censoring rates. Performance is assessed using Relative Average Absolute Error (RAAE), Average Absolute Loss (AALS), Average Absolute Shrinkage (AAS), and other evaluation criteria described in the accompanying paper.

This repository also provides a practical workflow for applying the proposed estimator to real survival datasets after verifying the assumptions of the Cox model.

---

# Methodology

The proposed methodology extends the Improved Aslam–Ahmad estimator to the Cox model by integrating shrinkage estimation with equality-constrained estimation.

Rather than relying solely on the unrestricted maximum partial likelihood estimator, the proposed approach incorporates structural information obtained from variable selection. Variables identified as inactive are treated as equality constraints, leading to a restricted estimator with reduced variance while maintaining satisfactory prediction performance.

The proposed estimation procedure consists of the following steps:

1. Fit the standard Cox model using all available covariates.

2. Verify the adequacy of the fitted Cox model through goodness-of-fit diagnostics.

3. Assess the functional form of continuous covariates and apply suitable transformations if necessary.

4. Perform variable selection using the LASSO Cox regression model.

5. Identify regression coefficients that are shrunk exactly to zero.

6. Construct the equality constraint matrix based on these inactive predictors.

7. Compute the restricted maximum partial likelihood estimator under the specified linear constraints.

8. Apply the Improved Aslam–Ahmad shrinkage estimator to obtain the final constrained estimator.

This methodology combines data-driven variable selection with constrained estimation, resulting in improved estimation efficiency, especially in the presence of multicollinearity and high-dimensional covariates.

---

# Limitations

Although the proposed estimator demonstrates favorable performance across a wide range of simulation settings, several limitations should be considered.

- The methodology assumes that the proportional hazards assumption of the Cox model is satisfied.

- Adequate model fit should be verified before applying the proposed estimator.

- Continuous covariates are assumed to have the appropriate functional form. Nonlinear effects should be modeled before estimation.

- The quality of the final estimator depends on the variable selection procedure. Incorrect selection by LASSO may produce inappropriate equality constraints.

- Equality constraints are treated as exact. If truly important variables are incorrectly constrained to zero, estimation bias may increase.

- The proposed estimator is designed for linear equality constraints and is not directly applicable to nonlinear constraint structures.

- Performance may deteriorate under extremely high censoring rates or very small sample sizes.

- The simulation study considers predefined correlation structures, sparsity levels, and sample sizes. Additional studies are required to evaluate performance under more complex data-generating mechanisms.

- The methodology assumes independent observations and does not directly address clustered survival data, recurrent events, or time-dependent covariates.

---

# Reproducing the Simulation Study

The simulation scripts are fully parameterized. To reproduce the results reported in the paper or conduct additional experiments, modify the following parameters.

## 1. Sample Size

```r
n
```

Typical values:

```text
50
100
200
```

---

## 2. Number of Covariates

```r
p
```

Example values:

```text
6
12
18
```

---

## 3. Sparsity Level

Specify the proportion of truly non-zero regression coefficients.

```r
cp
```

Examples:

```text
0.30
0.60
```

---

## 4. Correlation Between Covariates

```r
rho
```

Typical values:

```text
0.50
0.70
0.95
```

---

## 5. Censoring Level

Modify

```r
tau
```

Different values generate different censoring percentages used throughout the simulation study.

---

## 6. Number of Monte Carlo Replications

Modify

```r
N.sim
```

to increase or decrease the number of simulation replications.

---

## 7. Regression Coefficients

Specify the true coefficient vector

```r
beta
```

Users may define any sparse or dense coefficient configuration.

---

## 8. Random Seed

For reproducibility,

```r
set.seed(...)
```

may be changed to generate independent simulation runs.

---

# Applying the Method to Real Data

Before applying the Improved Aslam–Ahmad estimator, the assumptions of the Cox proportional hazards model should be carefully evaluated.

## Step 1. Fit the Standard Cox Model

Fit the ordinary Cox proportional hazards model using all candidate predictors.

---

## Step 2. Assess Overall Model Adequacy

Evaluate whether the Cox model adequately represents the observed survival data.

Recommended diagnostic methods include:

- Cox–Snell residual analysis
- Deviance residual plots
- Martingale residual diagnostics
- Other appropriate goodness-of-fit procedures

Proceed only if the Cox model demonstrates an acceptable fit.

---

## Step 3. Check the Functional Form of Continuous Covariates

Verify the linearity assumption for continuous predictors.

Common approaches include:

- Martingale residual plots
- Restricted cubic splines
- Fractional polynomial models
- Component-plus-residual diagnostics

If nonlinear relationships are detected, appropriate transformations or spline functions should be incorporated before variable selection.

---

## Step 4. Variable Selection

After confirming model adequacy, perform variable selection using the LASSO Cox regression model.

```text
LASSO Cox Regression
```

Predictors with non-zero estimated coefficients are retained.

Predictors whose coefficients shrink exactly to zero are regarded as inactive variables.

---

## Step 5. Construct the Equality Constraint Matrix

Use the LASSO results to construct the restriction matrix.

- Variables selected by LASSO remain unrestricted.
- Variables with zero coefficients are treated as equality constraints.

These constraints define the restricted parameter space for the proposed estimator.

---

## Step 6. Estimate the Improved Aslam–Ahmad Estimator

Fit the proposed estimator using

- the survival outcome,
- the selected predictors,
- and the equality constraint matrix.

The resulting estimator combines restricted maximum partial likelihood estimation with adaptive shrinkage.

---

## Step 7. Evaluate the Final Model

Report

- regression coefficient estimates,
- hazard ratios,
- confidence intervals,
- prediction accuracy,
- and all performance measures considered in the accompanying paper.

Comparisons may be performed against

- Standard Cox Regression
- Ridge Cox Regression
- LASSO Cox Regression
- Adaptive LASSO
- Elastic Net Cox Regression
- Other competing estimators.

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
Construct Equality Constraints
        │
        ▼
Restricted Maximum Partial Likelihood Estimation
        │
        ▼
Improved Aslam–Ahmad Estimation
        │
        ▼
Estimate Hazard Ratios and Evaluate Performance
```

This workflow ensures that the Improved Aslam–Ahmad estimator is applied only after validating the Cox proportional hazards model assumptions and constructing reliable equality constraints from the variable selection stage.
