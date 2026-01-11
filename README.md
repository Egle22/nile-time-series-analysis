# nile-time-series-analysis
State space analysis of Nile River flow using Kalman filtering, smoothing, and maximum likelihood estimation in R
# Nile River Time Series Analysis

A comprehensive time series analysis of the Nile River annual flow dataset using State Space Models, Kalman Filtering, Smoothing, and Maximum Likelihood Estimation.

## Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Methodology](#methodology)
- [Key Results](#key-results)
- [Project Structure](#project-structure)
- [Setup Instructions](#setup-instructions)
- [Usage](#usage)
- [Theoretical Background](#theoretical-background)


## Overview

This project implements a Local Level Model (LLM) to analyze 100 years of annual flow measurements of the Nile River at Aswan (1871-1970). The analysis demonstrates:

- Kalman Filtering for real-time state estimation
- Kalman Smoothing for optimal retrospective estimation
- Parameter Estimation via Maximum Likelihood (BFGS optimization)
- Simulation of conditional and unconditional samples
- Model Diagnostics through disturbance analysis

The Local Level Model represents the time series as a time-varying mean (state) plus observation noise, where the state evolves as a random walk.

## Dataset

Source: Built-in R dataset `Nile`

Description: Annual flow volume of the Nile River at Aswan, Egypt (1871-1970)

Sample Size: 100 observations

Descriptive Statistics:
- Mean: 919.4
- Standard Deviation: 169.23
- Minimum: 456.0
- Maximum: 1370.0

Notable features include an unusually low flow around 1910 and an exceptionally high flow around 1880.

## Methodology

### 1. Local Level Model Specification

The Local Level Model (LLM) is formulated as:

```
yₜ = αₜ + εₜ,  εₜ ~ N(0, σ²_ε)
αₜ₊₁ = αₜ + ηₜ,  ηₜ ~ N(0, σ²_η)
```

Where:
- yₜ = observed flow at time t
- αₜ = unobserved level (latent state)
- εₜ = observation disturbance (measurement error)
- ηₜ = state disturbance (level innovation)

### 2. Kalman Filter

Implements the recursive estimation procedure:

```
vₜ = yₜ - aₜ                    (prediction error)
Fₜ = Pₜ + σ²_ε                   (prediction variance)
Kₜ = Pₜ / Fₜ                     (Kalman gain)
aₜ|ₜ = aₜ + Kₜvₜ                 (filtered state)
Pₜ|ₜ = Pₜ(1 - Kₜ)               (filtered variance)
aₜ₊₁ = aₜ + Kₜvₜ                 (predicted state)
Pₜ₊₁ = Pₜ(1 - Kₜ) + σ²_η         (predicted variance)
```

### 3. Kalman Smoothing

Computes optimal estimates using all available data:

```
αₜ|Yₙ ~ N(α̂ₜ, Vₜ)
```

Where α̂ₜ = E[αₜ|Yₙ] uses data up to time n (end of series).

### 4. Parameter Estimation

Uses Concentrated Diffuse Log-Likelihood with BFGS optimization:

```
log Ldc ∝ -(n-1)/2 · log(σ̂²_ε) - 1/2 · Σ log(F*ₜ)
```

Reparameterization: q = σ²_η / σ²_ε (signal-to-noise ratio)

Optimization: BFGS method with ψ = log(q) for unconstrained optimization

### 5. Simulation

Two types of simulations:

Unconditional Simulation:
```
y⁺ₜ = α⁺ₜ + ε⁺ₜ
α⁺ₜ₊₁ = α⁺ₜ + η⁺ₜ
```

Conditional Simulation (given observed data):
```
ε̃ₜ = ε⁺ₜ + ε̂⁺ₜ - ε̂ₜ
α̃ₜ = yₜ - ε̃ₜ
η̃ₜ = α̃ₜ₊₁ - α̃ₜ
```

## Key Results

### Estimated Parameters

| Parameter | Symbol | Estimated Value |
|-----------|--------|----------------|
| Observation variance | σ²_ε | 15,099.0 |
| State variance | σ²_η | 1,469.1 |
| Signal-to-noise ratio | q | 0.0973 |

### Convergence

- BFGS Optimization: 19 iterations
- Log-Likelihood: -492.07
- Optimal ψ: -2.33
- Optimal q: 0.097

### Model Diagnostics

- Filtered state variance converges to steady state (P̄ ≈ 4,031.9) after approximately 10 iterations
- Prediction variance stabilizes quickly, indicating filter convergence
- Smoothed estimates show reduced uncertainty compared to filtered estimates

## Project Structure

```
nile-time-series-analysis/
├── README.md                    
├── code/
│   └── nile_analysis.R         
├── report/
│   └── Home_work_adv_ts.pdf    
└── .gitignore
```

## Setup Instructions

### Prerequisites

1. R (version 4.0 or higher)
2. RStudio (recommended)

### Installation

Clone the repository:
```bash
git clone https://github.com/Egle22/nile-time-series-analysis.git
cd nile-time-series-analysis
```

Install required R packages:
```r
install.packages(c(
  "TSA",
  "dlm",
  "statespacer",
  "graphics",
  "stats",
  "kableExtra",
  "dplyr",
  "ggplot2",
  "hrbrthemes",
  "float",
  "tidyverse"
))
```

## Usage

Open RStudio and set working directory:
```r
setwd("/path/to/nile-time-series-analysis")
```

Run the script:
```r
source("code/nile_analysis.R")
```

### Script Sections

The analysis is organized in numbered sections:

1. Libraries - Load required packages
2. Data Loading - Load and summarize Nile dataset
3. Data Visualization - Plot time series
4. Kalman Filter - Fit Local Level Model
5. Kalman Filter Visualization - Plot filtered quantities
6. Kalman Smoothing - Compute smoothed estimates
7. Smoothing Visualization - Plot smoothed disturbances
8. Comparison - Filtered vs Smoothed states
9. Simulation - Unconditional and conditional sampling
10. Simulation Visualization - Compare simulated vs actual
11. Parameter Estimation - CDLL and BFGS optimization
12. Optimization Results - Display iteration table

## Theoretical Background

### State Space Representation

The Local Level Model is an Unobserved Component Model cast in State Space Form (SSF):

Observation Equation:
```
yₜ = Zₜαₜ + εₜ
```

State Equation:
```
αₜ₊₁ = Tₜαₜ + Rₜηₜ
```

For the LLM:
- Zₜ = 1 (scalar)
- Tₜ = 1 (random walk transition)
- Rₜ = 1 (state noise multiplier)

### Diffuse Initialization

The initial state α₁ is treated with a diffuse prior:
- a₁ fixed at arbitrary value (e.g., y₁)
- P₁ → ∞ (infinite initial uncertainty)

This approach converges to the same MLE as finite initialization.

### BFGS Algorithm

Broyden–Fletcher–Goldfarb–Shanno (BFGS) is a quasi-Newton optimization method that approximates the inverse Hessian matrix and updates search direction iteratively. It is robust for unconstrained nonlinear optimization with low memory requirements.

