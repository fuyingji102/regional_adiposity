
# Regional Adiposity Shapes Brain and Cognition in Adults 

This repository contains the code used for the analyses in **MATLAB (R2017b; The MathWorks, Inc., USA), Python (version 3.6), and R (version 4.3.1)**.    

## Overview

The study investigated the associations between regional adiposity measures and brain morphology, function connectivity, white matter microstructure and cognition performance, using:  
1. **Linear regression models**
2. **Mediation analysis**

## Repository Structure

<pre>
├── scripts/
│   ├── 01_linear_regression.m  
│   ├── 02_mediation_analysis.R  
├── data/
│   └── example_data.csv
├── figures/
│   ├── regression_forest_plot.png  
│   └── mediation_path_diagram.png
└── README.md    
</pre>

## Analysis Pipeline

### 1. Data Preparation
- Import and clean raw dataset.  
- Calculate adiposity measures: trunk fat, leg fat, arm fat, visceral adipose tissue (VAT), BMI.
- Merge demographic covariates (age, sex, ethnicity, smoking, alcohol use, etc.).  

### 2. Linear Regression
- Models run separately for each adiposity measure.  
- Adjusted for covariates.
- Outputs: β estimates, 95% confidence intervals, p-values.  

Example:
```r
model <- lm(outcome ~ adiposity_measure + covariates, data = df)  
summary(model)


