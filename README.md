# Household Expenditure Analysis Across UK Regions 💷

## Description
This project investigates how demographic and economic factors influence household expenditure patterns across various regions of the UK. By employing **Bayesian and Hierarchical Bayesian models**, this analysis compares their performance against frequentist linear regression to understand the most effective approach for capturing regional variations in spending behavior.

The study provides valuable insights for policymakers, helping them design targeted fiscal and social policies to address economic challenges such as income inequality, inflation, and regional disparities.

---

## Features
- **Demographic and Economic Analysis**:
  - Factors considered: household income, household tenure, number of children, and economic position of the household reference person.
- **Comparison of Modeling Techniques**:
  - Frequentist Linear Regression.
  - Bayesian Regression (Nimble and STAN).
  - Hierarchical Bayesian Models (Nimble and STAN) with regional hierarchy.
- **Key Insights**:
  - Identified which variables have consistent impacts on spending across UK regions.
  - Evaluated the predictive accuracy and complexity of different models.
- **Policy-Oriented Results**:
  - Provides actionable insights for tailoring economic and fiscal policies.

---

## My Contributions
1. **Model Development**:
   - Built and implemented linear regression, Bayesian regression, and Hierarchical Bayesian models using Nimble and STAN.
2. **Variable Integration**:
   - Incorporated demographic and economic variables into models for robust analysis.
3. **Methodological Comparison**:
   - Compared frequentist and Bayesian approaches to validate findings and ensure reliability.
4. **Policy Recommendations**:
   - Highlighted the implications of expenditure patterns for government policy.

---

## Tech Stack
- **Programming Language**: R
- **Key Libraries**:
  - `Nimble` for Bayesian modeling
  - `rstan` for STAN-based modeling
  - `ggplot2` for visualizations
  - `dplyr` for data manipulation
- **Tools**: RStudio, Markdown

---

## Key Findings
1. **Bayesian Hierarchical Models**:
   - Achieved the best balance of predictive accuracy and model complexity.
   - Incorporating a random intercept for regional hierarchy provided robust insights into regional expenditure variations.
2. **Demographic Impact**:
   - Household income and the economic position of the reference person had significant effects on spending.
   - The number of children in the household showed varying impacts across regions.
3. **Model Comparison**:
   - Bayesian approaches outperformed linear regression in capturing regional nuances.
   - STAN and Nimble implementations yielded comparable results with slight differences in computational efficiency.

---

