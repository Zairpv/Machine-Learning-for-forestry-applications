# Machine Learning for Forestry Applications 🌲
## Estimations of Forest floor C stocks (2013 - 2023)

This repository contains a portfolio of scripts that apply machine learning (ML) techniques to model forest floor. 
It aims to support reproducible modeling workflows for analyzing environmental issues.

---

### 📄 Technical Report

The full technical report detailing the data processing and interpretation of results is available at the following link:

👉 [View the Technical Report](https://zairpv.github.io/Machine-Learning-for-forestry-applications/)

This report summarizes the analyses and methodologies implemented in the scripts and datasets contained within this repository, following a reproducible R Markdown workflow.


---

### 🔍 Objectives

- Explore and document the use of ML algorithms in forest science.
- Share reproducible workflows for carbon and biodiversity modeling.

---

### 📝 List of Scripts and details 

**1_Exploratory analyses for RF.R:**

- Estimate the descriptive statistics for target and predictor variables.  
- Analyze data distribution and normality.  
- Explore correlation between target and predictor variables.  
- Explore potential multicollinearity.  
- Split data into train and test subsets.  

---

**2_Recursive Feature Selection.R**
This script performs the first step of any machine learning process: the selection of best predictors.
It includes: 

- Reading training and testing data (Only 2013 data)
- Performs the Boruta algorithm 
- Filter "confirmed" predictors for subsequent random forest modeling
- Data visualization

---

**3_RandomForest_Example.R**
This script performs the second step: RF modeling. It includes:

- Fit RF modeling for C stocks data (2013)
- Tuning hyperparameters 
  - The selection of the best combination were evaluated as iterations increases in order to detect error stability
  - We tested 5, 10, 15, 20, 25, 50, 75 and 100 cross validation repetitions
  - Metrics of error stability were assesed by Non parametric friedman test
  - This script is an optimization framework for RF modeling 
- Assesing model uncertainty
  - Definition of observed vs predicted values
- Model validation
  - Cross validation metrics 
  - Error stabilization assessment
- Spatial predictions
  - Map creation
  - SD maps 
- Variable importance 
  - Comparison of Variable impiortance
  - Partial plots and ecological analysis
  - 3D Plotly graphs - Non linear relationships among dependant and predictor variable


**4. Final graphs - all models.R**
It contains the visualization of results:
  - ggplots of error stabilization
  - variable importance comparisons
  - spatial predictions assessments
  
---

### 🚪 Getting Started

- Clone or download the repository.
- Open the .Rproj file in RStudio.
- Review scripts under the scripts/ folder.
- Use the simulated dataset in data/ to reproduce the workflow.

---
