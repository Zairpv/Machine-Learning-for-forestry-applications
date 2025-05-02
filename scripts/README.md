# 📂 Folder `scripts/`

This folder contains R scripts developed to apply machine learning tools to ecological and forestry data. Each script is numbered and documented for easy review and execution.

---

# 📝 Script list

## `1_Exploratory analyses for RF.R`
This script performs an exploratory analysis prior to implementing the **Random Forest** model. It includes:

- Reading and preparing input data.
- Analysis of data distribution
- Correlation analysis
- Split data: training & testing datasets
- VIF determination
- data visualization

---

## `2_Recursive Feature Selection.R`
This script performs the first step of any machine learning process: the selection of best predictors.
It includes: 

- Reading training and testing data (Only 2013 data)
- Performs the Boruta algorithm 
- Filter "confirmed" predictors for subsequent random forest modeling
- Data visualization

---


## `3_RandomForest_Example.R`
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
  