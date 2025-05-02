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
