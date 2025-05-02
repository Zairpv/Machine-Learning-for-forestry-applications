################################################################################
# Script: 3_RandomForest_Example.R
# Author: Dra. Zaira Rosario Pérez-Vázquez
# Contact: zairpv@gmail.com
#
# Description:
# This script performs the steps for selecting appropiate predictors by employing
# several machine learning tools. For example, the algorithm Boruta which is implementing 
# for Random Forest (RF) model.
#
# Note:
# This is part of a larger portfolio on Machine Learning for Forestry Applications.
################################################################################

# 0. LIBRARIES -----------------------------------------------------------------
library(ggpubr)
library(tools)
library(sf)
library(randomForest)
library(randomForestSRC)
library(ggRandomForests)
library(ggplot2)
library(RColorBrewer)
library(plot3D)
library(dplyr)
library(parallel)
library(raster)
library(ModelMap)
library(rgdal)
library(sp)
library(extrafont)
library(caret)
library(e1071)
library(rfUtilities)
library(plotmo)
library(data.table)
library(ggpubr)
library(here)
library(tidyverse)
library(doParallel)
library(foreach)
library(CAST)
library(cowplot)
library(grid)
library(gridExtra)
library(forcats)
library(magick)
library(biscale)
library(gt)
library(pdp)
library(DescTools)
library(psych)
library(purrr)
library(PerformanceAnalytics)
library(coin)
library(emmeans)
library(stats)
library(agricolae)
library(emmeans)
library(coin)
library(Kendall)
library(akima)
library(vivid)
library(plotly)

# 1.DATA AND CONSOLE PROPERTIES ------------------------------------------------
options(max.print=1000000)
windowsFonts(A = windowsFont("Times New Roman"))
getwd()


# Train and test dataset 2013
data_train2013=read_csv(here("data","training_2013_85.csv"))
data_test2013=read_csv(here("data","testing_2013_15.csv"))

# Define as factors
View(data_train2013)
data_train2013$YEAR=as.factor(data_train2013$YEAR)
data_train2013$LAYER=as.factor(data_train2013$LAYER)
data_train2013$CONDITION=as.factor(data_train2013$CONDITION)
data_test2013$YEAR=as.factor(data_test2013$YEAR)
data_test2013$LAYER=as.factor(data_test2013$LAYER)
data_test2013$CONDITION=as.factor(data_test2013$CONDITION)

#Select covariates
names(data_train2013)
data2013_formodels=data.frame(data_train2013[c(1,4,5,6:18)])
names(data2013_formodels)=c("YEAR","FFLAYER","C",
                            "UTMX","UTMY",
                            "AGE","BA","DTH","RICH",
                            "H","GAP","COV","ELE","SLO","ASP","MAN")

# 2. RANDOM FOREST TRAINING ----------------------------------------------------
## 2.1. Expressions for plots --------------------------------------------------
Mgha_expression=expression("C stocks (Mg ha"^"-1"~")")
Mgha_expression2=expression(bold("a) C models (Mg ha"^"-1"~")"))
Mgha_expression3=expression("Mg ha"^"-1")
initialmodel_expresion=expression(paste("RF"[0]," (Default values)"))
finalmodel_expression=expression(paste("RF"[1], " (Optimized values)"))

## 2.2. Tuning hyperparameters -------------------------------------------------

#Create the RF function to find the best hyperparameters
customRF=list(type="Regression",
              library="randomForest",
              loop=NULL)
customRF$parameters=data.frame(parameter=c("ntree","mtry"),
                               class=rep("numeric",2),
                               label=c("ntree","mtry"))
customRF$grid=function(x,y,len=NULL,search="grid"){}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}

customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes


# Outline the grid of parameters
set.seed(12345)

tunegrid_all <- expand.grid(.ntree=c(500,600,700,800,900,1000),.mtry=c(1:12))

# Set grid search parameters of ten-fold cross validation with repetition
#Remember:The goal of RKCV is to reduce the variance in your performance estimates.
# More repetitions generally lead to lower variance but at a higher computational cost.

controltenfold_5rep <- trainControl(method="repeatedcv", number=10, repeats=5, search='grid')

# Outline the grid of parameters
set.seed(12345)
names(data_train2013)
all_Predictors2013_train <- data_train2013[c(4,8:18)]
all_Predictors2013_test <- data_test2013[c(4,8:18)]



## 2.3. Random forest application ----------------------------------------------
best_hyperpars_C2013_5rep<- train(x=all_Predictors2013_train, 
                                 y=data_train2013$C_STOCKS, 
                                 method=customRF, 
                                metric="RMSE", 
                                 tuneGrid=tunegrid_all, 
                                 trControl=controltenfold_5rep)
best_hyperpars_C2013_5rep
best_hyperpars_C2013_5rep$finalModel
View(best_hyperpars_C2013_5rep$results)
write.csv(best_hyperpars_C2013_5rep$results,here("outputs", "best_hyperpars_C2013_CV5rep.csv"),row.names = TRUE)

## 2.4. Observed vs predicted --------------------------------------------------
pred_C2013_train_5rep=predict(best_hyperpars_C2013_5rep$finalModel,data_train2013)
obs_C2013_train_5rep=data_train2013$C_STOCKS
Obs_Pre_C2013_train_5rep=data.frame(obs=obs_C2013_train_5rep,
                                    pred=pred_C2013_train_5rep)


#Metric errors CV - independent
sqrt(mean((obs_C2013_train_5rep - pred_C2013_train_5rep)^2))  #RMSE
1 - sum((obs_C2013_train_5rep - pred_C2013_train_5rep)^2) / sum((train_data$C_STOCKS - mean(train_data$C_STOCKS))^2) #R2
mean(abs(obs_C2013_train_5rep - pred_C2013_train_5rep)) #MAE
#How correlated are obs vs pred from cross validation method
cor(Obs_Pre_C2013_train_5rep$obs,Obs_Pre_C2013_train_5rep$pred,method = "pearson")
cor(Obs_Pre_C2013_train_5rep$obs,Obs_Pre_C2013_train_5rep$pred,method = "spearman")


# Test data:
pred_C2013_test_5rep=predict(best_hyperpars_C2013_5rep$finalModel,data_test2013)
obs_C2013_test_5rep=data_test2013$C_STOCKS
Obs_Pre_C2013_test_5rep=data.frame(obs=obs_C2013_test_5rep,
                                   pred=pred_C2013_test_5rep)

#Metric errors CV 
sqrt(mean((obs_C2013_test_5rep - pred_C2013_test_5rep)^2))  #RMSE
1 - sum((obs_C2013_test_5rep - pred_C2013_test_5rep)^2) / sum((test_data$C_STOCKS - mean(test_data$C_STOCKS))^2) #R2
mean(abs(obs_C2013_test_5rep - pred_C2013_test_5rep)) #MAE
#How correlated are obs vs pred from cross validation method
cor(Obs_Pre_C2013_test_5rep$obs,Obs_Pre_C2013_test_5rep$pred,method = "pearson")
cor(Obs_Pre_C2013_test_5rep$obs,Obs_Pre_C2013_test_5rep$pred,method = "spearman")

## 3.2. With 10 repetitions  ---------------------------------------------------
### 3.2.1 Model training -------------------------------------------------------
#best_hyperpars_C2013_10rep<- train(x=all_Predictors2013_train, 
#                                 y=data_train2013$C_STOCKS, 
#                                 method=customRF, 
#                                 metric="RMSE", 
#                                 tuneGrid=tunegrid_all, 
#                                 trControl=controltenfold_10rep)
best_hyperpars_C2013_10rep
best_hyperpars_C2013_10rep$finalModel
View(best_hyperpars_C2013_10rep$results)
write.csv(best_hyperpars_C2013_10rep$results,here("Final data", "best_hyperpars_C2013_CV10rep.csv"),row.names = TRUE)









## 3.3. With 15 repetitions  ---------------------------------------------------
### 3.3.1 Model training -------------------------------------------------------
#best_hyperpars_C2013_15rep<- train(x=all_Predictors2013_train, 
#                                 y=data_train2013$C_STOCKS, 
#                                 method=customRF, 
#                                 metric="RMSE", 
#                                 tuneGrid=tunegrid_all, 
#                                 trControl=controltenfold_15rep)
best_hyperpars_C2013_15rep
best_hyperpars_C2013_15rep$finalModel
View(best_hyperpars_C2013_15rep$results)
write.csv(best_hyperpars_C2013_15rep$results,here("Final data", "best_hyperpars_C2013_CV15rep.csv"),row.names = TRUE)


## 3.4. With 20 repetitions  ---------------------------------------------------
### 3.4.1 Model training -------------------------------------------------------
#best_hyperpars_C2013_20rep<- train(x=all_Predictors2013_train, 
#                                 y=data_train2013$C_STOCKS, 
#                                 method=customRF, 
#                                 metric="RMSE", 
#                                 tuneGrid=tunegrid_all, 
#                                 trControl=controltenfold_20rep)
best_hyperpars_C2013_20rep
best_hyperpars_C2013_20rep$finalModel
View(best_hyperpars_C2013_20rep$results)
write.csv(best_hyperpars_C2013_20rep$results,here("Final data","best_hyperpars_C2013_CV20rep.csv"),row.names = TRUE)





## 3.5. With 25 repetitions  ---------------------------------------------------
### 3.5.1 Model training -------------------------------------------------------
#best_hyperpars_C2013_25rep<- train(x=all_Predictors2013_train, 
#                                 y=data_train2013$C_STOCKS, 
#                                 method=customRF, 
#                                 metric="RMSE", 
#                                 tuneGrid=tunegrid_all, 
#                                 trControl=controltenfold_25rep)
best_hyperpars_C2013_25rep
best_hyperpars_C2013_25rep$finalModel
View(best_hyperpars_C2013_25rep$results)
write.csv(best_hyperpars_C2013_25rep$results,here("Final data", "best_hyperpars_C2013_CV25rep.csv"),row.names = TRUE)


#-----------------------------------------------------------------------------------
# 4. RUN DEFAULT MODELS WITH ALL DATA - !REPORTED!  --------------------------------------------
tunegrid_default= expand.grid(.ntree=c(500),.mtry=c(4))
alldata=rbind(data_train2013,data_test2013)
View(alldata)
## 4.1. With 5 repetitions  ---------------------------------------------------
### 4.1.1 Model training -------------------------------------------------------
default_alldata_C2013_5rep<- train(x=alldata[c(5,12:22)], 
                                   y=alldata$C_STOCKS, 
                                   method=customRF, 
                                   metric="RMSE", 
                                   tuneGrid=tunegrid_default, 
                                   trControl=controltenfold_5rep)
default_alldata_C2013_5rep
View(default_alldata_C2013_5rep$results)
default_alldata_C2013_5rep$finalModel


## 4.2. With 10 repetitions  ---------------------------------------------------
### 4.2.1 Model training -------------------------------------------------------
default_alldata_C2013_10rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_10rep)
default_alldata_C2013_10rep
View(default_alldata_C2013_10rep$results)
default_alldata_C2013_10rep$finalModel


## 4.3. With 15 repetitions  ---------------------------------------------------
### 4.3.1 Model training -------------------------------------------------------
default_alldata_C2013_15rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_15rep)
default_alldata_C2013_15rep
View(default_alldata_C2013_15rep$results)
default_alldata_C2013_15rep$finalModel

## 4.4. With 20 repetitions  ---------------------------------------------------
### 4.4.1 Model training -------------------------------------------------------
default_alldata_C2013_20rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_20rep)
default_alldata_C2013_20rep
View(default_alldata_C2013_20rep$results)
default_alldata_C2013_20rep$finalModel

## 4.5. With 25 repetitions  ---------------------------------------------------
### 4.5.1 Model training -------------------------------------------------------
default_alldata_C2013_25rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_25rep)
default_alldata_C2013_25rep
View(default_alldata_C2013_25rep$results)
default_alldata_C2013_25rep$finalModel

## 4.6. With 50 repetitions - SELECTED RF0 ---------------------------------------------------
### 4.6.1 Model training -------------------------------------------------------
default_alldata_C2013_50rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_50rep)
default_alldata_C2013_50rep
View(default_alldata_C2013_50rep$results)
default_alldata_C2013_50rep$finalModel
### 4.6.2. Observed vs predicted model  ------------------------------------
# Test data:
pred_C2013_DEFAULT_test=predict(default_alldata_C2013_50rep$finalModel,data_test2013)
obs_C2013_DEFAULT_test=data_test2013$C_STOCKS
obspred_C2013_DEFAULT_test=data.frame(obs=obs_C2013_DEFAULT_test,
                                      pred=pred_C2013_DEFAULT_test)
#Metric errors CV - independent
sqrt(mean((obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)^2))  #RMSE
1 - sum((obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)^2) / sum((obs_C2013_DEFAULT_test - mean(obs_C2013_DEFAULT_test))^2) #R2
mean(abs(obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_DEFAULT_test$obs,obspred_C2013_DEFAULT_test$pred,method = "spearman")
cor(obspred_C2013_DEFAULT_test$obs,obspred_C2013_DEFAULT_test$pred,method = "pearson")


# train data:
pred_C2013_DEFAULT_train=predict(default_alldata_C2013_50rep$finalModel,data_train2013)
obs_C2013_DEFAULT_train=data_train2013$C_STOCKS
obspred_C2013_DEFAULT_train=data.frame(obs=obs_C2013_DEFAULT_train,
                                       pred=pred_C2013_DEFAULT_train)
#Metric errors CV - independent
sqrt(mean((obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)^2))  #RMSE
1 - sum((obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)^2) / sum((obs_C2013_DEFAULT_train - mean(obs_C2013_DEFAULT_train))^2) #R2
mean(abs(obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_DEFAULT_train$obs,obspred_C2013_DEFAULT_train$pred,method = "spearman")
cor(obspred_C2013_DEFAULT_train$obs,obspred_C2013_DEFAULT_train$pred,method = "pearson")

### 4.6.3. Dataframe ObsPred  ------------------------------------

ObsPred_DEFAULT_BOTH=rbind(data.frame(variable="C stocks",year=2013,model="RF0",data="1. Train",obspred_C2013_DEFAULT_train),
                           data.frame(variable="C stocks",year=2013,model="RF0",data="2. Test",obspred_C2013_DEFAULT_test))
### 4.6.4. ntree RMSE data ----------------------------------
error_ntree_2013C_default=as.data.table(plot(default_alldata_C2013_50rep$finalModel))
error_ntree_2013C_default[, trees := .I]
error_ntree_2013C_default2 = melt(error_ntree_2013C_default, id.vars = "trees")
setnames(error_ntree_2013C_default2,"trees","error (MSE)")
error_ntree_2013C_default2=as.data.frame(error_ntree_2013C_default2)
names(error_ntree_2013C_default2)=c("trees","variable","error (MSE)")
max(sqrt(error_ntree_2013C_default2$`error (MSE)`)) #ylims
min(sqrt(error_ntree_2013C_default2$`error (MSE)`)) #ylims
ntreePLOT_C13_default_df=data.frame(error_ntree_2013C_default2,sqrt(error_ntree_2013C_default2$`error (MSE)`))
names(ntreePLOT_C13_default_df)=c("trees","variable","MSE","RMSE")
names(ntreePLOT_C13_default_df)

## 4.7. With 75 repetitions  ---------------------------------------------------
### 4.7.1 Model training -------------------------------------------------------
default_alldata_C2013_75rep<- train(x=alldata[c(5,12:22)], 
                                    y=alldata$C_STOCKS, 
                                    method=customRF, 
                                    metric="RMSE", 
                                    tuneGrid=tunegrid_default, 
                                    trControl=controltenfold_75rep)
default_alldata_C2013_75rep
View(default_alldata_C2013_75rep$results)
default_alldata_C2013_75rep$finalModel

## 4.8. With 100 repetitions  ---------------------------------------------------
### 4.8.1 Model training -------------------------------------------------------
default_alldata_C2013_100rep<- train(x=alldata[c(5,12:22)], 
                                     y=alldata$C_STOCKS, 
                                     method=customRF, 
                                     metric="RMSE", 
                                     tuneGrid=tunegrid_default, 
                                     trControl=controltenfold_100rep)
default_alldata_C2013_100rep
View(default_alldata_C2013_100rep$results)
default_alldata_C2013_100rep$finalModel
## 4.9. JOIN ALL RUNS  -------------------------------------------------------
default_models_CVREPS=rbind(data.frame(variable="C stocks",year=2013,model="RF0",CVreps=5,default_alldata_C2013_5rep$results,Modelselected="Descarted"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=10,default_alldata_C2013_10rep$results,Modelselected="Descarted"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=15,default_alldata_C2013_15rep$results,Modelselected="Descarted"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=20,default_alldata_C2013_20rep$results,Modelselected="Descarted"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=25,default_alldata_C2013_25rep$results,Modelselected="Descarted"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=50,default_alldata_C2013_50rep$results,Modelselected="Reported RF0"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=75,default_alldata_C2013_75rep$results,Modelselected="Reported RF0"),
                            data.frame(variable="C stocks",year=2013,model="RF0",CVreps=100,default_alldata_C2013_100rep$results,Modelselected="Descarted"))
View(default_models_CVREPS)
# ------------------------------------------------------------------------------
# 5. RUN DEFAULT MODELS WITH TRAIN DATA  --------------------------------------------
## 5.1. With 5 repetitions  ---------------------------------------------------
### 5.1.1 Model training -------------------------------------------------------
names(data_train2013)
default_C2013_5rep<- train(x=data_train2013[c(5,12:22)], 
                           y=data_train2013$C_STOCKS, 
                           method=customRF, 
                           metric="RMSE", 
                           tuneGrid=tunegrid_default, 
                           trControl=controltenfold_5rep)
default_C2013_5rep
View(default_C2013_5rep$results)
default_C2013_5rep$finalModel


## 5.2. With 10 repetitions  ---------------------------------------------------
### 5.2.1 Model training -------------------------------------------------------
names(data_train2013)
default_C2013_10rep<- train(x=data_train2013[c(5,12:22)], 
                            y=data_train2013$C_STOCKS, 
                            method=customRF, 
                            metric="RMSE", 
                            tuneGrid=tunegrid_default, 
                            trControl=controltenfold_10rep)
default_C2013_10rep
View(default_C2013_10rep$results)
default_C2013_10rep$finalModel


## 5.3. With 15 repetitions  ---------------------------------------------------
### 5.3.1 Model training -------------------------------------------------------
names(data_train2013)
default_C2013_15rep<- train(x=data_train2013[c(5,12:22)], 
                            y=data_train2013$C_STOCKS, 
                            method=customRF, 
                            metric="RMSE", 
                            tuneGrid=tunegrid_default, 
                            trControl=controltenfold_15rep)
default_C2013_15rep
View(default_C2013_15rep$results)
default_C2013_15rep$finalModel

## 5.4. With 20 repetitions  ---------------------------------------------------
### 5.4.1 Model training -------------------------------------------------------
names(data_train2013)
default_C2013_20rep<- train(x=data_train2013[c(5,12:22)], 
                            y=data_train2013$C_STOCKS, 
                            method=customRF, 
                            metric="RMSE", 
                            tuneGrid=tunegrid_default, 
                            trControl=controltenfold_20rep)
default_C2013_20rep
View(default_C2013_20rep$results)
default_C2013_20rep$finalModel

### 5.4.2. Observed vs predicted model  ------------------------------------
# Test data:
pred_C2013_DEFAULT_test=predict(default_C2013_20rep$finalModel,data_test2013)
obs_C2013_DEFAULT_test=data_test2013$C_STOCKS
obspred_C2013_DEFAULT_test=data.frame(obs=obs_C2013_DEFAULT_test,
                                      pred=pred_C2013_DEFAULT_test)
#Metric errors CV - independent
sqrt(mean((obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)^2))  #RMSE
1 - sum((obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)^2) / sum((obs_C2013_DEFAULT_test - mean(obs_C2013_DEFAULT_test))^2) #R2
mean(abs(obs_C2013_DEFAULT_test - pred_C2013_DEFAULT_test)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_DEFAULT_test$obs,obspred_C2013_DEFAULT_test$pred,method = "spearman")
cor(obspred_C2013_DEFAULT_test$obs,obspred_C2013_DEFAULT_test$pred,method = "pearson")


# train data:
pred_C2013_DEFAULT_train=predict(default_C2013_20rep$finalModel,data_train2013)
obs_C2013_DEFAULT_train=data_train2013$C_STOCKS
obspred_C2013_DEFAULT_train=data.frame(obs=obs_C2013_DEFAULT_train,
                                       pred=pred_C2013_DEFAULT_train)
#Metric errors CV - independent
sqrt(mean((obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)^2))  #RMSE
1 - sum((obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)^2) / sum((obs_C2013_DEFAULT_train - mean(obs_C2013_DEFAULT_train))^2) #R2
mean(abs(obs_C2013_DEFAULT_train - pred_C2013_DEFAULT_train)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_DEFAULT_train$obs,obspred_C2013_DEFAULT_train$pred,method = "spearman")
cor(obspred_C2013_DEFAULT_train$obs,obspred_C2013_DEFAULT_train$pred,method = "pearson")

### 5.4.3. Dataframe ObsPred  ------------------------------------

ObsPred_DEFAULT_BOTH=rbind(data.frame(variable="C stocks",year=2013,model="RF0",data="1. Train",obspred_C2013_DEFAULT_train),
                           data.frame(variable="C stocks",year=2013,model="RF0",data="2. Test",obspred_C2013_DEFAULT_test))
## 5.5. With 25 repetitions  ---------------------------------------------------
### 5.5.1 Model training -------------------------------------------------------
default_C2013_25rep<- train(x=data_train2013[c(5,12:22)], 
                            y=data_train2013$C_STOCKS, 
                            method=customRF, 
                            metric="RMSE", 
                            tuneGrid=tunegrid_default, 
                            trControl=controltenfold_25rep)
default_C2013_25rep
View(default_C2013_25rep$results)
default_C2013_25rep$finalModel

## 5.6. With 50 repetitions  ---------------------------------------------------
### 5.6.1 Model training -------------------------------------------------------
default_C2013_50rep<- train(x=data_train2013[c(5,12:22)], 
                            y=data_train2013$C_STOCKS, 
                            method=customRF, 
                            metric="RMSE", 
                            tuneGrid=tunegrid_default, 
                            trControl=controltenfold_50rep)
default_C2013_50rep
View(default_C2013_50rep$results)
default_C2013_50rep$finalModel
## 5.7. With 100 repetitions  ---------------------------------------------------
### 5.7.1 Model training -------------------------------------------------------
default_C2013_100rep<- train(x=data_train2013[c(5,12:22)], 
                             y=data_train2013$C_STOCKS, 
                             method=customRF, 
                             metric="RMSE", 
                             tuneGrid=tunegrid_default, 
                             trControl=controltenfold_100rep)
default_C2013_100rep
View(default_C2013_100rep$results)
default_C2013_100rep$finalModel
# ----------------------------------------------------------------------------------
# 6. RUN final MODELS ALL DATA --------------------------------------------
tunegrid_final= expand.grid(.ntree=c(600),.mtry=c(12))
View(alldata)
## 6.1. With 5 repetitions  ---------------------------------------------------
### 6.1.1 Model training -------------------------------------------------------
final_alldata_C2013_5rep<- train(x=alldata[c(5,12:22)], 
                                 y=alldata$C_STOCKS, 
                                 method=customRF, 
                                 metric="RMSE", 
                                 tuneGrid=tunegrid_final, 
                                 trControl=controltenfold_5rep)
final_alldata_C2013_5rep
View(final_alldata_C2013_5rep$results)
final_alldata_C2013_5rep$finalModel


## 6.2. With 10 repetitions  ---------------------------------------------------
### 6.2.1 Model training -------------------------------------------------------
final_alldata_C2013_10rep<- train(x=alldata[c(5,12:22)], 
                                  y=alldata$C_STOCKS, 
                                  method=customRF, 
                                  metric="RMSE", 
                                  tuneGrid=tunegrid_final, 
                                  trControl=controltenfold_10rep)
final_alldata_C2013_10rep
View(final_alldata_C2013_10rep$results)
final_alldata_C2013_10rep$finalModel


## 6.3. With 15 repetitions  ---------------------------------------------------
### 6.3.1 Model training -------------------------------------------------------
final_alldata_C2013_15rep<- train(x=alldata[c(5,12:22)], 
                                  y=alldata$C_STOCKS, 
                                  method=customRF, 
                                  metric="RMSE", 
                                  tuneGrid=tunegrid_final, 
                                  trControl=controltenfold_15rep)
final_alldata_C2013_15rep
View(final_alldata_C2013_15rep$results)
final_alldata_C2013_15rep$finalModel

## 6.4. With 20 repetitions  ---------------------------------------------------
### 6.4.1 Model training -------------------------------------------------------
final_alldata_C2013_20rep<- train(x=alldata[c(5,12:22)], 
                                  y=alldata$C_STOCKS, 
                                  method=customRF, 
                                  metric="RMSE", 
                                  tuneGrid=tunegrid_final, 
                                  trControl=controltenfold_20rep)
final_alldata_C2013_20rep
View(final_alldata_C2013_20rep$results)
final_alldata_C2013_20rep$finalModel

## 6.5. With 25 repetitions  ---------------------------------------------------
### 6.5.1 Model training -------------------------------------------------------
final_alldata_C2013_25rep<- train(x=alldata[c(5,12:22)], 
                                  y=alldata$C_STOCKS, 
                                  method=customRF, 
                                  metric="RMSE", 
                                  tuneGrid=tunegrid_final, 
                                  trControl=controltenfold_25rep)
final_alldata_C2013_25rep
View(final_alldata_C2013_25rep$results)
final_alldata_C2013_25rep$finalModel

## 6.6. With 50 repetitions  ---------------------------------------------------
### 6.6.1 Model training -------------------------------------------------------
final_alldata_C2013_50rep<- train(x=alldata[c(5,12:22)], 
                                  y=alldata$C_STOCKS, 
                                  method=customRF, 
                                  metric="RMSE", 
                                  tuneGrid=tunegrid_final, 
                                  trControl=controltenfold_50rep)
final_alldata_C2013_50rep
View(final_alldata_C2013_50rep$results)
final_alldata_C2013_50rep$finalModel
## 6.8. With 100 repetitions  ---------------------------------------------------
### 6.8.1 Model training -------------------------------------------------------
final_alldata_C2013_100rep<- train(x=alldata[c(5,12:22)], 
                                   y=alldata$C_STOCKS, 
                                   method=customRF, 
                                   metric="RMSE", 
                                   tuneGrid=tunegrid_final, 
                                   trControl=controltenfold_100rep)
final_alldata_C2013_100rep
View(final_alldata_C2013_100rep$results)
final_alldata_C2013_100rep$finalModel


# ----------------------------------------------------------------------------------
# 7. RUN final MODELS TRAIN DATA  ! REPORTED ¡ --------------------------------------------
tunegrid_final= expand.grid(.ntree=c(900),.mtry=c(9))
View(data_train2013)
## 7.1. With 5 repetitions  ---------------------------------------------------
### 7.1.1 Model training -------------------------------------------------------
final_C2013_5rep<- train(x=data_train2013[c(5,12:22)], 
                         y=data_train2013$C_STOCKS, 
                         method=customRF, 
                         metric="RMSE", 
                         tuneGrid=tunegrid_final, 
                         trControl=controltenfold_5rep)
final_C2013_5rep
View(final_C2013_5rep$results)
final_C2013_5rep$finalModel


## 7.2. With 10 repetitions  ---------------------------------------------------
### 7.2.1 Model training -------------------------------------------------------
final_C2013_10rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_final, 
                          trControl=controltenfold_10rep)
final_C2013_10rep
View(final_C2013_10rep$results)
final_C2013_10rep$finalModel


## 7.3. With 15 repetitions  ---------------------------------------------------
### 7.3.1 Model training -------------------------------------------------------
final_C2013_15rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_final, 
                          trControl=controltenfold_15rep)
final_C2013_15rep
View(final_C2013_15rep$results)
final_C2013_15rep$finalModel

## 7.4. With 20 repetitions  ---------------------------------------------------
### 7.4.1 Model training -------------------------------------------------------
final_C2013_20rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_final, 
                          trControl=controltenfold_20rep)
final_C2013_20rep
View(final_C2013_20rep$results)
final_C2013_20rep$finalModel

## 7.5. With 25 repetitions  ---------------------------------------------------
### 7.5.1 Model training -------------------------------------------------------
final_C2013_25rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_final, 
                          trControl=controltenfold_25rep)
final_C2013_25rep
View(final_C2013_25rep$results)
final_C2013_25rep$finalModel

## 7.6. With 50 repetitions  ---------------------------------------------------
### 7.6.1 Model training -------------------------------------------------------
names(data_train2013)
final_C2013_50rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_final, 
                          trControl=controltenfold_50rep)
final_C2013_50rep
View(final_C2013_50rep$results)
final_C2013_50rep$finalModel
## 7.7. With 75 repetitions  ---------------------------------------------------
### 7.7.1 Model training -------------------------------------------------------
final_C2013_75rep<- train(x=data_train2013[c(5,12:22)], 
                          y=data_train2013$C_STOCKS, 
                          method=customRF, 
                          metric="RMSE", 
                          tuneGrid=tunegrid_default, 
                          trControl=controltenfold_75rep)
final_C2013_75rep
View(final_C2013_75rep$results)
final_C2013_75rep$finalModel

## 7.8. With 100 repetitions - SELECTED RF1  ---------------------------------------------------
### 7.8.1 Model training -------------------------------------------------------
names(data_train2013)
final_C2013_100rep<- train(x=data_train2013[c(5,12:22)], 
                           y=data_train2013$C_STOCKS, 
                           method=customRF, 
                           metric="RMSE", 
                           tuneGrid=tunegrid_final, 
                           trControl=controltenfold_100rep)
final_C2013_100rep
View(final_C2013_100rep$results)
final_C2013_100rep$finalModel


### 7.8.2. Observed vs predicted model  ------------------------------------
# Test data:
pred_C2013_FINAL_OPTIM_test=predict(final_C2013_100rep$finalModel,data_test2013)
obs_C2013_FINAL_OPTIM_test=data_test2013$C_STOCKS
obspred_C2013_FINAL_OPTIM_test=data.frame(obs=obs_C2013_FINAL_OPTIM_test,
                                          pred=pred_C2013_FINAL_OPTIM_test)
#Metric errors CV - independent
sqrt(mean((obs_C2013_FINAL_OPTIM_test - pred_C2013_FINAL_OPTIM_test)^2))  #RMSE
1 - sum((obs_C2013_FINAL_OPTIM_test - pred_C2013_FINAL_OPTIM_test)^2) / sum((obs_C2013_FINAL_OPTIM_test - mean(obs_C2013_FINAL_OPTIM_test))^2) #R2
mean(abs(obs_C2013_FINAL_OPTIM_test - pred_C2013_FINAL_OPTIM_test)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_FINAL_OPTIM_test$obs,obspred_C2013_FINAL_OPTIM_test$pred,method = "spearman")
cor(obspred_C2013_FINAL_OPTIM_test$obs,obspred_C2013_FINAL_OPTIM_test$pred,method = "pearson")


# train data:
pred_C2013_FINAL_OPTIM_train=predict(final_C2013_100rep$finalModel,data_train2013)
obs_C2013_FINAL_OPTIM_train=data_train2013$C_STOCKS
obspred_C2013_FINAL_OPTIM_train=data.frame(obs=obs_C2013_FINAL_OPTIM_train,
                                           pred=pred_C2013_FINAL_OPTIM_train)
#Metric errors CV - independent
sqrt(mean((obs_C2013_FINAL_OPTIM_train - pred_C2013_FINAL_OPTIM_train)^2))  #RMSE
1 - sum((obs_C2013_FINAL_OPTIM_train - pred_C2013_FINAL_OPTIM_train)^2) / sum((obs_C2013_FINAL_OPTIM_train - mean(obs_C2013_FINAL_OPTIM_train))^2) #R2
mean(abs(obs_C2013_FINAL_OPTIM_train - pred_C2013_FINAL_OPTIM_train)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_FINAL_OPTIM_train$obs,obspred_C2013_FINAL_OPTIM_train$pred,method = "spearman")
cor(obspred_C2013_FINAL_OPTIM_train$obs,obspred_C2013_FINAL_OPTIM_train$pred,method = "pearson")

### 7.8.3. ntree RMSE data -------------------------------

error_ntree_2013C_optimized=as.data.table(plot(final_C2013_100rep$finalModel))
error_ntree_2013C_optimized[, trees := .I]
error_ntree_2013C_optimized2 = melt(error_ntree_2013C_optimized, id.vars = "trees")
setnames(error_ntree_2013C_optimized2,"trees","error (MSE)")
error_ntree_2013C_optimized2=as.data.frame(error_ntree_2013C_optimized2)
names(error_ntree_2013C_optimized2)=c("trees","variable","error (MSE)")
max(sqrt(error_ntree_2013C_optimized2$`error (MSE)`)) #ylims
min(sqrt(error_ntree_2013C_optimized2$`error (MSE)`)) #ylims
ntreePLOT_C13_optimized_df=data.frame(error_ntree_2013C_optimized2,sqrt(error_ntree_2013C_optimized2$`error (MSE)`))
names(ntreePLOT_C13_optimized_df)=c("trees","variable","MSE","RMSE")
names(ntreePLOT_C13_optimized_df)
#write.csv(ntreePLOT_C13_optimized_df,here("Outputs", "ntreePLOT_C13_optimized.csv"),row.names = TRUE)

## 7.9 JOIN ALL RUNS -----------------------------------------------------------
final_C2013_5rep$results
final_models_CVREPS=rbind(data.frame(variable="C stocks",year=2013,model="RF1",CVreps=5,final_C2013_5rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=10,final_C2013_10rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=15,final_C2013_15rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=20,final_C2013_20rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=25,final_C2013_25rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=50,final_C2013_50rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=75,final_C2013_75rep$results,Modelselected="Descarted"),
                          data.frame(variable="C stocks",year=2013,model="RF1",CVreps=100,final_C2013_100rep$results,Modelselected="Reported RF1"))
View(final_models_CVREPS)
# ----------------------------------------------------------------------------------
# 8. DATA FOR CV REPETITIONS IMPACT ---------------------------------------------------
selected_models_CVREPS=rbind(default_models_CVREPS,final_models_CVREPS)
write.csv(selected_models_CVREPS,here("Final data", "selected_models_CVREPS_2013.csv"),row.names = TRUE)

# ----------------------------------------------------------------------------------
# 9. RUN FINAL MODEL WITH randomForest function ------------------------------------------------------
## 9.1. Model iteration  ------------------------------------------------------
#models_C13_1000it <- map(1:1000, ~ {print(.x); randomForest(C_STOCKS~LAYER+STAND_AGE+
#                                                           BASAL_AREA+DOM_HEIGHT+SPECIES_RICHNESS+
#                                                            SHANNON_INDEX+GAP_FRACTION+CANOPY_COVER+
#                                                            ELEVATION+SLOPE+ASPECT+CONDITION,
#                                                          data = data_train2013,
#                                                          mtry=12,ntree=600,importance=TRUE)})
# Capture the model outputs:
#model_output_OPTIM_C2013_1000R<- capture.output(models_C13_1000it)
# Combine the results in a character vector:
#model_output_OPTIM_C2013_1000R<- unlist(model_output_OPTIM_C2013_1000R)
# Write the results within a text file:
#writeLines(model_output_OPTIM_C2013_1000R,here("Outputs","OPTIM_C2013_1000rep.txt"))
#Select the best model within the repetitions:
(FINAL_RF1_C2013_1000rep=models_C13_1000it[[844]])


## 9.2. Observed vs predicted model  ------------------------------------
# Test data:
pred_C2013_FINAL_ALONE_test=predict(FINAL_RF1_C2013_1000rep,data_test2013)
obs_C2013_FINAL_ALONE_test=data_test2013$C_STOCKS
obspred_C2013_FINAL_ALONE_test=data.frame(obs=obs_C2013_FINAL_ALONE_test,
                                          pred=pred_C2013_FINAL_ALONE_test)
#Metric errors CV - independent
sqrt(mean((obs_C2013_FINAL_ALONE_test - pred_C2013_FINAL_ALONE_test)^2))  #RMSE
1 - sum((obs_C2013_FINAL_ALONE_test - pred_C2013_FINAL_ALONE_test)^2) / sum((obs_C2013_FINAL_ALONE_test - mean(obs_C2013_FINAL_ALONE_test))^2) #R2
mean(abs(obs_C2013_FINAL_ALONE_test - pred_C2013_FINAL_ALONE_test)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_FINAL_ALONE_test$obs,obspred_C2013_FINAL_ALONE_test$pred,method = "spearman")
cor(obspred_C2013_FINAL_ALONE_test$obs,obspred_C2013_FINAL_ALONE_test$pred,method = "pearson")


# train data:
pred_C2013_FINAL_ALONE_train=predict(FINAL_RF1_C2013_1000rep,data_train2013)
obs_C2013_FINAL_ALONE_train=data_train2013$C_STOCKS
obspred_C2013_FINAL_ALONE_train=data.frame(obs=obs_C2013_FINAL_ALONE_train,
                                           pred=pred_C2013_FINAL_ALONE_train)
#Metric errors CV - independent
sqrt(mean((obs_C2013_FINAL_ALONE_train - pred_C2013_FINAL_ALONE_train)^2))  #RMSE
1 - sum((obs_C2013_FINAL_ALONE_train - pred_C2013_FINAL_ALONE_train)^2) / sum((obs_C2013_FINAL_ALONE_train - mean(obs_C2013_FINAL_ALONE_train))^2) #R2
mean(abs(obs_C2013_FINAL_ALONE_train - pred_C2013_FINAL_ALONE_train)) #MAE
#How correlated are obs vs pred from cross validation method
cor(obspred_C2013_FINAL_ALONE_train$obs,obspred_C2013_FINAL_ALONE_train$pred,method = "spearman")
cor(obspred_C2013_FINAL_ALONE_train$obs,obspred_C2013_FINAL_ALONE_train$pred,method = "pearson")


## 9.3. Dataframe ObsPred  ------------------------------------
ObsPred_FINAL_BOTH=rbind(data.frame(variable="C stocks",year=2013,model="RF1",data="1. Train",obspred_C2013_FINAL_ALONE_train),
                         data.frame(variable="C stocks",year=2013,model="RF1",data="2. Test",obspred_C2013_FINAL_ALONE_test))

## 9.4. Ntree RMSE data --------------------------------

error_ntree_2013C_final=as.data.table(plot(FINAL_RF1_C2013_1000rep))
error_ntree_2013C_final[, trees := .I]
error_ntree_2013C_final2 = melt(error_ntree_2013C_final, id.vars = "trees")
setnames(error_ntree_2013C_final2,"trees","error (MSE)")
error_ntree_2013C_final2=as.data.frame(error_ntree_2013C_final2)
names(error_ntree_2013C_final2)=c("trees","variable","error (MSE)")
max(sqrt(error_ntree_2013C_final2$`error (MSE)`)) #ylims
min(sqrt(error_ntree_2013C_final2$`error (MSE)`)) #ylims
ntreePLOT_C13_final_df=data.frame(error_ntree_2013C_final2,sqrt(error_ntree_2013C_final2$`error (MSE)`))
names(ntreePLOT_C13_final_df)=c("trees","variable","MSE","RMSE")
names(ntreePLOT_C13_final_df)
#write.csv(ntreePLOT_C13_final_df,here("Outputs", "ntreePLOT_C13_final.csv"),row.names = TRUE)


## 9.5. Vimp data IncMSE criteria -----------------------------------------------------

abbreviations_predictors=c("FFLAYER","AGE","BA","DTH","RICH","H","GAP","COV",
                           "ELE","SLO","ASP","MAN")
categ_predictors=c("Forest floor stage","Canopy structure","Canopy structure",
                   "Canopy structure","Forest diversity","Forest diversity",
                   "Canopy structure","Canopy structure",
                   "Topographic features","Topographic features","Topographic features",
                   "Management condition")

VIMP_C2013_FINAL=data.frame(categ_predictors,abbreviations_predictors,
                            Variable="C stocks",Year=2013,Model="RF1",
                            importance(FINAL_RF1_C2013_1000rep))
View(VIMP_C2013_FINAL)
write.csv(VIMP_C2013_FINAL,here("Final data", "VIMP_C13_FINAL.csv"),row.names = TRUE)
rownames(VIMP_C2013_FINAL)=abbreviations_predictors
VIMP_C2013_FINAL$Year=as.factor(VIMP_C2013_FINAL$Year)

VIMP_C2013_FINAL$Plotorder=1
VIMP_C2013_FINAL$COVARS_ORDER=1
varImpPlot(FINAL_RF1_C2013_1000rep)
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "BA"]<- 11
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "AGE"]<- 10
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "COV"]<- 9
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "DTH"]<- 8
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "RICH"]<-7
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "H"]<- 6
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "MAN"]<- 5
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "ELE"]<- 4
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "GAP"]<- 3
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "SLO"]<- 2
VIMP_C2013_FINAL$COVARS_ORDER[VIMP_C2013_FINAL$abbreviations_predictors  == "ASP"]<- 1


ggplot(VIMP_C2013_FINAL, aes(x=reorder(abbreviations_predictors, COVARS_ORDER), y=X.IncMSE,color=Year)) + 
  geom_segment(aes(x=reorder(abbreviations_predictors, COVARS_ORDER), 
                   xend=abbreviations_predictors, 
                   color=Year,
                   y=0, 
                   yend=X.IncMSE)) + 
  geom_point(size = 4,pch=20,alpha=0.8) +
  coord_flip()+
  ylab("%IncMSE")+
  xlab("Predictor variables")+
  theme_bw()+
  scale_color_manual(values=c("#FF004D","#184e77","#76c893"),name="Year:")+
  facet_wrap(~reorder(Variable, Plotorder,FUN=max),ncol = 3,nrow = 1)+
  theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
        axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
        text = element_text(size=12,family = "A"),
        strip.text = element_text(size = 12,face = "bold"),
        legend.position="bottom",
        legend.title = element_text(size = 11, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
        panel.grid.major.y = element_line(colour = 'gray', linetype = 'dotted'),
        axis.text = element_text(size = 12))




# ----------------------------------------------------------------------------------
# 10. DATA FOR OBSERVED VS PREDICTED COMPARISON --------------------------------------
OBS_PRED_RFMODELS_ALL=rbind(ObsPred_FINAL_BOTH,ObsPred_DEFAULT_BOTH)
write.csv(OBS_PRED_RFMODELS_ALL,here("Final data", "OBS_PRED_RFMODELS_ALL_2013.csv"),row.names = TRUE)
names(OBS_PRED_RFMODELS_ALL)

# Kendall coefficient by data and model
kendall_results <- OBS_PRED_RFMODELS_ALL %>%
  group_by(data, model) %>%
  summarize(kendall_tau = cor.test(obs, pred, method = "kendall")$estimate)
print(kendall_results)

# 11. DATA FOR NTREE ERROR ESTABILIZATION --------------------------------------
ntreePLOT_C13_default_df
ntreePLOT_C13_optimized_df
ntreePLOT_C13_final_df

# Plot all models
all_ntrees=rbind(data.frame(ntreePLOT_C13_default_df,model="1. Initial model",Targetvariable="C stocks",Year=2013),
                 data.frame(ntreePLOT_C13_optimized_df,model="3. Optimized model",Targetvariable="C stocks",Year=2013),
                 data.frame(ntreePLOT_C13_final_df,model="2. Final model",Targetvariable="C stocks",Year=2013))
View(all_ntrees)

ggplot(data = all_ntrees, aes(x = trees, y = RMSE,color=model)) +
  geom_line(linewidth=0.52)+
  theme_bw()+
  xlim(0,1000)+
  scale_y_continuous(limits = c(3.1, 5.7),breaks = c(3.5,4,4.5,5,5.5)) +
  ggtitle("2013\n",subtitle = "a) C stocks")+
  ylab(Mgha_expression3)+
  xlab("")+
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 200)) +
  annotate("point",y=3.149019,x=600,color="blue",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
  annotate("point",y=3.296195,x=500,color="red",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
  annotate("point",y=3.240893,x=600,color="green",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
  scale_color_manual(values = c("1. Initial model"="red",
                                "2. Final model"="blue",
                                "3. Optimized model"="green"),
                     name="Random forest models:",
                     labels=c(initialmodel_expresion,
                              finalmodel_expression,
                              "Optimized mode")) +
  theme(text = element_text(size=11,family = "A"),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
        panel.background = element_rect(fill='white', colour='black'))

# Save the data for final and default model
ntree_data_bothmodels=rbind(data.frame(variable="C stocks",year=2013,model="RF0",ntreePLOT_C13_default_df[c(1,3:4)]),
                            data.frame(variable="C stocks",year=2013,model="RF1",ntreePLOT_C13_final_df[c(1,3:4)]))
write.csv(ntree_data_bothmodels,here("Final data", "ntree_data_bothmodels_2013.csv"),row.names = TRUE)
# 12. ANOVA TWO-WAY -------------------------------------
df_finalmodels_C2013 <- selected_models_CVREPS %>%
  mutate(CVreps = as.numeric(CVreps),
         ntree = as.factor(ntree))
anova_results_C2013 <- aov(RMSE ~ CVreps * model, data = df_finalmodels_C2013)
summary(anova_results_C2013)
emmeans(anova_results_C2013, pairwise ~ model)
shapiro.test(residuals(anova_results_C2013)) #Normal residuals 

# 13. FRIEDMAN TEST -----------------------------------
df_finalmodels_C2013$model <- as.factor(df_finalmodels_C2013$model)
df_finalmodels_C2013$CVreps <- as.factor(df_finalmodels_C2013$CVreps)

friedman_C2013=friedman_test(RMSE ~ model | CVreps, data = df_finalmodels_C2013)
print(friedman_C2013)
kendalls_w <- statistic(friedman_C2013) / (nrow(df_finalmodels_C2013) * (nlevels(df_finalmodels_C2013$model) - 1))
kendalls_w # Calculate Kendall's W

friedman2_C2013=friedman_test(MAE ~ model | CVreps, data = df_finalmodels_C2013)
print(friedman2_C2013)
kendalls_w2 <- statistic(friedman2_C2013) / (nrow(df_finalmodels_C2013) * (nlevels(df_finalmodels_C2013$model) - 1))
kendalls_w2 # Calculate Kendall's W

friedman3_C2013=friedman_test(Rsquared ~ model | CVreps, data = df_finalmodels_C2013)
print(friedman2_C2013)
kendalls_w3 <- statistic(friedman3_C2013) / (nrow(df_finalmodels_C2013) * (nlevels(df_finalmodels_C2013$model) - 1))
kendalls_w3 # Calculate Kendall's W

# 14. POSTHOC TESTS
for (i in unique(df_finalmodels_C2013$CVreps)) {
  print(paste("Comparaciones para la repetición", i))
  print(pairwise.wilcox.test(df_finalmodels_C2013$RMSE[df_finalmodels_C2013$CVreps == i], 
                             df_finalmodels_C2013$model[df_finalmodels_C2013$CVreps == i], p.adjust.method = "bonferroni"))
}

# post-hoc comparisons between repetitions within each model RFO and RF1
for (i in unique(df_finalmodels_C2013$model)) {
  print(paste("Comparaciones para el modelo", i))
  print(pairwise.wilcox.test(df_finalmodels_C2013$RMSE[df_finalmodels_C2013$model == i], 
                             df_finalmodels_C2013$CVreps[df_finalmodels_C2013$model == i], p.adjust.method = "bonferroni"))
}

wilcox.test(df_finalmodels_C2013$RMSE[df_finalmodels_C2013$model == "RF0"], 
            df_finalmodels_C2013$RMSE[df_finalmodels_C2013$model == "RF1"])

# -------------------------------------CHAPTER 2---------------------------------------------
# 14. PARTIAL PLOTS -----------------------------------------------------
?partialPlot()

## 14.1. INDIVIDUAL EFFECT --------------------------------------------------------------------

#Layer
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "LAYER", main = "a",ylab = "Mg/ha",xlab = "Forest floor layer")
partial(FINAL_RF1_C2013_1000rep,pred.var = "LAYER",plot=TRUE)
(LAYER_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "LAYER"))
LAYER_C13_df=data.frame(LAYER_C2013,Year=("2013"),Variable="C stocks")
write.csv(LAYER_C13_df,here("Partialplots", "PARTIAL_LAYER_C13.csv"),row.names = TRUE)

#Stand age
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "STAND_AGE", main = "a",ylab = "Mg/ha",xlab = "Stand age")
partial(FINAL_RF1_C2013_1000rep,pred.var = "STAND_AGE",plot=TRUE)
(STAND_AGE_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "STAND_AGE"))
STAND_AGE_C13_df=data.frame(STAND_AGE_C2013,Year=("2013"),Variable="C stocks")
write.csv(STAND_AGE_C13_df,here("Partialplots", "PARTIAL_STAND_AGE_C13.csv"),row.names = TRUE)

#Basal area
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "BASAL_AREA", main = "a",ylab = "Mg/ha",xlab = "Basal area")
partial(FINAL_RF1_C2013_1000rep,pred.var = "BASAL_AREA",plot=TRUE)
(BASAL_AREA_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "BASAL_AREA"))
BASAL_AREA_C13_df=data.frame(BASAL_AREA_C2013,Year=("2013"),Variable="C stocks")
write.csv(BASAL_AREA_C13_df,here("Partialplots", "PARTIAL_BASAL_AREA_C13.csv"),row.names = TRUE)

#Dominant tree height
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "DOM_HEIGHT", main = "a",ylab = "Mg/ha",xlab = "Dominant height")
partial(FINAL_RF1_C2013_1000rep,pred.var = "DOM_HEIGHT",plot=TRUE)
(DOM_HEIGHT_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "DOM_HEIGHT"))
DOM_HEIGHT_C13_df=data.frame(DOM_HEIGHT_C2013,Year=("2013"),Variable="C stocks")
write.csv(DOM_HEIGHT_C13_df,here("Partialplots", "PARTIAL_DOM_HEIGHT_C13.csv"),row.names = TRUE)

#Species richness
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "SPECIES_RICHNESS", main = "a",ylab = "Mg/ha",xlab = "Species richness")
partial(FINAL_RF1_C2013_1000rep,pred.var = "SPECIES_RICHNESS",plot=TRUE)
(SPECIES_RICHNESS_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "SPECIES_RICHNESS"))
SPECIES_RICHNESS_C13_df=data.frame(SPECIES_RICHNESS_C2013,Year=("2013"),Variable="C stocks")
write.csv(SPECIES_RICHNESS_C13_df,here("Partialplots", "PARTIAL_SPECIES_RICHNESS_C13.csv"),row.names = TRUE)

#Shannon
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "SHANNON_INDEX", main = "a",ylab = "Mg/ha",xlab = "Shannon's diversity index")
partial(FINAL_RF1_C2013_1000rep,pred.var = "SHANNON_INDEX",plot=TRUE)
(SHANNON_INDEX_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "SHANNON_INDEX"))
SHANNON_INDEX_C13_df=data.frame(SHANNON_INDEX_C2013,Year=("2013"),Variable="C stocks")
write.csv(SHANNON_INDEX_C13_df,here("Partialplots", "PARTIAL_SHANNON_INDEX_C13.csv"),row.names = TRUE)

#Gap fraction
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "GAP_FRACTION", main = "a",ylab = "Mg/ha",xlab = "Gap fraction")
partial(FINAL_RF1_C2013_1000rep,pred.var = "GAP_FRACTION",plot=TRUE)
(GAP_FRACTION_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "GAP_FRACTION"))
GAP_FRACTION_C13_df=data.frame(GAP_FRACTION_C2013,Year=("2013"),Variable="C stocks")
write.csv(GAP_FRACTION_C13_df,here("Partialplots", "PARTIAL_GAP_FRACTION_C13.csv"),row.names = TRUE)

#Canopy cover
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "CANOPY_COVER", main = "a",ylab = "Mg/ha",xlab = "Canopy cover")
partial(FINAL_RF1_C2013_1000rep,pred.var = "CANOPY_COVER",plot=TRUE)
(CANOPY_COVER_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "CANOPY_COVER"))
CANOPY_COVER_C13_df=data.frame(CANOPY_COVER_C2013,Year=("2013"),Variable="C stocks")
write.csv(CANOPY_COVER_C13_df,here("Partialplots", "PARTIAL_CANOPY_COVER_C13.csv"),row.names = TRUE)

#Elevation
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "ELEVATION", main = "a",ylab = "Mg/ha",xlab = "Elevation")
partial(FINAL_RF1_C2013_1000rep,pred.var = "ELEVATION",plot=TRUE)
(ELEVATION_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "ELEVATION"))
ELEVATION_C13_df=data.frame(ELEVATION_C2013,Year=("2013"),Variable="C stocks")
write.csv(ELEVATION_C13_df,here("Partialplots", "PARTIAL_ELEVATION_C13.csv"),row.names = TRUE)

#Slope
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "SLOPE", main = "a",ylab = "Mg/ha",xlab = "Slope")
partial(FINAL_RF1_C2013_1000rep,pred.var = "SLOPE",plot=TRUE)
(SLOPE_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "SLOPE"))
SLOPE_C13_df=data.frame(SLOPE_C2013,Year=("2013"),Variable="C stocks")
write.csv(SLOPE_C13_df,here("Partialplots", "PARTIAL_SLOPE_C13.csv"),row.names = TRUE)

#Aspect
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "ASPECT", main = "a",ylab = "Mg/ha",xlab = "Aspect")
partial(FINAL_RF1_C2013_1000rep,pred.var = "ASPECT",plot=TRUE)
(ASPECT_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "ASPECT"))
ASPECT_C13_df=data.frame(ASPECT_C2013,Year=("2013"),Variable="C stocks")
write.csv(ASPECT_C13_df,here("Partialplots", "PARTIAL_ASPECT_C13.csv"),row.names = TRUE)

#Condition
partialPlot(FINAL_RF1_C2013_1000rep, data_train2013, x.var = "CONDITION", main = "a",ylab = "Mg/ha",xlab = "Management condition")
partial(FINAL_RF1_C2013_1000rep,pred.var = "CONDITION",plot=TRUE)
(CONDITION_C2013=partial(FINAL_RF1_C2013_1000rep,pred.var = "CONDITION"))
CONDITION_C13_df=data.frame(CONDITION_C2013,Year=("2013"),Variable="C stocks")
write.csv(CONDITION_C13_df,here("Partialplots", "PARTIAL_CONDITION_C13.csv"),row.names = TRUE)


## 14.2. EFFECT BY LAYER  --------------------------------------------------------------------

#Layer and Stand age
partial(FINAL_RF1_C2013_1000rep,pred.var = c("STAND_AGE","LAYER"),plot = TRUE)
SA_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("STAND_AGE","LAYER"))
SA_LAYER_C13_df=data.frame(SA_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(SA_LAYER_C13_df,here("Partialplots", "PARTIAL_SA_LAYER_C13.csv"),row.names = TRUE)

#Layer and Basal area
partial(FINAL_RF1_C2013_1000rep,pred.var = c("BASAL_AREA","LAYER"),plot = TRUE)
BA_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("BASAL_AREA","LAYER"))
BA_LAYER_C13_df=data.frame(BA_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(BA_LAYER_C13_df,here("Partialplots", "PARTIAL_BA_LAYER_C13.csv"),row.names = TRUE)

#Layer and Dom height
partial(FINAL_RF1_C2013_1000rep,pred.var = c("DOM_HEIGHT","LAYER"),plot = TRUE)
DTH_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("DOM_HEIGHT","LAYER"))
DHT_LAYER_C13_df=data.frame(DTH_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(DTH_LAYER_C13_df,here("Partialplots", "PARTIAL_DHT_LAYER_C13.csv"),row.names = TRUE)

#Layer and species richness
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SPECIES_RICHNESS","LAYER"),plot = TRUE)
RICH_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SPECIES_RICHNESS","LAYER"))
RICH_LAYER_C13_df=data.frame(RICH_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(RICH_LAYER_C13_df,here("Partialplots", "PARTIAL_RICH_LAYER_C13.csv"),row.names = TRUE)

#Layer and shannon
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SHANNON_INDEX","LAYER"),plot = TRUE)
SHANNON_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SHANNON_INDEX","LAYER"))
SHANNON_LAYER_C13_df=data.frame(SHANNON_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(SHANNON_LAYER_C13_df,here("Partialplots", "PARTIAL_SHANNON_LAYER_C13.csv"),row.names = TRUE)

#Layer and gap fraction
partial(FINAL_RF1_C2013_1000rep,pred.var = c("GAP_FRACTION","LAYER"),plot = TRUE)
GAP_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("GAP_FRACTION","LAYER"))
GAP_LAYER_C13_df=data.frame(GAP_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(GAP_LAYER_C13_df,here("Partialplots", "PARTIAL_GAP_LAYER_C13.csv"),row.names = TRUE)

#Layer and canopy cover
partial(FINAL_RF1_C2013_1000rep,pred.var = c("CANOPY_COVER","LAYER"),plot = TRUE)
COV_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("CANOPY_COVER","LAYER"))
COV_LAYER_C13_df=data.frame(COV_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(COV_LAYER_C13_df,here("Partialplots", "PARTIAL_COV_LAYER_C13.csv"),row.names = TRUE)

#Layer and elevation
partial(FINAL_RF1_C2013_1000rep,pred.var = c("ELEVATION","LAYER"),plot = TRUE)
ELE_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("ELEVATION","LAYER"))
ELE_LAYER_C13_df=data.frame(ELE_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(ELE_LAYER_C13_df,here("Partialplots", "PARTIAL_ELE_LAYER_C13.csv"),row.names = TRUE)

#Layer and aspect
partial(FINAL_RF1_C2013_1000rep,pred.var = c("ASPECT","LAYER"),plot = TRUE)
ASP_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("ASPECT","LAYER"))
ASP_LAYER_C13_df=data.frame(ASP_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(ASP_LAYER_C13_df,here("Partialplots", "PARTIAL_ASP_LAYER_C13.csv"),row.names = TRUE)

#Layer and slope
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SLOPE","LAYER"),plot = TRUE)
SLO_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SLOPE","LAYER"))
SLO_LAYER_C13_df=data.frame(SLO_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(SLO_LAYER_C13_df,here("Partialplots", "PARTIAL_SLO_LAYER_C13.csv"),row.names = TRUE)

#Layer and condition
partial(FINAL_RF1_C2013_1000rep,pred.var = c("CONDITION","LAYER"),plot = TRUE)
MCON_LAYER_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("CONDITION","LAYER"))
MCON_LAYER_C13_df=data.frame(MCON_LAYER_C13,Year=("2013"),Variable="C stocks")
write.csv(MCON_LAYER_C13_df,here("Partialplots", "PARTIAL_MCOND_LAYER_C13.csv"),row.names = TRUE)


## 14.3. EFFECT BY STAND AGE   --------------------------------------------------------------------

#Stand age and Basal area
partial(FINAL_RF1_C2013_1000rep,pred.var = c("BASAL_AREA","STAND_AGE"),plot = TRUE)
BA_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("BASAL_AREA","STAND_AGE"))
BA_SA_C13_df=data.frame(BA_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(BA_SA_C13_df,here("Partialplots", "PARTIAL_BA_SA_C13.csv"),row.names = TRUE)

#Stand age and Dom height
partial(FINAL_RF1_C2013_1000rep,pred.var = c("DOM_HEIGHT","STAND_AGE"),plot = TRUE)
DTH_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("DOM_HEIGHT","STAND_AGE"))
DTH_SA_C13_df=data.frame(DTH_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(DTH_SA_C13_df,here("Partialplots", "PARTIAL_DHT_SA_C13.csv"),row.names = TRUE)

#Stand age and species richness
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SPECIES_RICHNESS","STAND_AGE"),plot = TRUE)
RICH_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SPECIES_RICHNESS","STAND_AGE"))
RICH_SA_C13_df=data.frame(RICH_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(RICH_SA_C13_df,here("Partialplots", "PARTIAL_RICH_SA_C13.csv"),row.names = TRUE)

#Stand age and shannon
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SHANNON_INDEX","STAND_AGE"),plot = TRUE)
SHANNON_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SHANNON_INDEX","STAND_AGE"))
SHANNON_SA_C13_df=data.frame(SHANNON_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(SHANNON_SA_C13_df,here("Partialplots", "PARTIAL_SHANNON_SA_C13.csv"),row.names = TRUE)

#Stand age and gap fraction
partial(FINAL_RF1_C2013_1000rep,pred.var = c("GAP_FRACTION","STAND_AGE"),plot = TRUE)
GAP_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("GAP_FRACTION","STAND_AGE"))
GAP_SA_C13_df=data.frame(GAP_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(GAP_SA_C13_df,here("Partialplots", "PARTIAL_GAP_SA_C13.csv"),row.names = TRUE)

#Stand age and canopy cover
partial(FINAL_RF1_C2013_1000rep,pred.var = c("CANOPY_COVER","STAND_AGE"),plot = TRUE)
COV_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("CANOPY_COVER","STAND_AGE"))
COV_SA_C13_df=data.frame(COV_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(COV_SA_C13_df,here("Partialplots", "PARTIAL_COV_SA_C13.csv"),row.names = TRUE)

#Stand age and elevation
partial(FINAL_RF1_C2013_1000rep,pred.var = c("ELEVATION","STAND_AGE"),plot = TRUE)
ELE_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("ELEVATION","STAND_AGE"))
ELE_SA_C13_df=data.frame(ELE_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(ELE_SA_C13_df,here("Partialplots", "PARTIAL_ELE_SA_C13.csv"),row.names = TRUE)

#Stand age and aspect
partial(FINAL_RF1_C2013_1000rep,pred.var = c("ASPECT","STAND_AGE"),plot = TRUE)
ASP_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("ASPECT","STAND_AGE"))
ASP_SA_C13_df=data.frame(ASP_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(ASP_SA_C13_df,here("Partialplots", "PARTIAL_ASP_SA_C13.csv"),row.names = TRUE)

#Stand age and slope
partial(FINAL_RF1_C2013_1000rep,pred.var = c("SLOPE","STAND_AGE"),plot = TRUE)
SLO_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("SLOPE","STAND_AGE"))
SLO_SA_C13_df=data.frame(SLO_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(SLO_SA_C13_df,here("Partialplots", "PARTIAL_SLO_SA_C13.csv"),row.names = TRUE)

#Stand age and condition
partial(FINAL_RF1_C2013_1000rep,pred.var = c("CONDITION","STAND_AGE"),plot = TRUE)
MCON_SA_C13=partial(FINAL_RF1_C2013_1000rep,pred.var = c("CONDITION","STAND_AGE"))
MCON_SA_C13_df=data.frame(MCON_SA_C13,Year=("2013"),Variable="C stocks")
write.csv(MCON_SA_C13_df,here("Partialplots", "PARTIAL_MCOND_SA_C13.csv"),row.names = TRUE)


## 14.4 plots -----------------------------------------------------

#Descarted
rwb <- colorRampPalette(c("red", "white", "blue"))
ggplot(BA_SA_C13, aes(x = STAND_AGE, y = BASAL_AREA, z = yhat, fill = yhat)) +
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.5) + 
  scale_fill_distiller(name = "Centered\nlogit", palette = "Spectral") + 
  theme_bw() + 
  #facet_grid(~LAYER)+
  geom_rug()


#3D plot

dens_BA_SA_C13 <- interp(x = BA_SA_C13$STAND_AGE, y = BA_SA_C13$BASAL_AREA, z = BA_SA_C13$yhat)

# 3D partial dependence plot with a coloring scale
BA_SA_C13 <- plot_ly(x = dens_BA_SA_C13$x, 
                     y = dens_BA_SA_C13$y, 
                     z = dens_BA_SA_C13$z,
                     colors = c("blue", "grey", "red"),
                     type = "surface")
# Add axis labels for 3D plots
BA_SA_C13 <- BA_SA_C13 %>% layout(scene = list(xaxis = list(title = "Stand age"),
                                               yaxis = list(title = "Basal area"),
                                               zaxis = list(title = "C stocks")))
# Show the plot
show(BA_SA_C13)


#Another method
plotmo(FINAL_RF1_C2013_1000rep, pmethod="partdep", all1=FALSE, all2=1)

