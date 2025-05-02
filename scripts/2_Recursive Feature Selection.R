################################################################################
# Script: 2_Recursive Feature Selection.R
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
library(Boruta) 
library(ggplot2)
library(dplyr)
library(here)
library(readxl)

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


# 2. BORUTA ALGORITH WITH SEVERAL ITERATIONS -----------------------------------
#MaxRuns = 100
borutaRF_traindf_C13_100runs=Boruta(C~UTMX+UTMY+ELE+SLO+ASP+AGE+DTH+BA+GAP+COV+RICH+H+MAN+FFLAYER,
                                    data=data2013_formodels,
                                    doTrace=3,pValue=0.5,ntree = 500, mtry=4)
print(borutaRF_traindf_C13_100runs$finalDecision)
(BORzscores_traindf_C13_100runs<-attStats(borutaRF_traindf_C13_100runs))

#MaxRuns = 500
borutaRF_traindf_C13_500runs=Boruta(C~UTMX+UTMY+ELE+SLO+ASP+AGE+DTH+BA+GAP+COV+RICH+H+MAN+FFLAYER,
                                    data=data2013_formodels,
                                    doTrace=3,pValue=0.5,ntree = 500, mtry=4,maxRuns=500)
print(borutaRF_traindf_C13_500runs$finalDecision)
(BORzscores_traindf_C13_500runs<-attStats(borutaRF_traindf_C13_500runs))
plot(borutaRF_traindf_C13_500runs)

#MaxRuns = 1000
borutaRF_traindf_C13_1000runs=Boruta(C~UTMX+UTMY+ELE+SLO+ASP+AGE+DTH+BA+GAP+COV+RICH+H+MAN+FFLAYER,
                                     data=data2013_formodels,
                                     doTrace=3,pValue=0.5,ntree = 500, mtry=4,maxRuns=1000)
print(borutaRF_traindf_C13_1000runs$finalDecision)
(BORzscores_traindf_C13_1000runs<-attStats(borutaRF_traindf_C13_1000runs))

