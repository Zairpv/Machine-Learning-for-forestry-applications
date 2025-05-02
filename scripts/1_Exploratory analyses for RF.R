################################################################################
# Script: 1_Exploratory analyses for RF.R
# Author: Dra. Zaira Rosario Pérez-Vázquez
# Contact: zairpv@gmail.com
#
# Description:
# This script performs exploratory data analysis (EDA) to prepare input variables
# for a Random Forest (RF) model applied to forest and ecological data.
##
# Note:
# This is part of a larger portfolio on Machine Learning for Forestry Applications.
################################################################################

# 0. LIBRARIES -----------------------------------------------------------------
library(readr)
library(psych)
library(here)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(grid)
library(car)
library(corrplot)
library(ggcorrplot)
library(Hmisc)
library(PerformanceAnalytics)
library(xtable)
library(moments)
library(stats)
library(cowplot)
library(reshape2)
library(tidyverse)
library(paletteer)
library(readxl)


# 1.DATA AND CONSOLE PROPERTIES ------------------------------------------------
options(max.print=1000000)
windowsFonts(A = windowsFont("Times New Roman"))
data_path <- "data/RF_data_Carbon_ForestFloor.csv"
datasetFF <- read_csv(data_path)
glimpse(datasetFF)
summary(datasetFF)

# 2. PREPARE DATA---------------------------------------------------------------
## 2.1. Defining variables as factors
datasetFF$YEAR=as.factor(datasetFF$YEAR)
datasetFF$LAYER=as.factor(datasetFF$LAYER)
datasetFF$CONDITION=as.factor(datasetFF$CONDITION)
datasetFF$LAYER_NUM=as.factor(datasetFF$LAYER_NUM)
datasetFF$CON_NUM=as.factor(datasetFF$CON_NUM)

## 2.2. Split by year
data2013=datasetFF[datasetFF$YEAR=="2013",]
data2018=datasetFF[datasetFF$YEAR=="2018",]  
data2023=datasetFF[datasetFF$YEAR=="2023",]

## 2.3. Split by layer
data2013_L=data2013[data2013$LAYER=="L",]
data2013_FH=data2013[data2013$LAYER=="FH",]
data2018_L=data2018[data2018$LAYER=="L",]
data2018_FH=data2018[data2018$LAYER=="FH",]
data2023_L=data2023[data2023$LAYER=="L",]
data2023_FH=data2023[data2023$LAYER=="FH",]

# 3. EXPLORATORY ANALYSES ------------------------------------------------------

## 3.1. Descriptive statistics 
names(datasetFF)
descriptives_mainvariables=describeBy(datasetFF$C_STOCKS,list(datasetFF$YEAR,datasetFF$LAYER)) #Target variables
descriptives_mainvariables
descriptives_mainvariables_output<- capture.output(descriptives_mainvariables)
descriptives_mainvariables_output<- unlist(descriptives_mainvariables_output)
writeLines(descriptives_mainvariables_output,here("outputs","descriptives_main_variables.txt"))

descriptives_predictors=describeBy(datasetFF[6:20],list(datasetFF$YEAR,datasetFF$LAYER)) #Predictors
descriptives_predictors_output<- capture.output(descriptives_predictors)
descriptives_predictors_output<- unlist(descriptives_predictors_output)
writeLines(descriptives_predictors_output,here("outputs","descriptives_predictors_variables.txt"))

## 3.2. Interquantile ranges
quantile(data2013_L$C_STOCKS)
quantile(data2013_FH$C_STOCKS)
quantile(data2018_L$C_STOCKS)
quantile(data2018_FH$C_STOCKS)
quantile(data2023_L$C_STOCKS)
quantile(data2023_FH$C_STOCKS)

## 3.3. Coef. of variation
CoefVar(data2013_L$C_STOCKS)
CoefVar(data2013_FH$C_STOCKS)
CoefVar(data2018_L$C_STOCKS)
CoefVar(data2018_FH$C_STOCKS)
CoefVar(data2023_L$C_STOCKS)
CoefVar(data2023_FH$C_STOCKS)

# 4. DATA DISTRIBUTION ---------------------------------------------------------
Mgha_expression=expression(bold("a) C stocks (Mg ha"^"-1"~")"))
hist(data2013_L$C_STOCKS,xlab = Mgha_expression,main = "2013 - L")
hist(data2013_FH$C_STOCKS,xlab = Mgha_expression,main = "2013 - FH")
hist(data2018_L$C_STOCKS,xlab = Mgha_expression,main = "2018 - L")
hist(data2018_FH$C_STOCKS,xlab = Mgha_expression,main = "2018 - FH")
hist(data2023_L$C_STOCKS,xlab = Mgha_expression,main = "2023 - L")
hist(data2023_FH$C_STOCKS,xlab = Mgha_expression,main = "2023 - FH")

## 4.1. Normality
Mgha_expression1=expression("Mg ha"^"-1")

### L - 2013 
normCL13=seq(min(data2013_L$C_STOCKS),max(data2013_L$C_STOCKS),length=length(data2013_L$C_STOCKS))
fun_normCL13=dnorm(normCL13,mean = mean(data2013_L$C_STOCKS),sd=sd(data2013_L$C_STOCKS))
hist(data2013_L$C_STOCKS,probability = T,xlab = "p-value=0.0689",main = "2013")
lines(density(data2013_L$C_STOCKS),col=4,lwd=2)
lines(normCL13,fun_normCL13,col=2,lwd=2)
ks.test(data2013_L$C_STOCKS,pnorm,mean(data2013_L$C_STOCKS),sd(data2013_L$C_STOCKS))

(HIST_CL13=ggplot(data2013_L, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 0.5, fill = "lightgray") +
    stat_function(aes(color = "theoretical normal distribution"), 
                  fun = dnorm, args = list(mean = mean(data2013_L$C_STOCKS),
                                           sd = sd(data2013_L$C_STOCKS)),  size = 0.5)+
    geom_density(aes(color = "density"), size = 0.5,position = "identity") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T), color = "mean"),   # Ignore NA values for mean
               linetype="dashed", size=0.8)+
    scale_color_manual(name = "Data distribution\n               L layer:", 
                       values = c("density" = "blue",
                                  "theoretical normal distribution" = "red",
                                  "mean" = "orange"),
                       labels = c("Observed density", 
                                  "Mean value",
                                  "Theoretical normal"))+
    labs(x = Mgha_expression1, y = "Density")+
    labs(title = "2013")+
    annotate("text", x=4.85, y=0.68, label= "p = 0.069",family = "A") +
    theme_bw()+
    theme(text = element_text(size=12,family = "A"),
          legend.position = "bottom",
          legend.title = element_text(hjust = 0.5,size = 11,face = "bold",family = "A"),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5,family = "A"),
          plot.subtitle = element_text(size = 13,face = "bold",hjust = 0,family = "A"),
          axis.text = element_text(size = 12)))

### L - 2018 
normCL18=seq(min(data2018_L$C_STOCKS),max(data2018_L$C_STOCKS),length=length(data2018_L$C_STOCKS))
fun_normCL18=dnorm(normCL18,mean = mean(data2018_L$C_STOCKS),sd=sd(data2018_L$C_STOCKS))
hist(data2018_L$C_STOCKS,probability = T,xlab = "p-value=0.0689",main = "2018")
lines(density(data2018_L$C_STOCKS),col=4,lwd=2)
lines(normCL18,fun_normCL18,col=2,lwd=2)
ks.test(data2018_L$C_STOCKS,pnorm,mean(data2018_L$C_STOCKS),sd(data2018_L$C_STOCKS))

(HIST_CL18=ggplot(data2018_L, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 0.5, fill = "lightgray") +
    stat_function(aes(color = "theoretical normal distribution"), 
                  fun = dnorm, args = list(mean = mean(data2018_L$C_STOCKS),
                                           sd = sd(data2018_L$C_STOCKS)),  size = 0.5)+
    geom_density(aes(color = "density"), size = 0.5,position = "identity") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T), color = "mean"),   # Ignore NA values for mean
               linetype="dashed", size=0.8)+
    scale_color_manual(name = "Data distribution\n               L layer:", 
                       values = c("density" = "blue",
                                  "theoretical normal distribution" = "red",
                                  "mean" = "orange"),
                       labels = c("Observed density", 
                                  "Mean value",
                                  "Theoretical normal"))+
    labs(x = Mgha_expression1, y = "Density")+
    labs(title = "2018")+
    annotate("text", x=4.85, y=0.68, label= "p = 0.069",family = "A") +
    theme_bw()+
    theme(text = element_text(size=12,family = "A"),
          legend.position = "bottom",
          legend.title = element_text(hjust = 0.5,size = 11,face = "bold",family = "A"),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5,family = "A"),
          plot.subtitle = element_text(size = 13,face = "bold",hjust = 0,family = "A"),
          axis.text = element_text(size = 12)))


### L - 2023
normCL23=seq(min(data2023_L$C_STOCKS),max(data2023_L$C_STOCKS),length=length(data2023_L$C_STOCKS))
fun_normCL23=dnorm(normCL23,mean = mean(data2023_L$C_STOCKS),sd=sd(data2023_L$C_STOCKS))
hist(data2023_L$C_STOCKS,probability = T,xlab = "p-value=0.0689",main = "2023")
lines(density(data2023_L$C_STOCKS),col=4,lwd=2)
lines(normCL23,fun_normCL23,col=2,lwd=2)
ks.test(data2023_L$C_STOCKS,pnorm,mean(data2023_L$C_STOCKS),sd(data2023_L$C_STOCKS))

(HIST_CL23=ggplot(data2023_L, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 0.5, fill = "lightgray") +
    stat_function(aes(color = "theoretical normal distribution"), 
                  fun = dnorm, args = list(mean = mean(data2023_L$C_STOCKS),
                                           sd = sd(data2023_L$C_STOCKS)),  size = 0.5)+
    geom_density(aes(color = "density"), size = 0.5,position = "identity") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T), color = "mean"),   # Ignore NA values for mean
               linetype="dashed", size=0.8)+
    scale_color_manual(name = "Data distribution\n               L layer:", 
                       values = c("density" = "blue",
                                  "theoretical normal distribution" = "red",
                                  "mean" = "orange"),
                       labels = c("Observed density", 
                                  "Mean value",
                                  "Theoretical normal"))+
    labs(x = Mgha_expression1, y = "Density")+
    labs(title = "2023")+
    annotate("text", x=4.85, y=0.68, label= "p = 0.069",family = "A") +
    theme_bw()+
    theme(text = element_text(size=12,family = "A"),
          legend.position = "bottom",
          legend.title = element_text(hjust = 0.5,size = 11,face = "bold",family = "A"),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5,family = "A"),
          plot.subtitle = element_text(size = 13,face = "bold",hjust = 0,family = "A"),
          axis.text = element_text(size = 12)))



jpeg(filename = here("outputs","Histograms - L layer.jpeg"),width = 180,height = 70,units = "mm",res = 1000)
Supp1=ggarrange(HIST_CL13,HIST_CL18,HIST_CL23,
                ncol = 3,nrow = 1,
                common.legend = TRUE,legend="bottom")
Supp1
dev.off()


#C stocks - FH 13
normCFH13=seq(min(data2013_FH$C_STOCKS),max(data2013_FH$C_STOCKS),length=length(data2013_FH$C_STOCKS))
fun_normCFH13=dnorm(normCFH13,mean = mean(data2013_FH$C_STOCKS),sd=sd(data2013_FH$C_STOCKS))
hist(data2013_FH$C_STOCKS,probability = T,xlab = Mgha_expression,main = "2013 - FH")
lines(density(data2013_FH$C_STOCKS),col=4,lwd=2)
lines(normCFH13,fun_normCFH13,col=2,lwd=2)
ks.test(data2013_FH$C_STOCKS,pnorm,mean(data2013_FH$C_STOCKS),sd(data2013_FH$C_STOCKS))

(HIST_CFH13=ggplot(data2013_FH, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 4, fill = "lightgray", color = "black") +
    geom_density(color = "blue", size = 0.5,position = "identity") +
    labs(x = Mgha_expression1, y = "Density", title = "a) C stocks")+
    annotate("text", x=23.1, y=0.09, label= "p = 0.394",family = "A") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T)),   # Ignore NA values for mean
               color="orange", linetype="dashed", size=0.8)+
    stat_function(fun = dnorm, args = list(mean = mean(data2013_FH$C_STOCKS), 
                                           sd = sd(data2013_FH$C_STOCKS)), color = "red",  size = 0.5)+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
          plot.subtitle = element_text(size = 11,face = "bold",hjust = 0),
          axis.text = element_text(size = 11))+
    labs(title = "2013"))


#C stocks - FH 18
normCFH18=seq(min(data2018_FH$C_STOCKS),max(data2018_FH$C_STOCKS),length=length(data2018_FH$C_STOCKS))
fun_normCFH18=dnorm(normCFH18,mean = mean(data2018_FH$C_STOCKS),sd=sd(data2018_FH$C_STOCKS))
hist(data2018_FH$C_STOCKS,probability = T,xlab = Mgha_expression,main = "2018 - FH")
lines(density(data2018_FH$C_STOCKS),col=4,lwd=2)
lines(normCFH18,fun_normCFH18,col=2,lwd=2)
ks.test(data2018_FH$C_STOCKS,pnorm,mean(data2018_FH$C_STOCKS),sd(data2018_FH$C_STOCKS))

(HIST_CFH18=ggplot(data2018_FH, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 1.5, fill = "lightgray", color = "black") +
    geom_density(color = "blue", size = 0.5,position = "identity") +
    labs(x = Mgha_expression1, y = "Density", title = "a) C stocks")+
    annotate("text", x=8.6, y=0.22, label= "p = 0.825",family = "A") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T)),   # Ignore NA values for mean
               color="orange", linetype="dashed", size=0.8)+
    stat_function(fun = dnorm, args = list(mean = mean(data2018_FH$C_STOCKS), 
                                           sd = sd(data2018_FH$C_STOCKS)), color = "red",  size = 0.5)+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
          plot.subtitle = element_text(size = 11,face = "bold",hjust = 0),
          axis.text = element_text(size = 11))+
    labs(title = "2018"))


#C stocks - FH 23
normCFH23=seq(min(data2023_FH$C_STOCKS),max(data2023_FH$C_STOCKS),length=length(data2023_FH$C_STOCKS))
fun_normCFH23=dnorm(normCFH23,mean = mean(data2023_FH$C_STOCKS),sd=sd(data2023_FH$C_STOCKS))
hist(data2023_FH$C_STOCKS,probability = T,xlab = Mgha_expression,main = "2023 - FH")
lines(density(data2023_FH$C_STOCKS),col=4,lwd=2)
lines(normCFH23,fun_normCFH23,col=2,lwd=2)
ks.test(data2023_FH$C_STOCKS,pnorm,mean(data2023_FH$C_STOCKS),sd(data2023_FH$C_STOCKS))

(HIST_CFH23=ggplot(data2023_FH, aes(x = C_STOCKS)) +
    geom_histogram(aes(y=..density..),alpha=0.3, size = 0.4, color = "darkgrey",binwidth = 2.5, fill = "lightgray", color = "black") +
    geom_density(color = "blue", size = 0.5,position = "identity") +
    labs(x = Mgha_expression1, y = "Density", title = "a) C stocks")+
    annotate("text", x=13, y=0.18, label= "p = 0.788",family = "A") +
    geom_vline(aes(xintercept=mean(C_STOCKS, na.rm=T)),   # Ignore NA values for mean
               color="orange", linetype="dashed", size=0.8)+
    stat_function(fun = dnorm, args = list(mean = mean(data2023_FH$C_STOCKS), 
                                           sd = sd(data2023_FH$C_STOCKS)), color = "red",  size = 0.5)+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
          plot.subtitle = element_text(size = 11,face = "bold",hjust = 0),
          axis.text = element_text(size = 11))+
    labs(title = "2023"))



jpeg(filename = here("outputs","Histograms - FH layer.jpeg"),width = 180,height = 70,units = "mm",res = 1000)
Supp2=ggarrange(HIST_CFH13,HIST_CFH18,HIST_CFH23,
                ncol = 3,nrow = 1,
                common.legend = TRUE,legend="bottom")
Supp2
dev.off()


# 5. SPEARMAN CORRELATIONS -----------------------------------------------------

##  2013
names(data2013)
correlations_2013=data.frame(data2013[c(1,4,5,6:17,19,20)])
names(correlations_2013)=c("YEAR","LAYER","C",
                           "UTMX","UTMY",
                           "AGE","BA","DTH","RICH","H'","GAP","COV",
                           "ELE","SLO","ASP","FFLAYER","MAN")
names(correlations_2013)
(corr13_spearman=rcorr(as.matrix(correlations_2013[c(3:17)]),type = "spearman"))

corrs13_output<- capture.output(corr13_spearman)
corrs13_output<- unlist(corrs13_output)
writeLines(corrs13_output,here("outputs","corrs13_Spearman.txt"))

##  2018
names(data2018)
correlations_2018=data.frame(data2018[c(1,4,5,6:17,19,20)])
names(correlations_2018)=c("YEAR","LAYER","C",
                           "UTMX","UTMY",
                           "AGE","BA","DTH","RICH","H'","GAP","COV",
                           "ELE","SLO","ASP","FFLAYER","MAN")
names(correlations_2018)
(corr18_spearman=rcorr(as.matrix(correlations_2018[c(3:17)]),type = "spearman"))

corrs18_output<- capture.output(corr18_spearman)
corrs18_output<- unlist(corrs18_output)
writeLines(corrs18_output,here("outputs","corrs18_Spearman.txt"))

##  2023
names(data2023)
correlations_2023=data.frame(data2023[c(1,4,5,6:17,19,20)])
names(correlations_2023)=c("YEAR","LAYER","C",
                           "UTMX","UTMY",
                           "AGE","BA","DTH","RICH","H'","GAP","COV",
                           "ELE","SLO","ASP","FFLAYER","MAN")
names(correlations_2023)
(corr23_spearman=rcorr(as.matrix(correlations_2023[c(3:17)]),type = "spearman"))

corrs23_output<- capture.output(corr23_spearman)
corrs23_output<- unlist(corrs23_output)
writeLines(corrs23_output,here("outputs","corrs23_Spearman.txt"))


# 6. VIF -----------------------------------------------------------------------

modelC_2013VIF = lm(correlations_2013$C ~ 
                      correlations_2013$UTMX + 
                      correlations_2013$UTMY + 
                      correlations_2013$AGE + 
                      correlations_2013$BA + 
                      correlations_2013$DTH + 
                      correlations_2013$RICH + 
                      correlations_2013$`H'` + 
                      correlations_2013$GAP + 
                      correlations_2013$ELE + 
                      correlations_2013$SLO + 
                      correlations_2013$ASP,
                    data = correlations_2013)

# Compute VIF

vif_values <- vif(modelC_2013VIF)
print(vif_values)

# Optional: Identify variables with VIF > 5
high_vif <- vif_values[vif_values > 5]
print(high_vif)
vif_df <- data.frame(Variable = names(vif_values), VIF = vif_values)


# 7. TRAIN - TEST DATA SPLIT ---------------------------------------------------
#model training  (85%) model validation (15%)
#NOT RUN
#Note: To run again RF models, just read the csv previously splitted. 
#      It will be used the same split for all models
set.seed(1234567)

## 7.1. 2013
datasetsize_2013<-floor(nrow(data2013)*0.85) #training data  
index_2013<-sample(1:nrow(data2013),size = datasetsize_2013)
training_2013<-data2013[index_2013,]
testing_2013<-data2013[-index_2013,]
View(training_2013)
View(testing_2013)
write.csv(training_2013,here("data", "training_2013_85.csv"),row.names = TRUE)
write.csv(testing_2013,here("data", "testing_2013_15.csv"),row.names = TRUE)

## 7.2. 2018
datasetsize_2018<-floor(nrow(data2018)*0.85) #training data  
index_2018<-sample(1:nrow(data2018),size = datasetsize_2018)
training_2018<-data2018[index_2018,]
testing_2018<-data2018[-index_2018,]
View(training_2018)
View(testing_2018)
write.csv(training_2018,here("data", "training_2018_85.csv"),row.names = TRUE)
write.csv(testing_2018,here("data", "testing_2018_15.csv"),row.names = TRUE)

## 7.3. 2023
datasetsize_2023<-floor(nrow(data2023)*0.85) #training data  
index_2023<-sample(1:nrow(data2023),size = datasetsize_2023)
training_2023<-data2023[index_2023,]
testing_2023<-data2023[-index_2023,]
View(training_2023)
View(testing_2023)
write.csv(training_2023,here("data", "training_2023_85.csv"),row.names = TRUE)
write.csv(testing_2023,here("data", "testing_2023_15.csv"),row.names = TRUE)


