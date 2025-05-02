# CHAPTER I - OPTIMIZING THE RANDOM FOREST ALGORITM FOR SPATIOTEMPORAL MODELING OF 
#             CARBON AND NITROGEN STOCKS IN FOREST FLOOR LAYERS
#_______________________________________________________________________________

#Ctrl+L 
options(max.print=100000000)
windowsFonts(A = windowsFont("Times New Roman"))
# 0. Libraries -------------------------------------------------------------------------------
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
library(coin)
library(emmeans)
library(plotrix)
library(openair)

#1. Expressions for plots ---------------------------------------------------
Mgha_expression=expression("Mg ha"^"-1")
initialmodel_expresion=expression(paste("RF"[0]," (Default values)"))
finalmodel_expression=expression(paste("RF"[1], " (Optimized values)"))
Predicted_Mgha_expression=expression(bold("Predicted C stocks (Mg ha"^"-1"~")"))
Observed_Mgha_expression=expression(bold("Observed C stocks (Mg ha"^"-1"~")"))
Predicted_kgha_expression=expression(bold("Predicted N stocks (kg ha"^"-1"~")"))
Observed_kgha_expression=expression(bold("Observed N stocks (kg ha"^"-1"~")"))
Observed_ratio_expresion=expression(bold("Observed C/N ratio"))
Predicted_ratio_expresion=expression(bold("Predicted C/N ratio"))
Predictedvalues_RF1=expression(bold("Predicted values (RF"[1]~")"))
Evaluatedvariables_RF1=expression(bold(bold("Evaluated variables in RF"[1]~"models:")))
Observedvalues_expresion=expression(bold("Observed values"))
RMSE_expresion=expression(bold("Cross-validated RMSE"))
mtry_expresion=expression(bold("mtry values"))
ntree_expresion=expression(bold("ntree values"))
CVreps_expresion=expression(bold("10-fold cross-validation repetitions"))
kgha_expression=expression("kg ha"^"-1")
predictors_expresion=expression(bold("Predictor variables"))
IncMSE_expresion=expression(bold("%IncMSE"))
normalizedsd=expression(bold("Standard deviation (normalized)"))
#Mgha_expression1=expression("C stocks (Mg ha"^"-1"~")")
#Mgha_expression2=expression(bold("a) C models (Mg ha"^"-1"~")"))
#RMSE_Mgha_expression=expression("RMSE (Mg ha"^"-1"~")")

#------------------------------------CHAPTER 1----------------------------------
# 3. MTRY SELECTION IMPACT  -------------------------------------------------------------------------------
## 3.1. C 2013 ----------------------------------------------------------------
### 3.1.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_C2013_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_C2013_CV5rep.csv"),header = TRUE)
Hyperparameters_C2013_CV_15reps=read.csv(here("Hyperparameters","best_hyperpars_C2013_CV15rep.csv"),header = TRUE)

Hyperparameters_C2013_CV_5reps=Hyperparameters_C2013_CV_5reps[-1]
Hyperparameters_C2013_CV_15reps=Hyperparameters_C2013_CV_15reps[-1]

Hyperparameters_C2013_CV_5reps$ntree=as.factor(Hyperparameters_C2013_CV_5reps$ntree)
Hyperparameters_C2013_CV_15reps$ntree=as.factor(Hyperparameters_C2013_CV_15reps$ntree)

Hyperparameters_C2013_CV_5reps=data.frame(Variable="C stocks",Year=2013,
                                          CVrepetitions="5 reps",Hyperparameters_C2013_CV_5reps,
                                          ntreeselection="Descarted")
Hyperparameters_C2013_CV_15reps=data.frame(Variable="C stocks",Year=2013,
                                           CVrepetitions="15 reps",Hyperparameters_C2013_CV_15reps,
                                           ntreeselection="Descarted")

### 3.1.2. Default models:  -------------------------------
Hyperparameters_C2013_CV_5reps$ntreeselection[Hyperparameters_C2013_CV_5reps$ntree == 500]<- "RF0"
Hyperparameters_C2013_CV_15reps$ntreeselection[Hyperparameters_C2013_CV_15reps$ntree == 500]<- "RF0"
### 3.1.3. Optimized models ---------------------------------
Hyperparameters_C2013_CV_5reps$ntreeselection[Hyperparameters_C2013_CV_5reps$ntree == 900]<- "RF1"
Hyperparameters_C2013_CV_15reps$ntreeselection[Hyperparameters_C2013_CV_15reps$ntree == 600]<- "RF1"

View(Hyperparameters_C2013_CV_5reps)
View(Hyperparameters_C2013_CV_15reps)

### 3.1.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_C2013_CV_5reps=rbind(Hyperparameters_C2013_CV_5reps[Hyperparameters_C2013_CV_5reps$ntreeselection== "RF0",],
                                     Hyperparameters_C2013_CV_5reps[Hyperparameters_C2013_CV_5reps$ntreeselection=="RF1",])
selected_ntrees_C2013_CV_15reps=rbind(Hyperparameters_C2013_CV_15reps[Hyperparameters_C2013_CV_15reps$ntreeselection== "RF0",],
                                      Hyperparameters_C2013_CV_15reps[Hyperparameters_C2013_CV_15reps$ntreeselection=="RF1",])

### 3.1.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 5
#      RF1: ntree = 600;  mtry = 12;   CVRep = 15

selected_ntrees_C2013_CV_PLOT=rbind(selected_ntrees_C2013_CV_5reps[selected_ntrees_C2013_CV_5reps$ntreeselection== "RF0",],
                                    selected_ntrees_C2013_CV_15reps[selected_ntrees_C2013_CV_15reps$ntreeselection=="RF1",])
View(selected_ntrees_C2013_CV_PLOT)

# Crear el gráfico ajustado
MTRY_C2013 <- ggplot(selected_ntrees_C2013_CV_PLOT, aes(x=mtry, y=RMSE, color=ntreeselection)) +
  stat_summary(fun=mean, geom="line", alpha=0.4, size=0.7) +
  stat_summary(fun=mean, geom="point", size=0.4, pch=16) +
  geom_hline(yintercept = 3.262240, linetype = "dashed", color = "black", size=0.4) + # Umbral del modelo predeterminado
  annotate("point", y=3.262240, x=4, color="#C80036", size=1.8) +  # Punto localizado de los hiperparámetros - MODELO PREDETERMINADO 
  annotate("point", y=3.122904, x=12, color="#050C9C", size=1.8) + # Punto localizado de los hiperparámetros - MODELO OPTIMIZADO
  ylab(Mgha_expression) +
  xlab("") +
  scale_color_manual(values = c("#C80036", "#050C9C"),
                     name="Random forest models:",
                     labels=c(initialmodel_expresion, finalmodel_expression)) +
  theme_bw() +
  ggtitle("2013", subtitle = "a) C stocks") +
  scale_x_continuous(limits = c(1, 12), breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(3.0, 4.5), breaks = c(3.30, 3.60, 3.90, 4.20, 4.5)) +
  theme(
    axis.text.x = element_text(vjust=0.6, family = "A", size = 11),
    axis.text.y = element_text(vjust=0.6, family = "A", size = 11),
    text = element_text(size=11, family = "A"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 13, hjust = 0, face = "bold"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 13, face = "bold"),
    plot.margin = unit(c(0, 0.1, 0, 0), "cm")
  )

# Mostrar el gráfico
print(MTRY_C2013)



## 3.2. C 2018 ----------------------------------------------------------------
### 3.2.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_C2018_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_C2018_CV5rep.csv"),header = TRUE)
Hyperparameters_C2018_CV_15reps=read.csv(here("Hyperparameters","best_hyperpars_C2018_CV15rep.csv"),header = TRUE)

Hyperparameters_C2018_CV_5reps=Hyperparameters_C2018_CV_5reps[-1]
Hyperparameters_C2018_CV_15reps=Hyperparameters_C2018_CV_15reps[-1]

Hyperparameters_C2018_CV_5reps$ntree=as.factor(Hyperparameters_C2018_CV_5reps$ntree)
Hyperparameters_C2018_CV_15reps$ntree=as.factor(Hyperparameters_C2018_CV_15reps$ntree)

Hyperparameters_C2018_CV_5reps=data.frame(Variable="C stocks",Year=2018,
                                          CVrepetitions="5 reps",Hyperparameters_C2018_CV_5reps,
                                          ntreeselection="Descarted")
Hyperparameters_C2018_CV_15reps=data.frame(Variable="C stocks",Year=2018,
                                           CVrepetitions="15 reps",Hyperparameters_C2018_CV_15reps,
                                           ntreeselection="Descarted")

### 3.2.2. Default models:  -------------------------------
Hyperparameters_C2018_CV_5reps$ntreeselection[Hyperparameters_C2018_CV_5reps$ntree == 500]<- "RF0"
Hyperparameters_C2018_CV_15reps$ntreeselection[Hyperparameters_C2018_CV_15reps$ntree == 500]<- "RF0"
### 3.2.3. Optimized models ---------------------------------
Hyperparameters_C2018_CV_5reps$ntreeselection[Hyperparameters_C2018_CV_5reps$ntree == 1000]<- "RF1"
Hyperparameters_C2018_CV_15reps$ntreeselection[Hyperparameters_C2018_CV_15reps$ntree == 800]<- "RF1"

View(Hyperparameters_C2018_CV_5reps)
View(Hyperparameters_C2018_CV_15reps)

### 3.2.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_C2018_CV_5reps=rbind(Hyperparameters_C2018_CV_5reps[Hyperparameters_C2018_CV_5reps$ntreeselection== "RF0",],
                                     Hyperparameters_C2018_CV_5reps[Hyperparameters_C2018_CV_5reps$ntreeselection=="RF1",])
selected_ntrees_C2018_CV_15reps=rbind(Hyperparameters_C2018_CV_15reps[Hyperparameters_C2018_CV_15reps$ntreeselection== "RF0",],
                                      Hyperparameters_C2018_CV_15reps[Hyperparameters_C2018_CV_15reps$ntreeselection=="RF1",])

### 3.2.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 15
#      RF1: ntree = 1000;  mtry = 5;   CVRep = 5

selected_ntrees_C2018_CV_PLOT=rbind(selected_ntrees_C2018_CV_15reps[selected_ntrees_C2018_CV_15reps$ntreeselection== "RF0",],
                                    selected_ntrees_C2018_CV_5reps[selected_ntrees_C2018_CV_5reps$ntreeselection=="RF1",])
View(selected_ntrees_C2018_CV_PLOT)

(MTRY_C2018=ggplot(selected_ntrees_C2018_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 1.388339 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=1.388339,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=1.368969,x=5,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("2018",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(1.36,1.57),breaks = c(1.40,1.45,1.50,1.55))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position="bottom",
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))


## 3.3. C 2023 ----------------------------------------------------------------
### 3.3.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_C2023_CV_20reps=read.csv(here("Hyperparameters","best_hyperpars_C2023_CV20rep.csv"),header = TRUE)
Hyperparameters_C2023_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_C2023_CV5rep.csv"),header = TRUE)

Hyperparameters_C2023_CV_20reps=Hyperparameters_C2023_CV_20reps[-1]
Hyperparameters_C2023_CV_5reps=Hyperparameters_C2023_CV_5reps[-1]

Hyperparameters_C2023_CV_20reps$ntree=as.factor(Hyperparameters_C2023_CV_20reps$ntree)
Hyperparameters_C2023_CV_5reps$ntree=as.factor(Hyperparameters_C2023_CV_5reps$ntree)

Hyperparameters_C2023_CV_20reps=data.frame(Variable="C stocks",Year=2023,
                                          CVrepetitions="20 reps",Hyperparameters_C2023_CV_20reps,
                                          ntreeselection="Descarted")
Hyperparameters_C2023_CV_5reps=data.frame(Variable="C stocks",Year=2023,
                                           CVrepetitions="5 reps",Hyperparameters_C2023_CV_5reps,
                                           ntreeselection="Descarted")

### 3.3.2. Default models:  -------------------------------
Hyperparameters_C2023_CV_20reps$ntreeselection[Hyperparameters_C2023_CV_20reps$ntree == 500]<- "RF0"

### 3.3.3. Optimized models ---------------------------------
Hyperparameters_C2023_CV_5reps$ntreeselection[Hyperparameters_C2023_CV_5reps$ntree == 500]<- "RF1"

View(Hyperparameters_C2023_CV_5reps)
View(Hyperparameters_C2023_CV_20reps)

### 3.3.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_C2023_CV_20reps=rbind(Hyperparameters_C2023_CV_20reps[Hyperparameters_C2023_CV_20reps$ntreeselection== "RF0",],
                                     Hyperparameters_C2023_CV_20reps[Hyperparameters_C2023_CV_20reps$ntreeselection=="RF1",])
selected_ntrees_C2023_CV_5reps=rbind(Hyperparameters_C2023_CV_5reps[Hyperparameters_C2023_CV_5reps$ntreeselection== "RF0",],
                                      Hyperparameters_C2023_CV_5reps[Hyperparameters_C2023_CV_5reps$ntreeselection=="RF1",])

### 3.3.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 20
#      RF1: ntree = 500;  mtry = 5;   CVRep = 5

selected_ntrees_C2023_CV_PLOT=rbind(selected_ntrees_C2023_CV_20reps[selected_ntrees_C2023_CV_20reps$ntreeselection== "RF0",],
                                    selected_ntrees_C2023_CV_5reps[selected_ntrees_C2023_CV_5reps$ntreeselection=="RF1",])
View(selected_ntrees_C2023_CV_PLOT)

(MTRY_C2023=ggplot(selected_ntrees_C2023_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 1.486780 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=1.486780,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=1.470789,x=5,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("2023",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(1.45,1.9),breaks = c(1.50,1.60,1.70,1.80,1.9))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))

## 3.4. N 2013 ----------------------------------------------------------------
### 3.4.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_N2013_CV_10reps=read.csv(here("Hyperparameters","best_hyperpars_N2013_CV10rep.csv"),header = TRUE)
Hyperparameters_N2013_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_N2013_CV5rep.csv"),header = TRUE)

Hyperparameters_N2013_CV_10reps=Hyperparameters_N2013_CV_10reps[-1]
Hyperparameters_N2013_CV_5reps=Hyperparameters_N2013_CV_5reps[-1]

Hyperparameters_N2013_CV_10reps$ntree=as.factor(Hyperparameters_N2013_CV_10reps$ntree)
Hyperparameters_N2013_CV_5reps$ntree=as.factor(Hyperparameters_N2013_CV_5reps$ntree)

Hyperparameters_N2013_CV_10reps=data.frame(Variable="N stocks",Year=2013,
                                          CVrepetitions="10 reps",Hyperparameters_N2013_CV_10reps,
                                          ntreeselection="Descarted")
Hyperparameters_N2013_CV_5reps=data.frame(Variable="N stocks",Year=2013,
                                           CVrepetitions="5 reps",Hyperparameters_N2013_CV_5reps,
                                           ntreeselection="Descarted")

### 3.4.2. Default models:  -------------------------------
Hyperparameters_N2013_CV_10reps$ntreeselection[Hyperparameters_N2013_CV_10reps$ntree == 500]<- "RF0"
### 3.4.3. Optimized models ---------------------------------
Hyperparameters_N2013_CV_5reps$ntreeselection[Hyperparameters_N2013_CV_5reps$ntree == 800]<- "RF1"

View(Hyperparameters_N2013_CV_10reps)
View(Hyperparameters_N2013_CV_5reps)

### 3.4.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_N2013_CV_10reps=rbind(Hyperparameters_N2013_CV_10reps[Hyperparameters_N2013_CV_10reps$ntreeselection== "RF0",],
                                     Hyperparameters_N2013_CV_10reps[Hyperparameters_N2013_CV_10reps$ntreeselection=="RF1",])
selected_ntrees_N2013_CV_5reps=rbind(Hyperparameters_N2013_CV_5reps[Hyperparameters_N2013_CV_5reps$ntreeselection== "RF0",],
                                      Hyperparameters_N2013_CV_5reps[Hyperparameters_N2013_CV_5reps$ntreeselection=="RF1",])

### 3.4.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 10
#      RF1: ntree = 800;  mtry = 8;   CVRep = 5

selected_ntrees_N2013_CV_PLOT=rbind(selected_ntrees_N2013_CV_10reps[selected_ntrees_N2013_CV_10reps$ntreeselection== "RF0",],
                                    selected_ntrees_N2013_CV_5reps[selected_ntrees_N2013_CV_5reps$ntreeselection=="RF1",])
View(selected_ntrees_N2013_CV_PLOT)

(MTRY_N2013=ggplot(selected_ntrees_N2013_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 62.02852 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=62.02852,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=60.22298,x=8,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab(kgha_expression)+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "b) N stocks")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(60,80),breaks = c(60,65,70,75,80))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))


## 3.5. N 2018 ----------------------------------------------------------------
### 3.5.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_N2018_CV_10reps=read.csv(here("Hyperparameters","best_hyperpars_N2018_CV10rep.csv"),header = TRUE)
Hyperparameters_N2018_CV_15reps=read.csv(here("Hyperparameters","best_hyperpars_N2018_CV15rep.csv"),header = TRUE)

Hyperparameters_N2018_CV_10reps=Hyperparameters_N2018_CV_10reps[-1]
Hyperparameters_N2018_CV_15reps=Hyperparameters_N2018_CV_15reps[-1]

Hyperparameters_N2018_CV_10reps$ntree=as.factor(Hyperparameters_N2018_CV_10reps$ntree)
Hyperparameters_N2018_CV_15reps$ntree=as.factor(Hyperparameters_N2018_CV_15reps$ntree)

Hyperparameters_N2018_CV_10reps=data.frame(Variable="N stocks",Year=2013,
                                           CVrepetitions="10 reps",Hyperparameters_N2018_CV_10reps,
                                           ntreeselection="Descarted")
Hyperparameters_N2018_CV_15reps=data.frame(Variable="N stocks",Year=2013,
                                          CVrepetitions="15 reps",Hyperparameters_N2018_CV_15reps,
                                          ntreeselection="Descarted")

### 3.5.2. Default models:  -------------------------------
Hyperparameters_N2018_CV_10reps$ntreeselection[Hyperparameters_N2018_CV_10reps$ntree == 500]<- "RF0"
### 3.5.3. Optimized models ---------------------------------
Hyperparameters_N2018_CV_15reps$ntreeselection[Hyperparameters_N2018_CV_15reps$ntree == 500]<- "RF1"

View(Hyperparameters_N2018_CV_10reps)
View(Hyperparameters_N2018_CV_15reps)

### 3.5.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_N2018_CV_10reps=rbind(Hyperparameters_N2018_CV_10reps[Hyperparameters_N2018_CV_10reps$ntreeselection== "RF0",],
                                      Hyperparameters_N2018_CV_10reps[Hyperparameters_N2018_CV_10reps$ntreeselection=="RF1",])
selected_ntrees_N2018_CV_15reps=rbind(Hyperparameters_N2018_CV_15reps[Hyperparameters_N2018_CV_15reps$ntreeselection== "RF0",],
                                     Hyperparameters_N2018_CV_15reps[Hyperparameters_N2018_CV_15reps$ntreeselection=="RF1",])

### 3.5.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 10
#      RF1: ntree = 500;  mtry = 6;   CVRep = 15

selected_ntrees_N2018_CV_PLOT=rbind(selected_ntrees_N2018_CV_10reps[selected_ntrees_N2018_CV_10reps$ntreeselection== "RF0",],
                                    selected_ntrees_N2018_CV_15reps[selected_ntrees_N2018_CV_15reps$ntreeselection=="RF1",])
View(selected_ntrees_N2018_CV_PLOT)

(MTRY_N2018=ggplot(selected_ntrees_N2018_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 39.64387 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=39.64387,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=39.21273,x=6,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(39,47),breaks = c(39,41,43,45,47))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))

## 3.6. N 2023 ----------------------------------------------------------------
### 3.6.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_N2023_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_N2023_CV5rep.csv"),header = TRUE)
Hyperparameters_N2023_CV_10reps=read.csv(here("Hyperparameters","best_hyperpars_N2023_CV10rep.csv"),header = TRUE)

Hyperparameters_N2023_CV_5reps=Hyperparameters_N2023_CV_5reps[-1]
Hyperparameters_N2023_CV_10reps=Hyperparameters_N2023_CV_10reps[-1]

Hyperparameters_N2023_CV_5reps$ntree=as.factor(Hyperparameters_N2023_CV_5reps$ntree)
Hyperparameters_N2023_CV_10reps$ntree=as.factor(Hyperparameters_N2023_CV_10reps$ntree)

Hyperparameters_N2023_CV_5reps=data.frame(Variable="N stocks",Year=2023,
                                           CVrepetitions="5 reps",Hyperparameters_N2023_CV_5reps,
                                           ntreeselection="Descarted")
Hyperparameters_N2023_CV_10reps=data.frame(Variable="N stocks",Year=2023,
                                           CVrepetitions="10 reps",Hyperparameters_N2023_CV_10reps,
                                           ntreeselection="Descarted")

### 3.6.2. Default models:  -------------------------------
Hyperparameters_N2023_CV_5reps$ntreeselection[Hyperparameters_N2023_CV_5reps$ntree == 500]<- "RF0"
### 3.6.3. Optimized models ---------------------------------
Hyperparameters_N2023_CV_10reps$ntreeselection[Hyperparameters_N2023_CV_10reps$ntree == 900]<- "RF1"

View(Hyperparameters_N2023_CV_5reps)
View(Hyperparameters_N2023_CV_10reps)

### 3.6.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_N2023_CV_5reps=rbind(Hyperparameters_N2023_CV_5reps[Hyperparameters_N2023_CV_5reps$ntreeselection== "RF0",],
                                      Hyperparameters_N2023_CV_5reps[Hyperparameters_N2023_CV_5reps$ntreeselection=="RF1",])
selected_ntrees_N2023_CV_10reps=rbind(Hyperparameters_N2023_CV_10reps[Hyperparameters_N2023_CV_10reps$ntreeselection== "RF0",],
                                      Hyperparameters_N2023_CV_10reps[Hyperparameters_N2023_CV_10reps$ntreeselection=="RF1",])

### 3.6.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 5
#      RF1: ntree = 900;  mtry = 11;   CVRep = 10

selected_ntrees_N2023_CV_PLOT=rbind(selected_ntrees_N2023_CV_5reps[selected_ntrees_N2023_CV_5reps$ntreeselection== "RF0",],
                                    selected_ntrees_N2023_CV_10reps[selected_ntrees_N2023_CV_10reps$ntreeselection=="RF1",])
View(selected_ntrees_N2023_CV_PLOT)

(MTRY_N2023=ggplot(selected_ntrees_N2023_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 36.01659 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=36.01659,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=35.22512,x=11,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(35,44),breaks = c(36,38,40,42,44))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))


## 3.7. CN RATIO 2013 ----------------------------------------------------------------
### 3.7.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_CNRATIO2013_CV_15reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2013_CV15rep.csv"),header = TRUE)
Hyperparameters_CNRATIO2013_CV_10reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2013_CV10rep.csv"),header = TRUE)

Hyperparameters_CNRATIO2013_CV_15reps=Hyperparameters_CNRATIO2013_CV_15reps[-1]
Hyperparameters_CNRATIO2013_CV_10reps=Hyperparameters_CNRATIO2013_CV_10reps[-1]

Hyperparameters_CNRATIO2013_CV_15reps$ntree=as.factor(Hyperparameters_CNRATIO2013_CV_15reps$ntree)
Hyperparameters_CNRATIO2013_CV_10reps$ntree=as.factor(Hyperparameters_CNRATIO2013_CV_10reps$ntree)

Hyperparameters_CNRATIO2013_CV_15reps=data.frame(Variable="C/N ratio",Year=2013,
                                           CVrepetitions="15 reps",Hyperparameters_CNRATIO2013_CV_15reps,
                                           ntreeselection="Descarted")
Hyperparameters_CNRATIO2013_CV_10reps=data.frame(Variable="C/N ratio",Year=2013,
                                          CVrepetitions="10 reps",Hyperparameters_CNRATIO2013_CV_10reps,
                                          ntreeselection="Descarted")

### 3.7.2. Default models:  -------------------------------
Hyperparameters_CNRATIO2013_CV_15reps$ntreeselection[Hyperparameters_CNRATIO2013_CV_15reps$ntree == 500]<- "RF0"
### 3.7.3. Optimized models ---------------------------------
Hyperparameters_CNRATIO2013_CV_10reps$ntreeselection[Hyperparameters_CNRATIO2013_CV_10reps$ntree == 800]<- "RF1"

View(Hyperparameters_CNRATIO2013_CV_15reps)
View(Hyperparameters_CNRATIO2013_CV_10reps)

### 3.7.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_CNRATIO2013_CV_15reps=rbind(Hyperparameters_CNRATIO2013_CV_15reps[Hyperparameters_CNRATIO2013_CV_15reps$ntreeselection== "RF0",],
                                      Hyperparameters_CNRATIO2013_CV_15reps[Hyperparameters_CNRATIO2013_CV_15reps$ntreeselection=="RF1",])
selected_ntrees_CNRATIO2013_CV_10reps=rbind(Hyperparameters_CNRATIO2013_CV_10reps[Hyperparameters_CNRATIO2013_CV_10reps$ntreeselection== "RF0",],
                                     Hyperparameters_CNRATIO2013_CV_10reps[Hyperparameters_CNRATIO2013_CV_10reps$ntreeselection=="RF1",])

### 3.7.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 15
#      RF1: ntree = 800;  mtry = 8;   CVRep = 10

selected_ntrees_CNRATIO2013_CV_PLOT=rbind(selected_ntrees_CNRATIO2013_CV_15reps[selected_ntrees_CNRATIO2013_CV_15reps$ntreeselection== "RF0",],
                                    selected_ntrees_CNRATIO2013_CV_10reps[selected_ntrees_CNRATIO2013_CV_10reps$ntreeselection=="RF1",])
View(selected_ntrees_CNRATIO2013_CV_PLOT)

(MTRY_CNRATIO2013=ggplot(selected_ntrees_CNRATIO2013_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 6.277482 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=6.277482,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=5.915066,x=8,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("ratio")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "c) C/N ratio")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(5.8,10),breaks = c(6,7,8,9,10))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))

## 3.8. CN RATIO 2018 ----------------------------------------------------------------
### 3.8.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_CNRATIO2018_CV_20reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2018_CV20rep.csv"),header = TRUE)
Hyperparameters_CNRATIO2018_CV_5reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2018_CV5rep.csv"),header = TRUE)

Hyperparameters_CNRATIO2018_CV_20reps=Hyperparameters_CNRATIO2018_CV_20reps[-1]
Hyperparameters_CNRATIO2018_CV_5reps=Hyperparameters_CNRATIO2018_CV_5reps[-1]

Hyperparameters_CNRATIO2018_CV_20reps$ntree=as.factor(Hyperparameters_CNRATIO2018_CV_20reps$ntree)
Hyperparameters_CNRATIO2018_CV_5reps$ntree=as.factor(Hyperparameters_CNRATIO2018_CV_5reps$ntree)

Hyperparameters_CNRATIO2018_CV_20reps=data.frame(Variable="C/N ratio",Year=2018,
                                                 CVrepetitions="20 reps",Hyperparameters_CNRATIO2018_CV_20reps,
                                                 ntreeselection="Descarted")
Hyperparameters_CNRATIO2018_CV_5reps=data.frame(Variable="C/N ratio",Year=2018,
                                                 CVrepetitions="5 reps",Hyperparameters_CNRATIO2018_CV_5reps,
                                                 ntreeselection="Descarted")

### 3.8.2. Default models:  -------------------------------
Hyperparameters_CNRATIO2018_CV_20reps$ntreeselection[Hyperparameters_CNRATIO2018_CV_20reps$ntree == 500]<- "RF0"
### 3.8.3. Optimized models ---------------------------------
Hyperparameters_CNRATIO2018_CV_5reps$ntreeselection[Hyperparameters_CNRATIO2018_CV_5reps$ntree == 500]<- "RF1"

View(Hyperparameters_CNRATIO2018_CV_20reps)
View(Hyperparameters_CNRATIO2018_CV_5reps)

### 3.8.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_CNRATIO2018_CV_20reps=rbind(Hyperparameters_CNRATIO2018_CV_20reps[Hyperparameters_CNRATIO2018_CV_20reps$ntreeselection== "RF0",],
                                            Hyperparameters_CNRATIO2018_CV_20reps[Hyperparameters_CNRATIO2018_CV_20reps$ntreeselection=="RF1",])
selected_ntrees_CNRATIO2018_CV_5reps=rbind(Hyperparameters_CNRATIO2018_CV_5reps[Hyperparameters_CNRATIO2018_CV_5reps$ntreeselection== "RF0",],
                                            Hyperparameters_CNRATIO2018_CV_5reps[Hyperparameters_CNRATIO2018_CV_5reps$ntreeselection=="RF1",])

### 3.8.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 20
#      RF1: ntree = 500;  mtry = 7;   CVRep = 5

selected_ntrees_CNRATIO2018_CV_PLOT=rbind(selected_ntrees_CNRATIO2018_CV_20reps[selected_ntrees_CNRATIO2018_CV_20reps$ntreeselection== "RF0",],
                                          selected_ntrees_CNRATIO2018_CV_5reps[selected_ntrees_CNRATIO2018_CV_5reps$ntreeselection=="RF1",])
View(selected_ntrees_CNRATIO2018_CV_PLOT)

(MTRY_CNRATIO2018=ggplot(selected_ntrees_CNRATIO2018_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 3.955996 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=3.955996,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=3.815766,x=7,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(3.8,5.6),breaks = c(4.0,4.5,5.0,5.5))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))


## 3.9. CN RATIO 2023 ----------------------------------------------------------------
### 3.9.1. Selected dataframes for plot mtry selection -----------------------

Hyperparameters_CNRATIO2023_CV_10reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2023_CV10rep.csv"),header = TRUE)
Hyperparameters_CNRATIO2023_CV_20reps=read.csv(here("Hyperparameters","best_hyperpars_CNratio2023_CV20rep.csv"),header = TRUE)

Hyperparameters_CNRATIO2023_CV_10reps=Hyperparameters_CNRATIO2023_CV_10reps[-1]
Hyperparameters_CNRATIO2023_CV_20reps=Hyperparameters_CNRATIO2023_CV_20reps[-1]

Hyperparameters_CNRATIO2023_CV_10reps$ntree=as.factor(Hyperparameters_CNRATIO2023_CV_10reps$ntree)
Hyperparameters_CNRATIO2023_CV_20reps$ntree=as.factor(Hyperparameters_CNRATIO2023_CV_20reps$ntree)

Hyperparameters_CNRATIO2023_CV_10reps=data.frame(Variable="C/N ratio",Year=2023,
                                                 CVrepetitions="10 reps",Hyperparameters_CNRATIO2023_CV_10reps,
                                                 ntreeselection="Descarted")
Hyperparameters_CNRATIO2023_CV_20reps=data.frame(Variable="C/N ratio",Year=2023,
                                                 CVrepetitions="20 reps",Hyperparameters_CNRATIO2023_CV_20reps,
                                                 ntreeselection="Descarted")

### 3.9.2. Default models:  -------------------------------
Hyperparameters_CNRATIO2023_CV_10reps$ntreeselection[Hyperparameters_CNRATIO2023_CV_10reps$ntree == 500]<- "RF0"
### 3.9.3. Optimized models ---------------------------------
Hyperparameters_CNRATIO2023_CV_20reps$ntreeselection[Hyperparameters_CNRATIO2023_CV_20reps$ntree == 500]<- "RF1"

View(Hyperparameters_CNRATIO2023_CV_10reps)
View(Hyperparameters_CNRATIO2023_CV_20reps)

### 3.9.4. Filter RF0 AND RF1 ----------------------------------------
selected_ntrees_CNRATIO2023_CV_10reps=rbind(Hyperparameters_CNRATIO2023_CV_10reps[Hyperparameters_CNRATIO2023_CV_10reps$ntreeselection== "RF0",],
                                            Hyperparameters_CNRATIO2023_CV_10reps[Hyperparameters_CNRATIO2023_CV_10reps$ntreeselection=="RF1",])
selected_ntrees_CNRATIO2023_CV_20reps=rbind(Hyperparameters_CNRATIO2023_CV_20reps[Hyperparameters_CNRATIO2023_CV_20reps$ntreeselection== "RF0",],
                                            Hyperparameters_CNRATIO2023_CV_20reps[Hyperparameters_CNRATIO2023_CV_20reps$ntreeselection=="RF1",])

### 3.9.5  PLOT  ----------------------

#Note: RF0: ntree = 500;  mtry = 4    CVRep = 10
#      RF1: ntree = 500;  mtry = 7;   CVRep = 20

selected_ntrees_CNRATIO2023_CV_PLOT=rbind(selected_ntrees_CNRATIO2023_CV_10reps[selected_ntrees_CNRATIO2023_CV_10reps$ntreeselection== "RF0",],
                                          selected_ntrees_CNRATIO2023_CV_20reps[selected_ntrees_CNRATIO2023_CV_20reps$ntreeselection=="RF1",])
View(selected_ntrees_CNRATIO2023_CV_PLOT)

(MTRY_CNRATIO2023=ggplot(selected_ntrees_CNRATIO2023_CV_PLOT,aes(x=mtry,y=RMSE,color=ntreeselection))+
    stat_summary(fun=mean,geom="line",alpha=0.4,size=0.7)+
    stat_summary(fun=mean,geom="point",size=0.4,pch=16)+
    geom_hline(yintercept = 6.191687 , linetype = "dashed", color = "black",size=0.4) + #Threshold of DEFAULT MODEL
    annotate("point",y=6.191687,x=4,color="#C80036",size=1.8)+  # Located point of hyperparameters - DEFAULT MODEL 
    annotate("point",y=5.933402,x=7,color="#050C9C",size=1.8)+ # Located point of hyperparameters - OPTIMIZED MODEL
    ylab("")+
    xlab("")+
    scale_color_manual(values = c("#C80036", "#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression))+
    theme_bw()+
    ggtitle("",subtitle = "")+
    scale_x_continuous(limits = c(1, 12),breaks = c(2,4,6,8,10,12) )+
    scale_y_continuous(limits = c(5.8,10),breaks = c(6,7,8,9,10))+
    theme(axis.text.x = element_text(vjust=0.6,family = "A",size = 11),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 11),
          text = element_text(size=11,family = "A"),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.position="bottom",
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size=13),
          legend.title = element_text(size = 13, face = "bold")))
## 3.10 Figures ---------------------------------------------------------------

jpeg(filename = here("Final figures", "1. mtrys_fullcolumnsize1000.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
mtry_error_plots=ggarrange(MTRY_C2013,MTRY_C2018,MTRY_C2023,
                           MTRY_N2013,MTRY_N2018,MTRY_N2023,
                           MTRY_CNRATIO2013,MTRY_CNRATIO2018,MTRY_CNRATIO2023,
                            ncol = 3,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
mtry_error_plots
annotate_figure(mtry_error_plots,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.2, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(mtry_expresion, vjust = -3.0,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()

# Shorter and edited in Paint
jpeg(filename = here("Final figures", "1.1 mtrys_fullcolumnsize1000_190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
mtry_error_plots=ggarrange(MTRY_C2013,MTRY_C2018,MTRY_C2023,
                           MTRY_N2013,MTRY_N2018,MTRY_N2023,
                           MTRY_CNRATIO2013,MTRY_CNRATIO2018,MTRY_CNRATIO2023,
                           ncol = 3,nrow = 3,
                           common.legend = TRUE,
                           legend = "bottom")
mtry_error_plots
annotate_figure(mtry_error_plots,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.2, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(mtry_expresion, vjust = -3,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()




#-------------------------------------------------------------------------------
# 4. CV IMPACT (TRAIN FUNCTION - ORIGINAL) -------------------------------------
## 4.1. Dataframes  C stocks ----------------------------------------------------------------
#C 2013
selected_models_CVREPS_2013=read.csv(here("Hyperparameters","selected_models_CVREPS_2013.csv"),header = TRUE)
selected_models_CVREPS_2013=selected_models_CVREPS_2013[-1]

#C 2018 
selected_models_CVREPS_2018=read.csv(here("Hyperparameters","selected_models_CVREPS_2018.csv"),header = TRUE)
selected_models_CVREPS_2018=selected_models_CVREPS_2018[-1]

#C 2023 
selected_models_CVREPS_2023=read.csv(here("Hyperparameters","selected_models_CVREPS_2023.csv"),header = TRUE)
selected_models_CVREPS_2023=selected_models_CVREPS_2023[-1]

# Combine all models
CVreps_all=rbind(selected_models_CVREPS_2013,
                 selected_models_CVREPS_2018,
                 selected_models_CVREPS_2023)

CVreps_all$year=as.factor(CVreps_all$year)
CVreps_all$variable=as.factor(CVreps_all$variable)


CVREPETS_Cstocks <- ggplot(CVreps_all, aes(x = CVreps, y = RMSE, color = year, linetype = model)) +
  geom_line(size = 0.4, alpha = 0.8) +
  geom_point(size = 0.8) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 50, 75, 100)) +
  ylab(Mgha_expression) +
  xlab("") +
  theme_bw() +
  ggtitle("a) C stocks") +
  scale_color_manual(values = c("#C80036", "#050C9C", "orange"),
                     name = "Year:",
                     labels = c("2013", "2018", "2023")) +
  scale_linetype_manual(values = c("longdash", "solid"),
                        name = "Random forest models:\n\n",
                        labels = c(initialmodel_expresion, finalmodel_expression)) +
  guides(color = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 2)) +  # Ajusta el número de filas según sea necesario
  labs(linetype = "") +  # Elimina la etiqueta de la leyenda para "linetype"
  # facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(limits = c(1.3, 3.5), breaks = c(1.5, 2, 2.5, 3, 3.5)) +
  theme(axis.text.x = element_text(vjust = 0.6, family = "A", size = 10),
        axis.text.y = element_text(vjust = 0.6, family = "A", size = 10),
        text = element_text(size = 10, family = "A"),
        plot.title = element_text(size = 14, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 13, hjust = 0, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11, face = "bold"))
CVREPETS_Cstocks


## 4.2. Dataframes  N stocks ----------------------------------------------------------------
#N 2013
selected_models_CVREPS_2013N=read.csv(here("Hyperparameters","selected_models_CVREPS_N2013.csv"),header = TRUE)
selected_models_CVREPS_2013N=selected_models_CVREPS_2013N[-1]

#N 2018
selected_models_CVREPS_2018N=read.csv(here("Hyperparameters","selected_models_CVREPS_N2018.csv"),header = TRUE)
selected_models_CVREPS_2018N=selected_models_CVREPS_2018N[-1]

#N 2023
selected_models_CVREPS_2023N=read.csv(here("Hyperparameters","selected_models_CVREPS_N2023.csv"),header = TRUE)
selected_models_CVREPS_2023N=selected_models_CVREPS_2023N[-1]


# Combine all models
CVreps_allN=rbind(selected_models_CVREPS_2013N,
                 selected_models_CVREPS_2018N,
                 selected_models_CVREPS_2023N)

CVreps_allN$year=as.factor(CVreps_allN$year)
CVreps_allN$variable=as.factor(CVreps_allN$variable)


CVREPETS_Nstocks <- ggplot(CVreps_allN, aes(x = CVreps, y = RMSE, color = year, linetype = model)) +
  geom_line(size = 0.4, alpha = 0.8) +
  geom_point(size = 0.8) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 50, 75, 100)) +
  ylab(kgha_expression) +
  xlab("") +
  theme_bw() +
  ggtitle("b) N stocks") +
  scale_color_manual(values = c("#C80036", "#050C9C", "orange"),
                     name = "Year:",
                     labels = c("2013", "2018", "2023")) +
  scale_linetype_manual(values = c("longdash", "solid"),
                        name = "Random forest models:\n\n",
                        labels = c(initialmodel_expresion, finalmodel_expression)) +
  guides(color = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 2)) +  # Ajusta el número de filas según sea necesario
  labs(linetype = "") +  # Elimina la etiqueta de la leyenda para "linetype"
  # facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(limits = c(30,70), breaks = c(30,40,50,60,70)) +
  theme(axis.text.x = element_text(vjust = 0.6, family = "A", size = 10),
        axis.text.y = element_text(vjust = 0.6, family = "A", size = 10),
        text = element_text(size = 10, family = "A"),
        plot.title = element_text(size = 14, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 13, hjust = 0, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11, face = "bold"))
CVREPETS_Nstocks

## 4.2. Dataframes  C/N ratio ----------------------------------------------------------------
#CNRATIO 2013
selected_models_CVREPS_2013CNR=read.csv(here("Hyperparameters","selected_models_CVREPS_CNRATIO2013.csv"),header = TRUE)
selected_models_CVREPS_2013CNR=selected_models_CVREPS_2013CNR[-1]

#CN RATIO 2018
selected_models_CVREPS_2018CNR=read.csv(here("Hyperparameters","selected_models_CVREPS_CNRATIO2018.csv"),header = TRUE)
selected_models_CVREPS_2018CNR=selected_models_CVREPS_2018CNR[-1]

#CNRATIO 2023
selected_models_CVREPS_2023CNR=read.csv(here("Hyperparameters","selected_models_CVREPS_CNRATIO2023.csv"),header = TRUE)
selected_models_CVREPS_2023CNR=selected_models_CVREPS_2023CNR[-1]


# Combine all models
names(selected_models_CVREPS_2018CNR)
CVreps_allCNR=rbind(selected_models_CVREPS_2013CNR,
                  selected_models_CVREPS_2018CNR,
                  selected_models_CVREPS_2023CNR)

CVreps_allCNR$year=as.factor(CVreps_allCNR$year)
CVreps_allCNR$variable=as.factor(CVreps_allCNR$variable)


CVREPETS_CNratio <- ggplot(CVreps_allCNR, aes(x = CVreps, y = RMSE, color = year, linetype = model)) +
  geom_line(size = 0.4, alpha = 0.8) +
  geom_point(size = 0.8) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25, 50, 75, 100)) +
  ylab("ratio") +
  xlab("") +
  theme_bw() +
  ggtitle("c) C/N ratio") +
  scale_color_manual(values = c("#C80036", "#050C9C", "orange"),
                     name = "Year:",
                     labels = c("2013", "2018", "2023")) +
  scale_linetype_manual(values = c("longdash", "solid"),
                        name = "Random forest models:\n\n",
                        labels = c(initialmodel_expresion, finalmodel_expression)) +
  guides(color = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 2)) +  # Ajusta el número de filas según sea necesario
  labs(linetype = "") +  # Elimina la etiqueta de la leyenda para "linetype"
  # facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(limits = c(2,6), breaks = c(2,3,4,5,6)) +
  theme(axis.text.x = element_text(vjust = 0.6, family = "A", size = 10),
        axis.text.y = element_text(vjust = 0.6, family = "A", size = 10),
        text = element_text(size = 10, family = "A"),
        plot.title = element_text(size = 14, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 13, hjust = 0, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11, face = "bold"))
CVREPETS_CNratio

## 4.4 Figure ---------------------------------------------------------------


jpeg(filename = here("Final figures", "2. CV_repetitions_RMSE_fullcolumnsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)

CVREPETS=ggarrange(CVREPETS_Cstocks,
                   CVREPETS_Nstocks,
                   CVREPETS_CNratio,
                   ncol = 1,nrow = 3,
                   common.legend = TRUE,
                   legend = "bottom")
CVREPETS
annotate_figure(CVREPETS,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.2, gp = gpar(cex = 1, fontfamily = "A")),
                bottom = textGrob(CVreps_expresion, vjust = -5.5,gp = gpar(cex = 1, fontfamily = "A")))

dev.off()


#Shorter and edited in paint
jpeg(filename = here("Final figures", "2.1 CV_repetitions_RMSE_fullcolumnsize190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)

CVREPETS=ggarrange(CVREPETS_Cstocks,
                   CVREPETS_Nstocks,
                   CVREPETS_CNratio,
                   ncol = 1,nrow = 3,
                   common.legend = TRUE,
                   legend = "bottom")
CVREPETS
annotate_figure(CVREPETS,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.2, gp = gpar(cex = 1, fontfamily = "A")),
                bottom = textGrob(CVreps_expresion, vjust = -5.5,gp = gpar(cex = 1, fontfamily = "A")))

dev.off()
#-------------------------------------------------------------------------------
# 5. OBSERVED VS PREDICTED BOTH MODELS
## 5.1. C 2013  ----------------------------------------------------------------
#Load data from final and default
ObsPred_C2013_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_2013.csv"),header = TRUE)
ObsPred_C2013_=ObsPred_C2013_[-1]
#RF0= Train (0.92) Test (0.8)
#RF1= Train (0.94) Test (0.88)
View(ObsPred_C2013_)

#Dataframe of R2's
Rsquared2013=rbind(data.frame(model="RF0",data="1. Train",x=2.2,y=24.8,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                   data.frame(model="RF0",data="2. Test",x=2.2,y=22.8,Rsquared=0.80,plotlabel=paste0("R^2",":0.81")),
                   data.frame(model="RF1",data="1. Train",x=2.2,y=24.8,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")),
                   data.frame(model="RF1",data="2. Test",x=2.2,y=22.8,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")))

#Plot both models
(OBS_PRED_C2013_bothmod=ggplot(ObsPred_C2013_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,25),breaks = c(0,5,10,15,20,25))+
    scale_x_continuous(limits = c(0,25),breaks = c(0,5,10,15,20,25))+
    xlab("")+
    ylab("")+
    ggtitle("a) 2013")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))




## 5.2. C 2018  ----------------------------------------------------------------
#Load data from final and default
ObsPred_C2018_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_2018.csv"),header = TRUE)
ObsPred_C2018_=ObsPred_C2018_[-1]
#RF0= Train (0.87) Test (0.65)
#RF1= Train (0.88) Test (0.67)
View(ObsPred_C2018_)

#Dataframe of R2's
Rsquared2018=rbind(data.frame(model="RF0",data="1. Train",x=0.8,y=9.8,Rsquared=0.87,plotlabel=paste0("R^2",":0.87")),
                   data.frame(model="RF0",data="2. Test",x=0.8,y=9,Rsquared=0.65,plotlabel=paste0("R^2",":0.65")),
                   data.frame(model="RF1",data="1. Train",x=0.8,y=9.8,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")),
                   data.frame(model="RF1",data="2. Test",x=0.8,y=9,Rsquared=0.67,plotlabel=paste0("R^2",":0.67")))

#Plot both models
(OBS_PRED_C2018_bothmod=ggplot(ObsPred_C2018_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))+
    scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))+
    xlab("")+
    ylab("")+
    ggtitle("b) 2018")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.3. C 2023  ----------------------------------------------------------------
#Load data from final and default
ObsPred_C2023_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_2023.csv"),header = TRUE)
ObsPred_C2023_=ObsPred_C2023_[-1]
#RF0= Train (0.91) Test (0.81)
#RF1= Train (0.95) Test (0.89)

#Dataframe of R2's
Rsquared2023=rbind(data.frame(model="RF0",data="1. Train",x=0.8,y=11.8,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                   data.frame(model="RF0",data="2. Test",x=0.8,y=10.8,Rsquared=0.81,plotlabel=paste0("R^2",":0.81")),
                   data.frame(model="RF1",data="1. Train",x=0.8,y=11.8,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                   data.frame(model="RF1",data="2. Test",x=0.8,y=10.8,Rsquared=0.89,plotlabel=paste0("R^2",":0.89")))

#Plot both models
(OBS_PRED_C2023_bothmod=ggplot(ObsPred_C2023_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,12),breaks = c(0,3,6,9,12))+
    scale_x_continuous(limits = c(0,12),breaks = c(0,3,6,9,12))+
    xlab("")+
    ylab("")+
    ggtitle("c) 2023")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.4. N 2013  ----------------------------------------------------------------
#Load data from final and default
ObsPred_N2013_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_N2013.csv"),header = TRUE)
ObsPred_N2013_=ObsPred_N2013_[-1]
#RF0= Train (0.91) Test (0.82)
#RF1= Train (0.92) Test (0.84)
View(ObsPred_N2013_)

#Dataframe of R2's
Rsquared2013N=rbind(data.frame(model="RF0",data="1. Train",x=35,y=490,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                    data.frame(model="RF0",data="2. Test",x=35,y=453,Rsquared=0.82,plotlabel=paste0("R^2",":0.82")),
                    data.frame(model="RF1",data="1. Train",x=35,y=490,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                    data.frame(model="RF1",data="2. Test",x=35,y=453,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")))

#Plot both models
(OBS_PRED_N2013_bothmod=ggplot(ObsPred_N2013_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013N,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,500),breaks = c(0,100,200,300,400,500))+
    scale_x_continuous(limits = c(0,500),breaks = c(0,100,200,300,400,500))+
    xlab("")+
    ylab("")+
    ggtitle("a) 2013")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.5. N 2018  ----------------------------------------------------------------
#Load data from final and default
ObsPred_N2018_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_N2018.csv"),header = TRUE)
ObsPred_N2018_=ObsPred_N2018_[-1]
#RF0= Train (0.88) Test (0.80)
#RF1= Train (0.90) Test (0.81)
View(ObsPred_N2018_)

#Dataframe of R2's
Rsquared2018N=rbind(data.frame(model="RF0",data="1. Train",x=17,y=245,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")),
                    data.frame(model="RF0",data="2. Test",x=13,y=227,Rsquared=0.80,plotlabel=paste0("R^2",":0.80")),
                    data.frame(model="RF1",data="1. Train",x=12.9,y=245,Rsquared=0.900,plotlabel=paste0("R^2",":0.90")),
                    data.frame(model="RF1",data="2. Test",x=16.9,y=227,Rsquared=0.81,plotlabel=paste0("R^2",":0.81")))

#Plot both models
(OBS_PRED_N2018_bothmod=ggplot(ObsPred_N2018_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018N,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,250),breaks = c(0,50,100,150,200,250))+
    scale_x_continuous(limits = c(0,250),breaks = c(0,50,100,150,200,250))+
    xlab("")+
    ylab("")+
    ggtitle("b) 2018")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))


## 5.6. N 2023  ----------------------------------------------------------------
#Load data from final and default
ObsPred_N2023_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_N2023.csv"),header = TRUE)
ObsPred_N2023_=ObsPred_N2023_[-1]
#RF0= Train (0.90) Test (0.84)
#RF1= Train (0.92) Test (0.87)
View(ObsPred_N2023_)

#Dataframe of R2's
Rsquared2023N=rbind(data.frame(model="RF0",data="1. Train",x=17,y=283,Rsquared=0.90,plotlabel=paste0("R^2",":0.90")),
                    data.frame(model="RF0",data="2. Test",x=19.6,y=259,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")),
                    data.frame(model="RF1",data="1. Train",x=18.1,y=283,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                    data.frame(model="RF1",data="2. Test",x=19.3,y=259,Rsquared=0.87,plotlabel=paste0("R^2",":0.87")))

#Plot both models
(OBS_PRED_N2023_bothmod=ggplot(ObsPred_N2023_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023N,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,300),breaks = c(0,60,120,180,240,300))+
    scale_x_continuous(limits = c(0,300),breaks = c(0,60,120,180,240,300))+
    xlab("")+
    ylab("")+
    ggtitle("c) 2023")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.7. C/N RATIO 2013  ----------------------------------------------------------------
#Load data from final and default
ObsPred_CNR2013_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_CNR2013.csv"),header = TRUE)
ObsPred_CNR2013_=ObsPred_CNR2013_[-1]
#RF0= Train (0.94) Test (0.91)
#RF1= Train (0.96) Test (0.94)
View(ObsPred_N2013_)

#Dataframe of R2's
Rsquared2013CNR=rbind(data.frame(model="RF0",data="1. Train",x=30.5,y=69,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")),
                      data.frame(model="RF0",data="2. Test",x=30.5,y=66,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                      data.frame(model="RF1",data="1. Train",x=30.5,y=69,Rsquared=0.96,plotlabel=paste0("R^2",":0.96")),
                      data.frame(model="RF1",data="2. Test",x=30.5,y=66,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")))

#Plot both models
(OBS_PRED_CNR2013_bothmod=ggplot(ObsPred_CNR2013_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013CNR,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(28,70),breaks = c(30,40,50,60,70))+
    scale_x_continuous(limits = c(28,70),breaks = c(30,40,50,60,70))+
    xlab("")+
    ylab("")+
    ggtitle("a) 2013")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.8. C/N RATIO 2018  ----------------------------------------------------------------
#Load data from final and default
ObsPred_CNR2018_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_CNR2018.csv"),header = TRUE)
ObsPred_CNR2018_=ObsPred_CNR2018_[-1]
#RF0= Train (0.92) Test (0.84)
#RF1= Train (0.95) Test (0.91)
View(ObsPred_N2018_)

#Dataframe of R2's
Rsquared2018CNR=rbind(data.frame(model="RF0",data="1. Train",x=31.8,y=54,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                      data.frame(model="RF0",data="2. Test",x=31.8,y=52,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")),
                      data.frame(model="RF1",data="1. Train",x=31.8,y=54,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                      data.frame(model="RF1",data="2. Test",x=31.8,y=52,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")))

#Plot both models
(OBS_PRED_CNR2018_bothmod=ggplot(ObsPred_CNR2018_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018CNR,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(30,55),breaks = c(30,35,40,45,50,55))+
    scale_x_continuous(limits = c(30,55),breaks = c(30,35,40,45,50,55))+
    xlab("")+
    ylab("")+
    ggtitle("b) 2018")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 5.9. C/N RATIO 2023  ----------------------------------------------------------------
#Load data from final and default
ObsPred_CNR2023_=read.csv(here("Observed_vs_predicted","OBS_PRED_RFMODELS_ALL_CNR2023.csv"),header = TRUE)
ObsPred_CNR2023_=ObsPred_CNR2023_[-1]
#RF0= Train (0.92) Test (0.91)
#RF1= Train (0.95) Test (0.92)

#Dataframe of R2's
Rsquared2023CNR=rbind(data.frame(model="RF0",data="1. Train",x=35,y=79,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                      data.frame(model="RF0",data="2. Test",x=35,y=75,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                      data.frame(model="RF1",data="1. Train",x=35,y=79,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                      data.frame(model="RF1",data="2. Test",x=35,y=75,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")))

#Plot both models
(OBS_PRED_CNR2023_bothmod=ggplot(ObsPred_CNR2023_, aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023CNR,aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(31,81),breaks = c(30,40,50,60,70,80))+
    scale_x_continuous(limits = c(31,81),breaks = c(30,40,50,60,70,80))+
    xlab("")+
    ylab("")+
    ggtitle("c) 2023")+
    theme_bw()+
    facet_wrap(~model)+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

#-------------------------------------------------------------------------------
# 6. OBSERVED VS PREDICTED FINAL MODEL
## 6.1. C 2013  ----------------------------------------------------------------
#RF0= Train (0.92) Test (0.8)
#RF1= Train (0.94) Test (0.88)
#Final model
Rsquared2013_2=rbind(data.frame(model="RF0",data="1. Train",x=3.1,y=24.8,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                   data.frame(model="RF0",data="2. Test",x=3.1,y=22.8,Rsquared=0.80,plotlabel=paste0("R^2",":0.81")),
                   data.frame(model="RF1",data="1. Train",x=3.1,y=24.8,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")),
                   data.frame(model="RF1",data="2. Test",x=3.1,y=22.8,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")))


(OBS_PRED_C2013_FINAL=ggplot(ObsPred_C2013_[ObsPred_C2013_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013_2[Rsquared2013_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,25),breaks = c(0,5,10,15,20,25))+
    scale_x_continuous(limits = c(0,25),breaks = c(0,5,10,15,20,25))+
    xlab("")+
    ylab(Mgha_expression)+
    ggtitle("2013",subtitle = "a) C stocks")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0,0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))


## 6.2. C 2018  ----------------------------------------------------------------
#RF0= Train (0.87) Test (0.65)
#RF1= Train (0.88) Test (0.67)

#Final model
Rsquared2018_2=rbind(data.frame(model="RF0",data="1. Train",x=1.2,y=9.8,Rsquared=0.87,plotlabel=paste0("R^2",":0.87")),
                   data.frame(model="RF0",data="2. Test",x=1.2,y=9,Rsquared=0.65,plotlabel=paste0("R^2",":0.65")),
                   data.frame(model="RF1",data="1. Train",x=1.2,y=9.8,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")),
                   data.frame(model="RF1",data="2. Test",x=1.2,y=9,Rsquared=0.67,plotlabel=paste0("R^2",":0.67")))

(OBS_PRED_C2018_FINAL=ggplot(ObsPred_C2018_[ObsPred_C2018_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018_2[Rsquared2018_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))+
    scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))+
    xlab(Mgha_expression)+
    ylab("")+
    ggtitle("2018",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.3. C 2023  ----------------------------------------------------------------
#RF0= Train (0.91) Test (0.81)
#RF1= Train (0.95) Test (0.89)

#Final model
Rsquared2023_2=rbind(data.frame(model="RF0",data="1. Train",x=1.4,y=11.8,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                   data.frame(model="RF0",data="2. Test",x=1.4,y=10.8,Rsquared=0.81,plotlabel=paste0("R^2",":0.81")),
                   data.frame(model="RF1",data="1. Train",x=1.4,y=11.8,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                   data.frame(model="RF1",data="2. Test",x=1.4,y=10.8,Rsquared=0.89,plotlabel=paste0("R^2",":0.89")))


(OBS_PRED_C2023_FINAL=ggplot(ObsPred_C2023_[ObsPred_C2023_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023_2[Rsquared2023_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,12),breaks = c(0,3,6,9,12))+
    scale_x_continuous(limits = c(0,12),breaks = c(0,3,6,9,12))+
    xlab("")+
    ylab("")+
    ggtitle("2023",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))




## 6.4. N 2013  ----------------------------------------------------------------
#RF0= Train (0.91) Test (0.82)
#RF1= Train (0.92) Test (0.84)

#Final model
Rsquared2013N_2=rbind(data.frame(model="RF0",data="1. Train",x=63,y=492,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                    data.frame(model="RF0",data="2. Test",x=63,y=457,Rsquared=0.82,plotlabel=paste0("R^2",":0.82")),
                    data.frame(model="RF1",data="1. Train",x=63,y=492,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                    data.frame(model="RF1",data="2. Test",x=63,y=457,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")))


(OBS_PRED_N2013_FINAL=ggplot(ObsPred_N2013_[ObsPred_N2013_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013N_2[Rsquared2013N_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,500),breaks = c(0,100,200,300,400,500))+
    scale_x_continuous(limits = c(0,500),breaks = c(0,100,200,300,400,500))+
    xlab("")+
    ylab(kgha_expression)+
    ggtitle("",subtitle = "b) N stocks")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.5. N 2018  ----------------------------------------------------------------
#RF0= Train (0.88) Test (0.80)
#RF1= Train (0.90) Test (0.81)

#Final model
Rsquared2018N_2=rbind(data.frame(model="RF0",data="1. Train",x=32,y=245,Rsquared=0.88,plotlabel=paste0("R^2",":0.88")),
                    data.frame(model="RF0",data="2. Test",x=27.8,y=226,Rsquared=0.80,plotlabel=paste0("R^2",":0.80")),
                    data.frame(model="RF1",data="1. Train",x=27.8,y=245,Rsquared=0.900,plotlabel=paste0("R^2",":0.90")),
                    data.frame(model="RF1",data="2. Test",x=32,y=226,Rsquared=0.81,plotlabel=paste0("R^2",":0.81")))

(OBS_PRED_N2018_FINAL=ggplot(ObsPred_N2018_[ObsPred_N2018_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018N_2[Rsquared2018N_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,250),breaks = c(0,50,100,150,200,250))+
    scale_x_continuous(limits = c(0,250),breaks = c(0,50,100,150,200,250))+
    xlab(kgha_expression)+
    ylab("")+
    ggtitle("",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.6. N 2023  ----------------------------------------------------------------
#RF0= Train (0.90) Test (0.84)
#RF1= Train (0.92) Test (0.87)

#Final model
#Dataframe of R2's
Rsquared2023N_2=rbind(data.frame(model="RF0",data="1. Train",x=35,y=290,Rsquared=0.90,plotlabel=paste0("R^2",":0.90")),
                    data.frame(model="RF0",data="2. Test",x=35.3,y=265,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")),
                    data.frame(model="RF1",data="1. Train",x=35.9,y=290,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                    data.frame(model="RF1",data="2. Test",x=36.3,y=265,Rsquared=0.87,plotlabel=paste0("R^2",":0.87")))


(OBS_PRED_N2023_FINAL=ggplot(ObsPred_N2023_[ObsPred_N2023_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023N_2[Rsquared2023N_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(0,300),breaks = c(0,60,120,180,240,300))+
    scale_x_continuous(limits = c(0,300),breaks = c(0,60,120,180,240,300))+
    xlab("")+
    ylab("")+
    ggtitle("",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.7. C/N RATIO 2013  ----------------------------------------------------------------
#RF0= Train (0.94) Test (0.91)
#RF1= Train (0.96) Test (0.94)

#Final model
Rsquared2013CNR_2=rbind(data.frame(model="RF0",data="1. Train",x=33.1,y=69,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")),
                      data.frame(model="RF0",data="2. Test",x=33.1,y=66,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                      data.frame(model="RF1",data="1. Train",x=33.1,y=69,Rsquared=0.96,plotlabel=paste0("R^2",":0.96")),
                      data.frame(model="RF1",data="2. Test",x=33.1,y=66,Rsquared=0.94,plotlabel=paste0("R^2",":0.94")))


(OBS_PRED_CNR2013_FINAL=ggplot(ObsPred_CNR2013_[ObsPred_CNR2013_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2013CNR_2[Rsquared2013CNR_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(28,70),breaks = c(30,40,50,60,70))+
    scale_x_continuous(limits = c(28,70),breaks = c(30,40,50,60,70))+
    xlab("")+
    ylab("ratio")+
    ggtitle("",subtitle = "c) C/N ratio")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.8. C/N RATIO 2018  ----------------------------------------------------------------
#RF0= Train (0.92) Test (0.84)
#RF1= Train (0.95) Test (0.91)

#Final model
Rsquared2018CNR_2=rbind(data.frame(model="RF0",data="1. Train",x=33.1,y=54.1,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                      data.frame(model="RF0",data="2. Test",x=33.1,y=52.1,Rsquared=0.84,plotlabel=paste0("R^2",":0.84")),
                      data.frame(model="RF1",data="1. Train",x=33.1,y=54.1,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                      data.frame(model="RF1",data="2. Test",x=33.1,y=52.1,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")))


(OBS_PRED_CNR2018_FINAL=ggplot(ObsPred_CNR2018_[ObsPred_CNR2018_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2018CNR_2[Rsquared2018CNR_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(30,55),breaks = c(30,35,40,45,50,55))+
    scale_x_continuous(limits = c(30,55),breaks = c(30,35,40,45,50,55))+
    xlab("ratio")+
    ylab("")+
    ggtitle("",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

## 6.9. C/N RATIO 2023  ----------------------------------------------------------------
#RF0= Train (0.92) Test (0.84)
#RF1= Train (0.95) Test (0.91)

#Final model
Rsquared2023CNR_2=rbind(data.frame(model="RF0",data="1. Train",x=37,y=80,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")),
                      data.frame(model="RF0",data="2. Test",x=37,y=76,Rsquared=0.91,plotlabel=paste0("R^2",":0.91")),
                      data.frame(model="RF1",data="1. Train",x=37,y=80,Rsquared=0.95,plotlabel=paste0("R^2",":0.95")),
                      data.frame(model="RF1",data="2. Test",x=37,y=76,Rsquared=0.92,plotlabel=paste0("R^2",":0.92")))


(OBS_PRED_CNR2023_FINAL=ggplot(ObsPred_CNR2023_[ObsPred_CNR2023_$model == "RF1", ], aes(x = obs, y = pred, color = data)) +
    geom_point(size=0.5,alpha=0.8) +
    geom_text(data = Rsquared2023CNR_2[Rsquared2023CNR_2$model=="RF1",],aes(x=x,y=y,label = plotlabel,color = data),
              family="A",show.legend = FALSE,parse = TRUE)+
    geom_abline(intercept = 0, slope = 1, color="#6C5F5B",size=.3,linetype = "dashed") +
    stat_smooth(method = "lm", se = FALSE, size = 0.4) +
    labs(x = "Observed", y = "Predicted",color = "Data split:") +
    scale_color_manual(values = c("#059212", "#003285"), 
                       labels = c("Train (85%)","Test (15%)")) +
    scale_y_continuous(limits = c(31,81),breaks = c(30,40,50,60,70,80))+
    scale_x_continuous(limits = c(31,81),breaks = c(30,40,50,60,70,80))+
    xlab("")+
    ylab("")+
    ggtitle("",subtitle = "")+
    theme_bw()+
    theme(text = element_text(size=11,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.x = element_text(vjust=0.6,family = "A",size = 10.2),
          axis.text.y = element_text(vjust=0.6,family = "A",size = 10.2),
          legend.position="bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0,size = 13,face = "bold")))

#-------------------------------------------------------------------------------
# 7. OBSERVED VS PREDICTED FIGURES
### 7.1. Both models ---------------------------------------
#### C stocks --------------------------------------------------
jpeg(filename = here("Final figures", "3. Obs_Pred_RF0_RF1_Cstocks_fullsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
Obs_Vs_Pred_plots=ggarrange(OBS_PRED_C2013_bothmod,
                            OBS_PRED_C2018_bothmod,
                            OBS_PRED_C2023_bothmod,
                            ncol = 1,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
Obs_Vs_Pred_plots
annotate_figure(Obs_Vs_Pred_plots,
                left = textGrob(Predicted_Mgha_expression, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_Mgha_expression, vjust = -3,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()

#Shorter and edited in Paint
jpeg(filename = here("Final figures", "3.1. Obs_Pred_RF0_RF1_Cstocks_fullsize190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
Obs_Vs_Pred_plots=ggarrange(OBS_PRED_C2013_bothmod,
                            OBS_PRED_C2018_bothmod,
                            OBS_PRED_C2023_bothmod,
                            ncol = 1,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
Obs_Vs_Pred_plots
annotate_figure(Obs_Vs_Pred_plots,
                left = textGrob(Predicted_Mgha_expression, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_Mgha_expression, vjust = -3,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()


#### N stocks --------------------------------------------------
jpeg(filename = here("Final figures", "4. Obs_Pred_RF0_RF1_Nstocks_fullsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
Obs_Vs_Pred_plotsN=ggarrange(OBS_PRED_N2013_bothmod,
                            OBS_PRED_N2018_bothmod,
                            OBS_PRED_N2023_bothmod,
                            ncol = 1,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
Obs_Vs_Pred_plotsN
annotate_figure(Obs_Vs_Pred_plotsN,
                left = textGrob(Predicted_kgha_expression, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_kgha_expression, vjust = -3,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()

#Shorter and edited in Paint
jpeg(filename = here("Final figures", "4.1, Obs_Pred_RF0_RF1_Nstocks_190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
Obs_Vs_Pred_plotsN=ggarrange(OBS_PRED_N2013_bothmod,
                             OBS_PRED_N2018_bothmod,
                             OBS_PRED_N2023_bothmod,
                             ncol = 1,nrow = 3,
                             common.legend = TRUE,
                             legend = "bottom")
Obs_Vs_Pred_plotsN
annotate_figure(Obs_Vs_Pred_plotsN,
                left = textGrob(Predicted_kgha_expression, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_kgha_expression, vjust = -3,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()

#### C/N ratio -------------------------------------------------------
jpeg(filename = here("Final figures", "5. Obs_Pred_RF0_RF1_CNratio_fullsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
Obs_Vs_Pred_plotsCNR=ggarrange(OBS_PRED_CNR2013_bothmod,
                             OBS_PRED_CNR2018_bothmod,
                             OBS_PRED_CNR2023_bothmod,
                             ncol = 1,nrow = 3,
                             common.legend = TRUE,
                             legend = "bottom")
Obs_Vs_Pred_plotsCNR
annotate_figure(Obs_Vs_Pred_plotsCNR,
                left = textGrob(Predicted_ratio_expresion, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_ratio_expresion, vjust = -4.6,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()


#Shorter and edited in paint
jpeg(filename = here("Final figures", "5.1 Obs_Pred_RF0_RF1_CNratio_fullsize190,220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
Obs_Vs_Pred_plotsCNR=ggarrange(OBS_PRED_CNR2013_bothmod,
                               OBS_PRED_CNR2018_bothmod,
                               OBS_PRED_CNR2023_bothmod,
                               ncol = 1,nrow = 3,
                               common.legend = TRUE,
                               legend = "bottom")
Obs_Vs_Pred_plotsCNR
annotate_figure(Obs_Vs_Pred_plotsCNR,
                left = textGrob(Predicted_ratio_expresion, rot = 90, vjust = 1, gp = gpar(cex = 0.9, fontfamily = "A")),
                bottom = textGrob(Observed_ratio_expresion, vjust = -4.6,gp = gpar(cex = 0.9, fontfamily = "A")))

dev.off()



### 7.2. RF1 models ----------------------------------
#Needs edition in pain to edit general x label axis and legend
jpeg(filename = here("Final figures", "6. Obs_vs_Pred_ALL_RF1_fullsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
Obs_Vs_Pred_plots_final=ggarrange(OBS_PRED_C2013_FINAL,
                                  OBS_PRED_C2018_FINAL,
                                  OBS_PRED_C2023_FINAL,
                                  OBS_PRED_N2013_FINAL,
                                  OBS_PRED_N2018_FINAL,
                                  OBS_PRED_N2023_FINAL,
                                  OBS_PRED_CNR2013_FINAL,
                                  OBS_PRED_CNR2018_FINAL,
                                  OBS_PRED_CNR2023_FINAL,
                            ncol = 3,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
Obs_Vs_Pred_plots_final
annotate_figure(Obs_Vs_Pred_plots_final,
                left = textGrob(Predictedvalues_RF1, rot = 90, vjust = 0.3, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(Observedvalues_expresion,,hjust = 2.7, vjust = -3.1,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()


#Shorter Needs edition in pain to edit general x label axis and legend
jpeg(filename = here("Final figures", "6.1. Obs_vs_Pred_ALL_RF1_fullsize190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
Obs_Vs_Pred_plots_final=ggarrange(OBS_PRED_C2013_FINAL,
                                  OBS_PRED_C2018_FINAL,
                                  OBS_PRED_C2023_FINAL,
                                  OBS_PRED_N2013_FINAL,
                                  OBS_PRED_N2018_FINAL,
                                  OBS_PRED_N2023_FINAL,
                                  OBS_PRED_CNR2013_FINAL,
                                  OBS_PRED_CNR2018_FINAL,
                                  OBS_PRED_CNR2023_FINAL,
                                  ncol = 3,nrow = 3,
                                  common.legend = TRUE,
                                  legend = "bottom")
Obs_Vs_Pred_plots_final
annotate_figure(Obs_Vs_Pred_plots_final,
                left = textGrob(Predictedvalues_RF1, rot = 90, vjust = 0.3, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(Observedvalues_expresion,,hjust = 2.7, vjust = -3.1,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()



#-------------------------------------------------------------------------------
# 8.NTREE ERROR 
# 8.1. C 2013 ------------------------------------------------------------------
ntree_C2013=read.csv(here("Hyperparameters","ntree_data_bothmodels_2013.csv"),header = TRUE)
View(ntree_C2013)
(NTREE_C2013=ggplot(data = ntree_C2013, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(3.1, 5),breaks = c(3.5,4,4.5,5)) +
    ggtitle("2013",subtitle = "a) C stocks")+
    ylab(Mgha_expression)+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=3.149019,x=600,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=3.476576,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.2. C 2018 ------------------------------------------------------------------
ntree_C2018=read.csv(here("Hyperparameters","ntree_data_bothmodels_2018.csv"),header = TRUE)
View(ntree_C2018)
(NTREE_C2018=ggplot(data = ntree_C2018, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(1.3,1.9),breaks = c(1.3,1.5,1.7,1.9)) +
    ggtitle("2018",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=1.378992,x=1000,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=1.401540,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))



# 8.3. C 2023 ------------------------------------------------------------------
ntree_C2023=read.csv(here("Hyperparameters","ntree_data_bothmodels_2023.csv"),header = TRUE)
View(ntree_C2023)
(NTREE_C2023=ggplot(data = ntree_C2023, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(1.45,2.5),breaks = c(1.5,1.75,2,2.25,2.5)) +
    ggtitle("2023",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=1.499112,x=500,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=1.540701,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))



# 8.4. N 2013 ------------------------------------------------------------------
ntree_N2013=read.csv(here("Hyperparameters","ntree_data_bothmodels_N2013.csv"),header = TRUE)
View(ntree_N2013)
(NTREE_N2013=ggplot(data = ntree_N2013, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(58, 105),breaks = c(60,70,80,90,100)) +
        ggtitle("",subtitle = "b) N stocks")+
    ylab(kgha_expression)+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=61.07979,x=800,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=68.49563,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.5. N 2018 ------------------------------------------------------------------
ntree_N2018=read.csv(here("Hyperparameters","ntree_data_bothmodels_N2018.csv"),header = TRUE)
View(ntree_N2018)
(NTREE_N2018=ggplot(data = ntree_N2018, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(38, 60),breaks = c(40,45,50,55,60)) +
    ggtitle("",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=39.35826,x=500,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=40.48926,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.6. N 2023 ------------------------------------------------------------------
ntree_N2023=read.csv(here("Hyperparameters","ntree_data_bothmodels_N2023.csv"),header = TRUE)
View(ntree_N2023)
(NTREE_N2023=ggplot(data = ntree_N2023, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(34, 56),breaks = c(35,40,45,50,55)) +
    ggtitle("",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=36.01226,x=900,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=37.60615,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.7. C/N ratio 2013 ------------------------------------------------------------------
ntree_CNR2013=read.csv(here("Hyperparameters","ntree_data_bothmodels_CNR2013.csv"),header = TRUE)
View(ntree_CNR2013)
(ntree_CNR2013_=ggplot(data = ntree_CNR2013, aes(x = trees, y = modRMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(0.8,5),breaks = c(1,2,3,4,5)) +
    ggtitle("",subtitle = "c) C/N ratio")+
    ylab(expression("ratio"))+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=0.927395,x=800,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=1.338218,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.8. C/N ratio 2018 ------------------------------------------------------------------
ntree_CNR2018=read.csv(here("Hyperparameters","ntree_data_bothmodels_CNR2018.csv"),header = TRUE)
View(ntree_CNR2018)
(ntree_CNR2018_=ggplot(data = ntree_CNR2018, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(3.6,8),breaks = c(4,5,6,7,8)) +
    ggtitle("",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=3.767067,x=500,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=4.015744,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.9. C/N ratio 2023 ------------------------------------------------------------------
ntree_CNR2023=read.csv(here("Hyperparameters","ntree_data_bothmodels_CNR2023.csv"),header = TRUE)
View(ntree_CNR2023)
(ntree_CNR2023_=ggplot(data = ntree_CNR2023, aes(x = trees, y = RMSE,color=model)) +
    geom_line(linewidth=0.3)+
    theme_bw()+
    scale_y_continuous(limits = c(6,10),breaks = c(6,7,8,9,10)) +
    ggtitle("",subtitle = "")+
    ylab("")+
    xlab("")+
    scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250)) +
    annotate("point",y=6.012532,x=500,color="#050C9C",size=1.5)+ # Located point of hyperparameters - OPTIMIZED MODEL
    annotate("point",y=6.338218,x=500,color="#C80036",size=1.5)+  # Located point of hyperparameters - DEFAULT MODEL
    scale_color_manual(values = c("RF0"="#C80036","RF1"="#050C9C"),
                       name="Random forest models:",
                       labels=c(initialmodel_expresion,
                                finalmodel_expression)) +
    theme(text = element_text(size=10.2,family = "A"),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          axis.text = element_text(size = 10.2),
          axis.text.y = element_text(size = 10.2),
          legend.title = element_text(hjust = 0.5,size = 12,face = "bold",family = "A"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          plot.subtitle = element_text(size = 13,hjust = 0,face = "bold"),
          panel.background = element_rect(fill='white', colour='black')))

# 8.10 Figure ---------------------------------------------------------------

jpeg(filename = here("Final figures", "7. ntree_error_plots_fullsize.jpeg"),
     width = 190,height = 240,units = "mm",res = 1000)
ntree_plots=ggarrange(NTREE_C2013,NTREE_C2018,NTREE_C2023,
                      NTREE_N2013,NTREE_N2018,NTREE_N2023,
                      ntree_CNR2013_,ntree_CNR2018_,ntree_CNR2023_,
                            ncol = 3,nrow = 3,
                            common.legend = TRUE,
                            legend = "bottom")
ntree_plots
annotate_figure(ntree_plots,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.3, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(ntree_expresion, vjust = -3.6,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()

#shorter and edited in Paint
jpeg(filename = here("Final figures", "7.1. ntree_error_plots_fullsize190x220.jpeg"),
     width = 190,height = 220,units = "mm",res = 1000)
ntree_plots=ggarrange(NTREE_C2013,NTREE_C2018,NTREE_C2023,
                      NTREE_N2013,NTREE_N2018,NTREE_N2023,
                      ntree_CNR2013_,ntree_CNR2018_,ntree_CNR2023_,
                      ncol = 3,nrow = 3,
                      common.legend = TRUE,
                      legend = "bottom")
ntree_plots
annotate_figure(ntree_plots,
                left = textGrob(RMSE_expresion, rot = 90, vjust = 0.3, gp = gpar(cex = 1.1, fontfamily = "A")),
                bottom = textGrob(ntree_expresion, vjust = -3.6,gp = gpar(cex = 1.1, fontfamily = "A")))

dev.off()



#-------------------------------------------------------------------------------
# 9. VIMP
# 9.1. C 2013 ------------------------------------------------------------------
Vimp_C2013_FINAL=read.csv(here("Vimp","VIMP_C13_FINAL.csv"),header = TRUE)
View(Vimp_C2013_FINAL)
Vimp_C2013_FINAL$Year=as.factor(Vimp_C2013_FINAL$Year)
Vimp_C2013_FINAL$Plotorder=1
Vimp_C2013_FINAL$COVARS_ORDER=1

Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_C2013_FINAL$COVARS_ORDER[Vimp_C2013_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.2. C 2018 ------------------------------------------------------------------
Vimp_C2018_FINAL=read.csv(here("Vimp","VIMP_C2018_FINAL.csv"),header = TRUE)
View(Vimp_C2018_FINAL)
Vimp_C2018_FINAL$Year=as.factor(Vimp_C2018_FINAL$Year)
Vimp_C2018_FINAL$Plotorder=1
Vimp_C2018_FINAL$COVARS_ORDER=1
# 9.3. C 2023 ------------------------------------------------------------------
Vimp_C2023_FINAL=read.csv(here("Vimp","VIMP_C23_FINAL.csv"),header = TRUE)
View(Vimp_C2023_FINAL)
Vimp_C2023_FINAL$Year=as.factor(Vimp_C2023_FINAL$Year)
Vimp_C2023_FINAL$Plotorder=1
Vimp_C2023_FINAL$COVARS_ORDER=1


# 9.4. N 2013 ------------------------------------------------------------------
Vimp_N2013_FINAL=read.csv(here("Vimp","VIMP_N13_FINAL.csv"),header = TRUE)
View(Vimp_N2013_FINAL)
Vimp_N2013_FINAL$Year=as.factor(Vimp_N2013_FINAL$Year)
Vimp_N2013_FINAL$Plotorder=2
Vimp_N2013_FINAL$COVARS_ORDER=1

Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_N2013_FINAL$COVARS_ORDER[Vimp_N2013_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.5. N 2018 ------------------------------------------------------------------
Vimp_N2018_FINAL=read.csv(here("Vimp","VIMP_N2018_FINAL.csv"),header = TRUE)
View(Vimp_N2018_FINAL)
Vimp_N2018_FINAL$Year=as.factor(Vimp_N2018_FINAL$Year)
Vimp_N2018_FINAL$Plotorder=2
Vimp_N2018_FINAL$COVARS_ORDER=1

Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_N2018_FINAL$COVARS_ORDER[Vimp_N2018_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.6. N 2023 ------------------------------------------------------------------
Vimp_N2023_FINAL=read.csv(here("Vimp","VIMP_N23_FINAL.csv"),header = TRUE)
View(Vimp_N2023_FINAL)
Vimp_N2023_FINAL$Year=as.factor(Vimp_N2023_FINAL$Year)
Vimp_N2023_FINAL$Plotorder=2
Vimp_N2023_FINAL$COVARS_ORDER=1

Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_N2023_FINAL$COVARS_ORDER[Vimp_N2023_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.7. CNRatio 2013 ------------------------------------------------------------------
Vimp_CNR2013_FINAL=read.csv(here("Vimp","VIMP_CNR2013_FINAL.csv"),header = TRUE)
View(Vimp_CNR2013_FINAL)
Vimp_CNR2013_FINAL$Year=as.factor(Vimp_CNR2013_FINAL$Year)
Vimp_CNR2013_FINAL$Plotorder=3
Vimp_CNR2013_FINAL$COVARS_ORDER=1

Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_CNR2013_FINAL$COVARS_ORDER[Vimp_CNR2013_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.7. CNRatio 2018 ------------------------------------------------------------------
Vimp_CNR2018_FINAL=read.csv(here("Vimp","VIMP_CNR2018_FINAL.csv"),header = TRUE)
View(Vimp_CNR2018_FINAL)
Vimp_CNR2018_FINAL$Year=as.factor(Vimp_CNR2018_FINAL$Year)
Vimp_CNR2018_FINAL$Plotorder=3
Vimp_CNR2018_FINAL$COVARS_ORDER=1

Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_CNR2018_FINAL$COVARS_ORDER[Vimp_CNR2018_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.7. CNRatio 2023 ------------------------------------------------------------------
Vimp_CNR2023_FINAL=read.csv(here("Vimp","VIMP_CNR2023_FINAL.csv"),header = TRUE)
View(Vimp_CNR2023_FINAL)
Vimp_CNR2023_FINAL$Year=as.factor(Vimp_CNR2023_FINAL$Year)
Vimp_CNR2023_FINAL$Plotorder=3
Vimp_CNR2023_FINAL$COVARS_ORDER=1

Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "FFLAYER"]<- 12
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "BA"]<- 11
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "AGE"]<- 10
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "COV"]<- 9
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "DTH"]<- 8
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "RICH"]<-7
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "H"]<- 6
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "MAN"]<- 5
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "ELE"]<- 4
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "GAP"]<- 3
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "SLO"]<- 2
Vimp_CNR2023_FINAL$COVARS_ORDER[Vimp_CNR2023_FINAL$abbreviations_predictors  == "ASP"]<- 1

# 9.10 Merge all data -Plots ----------------------------------------------------------
VIMP_FINAL=rbind(Vimp_C2013_FINAL,
                 Vimp_C2018_FINAL,
                 Vimp_C2023_FINAL,
                 Vimp_N2013_FINAL,
                 Vimp_N2018_FINAL,
                 Vimp_N2023_FINAL,
                 Vimp_CNR2013_FINAL,
                 Vimp_CNR2018_FINAL,
                 Vimp_CNR2023_FINAL)

View(VIMP_FINAL)
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "FFLAYER"]<- "FFLAYER"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "BA"]<- "BA"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "AGE"]<- "AGE"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "COV"]<- "COV"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "DTH"]<- "DTH"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "RICH"]<- "RICH"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "H"]<- "H"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "MAN"]<- "MCOND"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "ELE"]<- "ELEV"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "GAP"]<- "GAP"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "SLO"]<- "SLO"
VIMP_FINAL$abbreviations_predictors[VIMP_FINAL$abbreviations_predictors  == "ASP"]<- "ASP"


(VIMP_PLOT=ggplot(VIMP_FINAL, aes(x=reorder(abbreviations_predictors, COVARS_ORDER), y=X.IncMSE,color=Year)) + 
    geom_segment(aes(x=reorder(abbreviations_predictors, COVARS_ORDER), 
                     xend=abbreviations_predictors, 
                     color=Year,
                     y=0, 
                     yend=X.IncMSE)) + 
    geom_point(size = 4,pch=20,alpha=0.8) +
    coord_flip()+
    ylab(IncMSE_expresion)+
    xlab(predictors_expresion)+
    theme_bw()+
    scale_y_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100))+
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
          axis.text = element_text(size = 12)))


#Figure 2 (option2)- 1 column size (180 x 140) - for research publication
jpeg(filename = here("Final figures", "8.Vimp_fullsize1000.jpeg"),
     width = 190,height = 100,units = "mm",res = 1000)
VIMP_PLOT
dev.off()


#-------------------------------------------------------------------------------
# 10. TAYLOR DIAGRAM ------------------------------------------------------
taylor.diagram()
RF1_C2013=ObsPred_C2013_[ObsPred_C2013_$model == "RF1", ]
RF1_C2018=ObsPred_C2018_[ObsPred_C2018_$model == "RF1", ]
RF1_C2023=ObsPred_C2023_[ObsPred_C2023_$model == "RF1", ]
RF1_N2013=ObsPred_N2013_[ObsPred_N2013_$model == "RF1", ]
RF1_N2018=ObsPred_N2018_[ObsPred_N2018_$model == "RF1", ]
RF1_N2023=ObsPred_N2023_[ObsPred_N2023_$model == "RF1", ]
RF1_CNR2013=ObsPred_CNR2013_[ObsPred_CNR2013_$model == "RF1", ]
RF1_CNR2018=ObsPred_CNR2018_[ObsPred_CNR2018_$model == "RF1", ]
RF1_CNR2023=ObsPred_CNR2023_[ObsPred_CNR2023_$model == "RF1", ]

allmodels_taylor=rbind(data.frame(RF1_C2013,Evaluated_model= "M1: C stocks - 2013",Labelvariable="a) C stocks"),
                       data.frame(RF1_C2018,Evaluated_model= "M2: C stocks - 2018",Labelvariable="a) C stocks"),
                       data.frame(RF1_C2023,Evaluated_model= "M3: C stocks - 2023",Labelvariable="a) C stocks"),
                       data.frame(RF1_N2013,Evaluated_model= "M4: N stocks - 2013",Labelvariable="b) N stocks"),
                      data.frame(RF1_N2018,Evaluated_model= "M5: N stocks - 2018",Labelvariable="b) N stocks"),
                      data.frame(RF1_N2023,Evaluated_model= "M6: N stocks - 2023",Labelvariable="b) N stocks"),
                      data.frame(RF1_CNR2013,Evaluated_model= "M7: C/N ratio - 2013",Labelvariable="c) C/N ratio"),
                      data.frame(RF1_CNR2018,Evaluated_model= "M8: C/N ratio - 2018",Labelvariable="c) C/N ratio"),
                      data.frame(RF1_CNR2023,Evaluated_model= "M9: C/N ratio - 2023",Labelvariable="c) C/N ratio"))
str(allmodels_taylor)

taylor.diagram(RF1_C2013$obs,RF1_C2013$pred,col = "blue",pch=19)
taylor.diagram(RF1_C2018$obs,RF1_C2018$pred,col = "blue",pch=19)
taylor.diagram(RF1_C2023$obs,RF1_C2023$pred,col = "blue",pch=19)
taylor.diagram(RF1_N2013$obs,RF1_N2013$pred,col = "blue",pch=19)
taylor.diagram(RF1_N2018$obs,RF1_N2018$pred,col = "blue",pch=19)
taylor.diagram(RF1_N2023$obs,RF1_N2023$pred,col = "blue",pch=19)
taylor.diagram(RF1_CNR2013$obs,RF1_CNR2013$pred,col = "blue",pch=19)
taylor.diagram(RF1_CNR2018$obs,RF1_CNR2018$pred,col = "blue",pch=19)
taylor.diagram(RF1_CNR2023$obs,RF1_CNR2023$pred,col = "blue",pch=19)

allmodels_taylor$year=as.factor(allmodels_taylor$year)
allmodels_taylor$variable=as.factor(allmodels_taylor$variable)
allmodels_taylor$Labelvariable=as.factor(allmodels_taylor$Labelvariable)


TaylorDiagram(allmodels_taylor[allmodels_taylor$variable == "C stocks", ],
              obs = "obs",
              mod = "pred",
              group = "Evaluated_model",
              pch=16,
              normalise = TRUE,
              cex=1.5)

TaylorDiagram(allmodels_taylor[allmodels_taylor$variable == "N stocks", ],
              obs = "obs",
              normalise = TRUE,
              mod = "pred",
              group = "Evaluated_model")

TaylorDiagram(allmodels_taylor[allmodels_taylor$variable == "C/N ratio", ],
              obs = "obs",
              mod = "pred",normalise = TRUE,
              group = "Evaluated_model")


# Crear una paleta de 9 colores personalizados
my_palette <- c("#713ABE", "#478CCF", "#FF7D29",
                "#BC5A94", "#FFA62F", "#74512D",
                "#EF5A6F", "#40A578", "#32012F")


jpeg(filename = here("Final figures", "9.Taylordiagram_fullsize1000.jpeg"),
     width = 190,height = 95,units = "mm",res = 1000)
TaylorDiagram(allmodels_taylor,
              obs = "obs",
              mod = "pred",
              group = c("Evaluated_model"),
              annotate = "centered\nRMS error",
              normalise = TRUE,
              type = "Labelvariable",
              key = TRUE,
              key.pos = "bottom",
              fontsize=9.5,
              fontfamily="A",
              xlab=normalizedsd,
              ylab=normalizedsd,
              par.settings=list(grid.pars=list(fontfamily="A")),
              cor.col="#C70039",
              rms.col="#191D88",
              key.columns = 3,
              auto.text = TRUE,
              key.title = Evaluatedvariables_RF1,
              pch=16,
              cex=0.8,
              col=my_palette)
dev.off()


#11. Flow diagram----------------------------------------------------------
# Instalar y cargar la librería DiagrammeR si no está instalada
  if (!require(DiagrammeR)) install.packages("DiagrammeR")
library(DiagrammeR)

# Crear el diagrama de flujo
DiagrammeR::grViz("
digraph flowchart {
  node [fontname = Helvetica, shape = box, style = filled, fillcolor = lightgray]
  
  A1 [label = 'Exploratory analysis of the target variables']
  A2 [label = 'Check Spearman’s correlation between FF properties and predictors and among predictors']
  A3 [label = 'Split data into train and test subsets']
  A4 [label = 'Run Boruta algorithm to verify if all predictors influence the target variable']
  A5 [label = 'Compare Boruta’s decision with Spearman’s correlation to evaluate similarities']
  A6 [label = 'Remove the variables rejected by the Boruta algorithm']
  A7 [label = 'Build RF0 model with default hyperparameters (ntree = 500, mtry = 4)']
  A8 [label = 'Check performance of RF0 with repeated ten-fold cross-validation']
  A9 [label = 'Perform cross-validation with 5, 10, 15, 25, 50, 75, 100 repetitions']
  A10 [label = 'Record RF0 error metrics of each repetition']
  A11 [label = 'Use train() function to find the best hyperparameter values']
  A12 [label = 'Generate optimized RF1 models for each dependent variable for 2013, 2018, 2023']
  A13 [label = 'Validate RF1 models using cross-validation with the same number of repetitions']
  A14 [label = 'Record results of RMSE, MAE, and R-squared of RF1 models in train and test subsets']
  A15 [label = 'Perform Freiman test to check RMSE decrease between RF0 and RF1']
  A16 [label = 'Check statistical significance of RMSE differences and cross-validation repetition variations']
  A17 [label = 'Plot RF0 and RF1 curves comparing RMSE vs ntree values']
  A18 [label = 'Plot RF0 and RF1 curves comparing RMSE vs mtry values']
  A19 [label = 'Obtain observed vs predicted plots for each model (C stocks, N stocks, C/N ratio) for 2013, 2018, 2023']
  A20 [label = 'Generate observed vs predicted plots comparing RF0 and RF1']
  A21 [label = 'Plot Taylor Diagram for C stocks, N stocks, and C/N ratio']
  A22 [label = 'Obtain importance variable plots (%incMSE) for C stocks, N stocks, and C/N ratio']
  A23 [label = 'Generate map predictions with modelmap library for each variable and evaluated year']
  A24 [label = 'Plot RMSE values reached in ten-fold cross-validation']
  
  subgraph cluster0 {
    label = 'Data Preparation';
    style=filled; color=lightyellow;
    A1 -> A2 -> A3 -> A4;
  }
  
  subgraph cluster1 {
    label = 'Model Building and Evaluation';
    style=filled; color=lightblue;
    A4 -> A5 -> A6 -> A7 -> A8 -> A9 -> A10 -> A11 -> A12 -> A13 -> A14 -> A15 -> A16;
  }
  
  subgraph cluster2 {
    label = 'Visualization and Interpretation';
    style=filled; color=lightgreen;
    A16 -> A17 -> A18 -> A19 -> A20 -> A21 -> A22 -> A23 -> A24;
  }

  # Conexiones adicionales
  A2 -> A4;
  A3 -> A4;
  A5 -> A7;
  A6 -> A7;
}
")

# 8. PARTIAL PLOTS AGE ----------------------------------------------------
AGE_C2013_df=read.csv(here("PartialPlots","PARTIAL_AGE_C13_FINAL.csv"),header = TRUE)
AGE_C2013_df=AGE_C2013_df[-1]
View(AGE_C2013_df)

#Merge the three years
AGE_C_trend=rbind(SA_C2013_df,SA_C2018_df,SA_C2023_df)
AGE_C_trend$Year=as.factor(SA_C_trend$Year)

#Merge the three years
SA_C_trend=rbind(SA_C2013_df,SA_C2018_df,SA_C2023_df)
SA_C_trend$Year=as.factor(SA_C_trend$Year)