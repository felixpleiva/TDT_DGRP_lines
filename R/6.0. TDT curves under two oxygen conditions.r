# -------------------------------------------------------------------------------------
# Script to analyse the effect of body size-related traits on survival time of 
# DGRP lines under two oxygen conditions
# -------------------------------------------------------------------------------------
# Data generated by F�lix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# Script created by F�lix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# ------------------------------------------------------------------------------
# Cite as:

# Leiva FP, Santos M, Rezende EL & Verberk WCEP. (2021). Paper data and code of
# manuscript: Intraspecific variation on heat tolerance in a model ectotherm:
# the role of oxygen, cell size, body size and duration. Zenodo.
# https://doi.org/10.5281/zenodo.5120028.
# ------------------------------------------------------------------------------
# Cleaning working space
rm(list=ls())
today<-format(Sys.Date(),"%Y%m%d")
# ------------------------------------------------------------------------------
# set the working directory to the folder containing this script:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()# check directory
# -------------------------------------------------------------------------------------
#Libraries
library(nlme)
library(car)
library(visreg)
library(MuMIn)
library(AICcmodavg)
library(MASS)
# ------------------------------------------------------------------------------
#Load data of hypoxia (6 lines)
all.lines.hypo<-read.csv("../Outputs/4.1.2. Survival time of Drosophila melanogaster and its associated size traits at two oxygen levels.csv")

10^(min(all.lines.hypo$surv.time))#lower range in minutes
10^(max(all.lines.hypo$surv.time))#upper range in minutes
10^(median(all.lines.hypo$surv.time))#median in minutes
# ------------------------------------------------------------------------------
#Step 1: aggregate data into median values
# ------------------------------------------------------------------------------
dhypo.median=aggregate(surv.time~test.temp+cell.area+
                       cell.number+stock+test.oxygen+iod+fw+ws+sex,all.lines.hypo,median)

#see structure of data
str(dhypo.median)

# transform structure of data
dhypo.median$test.temp<-as.numeric(dhypo.median$test.temp)
dhypo.median$test.oxygen<-as.numeric(dhypo.median$test.oxygen)
dhypo.median$stock<-as.factor(dhypo.median$stock)
dhypo.median$sex<-as.factor(dhypo.median$sex)

# ------------------------------------------------------------------------------
# Step 2: Correct traits for effect of sex so that the analysis works with flies 
# that have trait values that are relatively large or small for their sex
# ------------------------------------------------------------------------------
# relative cell area
dhypo.median$cell.area.rel<-as.numeric(scale(resid(lm(cell.area~sex,data=dhypo.median))))

# relative cell number
dhypo.median$cell.number.rel<-as.numeric(scale(resid(lm(cell.number~sex,data=dhypo.median))))

# relative fresh mass
dhypo.median$fw.rel<-as.numeric(scale(resid(lm(fw~sex,data=dhypo.median))))

# relative interocular distance
dhypo.median$iod.rel<-as.numeric(scale(resid(lm(iod~sex,data=dhypo.median))))

# relative wing size
dhypo.median$ws.rel<-as.numeric(scale(resid(lm(ws~sex,data=dhypo.median))))

# ------------------------------------------------------------------------------
#Step 3: Start with the full model and calculate VIFs
# ------------------------------------------------------------------------------
fhypomedian<-gls(surv.time~
                       scale(test.temp)
               +scale(test.temp)*scale(test.oxygen)
               +scale(test.temp)*cell.area.rel*scale(test.oxygen)
               +scale(test.temp)*cell.number.rel*scale(test.oxygen)
               +scale(test.temp)*ws.rel*scale(test.oxygen)
               +scale(test.temp)*fw.rel*scale(test.oxygen)
               +scale(test.temp)*iod.rel*scale(test.oxygen)
               +scale(test.temp)*sex*scale(test.oxygen)
               ,data=dhypo.median,na.action = na.omit,method = "ML")

vif(fhypomedian)# high VIFs

fhypomedian<-gls(surv.time~
                 scale(test.temp)
               +scale(test.oxygen)
               +scale(test.oxygen)*scale(test.temp)
               +scale(test.temp)*cell.area.rel*scale(test.oxygen)
               +scale(test.temp)*cell.number.rel*scale(test.oxygen)
               #+scale(test.temp)*ws.rel*scale(test.oxygen)
               +scale(test.temp)*fw.rel*scale(test.oxygen)
               +scale(test.temp)*iod.rel*scale(test.oxygen)
               +scale(test.temp)*sex*scale(test.oxygen)
               ,data=dhypo.median,na.action = na.omit,method = "ML")

vif(fhypomedian)# VIFs are now 0K!

anova(fhypomedian)
drop1(fhypomedian)
stepAIC(fhypomedian)
AICc(fhypomedian)#-77.73107
# ------------------------------------------------------------------------------
# Step 4: Analyse a list of candidate models (except for wing size) and compare 
# their AICs
# ------------------------------------------------------------------------------
m0<-gls(surv.time~
            scale(test.temp)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m0.1<-gls(surv.time~
          scale(test.temp)*scale(test.oxygen)
        ,data=dhypo.median,na.action = na.omit,method = "ML")

m1.1<-gls(surv.time~
                scale(test.temp)*scale(test.oxygen)
                +scale(test.temp)*sex*scale(test.oxygen)
               ,data=dhypo.median,na.action = na.omit,method = "ML")

m1.2<-gls(surv.time~
                 scale(test.temp)*scale(test.oxygen)
               +scale(test.temp)*cell.area.rel*scale(test.oxygen)
               ,data=dhypo.median,na.action = na.omit,method = "ML")

m1.3<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*cell.number.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m1.4<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*fw.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m1.5<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*iod.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

fit.list.step4 <- list(m0,m0.1,m1.1,m1.2,m1.3,m1.4,m1.5)

fit.names.step4 <-c("ST",
                 "ST * O2",  
                 "ST * O2 + ST * sex * O2",
                 "ST * O2 + ST * cell area * O2",
                 "ST * O2 + ST * cell number * O2",
                 "ST * O2 + ST * FM * O2",
                 "ST * O2 + ST * IOD * O2")

#compare by using AICc
fit.step4<-aictab(fit.list.step4,fit.names.step4, second.ord = T,sort = TRUE, digits = 3, LL=TRUE)
fit.step4

#                                 K    AICc Delta_AICc AICcWt Cum.Wt    LL
# ST * O2 + ST * sex * O2         9 -102.26       0.00   0.50   0.50 61.23
# ST * O2 + ST * cell area * O2   9 -101.85       0.41   0.41   0.91 61.02
# ST * O2                         5  -98.32       3.94   0.07   0.98 54.51
# ST * O2 + ST * IOD * O2         9  -93.80       8.46   0.01   0.99 57.00
# ST * O2 + ST * FM * O2          9  -92.58       9.68   0.00   1.00 56.39
# ST * O2 + ST * cell number * O2 9  -92.47       9.80   0.00   1.00 56.33
# ST                              3  -83.93      18.33   0.00   1.00 45.10

# sex and test temperature and oxygen interaction improve model fit. Important to
# highlight is the effect of test temperature and oxygen interaction and cell area. 
# It seems quite important also!!!
# ------------------------------------------------------------------------------
# Step 5: Adding to the best model from step 4 the interactive effect of traits 
# and test temperature and then compare their AICs
# ------------------------------------------------------------------------------
m2.1<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m2.2<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.number.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m2.3<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*fw.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m2.4<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*iod.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

fit.list.step5 <- list(m1.1,m2.1,m2.2,m2.3,m2.4)

fit.names.step5 <-c("ST * O2 + ST * sex * O2",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Cell area",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Cell number",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Fresh mass",
                    "ST * O2 + ST * sex * O2 + ST * O2 * IOD"
                    )

#compare by using AICc
fit.step5<-aictab(fit.list.step5,fit.names.step5, second.ord = T,sort = TRUE, digits = 3, LL=TRUE)
fit.step5

#                                                  K    AICc Delta_AICc AICcWt Cum.Wt    LL
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area   13 -107.07       0.00   0.90   0.90 68.87
# ST * O2 + ST * sex * O2                          9 -102.26       4.81   0.08   0.99 61.23
# ST * O2 + ST * sex * O2 + ST * O2 * IOD         13  -97.55       9.52   0.01   0.99 64.11
# ST * O2 + ST * sex * O2 + ST * O2 * Fresh mass  13  -96.25      10.82   0.00   1.00 63.46
# ST * O2 + ST * sex * O2 + ST * O2 * Cell number 13  -95.62      11.45   0.00   1.00 63.14

# ------------------------------------------------------------------------------
# Step 6: Adding to the best model from step 5 the interactive effect of traits 
# and test temperature and then compare their AICs
# ------------------------------------------------------------------------------
m3.1<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*fw.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m3.2<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*cell.number.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

m3.3<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*iod.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit,method = "ML")

fit.list.step6 <- list(m2.1,m3.1,m3.2,m3.3)


fit.names.step6 <-c("ST * O2 + ST * sex * O2 + ST * O2 * Cell area",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * FM",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * CN",
                    "ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * IOD")

#compare by using AICc
fit.step6<-aictab(fit.list.step6,fit.names.step6, second.ord = T,sort = TRUE, digits = 3, LL=TRUE)
fit.step6

#                                                               K    AICc Delta_AICc AICcWt Cum.Wt    LL
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area                 13 -107.07       0.00   0.86   0.86 68.87
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * IOD 17 -103.09       3.98   0.12   0.98 72.68
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * FM  17  -98.28       8.79   0.01   0.99 70.27
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area + ST * O2 * CN  17  -97.32       9.75   0.01   1.00 69.79

# model including sex and relative fresh mass works better, though adding 
# test temp and interocular distance or test temp and cell area is also important
# (DeltaAICc lower than 2).   
# ------------------------------------------------------------------------------
# Lets drop the best model

# (1) ST * O2 + ST * sex * O2 + ST * O2 * Cell area 
drop1(m2.1,test = "Chisq")
summary(m2.1)
anova(m2.1)

# Cell size interacts with stress temp and also with oxygen. The most significant
# effect os due to test.temp and test.temp and oxygen interaction. Similarly as in
# the experiment under normoxia , there is a significant effects of the interaction
# between sex and test temp. 
# ------------------------------------------------------------------------------
# Step 7: Lets now exclude the non significant variables
# ------------------------------------------------------------------------------
#exclude the three-ways interaction

m4.1<-gls(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex
          +scale(test.temp)*cell.area.rel
          +scale(test.oxygen)*cell.area.rel
          ,data=dhypo.median,na.action = na.omit,method = "ML")

drop1(m4.1,test = "Chisq")
summary(m4.1)
anova(m4.1)

# account only for the interactive effect
m4.2<-gls(surv.time~
            scale(test.temp):scale(test.oxygen)
          +scale(test.temp)*sex
          +scale(test.temp):cell.area.rel
          +scale(test.oxygen):cell.area.rel
          ,data=dhypo.median,na.action = na.omit,method = "ML")
drop1(m4.2,test = "Chisq")
summary(m4.2)
anova(m4.2)

fit.list.step7 <- list(m2.1,m4.1,m4.2)

fit.names.step7 <-c("ST * O2 + ST * sex * O2 + ST * O2 * Cell area",
                    "ST * O2 + ST * Sex + ST * Cell area + O2 * Cell area",
                    "ST : O2 + ST * Sex + ST : Cell area + O2 : Cell area")

#compare by using AICc
fit.step7<-aictab(fit.list.step7,fit.names.step7, second.ord = T,sort = TRUE, digits = 3, LL=TRUE)
fit.step7

#                                                       K    AICc Delta_AICc AICcWt Cum.Wt    LL
# ST : O2 + ST * Sex + ST : Cell area + O2 : Cell area  8 -113.42       0.00   0.58   0.58 65.58
# ST * O2 + ST * Sex + ST * Cell area + O2 * Cell area 10 -112.67       0.75   0.40   0.98 67.69
# ST * O2 + ST * sex * O2 + ST * O2 * Cell area        13 -107.07       6.35   0.02   1.00 68.87
#-------------------------------------------------------------------------------
# TABLE S10 SUPPLEMENTARY INFORMATION
#-------------------------------------------------------------------------------
m0<-lm(surv.time~
            scale(test.temp)
          ,data=dhypo.median,na.action = na.omit)
summary(m0)
m0.1<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m0.1)
  
m1.1<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*sex*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m1.1)  

m1.2<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*cell.area.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m1.2)
  
m1.3<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*cell.number.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m1.3)

m1.4<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*fw.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m1.4)  

m1.5<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*iod.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m1.5)  

m2.1<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*sex*scale(test.oxygen)
            +scale(test.temp)*cell.area.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m2.1)  

m2.2<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*sex*scale(test.oxygen)
            +scale(test.temp)*cell.number.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m2.2)  

m2.3<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*sex*scale(test.oxygen)
            +scale(test.temp)*fw.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m2.3)  

m2.4<-lm(surv.time~
              scale(test.temp)*scale(test.oxygen)
            +scale(test.temp)*sex*scale(test.oxygen)
            +scale(test.temp)*iod.rel*scale(test.oxygen)
            ,data=dhypo.median,na.action = na.omit)
summary(m2.4)  

m3.1<-lm(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*fw.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit)
summary(m3.1)

m3.2<-lm(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*cell.number.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit)
summary(m3.2)

m3.3<-lm(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex*scale(test.oxygen)
          +scale(test.temp)*cell.area.rel*scale(test.oxygen)
          +scale(test.temp)*iod.rel*scale(test.oxygen)
          ,data=dhypo.median,na.action = na.omit)
summary(m3.3)

m4.1<-lm(surv.time~
            scale(test.temp)*scale(test.oxygen)
          +scale(test.temp)*sex
          +scale(test.temp)*cell.area.rel
          +scale(test.oxygen)*cell.area.rel
          ,data=dhypo.median,na.action = na.omit)


fit.list.table.S10 <- list(m0,m0.1,m1.1,m1.2,m1.3,m1.4,m1.5,
                         m2.1,m2.2,m2.3,m2.4,
                         m3.1,m3.2,m3.3,
                         m4.1)

fit.names.table.S10 <-c("TT",
                      "TT * O2",  
                      "TT * O2 + TT * sex * O2",
                      "TT * O2 + TT * CA * O2",
                      "TT * O2 + TT * CN * O2",
                      "TT * O2 + TT * FM * O2",
                      "TT * O2 + TT * IOD * O2",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA",
                      "TT * O2 + TT * sex * O2 + TT * O2 * IOD",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * FM",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * CN",
                      "TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * IOD",
                      "TT * O2 + TT * sex + TT * CA + O2 * CA"
                      )

#compare by using AICc
fit.table.S10<-aictab(fit.list.table.S10,fit.names.table.S10, second.ord = T,sort = TRUE, digits = 3, LL=TRUE)
fit.table.S10
# Model selection based on AICc:
#   
#                                                         K    AICc Delta_AICc AICcWt Cum.Wt    LL
# TT * O2 + TT * sex + TT * CA + O2 * CA                 10 -112.67       0.00   0.92   0.92 67.69
# TT * O2 + TT * sex * O2 + TT * O2 * CA                 13 -107.07       5.60   0.06   0.98 68.87
# TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * IOD 17 -103.33       9.34   0.01   0.99 72.80
# TT * O2 + TT * sex * O2                                 9 -102.26      10.41   0.01   0.99 61.23
# TT * O2 + TT * CA * O2                                  9 -101.85      10.81   0.00   1.00 61.02
# TT * O2                                                 5  -98.32      14.35   0.00   1.00 54.51
# TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * FM  17  -98.28      14.39   0.00   1.00 70.27
# TT * O2 + TT * sex * O2 + TT * O2 * IOD                13  -97.87      14.79   0.00   1.00 64.27
# TT * O2 + TT * sex * O2 + TT * O2 * CA + TT * O2 * CN  17  -97.32      15.35   0.00   1.00 69.79
# TT * O2 + TT * sex * O2 + TT * O2 * CA                 13  -96.25      16.42   0.00   1.00 63.46
# TT * O2 + TT * sex * O2 + TT * O2 * CA                 13  -95.62      17.05   0.00   1.00 63.14
# TT * O2 + TT * IOD * O2                                 9  -94.08      18.59   0.00   1.00 57.14
# TT * O2 + TT * FM * O2                                  9  -92.58      20.08   0.00   1.00 56.39
# TT * O2 + TT * CN * O2                                  9  -92.47      20.20   0.00   1.00 56.33
# TT                                                      3  -83.93      28.73   0.00   1.00 45.10

#export table
write.csv(fit.table.S10,"../Outputs/6.1.1. Table S10 Model comparison of survival time at two kPa.csv",row.names = FALSE)

#-------------------------------------------------------------------------------
# TABLE 3 MAIN MANUSCRIPT
#-------------------------------------------------------------------------------
hyp_0 <- lm(surv.time ~ scale(test.temp), 
            data = dhypo.median, na.action = na.omit)
summary(hyp_0)

hyp_1 <- lm(surv.time ~ scale(test.temp) + 
              scale(test.temp) * scale(test.oxygen),
            data = dhypo.median, na.action = na.omit)
summary(hyp_1)


hyp_2 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * sex, 
            data = dhypo.median, na.action = na.omit)
summary(hyp_2)

hyp_3 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * cell.area.rel,
            data = dhypo.median, na.action = na.omit)
summary(hyp_3)

hyp_4 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * fw.rel,
            data = dhypo.median, na.action = na.omit)
summary(hyp_4)


hyp_5 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * sex +
              scale(test.temp) * scale(test.oxygen) * cell.area.rel,
            data = dhypo.median, na.action = na.omit)
summary(hyp_5)


hyp_6 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * sex +
              scale(test.temp) * scale(test.oxygen) * fw.rel,
            data = dhypo.median, na.action = na.omit)
summary(hyp_6)

hyp_7 <- lm(surv.time ~ scale(test.temp) +
              scale(test.temp) * scale(test.oxygen) * fw.rel +
              scale(test.temp) * scale(test.oxygen) * cell.area.rel,
            data = dhypo.median, na.action = na.omit)
summary(hyp_7)


hyp_full <- lm(surv.time ~ scale(test.temp) +
                 scale(test.temp)* scale(test.oxygen) * sex +
                 scale(test.temp) * fw.rel +
                 scale(test.temp) * cell.area.rel,
               data = dhypo.median, na.action = na.omit)
summary(hyp_full)


fit.list.table.3 <- list(hyp_0,
                         hyp_1,
                         hyp_2,
                         hyp_3,
                         hyp_4,
                         hyp_5,
                         hyp_6,
                         hyp_7,
                         hyp_full
)

fit.names.table.3 <-c(
  "TT",
  "TT * O2",
  "TT + TT * O2 * Sex",
  "TT + TT * O2 * CA",
  "TT + TT * O2 * FM",
  "TT + TT * O2 * Sex + TT * O2 * CA",
  "TT + TT * O2 * Sex + TT * O2 * FM",
  "TT + TT * O2 * CA  + TT * O2 * FM",
  "TT + TT * O2 * Sex + TT * O2 * CA + TT * O2 * FM"
)

#compare by using AICc
fit.table.3<-aictab(fit.list.table.3,fit.names.table.3, second.ord = T,sort = FALSE, digits = 3, LL=TRUE)
fit.table.3

#export table
write.csv(fit.table.3,"../Outputs/6.1.2. Table 3 Model comparison of survival time at two kPa.csv",row.names = FALSE)

#-------------------------------------------------------------------------------
# TABLE S12 SUPPLEMENTARY INFORMATION
#-------------------------------------------------------------------------------
summary(hyp_5)
#                                                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                        1.590748   0.018136  87.711  < 2e-16 ***
# scale(test.temp)                                  -0.404354   0.018269 -22.134  < 2e-16 ***
# scale(test.oxygen)                                 0.006814   0.018238   0.374  0.70967    
# sexmale                                           -0.054366   0.025616  -2.122  0.03690 *  
# cell.area.rel                                     -0.013113   0.012971  -1.011  0.31508    
# scale(test.temp):scale(test.oxygen)               -0.049758   0.018376  -2.708  0.00828 ** 
# scale(test.temp):sexmale                           0.078367   0.025756   3.043  0.00317 ** 
# scale(test.oxygen):sexmale                         0.029266   0.025756   1.136  0.25925    
# scale(test.temp):cell.area.rel                     0.037239   0.013372   2.785  0.00668 ** 
# scale(test.oxygen):cell.area.rel                   0.035293   0.013063   2.702  0.00842 ** 
# scale(test.temp):scale(test.oxygen):sexmale       -0.010193   0.025901  -0.394  0.69497    
# scale(test.temp):scale(test.oxygen):cell.area.rel -0.010710   0.013497  -0.794  0.42981 

anova(hyp_5)
#                                                   Df  Sum Sq Mean Sq  F value    Pr(>F)    
# scale(test.temp)                                   1 11.9232 11.9232 791.3415 < 2.2e-16 ***
# scale(test.oxygen)                                 1  0.0510  0.0510   3.3873  0.069407 .  
# sex                                                1  0.0666  0.0666   4.4225  0.038611 *  
# cell.area.rel                                      1  0.0132  0.0132   0.8731  0.352920    
# scale(test.temp):scale(test.oxygen)                1  0.3178  0.3178  21.0910 1.608e-05 ***
# scale(test.temp):sex                               1  0.1356  0.1356   8.9982  0.003602 ** 
# scale(test.oxygen):sex                             1  0.0193  0.0193   1.2823  0.260846    
# scale(test.temp):cell.area.rel                     1  0.0936  0.0936   6.2130  0.014747 *  
# scale(test.oxygen):cell.area.rel                   1  0.1063  0.1063   7.0544  0.009541 ** 
# scale(test.temp):scale(test.oxygen):sex            1  0.0023  0.0023   0.1556  0.694298    
# scale(test.temp):scale(test.oxygen):cell.area.rel  1  0.0095  0.0095   0.6297  0.429815    
# Residuals                                         80  1.2054  0.0151
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("../Outputs/6.1.3. TDT curves under two oxygen conditions_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################