# ------------------------------------------------------------------------------
# Script to predict interocular distance for remaining DGRP lines
# ------------------------------------------------------------------------------
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

# check directory
getwd()
# ------------------------------------------------------------------------------
#Libraries
library(lme4)
library(PerformanceAnalytics)
#-------------------------------------------------------------------------------
# load data
datos.20  <- read.csv("../Outputs/1.1.1. Individual body traits (no iod) for all DGRP lines.csv")
datos.14  <- read.csv("../Outputs/1.1.2. Individual body traits for 14 DGRP lines.csv")
datos.14  <- subset(datos.14,datos.14$stock!=25203) # problem with the scale of the stereomicroscope
#-------------------------------------------------------------------------------
# Select columns of interest

names(datos.20)
datos.20 <- datos.20[c(1:8,11,12)]

names(datos.14)
datos.14 <- datos.14[c(1:8,11:13)]
#-------------------------------------------------------------------------------
# Transforming variables wrongly assigned

str(datos.20)

datos.20$stock     <- as.factor(datos.20$stock)
datos.20$genotype  <- as.factor(datos.20$genotype)
datos.20$age       <- as.numeric(datos.20$age)
datos.20$plate     <- as.factor(datos.20$plate)
datos.20$sex       <- as.factor(datos.20$sex)
datos.20$well      <- as.factor(datos.20$well)

str(datos.14)

datos.14$stock     <- as.factor(datos.14$stock)
datos.14$genotype  <- as.factor(datos.14$genotype)
datos.14$age       <- as.numeric(datos.14$age)
datos.14$plate     <- as.factor(datos.14$plate)
datos.14$well      <- as.factor(datos.14$well)
datos.14$sex       <- as.factor(datos.14$sex)

# check again
str(datos.20)

str(datos.14)
#-------------------------------------------------------------------------------
# Predicting IOD for the remaining seven stocks. For this, I will explore what traits
# shows high correlations with the observed data of iod

# subset by sex

male.20    <- subset(datos.20,datos.14$sex=="male")
female.20  <- subset(datos.20,datos.14$sex!="male")

male.14    <- subset(datos.14,datos.14$sex=="male")
female.14  <- subset(datos.14,datos.14$sex!="male")

# subtract data frames for males
male.rem   <- male.20[ !(male.20$stock %in% male.14$stock), ]

# subtract data frames for females
female.rem <- female.20[ !(female.20$stock %in% female.14$stock), ]

# Plotting correlations for males
x11()
chart.Correlation(male.14[7:11])

# Plotting correlations for females
x11()
chart.Correlation(female.14[7:11])

#-------------------------------------------------------------------------------
# I will use fresh mass to predict interocular distance because of the high
# correlation coefficient with IOD

# Modelling with random effects (stock)
 
set.seed(6955)

# females
fit.female   <- lmer(iod ~ fw + (1|stock), data = female.14)
coef(summary(fit.female))[ , "Estimate"]
# (Intercept)         fw 
# 421.23355           48.48572

# males
fit.male     <- lmer(iod ~ fw + (1|stock), data = male.14)
coef(summary(fit.male))[ , "Estimate"]
# (Intercept)         fw 
# 387.88406           65.95476

# Predict interocular distance (IOD) for females
names(female.rem)
female.rem$iod<-predict(fit.female,newdata=female.rem,re.form=NA)

# Predict iod for males
names(male.rem)
male.rem$iod<-predict(fit.male,newdata=male.rem,re.form=NA)

#merge data frames
all_traits  <- rbind(datos.14, female.rem, male.rem)
xtabs(formula = ~ stock + sex, data = all_traits)

#-------------------------------------------------------------------------------
# Rescale all variables (z-scores, just to visualize in a better way the data) by
# sex

female  <- subset(all_traits, all_traits$sex == "female") # females
male    <- subset(all_traits, all_traits$sex == "male") # males

#female
female$fresh.mass.z   <-scale(female$fw)
female$iod.z          <-scale(female$iod)
female$wing.area.z    <-scale(female$ws)
female$cell.number.z  <-scale(female$cell.number)
female$cell.area.z    <-scale(female$cell.area)

# male
male$fresh.mass.z     <-scale(male$fw)
male$iod.z            <-scale(male$iod)
male$wing.area.z      <-scale(male$ws)
male$cell.number.z    <-scale(male$cell.number)
male$cell.area.z      <-scale(male$cell.area)
#-------------------------------------------------------------------------------
# merge females and males again
all_traits <- rbind(female, male)

# Export data frame
write.csv(all_traits, "../Outputs/2.1.1. Phenotypic traits data of DGRP lines.csv",row.names = FALSE)
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes

sink("../Outputs/2.1.2. Predictions to calculate IOD of remaining flies_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################