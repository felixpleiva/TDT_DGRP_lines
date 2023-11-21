# ------------------------------------------------------------------------------
# Data generated by Felix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# Script created by Felix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# ------------------------------------------------------------------------------
# Cite as:

# Leiva FP, Santos M, Rezende EL & Verberk WCEP. (2021). Paper data and code of
# manuscript: Intraspecific variation on heat tolerance in a model ectotherm:
# the role of oxygen, cell size and body size. Zenodo.
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
#-------------------------------------------------------------------------------
# Libraries
library(xlsx)
#-------------------------------------------------------------------------------
# create a folder for the outputs produced by running this script (if it doesn't
# already exist)
if (!file.exists("../Outputs")) dir.create("../Outputs")
# '../' goes up one level from the current working directory, so this creates
# the 'Outputs' folder just outside the 'R' folder
# ------------------------------------------------------------------------------
# Load data of body size related traits measured on flies
datos<-read.xlsx("../Data/results on CS available.xlsx", sheetName = "results on CS available")
#-------------------------------------------------------------------------------
# Select columns of interest
names(datos)
datos <- datos[c(1:5,7:10,12,13)]
#-------------------------------------------------------------------------------
#Rename columns
names(datos)
names(datos)[c(1:11)] = c("stock", "genotype", "sex", "fw", "ws", "well",
                          "area", "trichomes", "iod", "age", "plate")
#-------------------------------------------------------------------------------
# Transform variables wrongly assigned
str(datos)
datos$stock     <- as.factor(datos$stock) 
datos$fw        <- as.numeric(datos$fw) 
datos$ws        <- as.numeric(datos$ws) 
datos$trichomes <- as.numeric(datos$trichomes) 
datos$area      <- as.numeric(datos$area) 
datos$plate     <- as.factor(datos$plate)
datos$sex       <- as.factor(datos$sex)
datos$genotype  <- as.factor(datos$genotype)
datos$well      <- as.factor(datos$well)

# Check again
str(datos)
#-------------------------------------------------------------------------------
# Average of the three measurements wing size
datos2 <- aggregate(cbind(ws, fw, trichomes, area) ~ 
                      stock + 
                      genotype + 
                      sex + 
                      age + 
                      plate + 
                      well, 
                    data=datos,mean)
# The data frame "datos2" has 20 stocks
length(unique(datos2$stock))
#-------------------------------------------------------------------------------
# Average of the three measurements for all the traits, including IOD
datos3 <- aggregate(cbind(ws, fw, trichomes, area, iod) ~ 
                      stock + 
                      genotype + 
                      sex + 
                      age + 
                      plate + 
                      well, 
                    data=datos,mean)

# The data frame datos3 has 14 stocks
length(unique(datos3$stock))
#-------------------------------------------------------------------------------
#Calculate cell area and cell number for each data frame

# datos2
datos2$cell.area   <- datos2$area/datos2$trichomes # cell area
datos2$cell.number <- datos2$ws/datos2$cell.area   # cell number

# datos3
datos3$cell.area   <- datos3$area/datos3$trichomes # cell area
datos3$cell.number <- datos3$ws/datos3$cell.area   # cell number

#number of data for sex and each stock 
xtabs(formula = ~ sex + stock, data = datos3)
#-------------------------------------------------------------------------------
# Export both data frames. I am doing this because there are six stocks without data
# of interocular distance. See scripts "2.0. Predictions to calculate IOD of
# remaining flies"
# 
setwd("../Outputs/")

write.csv(datos2, "1.1.1. Individual body traits (no iod) for all DGRP lines.csv",row.names = FALSE)
write.csv(datos3, "1.1.2. Individual body traits for 14 DGRP lines.csv",row.names = FALSE)
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("../Outputs/1.1.3. Cell area for each DGRP line_R session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
