# ------------------------------------------------------------
# Script to explore and calculate survival times of DGRP lines
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

getwd()# check directory
#-------------------------------------------------------------------------------
#Libraries
library(xlsx)
library(dplyr)
# ------------------------------------------------------------
# read data of survival time
datos <- read.xlsx("../Data/results on TDT available.xlsx", sheetName = "results on TDT available")

# Check data structure
str(datos)
# ------------------------------------------------------------
# convert variables
datos$NUMBER_STOCK        <- as.factor(datos$NUMBER_STOCK) 
datos$FLASK_NUMBER        <- as.factor(datos$FLASK_NUMBER) 
datos$STRESS_TEMP         <- as.numeric(datos$STRESS_TEMP) 
datos$RUN_NUMBER          <- as.factor(datos$RUN_NUMBER) 
datos$DATE_REV            <- as.factor(datos$DATE_REV)
datos$TIME_SURV_MIN       <- as.numeric(datos$TIME_SURV_MIN)
datos$TIME_SURV_SEC       <- as.numeric(datos$TIME_SURV_SEC)
datos$VIDEO_DURATION_MIN  <- as.numeric(datos$VIDEO_DURATION_MIN)
datos$VIDEO_DURATION_SEC  <- as.numeric(datos$VIDEO_DURATION_SEC)
datos$ACCLIM_TIME         <- as.numeric(datos$ACCLIM_TIME)
datos$TEST_OXYGEN         <- as.numeric(datos$TEST_OXYGEN)
datos$PLATE_NAME          <- as.factor(datos$PLATE_NAME)
datos$WELL                <- as.factor(datos$WELL)
datos$READER              <- as.factor(datos$READER)
datos$GENOTYPE            <- as.factor(datos$GENOTYPE)
datos$SEX                 <- as.factor(datos$SEX)
datos$VIDEO_NAME          <- as.factor(datos$VIDEO_NAME)
datos$DATE_REV            <- as.factor(datos$DATE_REV)
# ------------------------------------------------------------
# check again
str(datos)
head(datos)
# ------------------------------------------------------------
# convert time (survival and video duration) to decimal of minutes
datos$surv.time  <- log10(datos$TIME_SURV_MIN + (datos$TIME_SURV_SEC/60))
datos$video.time <- log10(datos$VIDEO_DURATION_MIN + (datos$VIDEO_DURATION_SEC/60))
# ------------------------------------------------------------
# select columns of interest
names(datos)
datos2           <- datos[c(2:7, 9, 10, 13, 18:28)]
# ------------------------------------------------------------
# Rename columns
names(datos2)
names(datos2)[c(1:20)] = c("video.name", "stock", "genotype", "sex", "test.temp",
                           "run", "date.rev", "id.ind", "criteria", "reader",
                           "acclim.temp", "date.exp", "age", "test.oxygen", 
                           "acclim.oxygen", "acclim.time", "plate.id", "well.id",
                           "surv.time", "video.time")
head(datos2)
# ------------------------------------------------------------
# Select only flies measured at 21 KPa
all.21<-filter(datos2, test.oxygen == 21)
# ------------------------------------------------------------
# Select only stocks used at two oxygen levels
stock25180 <- filter(datos2, stock == "25180")
stock25182 <- filter(datos2, stock == "25182")
stock28247 <- filter(datos2, stock == "28247")
stock28196 <- filter(datos2, stock == "28196")
stock25203 <- filter(datos2, stock == "25203")
stock25201 <- filter(datos2, stock == "25201")
# ------------------------------------------------------------
# merge the six stocks
all.hypo <- data.frame(rbind(stock25180, 
                             stock25182, 
                             stock28247, 
                             stock28196, 
                             stock25203, 
                             stock25201))
# ------------------------------------------------------------
# Transform variables wrongly assigned (21 kPa)
str(all.21)
all.21$date.exp <- as.factor(all.21$date.exp)
# ------------------------------------------------------------
# Transform variables wrongly assigned (10 kPa)
str(all.hypo)
all.hypo$date.rev  <- as.factor(all.hypo$date.rev)
all.hypo$date.exp  <- as.factor(all.hypo$date.exp)
all.hypo$age       <- as.factor(all.hypo$age)
# ------------------------------------------------------------
#median per reader at 21 kPa
all.21    <-      aggregate(cbind(surv.time, video.time) ~ 
                                stock + 
                                genotype +
                                video.name + 
                                sex +
                                test.temp + 
                                run + 
                                id.ind + 
                                acclim.temp +
                                date.exp +
                                age + 
                                test.oxygen +
                                acclim.oxygen +
                                acclim.time,
                  data=all.21, median)
#-------------------------------------------------------------
#median per reader at 10 kPa
all.hypo    <-    aggregate(cbind(surv.time, video.time) ~ 
                              stock +
                              genotype +
                              video.name + 
                              sex +
                              test.temp +
                              run +
                              id.ind +
                              acclim.temp +
                              date.exp +
                              age +
                              test.oxygen +
                              acclim.oxygen +
                              acclim.time,
                    data=all.hypo, median)
#-------------------------------------------------------------
#load data with mean of size traits
traits<-read.csv("../Outputs/2.1.1. Phenotypic traits data of DGRP lines.csv")
#average of size traits
traits.mean<-aggregate(cbind(ws,
                             fw,
                             cell.area,
                             cell.number,
                             iod,wing.area.z,
                             fresh.mass.z,
                             cell.area.z,
                             cell.number.z,
                             iod.z) ~ stock +
                         sex +
                         genotype,
                       data=traits, median)
#------------------------------------------------------------
#merge both databases for 21 kPa
all.lines.21 <- merge(all.21, traits.mean, by = c("stock","genotype","sex"), all.y = TRUE)
#------------------------------------------------------------
#Selecting stocks used at two oxygen levels
stock25180.hypo   <- filter(traits.mean, stock == "25180")
stock25182.hypo   <- filter(traits.mean, stock == "25182")
stock28247.hypo   <- filter(traits.mean, stock == "28247")
stock28196.hypo   <- filter(traits.mean, stock == "28196")
stock25203.hypo   <- filter(traits.mean, stock == "25203")
stock25201.hypo   <- filter(traits.mean, stock == "25201")

#merge data frames
traits.mean.hypo  <- data.frame(rbind(stock25180.hypo,
                                      stock25182.hypo,
                                      stock28247.hypo,
                                      stock28196.hypo,
                                      stock25203.hypo,
                                      stock25201.hypo))
#add traits data
all.lines.10       <- merge(all.hypo, traits.mean.hypo, by = c("stock", "sex"), all.y = TRUE)
#------------------------------------------------------------
#export data frames for future analyses
write.csv(all.lines.21, "../Outputs/4.1.1. Survival time of Drosophila melanogaster and its associated size traits at 21 kPa.csv",row.names=FALSE)
write.csv(all.lines.10, "../Outputs/4.1.2. Survival time of Drosophila melanogaster and its associated size traits at two oxygen levels.csv",row.names=FALSE)
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("../Outputs/4.1.3. Exploratory analysis TDT curves_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################