# ------------------------------------------------------------
# Script to explore and calculate survival time of DGRP lines 
# Data generated by Felix P Leiva, Radboud University
# Created by Felix P Leiva
# ------------------------------------------------------------
# Cleaning working space
rm(list=ls()) 
today<-format(Sys.Date(),"%Y%m%d")
# ------------------------------------------------------------
setwd("C:/Users/Invunche/Dropbox/Radboud University/publicaciones/Thesis/3. Thermal tolerance on big and small cells/datos/TDT")
getwd()# check directory
# ------------------------------------------------------------
#Libraries
library(xlsx)
library(dplyr)
# ------------------------------------------------------------
#read data of survival time
dat<-read.xlsx("results on TDT available.xlsx", sheetName = "results on TDT available")
datos<-filter(dat,INCLUDE!="0")#Exclude flies died for another reasons (e.g. water leaked into vials during the trials)

#Check
str(datos)
# ------------------------------------------------------------
# convert variables
datos$NUMBER_STOCK<-as.factor(datos$NUMBER_STOCK) 
datos$FLASK_NUMBER<-as.factor(datos$FLASK_NUMBER) 
datos$STRESS_TEMP<-as.numeric(datos$STRESS_TEMP) 
datos$RUN_NUMBER<-as.factor(datos$RUN_NUMBER) 
datos$DATE_REV<-as.factor(datos$DATE_REV)
datos$TIME_SURV_MIN<-as.numeric(datos$TIME_SURV_MIN)
datos$TIME_SURV_SEC<-as.numeric(datos$TIME_SURV_SEC)
datos$VIDEO_DURATION_MIN<-as.numeric(datos$VIDEO_DURATION_MIN)
datos$VIDEO_DURATION_SEC<-as.numeric(datos$VIDEO_DURATION_SEC)
datos$ACCLIM_TIME<-as.numeric(datos$ACCLIM_TIME)
datos$TEST_OXYGEN<-as.numeric(datos$TEST_OXYGEN)
# ------------------------------------------------------------
# check again
str(datos)
head(datos)
# ------------------------------------------------------------
# convert survival time to decimal of minutes
datos$surv.time<-log10(datos$TIME_SURV_MIN+(datos$TIME_SURV_SEC/60))
datos$video.time<-log10(datos$VIDEO_DURATION_MIN+(datos$VIDEO_DURATION_SEC/60))
# ------------------------------------------------------------
#select columns of interest
names(datos)
datos2<-datos[c(2:7,9,10,13,18:28)]
# ------------------------------------------------------------
#Rename columns
names(datos2)
names(datos2)[c(1:20)] = c("video.name","stock","genotype","sex","test.temp","run","date.rev",
                           "id.ind","criteria","reader","acclim.temp","date.exp",
                           "age","test.oxygen","acclim.oxygen","acclim.time","plate.id","well.id","surv.time","video.time")
head(datos2)
# ------------------------------------------------------------
#Select only flies measured at 21 KPa
all.21<-filter(datos2,test.oxygen==21)
# ------------------------------------------------------------
#Selecting only stocks used at two oxygen levels
stock25180<-filter(datos2,stock=="25180")
stock25182<-filter(datos2,stock=="25182")
stock28247<-filter(datos2,stock=="28247")
stock28196<-filter(datos2,stock=="28196")
stock25203<-filter(datos2,stock=="25203")
stock25201<-filter(datos2,stock=="25201")
# ------------------------------------------------------------
#merging the six stocks
all.hypo<-data.frame(rbind(stock25180,stock25182,stock28247,stock28196,stock25203,stock25201))
# ------------------------------------------------------------
# Transform variables wrongly assigned (21 kPa)
all.21$stock<-as.factor(all.21$stock) 
all.21$test.temp<-as.numeric(all.21$test.temp) 
all.21$run<-as.factor(all.21$run) 
all.21$date.rev<-as.factor(all.21$date.rev)
all.21$id.ind<-as.factor(all.21$id.ind)
all.21$criteria<-as.factor(all.21$criteria)
all.21$acclim.temp<-as.factor(all.21$acclim.temp)
all.21$date.exp<-as.factor(all.21$date.exp)
all.21$age<-as.factor(all.21$age)
all.21$test.oxygen<-as.factor(all.21$test.oxygen)
all.21$acclim.oxygen<-as.factor(all.21$acclim.oxygen)
all.21$acclim.time<-as.factor(all.21$acclim.time)
# ------------------------------------------------------------
# Transform variables wrongly assigned (10 kPa)
all.hypo$stock<-as.factor(all.hypo$stock)
all.hypo$test.temp<-as.numeric(all.hypo$test.temp) 
all.hypo$run<-as.factor(all.hypo$run) 
all.hypo$date.rev<-as.factor(all.hypo$date.rev)
all.hypo$id.ind<-as.factor(all.hypo$id.ind)
all.hypo$criteria<-as.factor(all.hypo$criteria)
all.hypo$acclim.temp<-as.factor(all.hypo$acclim.temp)
all.hypo$date.exp<-as.factor(all.hypo$date.exp)
all.hypo$age<-as.factor(all.hypo$age)
all.hypo$test.oxygen<-as.factor(all.hypo$test.oxygen)
all.hypo$acclim.oxygen<-as.factor(all.hypo$acclim.oxygen)
all.hypo$acclim.time<-as.factor(all.hypo$acclim.time)
# ------------------------------------------------------------
#average per reader at 21 kPa
all.21<-aggregate(cbind(surv.time,video.time)~stock+genotype+
                    video.name+sex+test.temp+run+id.ind+acclim.temp+
                    date.exp+age+test.oxygen+acclim.oxygen+acclim.time,
                  data=all.21,mean)#na.rm=TRUE calcula medias excluyendo los "NA"
#-------------------------------------------------------------
#average per reader at 10 kPa
all.hypo<-aggregate(cbind(surv.time,video.time)~stock+genotype+
                      video.name+sex+test.temp+run+id.ind+acclim.temp+
                      date.exp+age+test.oxygen+acclim.oxygen+acclim.time,
                    data=all.hypo,mean)#na.rm=TRUE calcula medias excluyendo los "NA"
#-------------------------------------------------------------
#load data with mean of size traits
setwd("C:/Users/Invunche/Dropbox/Radboud University/publicaciones/Thesis/3. Thermal tolerance on big and small cells/manuscritos/Submission FE")
traits<-read.csv("Phenotypic traits data of DGRP lines.csv")
#average of size traits
traits.mean<-aggregate(cbind(ws,fw,cell.area,cell.number,iod,wing.area.z,fresh.mass.z,cell.area.z,
                             cell.number.z,iod.z)~stock+sex+genotype,data=traits,median)
#------------------------------------------------------------
#merge both databases for 21 kPa
all.lines.21<-merge(all.21,traits.mean,by=c("stock","genotype","sex"),all.y = TRUE)
#------------------------------------------------------------

#Selecting stocks used at two oxygen levels
stock25180.hypo<-filter(traits.mean,stock=="25180")
stock25182.hypo<-filter(traits.mean,stock=="25182")
stock28247.hypo<-filter(traits.mean,stock=="28247")
stock28196.hypo<-filter(traits.mean,stock=="28196")
stock25203.hypo<-filter(traits.mean,stock=="25203")
stock25201.hypo<-filter(traits.mean,stock=="25201")

#merge data frames
traits.mean.hypo<-data.frame(rbind(stock25180.hypo,stock25182.hypo,stock28247.hypo,
                                   stock28196.hypo,stock25203.hypo,stock25201.hypo))
#add traits data
all.lines.10<-merge(all.hypo,traits.mean.hypo,by=c("stock","sex"),all.y = TRUE)
#------------------------------------------------------------
write.csv(all.lines.21, "Survival time of Drosophila melanogaster and its associated size traits at 21 kPa.csv",row.names=FALSE)
write.csv(all.lines.10, "Survival time of Drosophila melanogaster and its associated size traits at two oxygen levels.csv",row.names=FALSE)
#------------------------------------------------------------
#############################################################
######################END OF SCRIPT##########################
#############################################################