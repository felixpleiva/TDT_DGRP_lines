# ------------------------------------------------------------------------------
# Data generated by F�lix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# Script created by F�lix P Leiva, Radboud University (e-mail: felixpleiva@gmail.com)
# ------------------------------------------------------------------------------
# Cite as:

# Leiva FP, Santos M, Rezende E, & Verberk WCEP. (2021, July 21). Draft version
# of paper data and code of manuscript: Intraspecific variation on heat
# tolerance in a model ectotherm: effects of body mass, cell size, oxygen and
# sex (1.0). Zenodo. https://doi.org/10.5281/zenodo.5120029
# ------------------------------------------------------------------------------
# Cleaning working space
rm(list=ls())
today<-format(Sys.Date(),"%Y%m%d")
#-------------------------------------------------------------------------------
#get working directory
setwd("C:/Users/Invunche/Dropbox/GitHub/TDT_DGRP_lines/Data")
getwd()# check directory
#-------------------------------------------------------------------------------
#Libraries
library(xlsx)
#-------------------------------------------------------------------------------
#Load data
datos<-read.xlsx("results on CS available.xlsx", sheetName = "results on CS available")
#-------------------------------------------------------------------------------
#Select columns of interest
names(datos)
datos<-datos[c(1:5,7:10,12,13)]
#-------------------------------------------------------------------------------
#Rename columns
names(datos)
names(datos)[c(1:11)] = c("stock","genotype","sex","fw","ws","well",
                          "area","trichomes","iod","age","plate")
#-------------------------------------------------------------------------------
#Transform variables wrongly assigned
str(datos)
datos$stock<-as.factor(datos$stock) 
datos$fw<-as.numeric(datos$fw) 
datos$ws<-as.numeric(datos$ws) 
datos$trichomes<-as.numeric(datos$trichomes) 
datos$area<-as.numeric(datos$area) 
datos$plate<-as.factor(datos$plate)
datos$sex<-as.factor(datos$sex)
datos$genotype<-as.factor(datos$genotype)
datos$well<-as.factor(datos$well)

#Check again
str(datos)
#-------------------------------------------------------------------------------
#Average of the three measurements wing size
datos2<-aggregate(cbind(ws,fw,trichomes,area)~stock+genotype+sex+age+plate+well,data=datos,mean)

#Including IOD
datos3<-aggregate(cbind(ws,fw,trichomes,area,iod)~stock+genotype+sex+age+plate+well,data=datos,mean)
#-------------------------------------------------------------------------------
#Calculate cell area and cell number for each data frame
# datos2
datos2$cell.area<-datos2$area/datos2$trichomes
datos2$cell.number<-datos2$ws/datos2$cell.area

# datos3
datos3$cell.area<-datos3$area/datos3$trichomes
datos3$cell.number<-datos3$ws/datos3$cell.area

#number of data for sex and each stock 
xtabs(formula = ~sex+stock,data=datos3)
#-------------------------------------------------------------------------------
#Export both data frames. I doing this because there are six stocks without data
#of interocular distance. See scripts "2.0. Predictions to calculate IOD of
#remaining flies"
setwd("C:/Users/Invunche/Dropbox/GitHub/TDT_DGRP_lines/Outputs/")

write.csv(datos2,"1.1.1. Individual body traits (no iod) for all DGRP lines.csv",row.names = FALSE)
write.csv(datos3,"1.1.2. Individual body traits for 14 DGRP lines.csv",row.names = FALSE)
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("C:/Users/Invunche/Dropbox/GitHub/TDT_DGRP_lines/Outputs/1.1.3. Cell area for each DGRP line_R session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################