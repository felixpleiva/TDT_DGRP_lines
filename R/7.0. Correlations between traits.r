# -------------------------------------------------------------------------------------
# Script to plot correlations between traits measured in twenty DGRP lines
# Data generated Felix P Leiva, Radboud University
# Created by F�lix P Leiva on 20201203 (YYYYMMDD)
# # -------------------------------------------------------------------------------------
# Clean working space
rm(list=ls())
today<-format(Sys.Date(),"%Y%m%d")
# -------------------------------------------------------------------------------------
setwd("C:/Users/Invunche/Dropbox/Radboud University/publicaciones/Thesis/3. Thermal tolerance on big and small cells/manuscritos/Submission FE")
getwd()# check directory
# -------------------------------------------------------------------------------------
#Libraries
library(smatr)
library(lme4)
library(ggplot2)
library(cowplot)
# -------------------------------------------------------------------------------------
# load data
traits<-read.csv("Phenotypic traits data of DGRP lines.csv")
head(traits)
str(traits)
traits$stock<-as.factor(traits$stock)
#reorder the columns
traits<-traits[,c(1:6,8,9,7,11,10)]
names(traits)
traits.mean<-aggregate(cbind(fw,iod,ws,cell.number,cell.area)~stock+genotype+sex,traits,mean)
# -------------------------------------------------------------------------------------
# FIGURE S1
# -------------------------------------------------------------------------------------
{
# pdf("Figure S1 Pairwise correlation between traits (individual values, females and males pooled).pdf",width = 7,height = 7)
# png("Figure S1 Pairwise correlation between traits (individual values, females and males pooled).png",width = 7,height = 7,units = "in",res = 300)
par(mfrow=c(1,1),tcl=-0.4, family="serif",
      omi=c(0,0,0,0))
my_cols <- c('#4DAC2660', '#A6761D60')
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y,cex=1.2, pch=16,col = my_cols[traits$sex])
}
# Create the plots
pairs(traits[,7:11],cex.labels = 1.4,gap = 0,
      labels = c("Fresh mass","IOD","Wing area","Cell number","Cell area"),
            lower.panel = panel.cor,
      upper.panel = upper.panel)
#dev.off()
}

# -------------------------------------------------------------------------------------
# FIGURE S2
# -------------------------------------------------------------------------------------
{
# pdf("Figure S2 Pairwise correlation between between traits (mean values for females and males pooled).pdf",width = 7,height = 7)
# png("Figure S2 Pairwise correlation between between traits (mean values for females and males pooled).png",width = 7,height = 7,units = "in",res = 300)
par(mfrow=c(1,1),tcl=-0.4, family="serif",
      omi=c(0,0,0,0))
my_cols <- c('#4DAC2660', '#A6761D60')
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y,cex=1.2, pch=16,col = my_cols[traits.mean$sex])
}
  # Create the plots
pairs(traits.mean[,4:8],cex.labels = 1.4,gap = 0,
      labels = c("Fresh mass","IOD","Wing area","Cell number","Cell area"),
      lower.panel = panel.cor,
      upper.panel = upper.panel)
#dev.off()
}

# -------------------------------------------------------------------------------------
# FIGURE S3
# -------------------------------------------------------------------------------------
female<-subset(traits.mean,sex=="female")
female$stock <- reorder(female$stock, female$cell.area)

levels(female$stock)#to check if order is 0K
# "25182" "25198" "25180" "29652" "28191" "28248" "25190" "25192" "28247" "28197" 
# "28141" "28135" "28258" "28153" "28198" "28173" "25203" "25201" "28180" "28196"

#Setting same colour for each DGRP line

col.female1=col=c("#CC3333","#993366","#666699","#339999","#339966","#669966","#666666","#996699","#CC6666","#FF6633",
                  "#FF9900","#FFCC33","#FFFF33","#CC9933","#996633","#c00000","#CC6699","#FF99CC","#CC9999","#999999")

{
# pdf("Figure S3 Pairwise correlation between between traits (mean values for females).pdf",width = 7,height = 7)
# png("Figure S3 Pairwise correlation between between traits (mean values for females).png",width = 7,height = 7,units = "in",res = 300)
par(mfrow=c(1,1),tcl=-0.4, family="serif",
    omi=c(0,0,0,0))

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y,cex=1.2, pch=16,col = col.female1[female$stock])
}
# Create the plots
pairs(female[,4:8],cex.labels = 1.4,gap = 0,
      labels = c("Fresh mass","IOD","Wing area","Cell number","Cell area"),
      lower.panel = panel.cor,
      upper.panel = upper.panel)
#dev.off()
}

# -------------------------------------------------------------------------------------
# FIGURE S4
# -------------------------------------------------------------------------------------
male<-subset(traits.mean,sex!="female")
male$stock <- reorder(male$stock, male$cell.area)

levels(male$stock)#to check if order is 0K
# "25182" "25180" "25190" "29652" "25198" "28191" "28198" "28141" "28248" "28258" 
# "25192" "28153" "28247" "28197" "28135" "28173" "25201" "28180" "25203" "28196"

#Setting same colour for each DGRP line

col.male1=col=c("#CC3333","#666699","#666666","#339999","#993366","#339966","#996633","#FF9900","#669966","#FFFF33",
                "#996699","#CC9933","#CC6666","#FF6633","#FFCC33","#c00000","#FF99CC","#CC9999","#CC6699","#999999")

{
# pdf("Figure S4 Pairwise correlation between between traits (mean values for males).pdf",width = 7,height = 7)
# png("Figure S4 Pairwise correlation between between traits (mean values for males).png",width = 7,height = 7,units = "in",res = 300)
par(mfrow=c(1,1),tcl=-0.4, family="serif",
    omi=c(0,0,0,0))
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y,cex=1.2, pch=16,col = col.male1[male$stock])
}
# Create the plots
pairs(male[,4:8],cex.labels = 1.4,gap = 0,
      labels = c("Fresh mass","IOD","Wing area","Cell number","Cell area"),
      lower.panel = panel.cor,
      upper.panel = upper.panel)
#dev.off()
}
# -------------------------------------------------------------------------------------
# Regression analyses between body mass and IOD, cell area, cell number and wing area
# -------------------------------------------------------------------------------------
# Log transform traits
traits$fw.log<-as.numeric(log10(traits$fw))
traits$iod.log<-as.numeric(log10(traits$iod))
traits$ws.log<-as.numeric(log10(traits$ws))
traits$cn.log<-as.numeric(log10(traits$cell.number))
traits$ca.log<-as.numeric(log10(traits$cell.area))

# subset by sex
data_m<-subset(traits,sex=="male") #males
data_f<-subset(traits,sex!="male") #females

# We are going to fit two models. One with the slope of the relationship between 
# fresh mass and interocular distance changing with stock, and one with just the intercept 
# changing with stock. We are then going to test whether the former has a 
# significant improvement on the latter.

# IOD vs. fresh mass in females
# mixed model with random variation in intercept and slope
f1<-lmer(iod.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_f)
summary(f1)

# mixed model with random variation in intercept
f2<-lmer(iod.log ~ fw.log + (1|stock), REML = TRUE, data=data_f)
summary(f2)

#Comparing the models using AIC and BIC:
AIC(f1,f2)
BIC(f1,f2)
#Model 2 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(f2))[ , "Estimate"]
# (Intercept)      fw.log 
# 2.67320413  0.08840174


data_f$pred <- predict(f2,re.form=NA)
data_f$pred1 <- predict(f2)
f2 <- ggplot(data_f,aes(fw.log,iod.log))+
  ylim(c(2.55,2.75))+xlim(c(-0.5,0.2))+
  geom_point(colour="#4DAC2660",size=2)
f2.plot=
  f2 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#4DAC26",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Interocular distance, IOD ","(", mu,"m",")"))) + 
  ggtitle("Females")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))+
  annotate('text', x = -0.15, y = 2.75, 
           label = "Log[10]~IOD== 2.673~+~0.088~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
f2.plot

# Wing area vs. fresh mass in females
# mixed model with random variation in intercept and slope
f3<-lmer(ws.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_f)
summary(f3)

# mixed model with random variation in intercept
f4<-lmer(ws.log ~ fw.log + (1|stock), REML = TRUE, data=data_f)
summary(f4)

#Comparing the models using AIC and BIC:
AIC(f3,f4)
BIC(f3,f4)
#Model 4 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(f4))[ , "Estimate"]
# (Intercept)      fw.log 
# 6.0181936   0.1044719
data_f$pred <- predict(f4,re.form=NA)
data_f$pred1 <- predict(f4)
f4 <- ggplot(data_f,aes(fw.log,ws.log))+
  geom_point(colour="#4DAC2660",size=2)+
  ylim(c(5.75,6.1))+xlim(c(-0.5,0.2))
f4.plot=
  f4 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#4DAC26",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Wing area, WA ","(", mu,"m"^"2",")")))+
  annotate('text', x = -0.15, y = 6.1, 
           label = "Log[10]~WA== 6.018~+~0.104~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
f4.plot
# Cell number vs. fresh mass in females
# mixed model with random variation in intercept and slope
f5<-lmer(cn.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_f)
summary(f5)

# mixed model with random variation in intercept
f6<-lmer(cn.log ~ fw.log + (1|stock), REML = TRUE, data=data_f)
summary(f6)

#Comparing the models using AIC and BIC:
AIC(f5,f6)
BIC(f5,f6)
#Model 6 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(f6))[ , "Estimate"]
# (Intercept)      fw.log 
# 3.79821278  0.08973608
data_f$pred <- predict(f6,re.form=NA)
data_f$pred1 <- predict(f6)
f6 <- ggplot(data_f,aes(fw.log,cn.log))+
  geom_point(colour="#4DAC2660",size=2)+
  ylim(c(3.65,3.90))+xlim(c(-0.5,0.2))
f6.plot=
  f6 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#4DAC26",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + ylab("Cell number") + 
  annotate('text', x = -0.15, y = 3.9, 
           label = "Log[10]~CN== 3.798~+~0.090~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()


# Cell area vs. fresh mass in females
# mixed model with random variation in intercept and slope
f7<-lmer(ca.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_f)
summary(f7)

# mixed model with random variation in intercept
f8<-lmer(ca.log ~ fw.log + (1|stock), REML = TRUE, data=data_f)
summary(f8)

#Comparing the models using AIC and BIC:
AIC(f7,f8)
BIC(f7,f8)
#Model 8 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(f8))[ , "Estimate"]
# (Intercept)      fw.log 
# 2.220520893 0.006606646
data_f$pred <- predict(f8,re.form=NA)
data_f$pred1 <- predict(f8)
f8 <- ggplot(data_f,aes(fw.log,ca.log))+
  geom_point(colour="#4DAC2660",size=2)+
  ylim(c(2.05,2.30))+xlim(c(-0.5,0.2))
f8.plot=
f8 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#4DAC26",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Cell area, CA ","(", mu,"m"^"2",")")))+ 
  annotate('text', x = -0.15, y = 2.3, 
           label = "Log[10]~CA== 2.221~+~0.007~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()

# IOD vs. fresh mass in males
# mixed model with random variation in intercept and slope
m1<-lmer(iod.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_m)
summary(m1)

# mixed model with random variation in intercept
m2<-lmer(iod.log ~ fw.log + (1|stock), REML = TRUE, data=data_m)
summary(m2)

#Comparing the models using AIC and BIC:
AIC(m1,m2)
BIC(m1,m2)
#Model 1 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(m1))[ , "Estimate"]
# (Intercept)      fw.log 
# 2.65325474  0.05571345
# 
data_m$pred <- predict(m1,re.form=NA)
data_m$pred1 <- predict(m1)
m1<- ggplot(data_m,aes(fw.log,iod.log))+
  geom_point(colour="#A6761D60",size=2)+
  ylim(c(2.55,2.75))+xlim(c(-0.5,0.2))
m1.plot=
m1 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#A6761D",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Interocular distance, IOD ","(", mu,"m",")"))) +
  ggtitle("Males")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))+
  annotate('text', x = -0.15, y = 2.75, 
           label = "Log[10]~IOD== 2.653~+~0.056~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
m1.plot
# Wing area vs. fresh mass in males

# mixed model with random variation in intercept and slope
m3<-lmer(ws.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_m)
summary(m3)

# mixed model with random variation in intercept
m4<-lmer(ws.log ~ fw.log + (1|stock), REML = TRUE, data=data_m)
summary(m4)

#Comparing the models using AIC and BIC:
AIC(m3,m4)
BIC(m3,m4)
#Model 3 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(m3))[ , "Estimate"]
# (Intercept)      fw.log 
#  5.9438629   0.1763224
# 
data_m$pred <- predict(m3,re.form=NA)
data_m$pred1 <- predict(m3)
m3<- ggplot(data_m,aes(fw.log,ws.log))+
  geom_point(colour="#A6761D60",size=2)+
  ylim(c(5.75,6.10))+xlim(c(-0.5,0.2))
m3.plot=
m3 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#A6761D",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Wing area, WA ","(", mu,"m"^"2",")")))+ 
  annotate('text', x = -0.15, y = 6.1, 
           label = "Log[10]~WA== 5.944~+~0.176~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
m3.plot
# Cell number vs. fresh mass in males
# 
# mixed model with random variation in intercept and slope
m5<-lmer(cn.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_m)
summary(m5)

# mixed model with random variation in intercept
m6<-lmer(cn.log ~ fw.log + (1|stock), REML = TRUE, data=data_m)
summary(m6)

#Comparing the models using AIC and BIC:
AIC(m5,m6)
BIC(m5,m6)
#Model 5 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(m5))[ , "Estimate"]
# (Intercept)      fw.log 
# 3.77642423  0.09811069
data_m$pred <- predict(m5,re.form=NA)
data_m$pred1 <- predict(m5)
m5<- ggplot(data_m,aes(fw.log,cn.log))+
  geom_point(colour="#A6761D60",size=2)+
  ylim(c(3.65,3.90))+xlim(c(-0.5,0.2))
m5.plot=
m5 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#A6761D",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + ylab("Cell number") + 
  annotate('text', x = -0.15, y = 3.9, 
           label = "Log[10]~CN== 3.776~+~0.098~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
m5.plot

# Cell area vs. fresh mass in males
# mixed model with random variation in intercept and slope
m7<-lmer(ca.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=data_m)
summary(m7)

# mixed model with random variation in intercept
m8<-lmer(ca.log ~ fw.log + (1|stock), REML = TRUE, data=data_m)
summary(m8)

#Comparing the models using AIC and BIC:
AIC(m7,m8)
BIC(m7,m8)
#Model 8 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(m8))[ , "Estimate"]
# (Intercept)      fw.log 
#  2.16620106  0.06477081
data_m$pred <- predict(m8,re.form=NA)
data_m$pred1 <- predict(m8)
m8<- ggplot(data_m,aes(fw.log,ca.log))+
  geom_point(colour="#A6761D60",size=2)+
  ylim(c(2.05,2.30))+xlim(c(-0.5,0.2))
m8.plot=
m8 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#A6761D",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Cell area, CA ","(", mu,"m"^"2",")")))+ 
  annotate('text', x = -0.15, y = 2.3, 
           label = "Log[10]~CA== 2.166~+~0.065~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
m8.plot

# IOD vs. fresh mass in all
# mixed model with random variation in intercept and slope
a1<-lmer(iod.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=traits)
summary(a1)

# mixed model with random variation in intercept
a2<-lmer(iod.log ~ fw.log + (1|stock), REML = TRUE, data=traits)
summary(a2)

#Comparing the models using AIC and BIC:
AIC(a1,a2)
BIC(a1,a2)
#Model 1 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(a1))[ , "Estimate"]
# (Intercept)      fw.log 
# 2.6662195   0.1615005
traits$pred <- predict(a1,re.form=NA)
traits$pred1 <- predict(a1)
a1<- ggplot(traits,aes(fw.log,iod.log))+
  geom_point(colour="#55555540",size=2)+
  ylim(c(2.55,2.75))+xlim(c(-0.5,0.2))
a1.plot=
a1 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#555555",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Interocular distance, IOD ","(", mu,"m",")"))) +
  ggtitle("Sex pooled")+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate('text', x = -0.15, y = 2.75, 
           label = "Log[10]~IOD== 2.666~+~0.162~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()
a1.plot

# Wing area vs. fresh mass in males

# mixed model with random variation in intercept and slope
a3<-lmer(ws.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=traits)
summary(a3)

# mixed model with random variation in intercept
a4<-lmer(ws.log ~ fw.log + (1|stock), REML = TRUE, data=traits)
summary(a4)

#Comparing the models using AIC and BIC:
AIC(a3,a4)
BIC(a3,a4)
#Model 3 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(a3))[ , "Estimate"]
# (Intercept)      fw.log 
#  5.9871142   0.4773002

traits$pred <- predict(a3,re.form=NA)
traits$pred1 <- predict(a3)
a3<- ggplot(traits,aes(fw.log,ws.log))+
  geom_point(colour="#55555540",size=2)
ylim(c(5.75,6.10))+xlim(c(-0.5,0.2))
a3.plot=
  a3 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#555555",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Wing area, WA ","(", mu,"m"^"2",")")))+ 
  annotate('text', x = -0.15, y = 6.1, 
           label = "Log[10]~WA== 5.987~+~0.477~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()

a3.plot

# Cell number vs. fresh mass in males
# 
# mixed model with random variation in intercept and slope
a5<-lmer(cn.log ~ fw.log + (1+fw|stock), REML = TRUE, data=traits)
summary(a5)

# mixed model with random variation in intercept
a6<-lmer(cn.log ~ fw.log + (1|stock), REML = TRUE, data=traits)
summary(a6)

#Comparing the models using AIC and BIC:
AIC(a5,a6)
BIC(a5,a6)
#Model 5 (variation in intercept and slope between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(a5))[ , "Estimate"]
# (Intercept)      fw.log 
# 3.7890481   0.1849083
traits$pred <- predict(a5,re.form=NA)
traits$pred1 <- predict(a5)
a5<- ggplot(traits,aes(fw.log,cn.log))+
  geom_point(colour="#55555540",size=2)+
  ylim(c(3.65,3.90))+xlim(c(-0.5,0.2))
a5.plot=
  a5 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#555555",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + ylab("Cell number") + 
  annotate('text', x = -0.15, y = 3.9, 
           label = "Log[10]~CN== 3.789~+~0.185~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()

a5.plot

# Cell area vs. fresh mass in all
# mixed model with random variation in intercept and slope
a7<-lmer(ca.log ~ fw.log + (1+fw.log|stock), REML = TRUE, data=traits)
summary(a7)

# mixed model with random variation in intercept
a8<-lmer(ca.log ~ fw.log + (1|stock), REML = TRUE, data=traits)
summary(a8)

#Comparing the models using AIC and BIC:
AIC(a7,a8)
BIC(a7,a8)
#Model 7 (variation in intercept between stock) is a better fit both by AIC and BIC

#PLOT
coef(summary(a7))[ , "Estimate"]
# (Intercept)      fw.log 
# 2.1980927   0.2883845
traits$pred <- predict(a7,re.form=NA)
traits$pred1 <- predict(a7)
a7<- ggplot(traits,aes(fw.log,ca.log))+
  geom_point(colour="#55555540",size=2)+
  ylim(c(2.05,2.30))+xlim(c(-0.5,0.2))
a7.plot=
  a7 +   geom_line(colour="gray",aes(y=pred1,group=stock)) +
  geom_line(colour="#555555",lwd=2,aes(y=pred,group=stock))+
  xlab("Fresh mass, FM (mg)") + 
  ylab(expression(paste("Cell area, CA ","(", mu,"m"^"2",")")))+
  annotate('text', x = -0.15, y = 2.3, 
           label = "Log[10]~CA== 2.198~+~0.288~Log[10]~FM",parse = TRUE,size=4)+
  theme_bw()

a7.plot


#Figure S5
Figure_S5<-plot_grid(f2.plot,m1.plot,a1.plot,
              f4.plot,m3.plot,a3.plot,
              f6.plot,m5.plot,a5.plot,
              f8.plot,m8.plot,a7.plot,
              labels=c("A","B","C",
                       "D","E","F",
                       "G","H","I",
                       "J","K","L")
              ,nrow = 4,ncol = 3,label_size=15)

#Store Plots
ggsave('Figure S5.pdf',Figure_S5,width=12,height=15)
ggsave('Figure S5.png',Figure_S5,width=14,height=17)
