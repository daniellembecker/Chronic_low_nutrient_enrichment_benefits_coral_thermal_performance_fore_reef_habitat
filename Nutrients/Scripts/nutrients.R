#Code to compare nutrients
#clear list 

rm(list=ls())

#load libraries 

library(ggpubr)
library(reshape)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(moments)
library(ghibli)
library(patchwork)
library(car)
library(pastecs)
library(psych)
library(graphics)
library(vegan)
library(here)
library(lme4)
library(lmerTest)

here()

##########################################################
#macroalgal CHN data 
#load calculated macroalgal N data
mydata <- read.csv("Nutrients/Data/macroalgal_data.csv")
View(mydata)

#need to organize data so it recognizes seperate sites

Enriched.site <- mydata %>%
  filter(treatment == "Enriched")


#need to organize data so it recognizes seperate sites

Control.site <- mydata %>%
  filter(treatment == "Control")

#analysis for site x N and graphs 
#check for normality, use normality plots

qqp(Control.site$N, "norm")


#check for normality, use normality plots


qqp(Enriched.site$N, "norm")


#check variance

SummaryByGroup <- mydata %>%
  group_by(treatment) %>%
  summarize(variance=var(N, na.rm=TRUE))
SummaryByGroup

boxplot(N~treatment, data=mydata, ylab= "N")

bartlett.test(N~treatment, data=mydata)

#p>0.05 means variances are equal, use students t test

data.summary<-mydata %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(N), se=sd(N)/sqrt(n())) #calculates mean 
data.summary

#use two-sample student t-test to test 

compare_means (N~treatment, data = mydata, method = "t.test")

a <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +  
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold("% Nitrogen Content"))) +
  theme(legend.position = "none") 

ggsave(filename = "Nutrients/Output/N.boxplot.png", device = "png", width = 10, height = 12)

###########################################################################
#water column nutrient data for just the in situ sites, not tanks
#load water column nutrient data
  watercol.dat <- read.csv("Nutrients/Data/watercolumn.nutrients.Oct.csv")
View(watercol.dat)


#need to organize data so it recognizes seperate sites

enrich.site <- watercol.dat %>%
  filter(treatment == "Enriched")


#need to organize data so it recognizes seperate sites

cont.site <- watercol.dat %>%
  filter(treatment == "Control")


#filter to just compare in situ values not tanks
watercol.dat.insitu <- watercol.dat %>%
  filter(treatment %in% c("Enriched", "Control"))

#analysis for site x N.N, P and NH4graphs 
#check for normality, use normality plots

qqp(cont.site$N.N, "norm")


#check for normality, use normality plots


qqp(enrich.site$N.N, "norm")


#analysis for site x N.N, P and NH4graphs 
#check for normality, use normality plots

qqp(cont.site$P, "norm")


#check for normality, use normality plots


qqp(enrich.site$P, "norm")


#analysis for site x N.N, P and NH4graphs 
#check for normality, use normality plots

qqp(cont.site$NH4, "norm")


#check for normality, use normality plots


qqp(enrich.site$NH4, "norm")


#check variance and make plot for N.N comparisons
SummaryByGroup <- watercol.dat.insitu %>%
  group_by(treatment) %>%
  summarize(variance=var(N.N, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.insitu %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(N.N), se=sd(N.N)/sqrt(n())) #calculates mean 
data.summary


boxplot(N.N~treatment, data=watercol.dat.insitu, ylab= "Nitrate")

bartlett.test(N.N~treatment, data=watercol.dat.insitu)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (N.N~treatment, data = watercol.dat.insitu, method = "t.test")

b <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Nitrate~(NO["3"]^~~{"-"})~+~Nitrite~(NO["2"]^~~{"-"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/N.N.boxplot.png", device = "png", width = 10, height = 12)

#check variance and make plot for P comparisons
SummaryByGroup <- watercol.dat.insitu %>%
  group_by(treatment) %>%
  summarize(variance=var(P, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.insitu %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(P), se=sd(P)/sqrt(n())) #calculates mean 
data.summary

boxplot(P~treatment, data=watercol.dat.insitu, ylab= "P")

bartlett.test(P~treatment, data=watercol.dat.insitu)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (P~treatment, data = watercol.dat.insitu, method = "t.test")

c <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Phosphate~(PO["4"]^~~{"3-"})~(mu*mol~L^{-1}))))  +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/P.boxplot.png", device = "png", width = 10, height = 12)

#check variance and make plot for NH4 comparisons
SummaryByGroup <- watercol.dat.insitu %>%
  group_by(treatment) %>%
  summarize(variance=var(NH4, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.insitu %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(NH4), se=sd(NH4)/sqrt(n())) #calculates mean 
data.summary

boxplot(NH4~treatment, data=watercol.dat.insitu, ylab= "NH4")

bartlett.test(NH4~treatment, data=watercol.dat.insitu)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (NH4~treatment, data = watercol.dat.insitu, method = "t.test")

d <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Ammonium~(NH["4"]^~~{"+"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/NH4.boxplot.png", device = "png", width = 10, height = 12)

#check variance and make plot for NH4 comparisons
SummaryByGroup <- watercol.dat.insitu %>%
  group_by(treatment) %>%
  summarize(variance=var(DIN.DIP, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.insitu %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(DIN.DIP), se=sd(DIN.DIP)/sqrt(n())) #calculates mean 
data.summary

boxplot(DIN.DIP~treatment, data=watercol.dat.insitu, ylab= "DIN:DIP")

bartlett.test(DIN.DIP~treatment, data=watercol.dat.insitu)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (DIN.DIP~treatment, data = watercol.dat.insitu, method = "t.test")

e <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=26), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold("DIN:DIP"))) +
  theme(legend.position = "none")

ggsave(filename = "Nutrients/Output/DIN:DIP.boxplot.png", device = "png", width = 10, height = 13)


##################################################################################
##################################################################################
#filter to just compare tank measurements 
watercol.dat.tanks <- watercol.dat %>%
  filter(treatment %in% c("Enriched Tank", "Control Tank"))

#check variance and make plot for N.N comparisons
SummaryByGroup <- watercol.dat.tanks %>%
  group_by(treatment) %>%
  summarize(variance=var(N.N, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.tanks %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(N.N), se=sd(N.N)/sqrt(n())) #calculates mean 
data.summary

boxplot(N.N~treatment, data=watercol.dat.tanks, ylab= "Nitrate")

bartlett.test(N.N~treatment, data=watercol.dat.tanks)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (N.N~treatment, data = watercol.dat.tanks, method = "t.test")

x <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Nitrate~(NO["3"]^~~{"-"})~+~Nitrite~(NO["2"]^~~{"-"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/N.N.tanks.png", device = "png", width = 10, height = 12)

#check variance and make plot for P comparisons
SummaryByGroup <- watercol.dat.tanks %>%
  group_by(treatment) %>%
  summarize(variance=var(P, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.tanks%>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(P), se=sd(P)/sqrt(n())) #calculates mean 
data.summary

boxplot(P~treatment, data=watercol.dat.tanks, ylab= "P")

bartlett.test(P~treatment, data=watercol.dat.tanks)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (P~treatment, data = watercol.dat.tanks, method = "t.test")

y <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Phosphate~(PO["4"]^~~{"3-"})~(mu*mol~L^{-1}))))  +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/P.tanks.png", device = "png", width = 10, height = 12)

#check variance and make plot for NH4 comparisons
SummaryByGroup <-watercol.dat.tanks %>%
  group_by(treatment) %>%
  summarize(variance=var(NH4, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.tanks %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(NH4), se=sd(NH4)/sqrt(n())) #calculates mean 
data.summary

boxplot(NH4~treatment, data=watercol.dat.tanks, ylab= "NH4")

bartlett.test(NH4~treatment, data=watercol.dat.tanks)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (NH4~treatment, data = watercol.dat.tanks, method = "t.test")

z <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Ammonium~(NH["4"]^~~{"+"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/NH4.tankst.png", device = "png", width = 10, height = 12)

#check variance and make plot for DIN:DIP comparisons
SummaryByGroup <-watercol.dat.tanks %>%
  group_by(treatment) %>%
  summarize(variance=var(DIN.DIP, na.rm=TRUE))
SummaryByGroup

data.summary<-watercol.dat.tanks %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(DIN.DIP), se=sd(DIN.DIP)/sqrt(n())) #calculates mean 
data.summary

boxplot(DIN.DIP~treatment, data=watercol.dat.tanks, ylab= "DIN:DIP")

bartlett.test(DIN.DIP~treatment, data=watercol.dat.tanks)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (DIN.DIP~treatment, data = watercol.dat.tanks, method = "t.test")

q <- ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold("DIN:DIP"))) +
  theme(legend.position = "none") 


ggsave(filename = "Nutrients/Output/DIN.DIP.tanks.png", device = "png", width = 10, height = 12)



#make plots for all nutrients in tank conditions

figure <- x + y + z + q +       #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 20, face = "bold"))   #edit the lettered text


figure


ggsave(filename = "Nutrients/Output/tank_nutrient_graphs.png", device = "png", width = 13, height = 13)


###########################################################################################################
#make plots for all nutrients on fore reef sites
#make plots for all nutrients in tank conditions

figure <- a + b + c + d + e  +         #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 20, face = "bold"))   #edit the lettered text


figure


ggsave(filename = "Nutrients/Output/nutrient_graphs.png", device = "png", width = 15, height = 13)

###########################################################################################################
#compare enriched tank to enriched water column nutrients

reef.tank.conditions <- read.csv("Nutrients/Data/watercolumn.nutrients.Oct.csv")
View(reef.tank.conditions)


data.summary<-reef.tank.conditions %>%
  group_by(treatment, condition) %>% #tells to group by treatment
  summarise(mean=mean(N.N), se=sd(N.N)/sqrt(n())) #calculates mean 
data.summary

ggplot(data.summary, aes(x=treatment, y=mean, color = condition)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Nitrate~(NO["3"]^~~{"-"})~+~Nitrite~(NO["2"]^~~{"-"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 

data.summary<-reef.tank.conditions %>%
  group_by(treatment, condition) %>% #tells to group by treatment
  summarise(mean=mean(NH4), se=sd(NH4)/sqrt(n())) #calculates mean 
data.summary

boxplot(NH4~treatment, data=reef.tank.conditions, ylab= "Nitrate")

bartlett.test(NH4~treatment, data=reef.tank.conditions)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (NH4~treatment, data = reef.tank.conditions, method = "t.test")

ggplot(data.summary, aes(x=treatment, y=mean, color = condition)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Ammonium~(NH["4"]^~~{"+"})~(mu*mol~L^{-1})))) +
  theme(legend.position = "none") 

data.summary<-reef.tank.conditions %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(P), se=sd(P)/sqrt(n())) #calculates mean 
data.summary

boxplot(P~treatment, data=reef.tank.conditions, ylab= "Nitrate")

bartlett.test(P~treatment, data=reef.tank.conditions)

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test 

compare_means (P~treatment, data = reef.tank.conditions, method = "t.test")

ggplot(data.summary, aes(x=treatment, y=mean)) + 
  geom_point(size=6) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=20), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=22), axis.title.y = element_text(face="bold", color="black", size=22),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(Phosphate~(PO["4"]^~~{"3-"})~(mu*mol~L^{-1}))))  +
  theme(legend.position = "none") 







  
