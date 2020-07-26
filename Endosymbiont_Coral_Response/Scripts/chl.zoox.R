#clear list 

rm(list=ls())

#load libraries 
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ghibli)
library(ggpubr)
library(patchwork)
library(wesanderson)
library(sjPlot)
library(MASS)
library(moments)
library(car)
library(lmerTest)
library(pastecs)
library(psych)
library(graphics)
library(lme4)
library(vegan)
library(here)

#set wd
here()

#load data

mydata <- read.csv("Endosymbiont_Coral_Response/Data/chl_zoox_sheet.csv")
View(mydata)

metadata <- read.csv("Endosymbiont_Coral_Response/Data/metadata.physio.csv")

mydata <- left_join(mydata, metadata )

#separate out by treatment to check normality, etc.

enriched <- mydata %>%
  filter(treatment=="enriched")

control <- mydata %>%
  filter(treatment=="control")


#rename control and enriched data in treatment column in data sheet to have capital letters for plots
mydata <- mydata %>% 
  mutate(treatment = recode(treatment, "enriched='Enriched'; control='Control'"))


#analysis for treatment x chlA.ugcm2 and graphs 
#check for normality, use normality plots

qqp(enriched$chlA.ugcm2, "norm")


#check for normality, use normality plots


qqp(control$chlA.ugcm2, "norm")


#check variance

SummaryByGroup <- mydata %>%
  group_by(treatment) %>%
  summarize(variance=var(chlA.ugcm2, na.rm=TRUE))
SummaryByGroup

boxplot(chlA.ugcm2~treatment, data=mydata, ylab= "chlA.ugcm2")

bartlett.test(chlA.ugcm2~treatment, data=mydata)

#use two-sample student t-test to test ug.chla, student t-test for pooled variance t-test
compare_means (chlA.ugcm2~treatment, data=mydata, var.equal = TRUE, method = "t.test")

#p>0.05 means variances are equal, use students t test

#mixed effects model with block as random factor
chl.mod <- lmer(chlA.ugcm2~treatment + (1|block), data = mydata)
anova(chl.mod)
tab_model(chl.mod)


#look at mean difference between treatments and total chlorophyll a content
data.summary<-mydata %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(chlA.ugcm2), se=sd(chlA.ugcm2)/sqrt(n())) #calculates mean 
data.summary


b <- ggplot(data.summary, aes(x=treatment, y=mean, col=treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Total Chlorophyll Content", "("*mu*g *~ cm^"-2"*")")))))  + #using quotations over numbers allow them to be bold
  theme(legend.position = "none")  

ggsave(filename = "Endosymbiont_Coral_Response/Output/chloro.png", device = "png", width = 5, height = 5)
 
#analysis for treatment x chlA.pg.cel and graphs 
#check for normality, use normality plots

qqp(enriched$chlA.pg.cell, "norm")


#check for normality, use normality plots


qqp(control$chlA.pg.cell, "norm")


#check variance

SummaryByGroup <- mydata %>%
  group_by(treatment) %>%
  summarize(variance=var(chlA.pg.cell, na.rm=TRUE))
SummaryByGroup

boxplot(chlA.pg.cell~treatment, data=mydata, ylab= "chlA.ugcm2")

bartlett.test(chlA.pg.cell~treatment, data=mydata)

#use two-sample student t-test to test ug.chla, student t-test for pooled variance t-test
compare_means (chlA.pg.cell~treatment, data=mydata, var.equal = TRUE, method = "t.test")

#p>0.05 means variances are equal, use students t test

#mixed effects model with block as random factor
chl.cell.mod <- lmer(chlA.pg.cell~treatment + (1|block), data = mydata)
anova(chl.cell.mod)
tab_model(chl.cell.mod)


#look at mean difference between treatments and total chlorophyll a content
data.summary<-mydata %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(chlA.pg.cell), se=sd(chlA.pg.cell)/sqrt(n())) #calculates mean 
data.summary


c <- ggplot(data.summary, aes(x=treatment, y=mean, col=treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Chlorophyll Content per Cell", "(" *pg*~ cell^"-1"*")")))))  + #using quotations over numbers allow them to be bold
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/chloro.cell.png", device = "png", width = 5, height = 5)

#stat_compare_means( method = "t.test") #to remove the space between x axis and bar graphs

#analysis for treatment x zoox.per.cm2 and graphs 
#check for normality, use normality plots
 
 qqp(enriched$zoox.per.cm2, "norm")
 
 
 #check for normality, use normality plots
 
 
 qqp(control$zoox.per.cm2, "norm")
 
 #check variance
 
 SummaryByGroup <- mydata %>%
   group_by(treatment) %>%
   summarize(variance=var(zoox.per.cm2, na.rm=TRUE))
 SummaryByGroup
 
 
 boxplot(zoox.per.cm2~treatment, data=mydata, ylab= "zoox.per.cm2")

 
 bartlett.test(zoox.per.cm2~treatment, data=mydata)
 
#p>0.05 means variances are equal, use students t test
#not equal variances, so use welchs 
 compare_means (zoox.per.cm2~treatment, data=mydata, na.rm=TRUE)
 
 #mixed effects model with block as random factor 
 zoox.mod <- lmer(zoox.per.cm2~treatment + (1|block), data =mydata)
 anova(zoox.mod)
 tab_model(zoox.mod)

 
#look at mean difference between treatments and edno densities
 data.summary<-mydata %>%
   group_by(treatment) %>% #tells to group by treatment
   summarise(mean=mean(zoox.per.cm2), se=sd(zoox.per.cm2)/sqrt(n())) #calculates mean 
 data.summary
 
a <- ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Endosymbiont Density", "(" *x*"10"^"6" *~cells *~ cm^"-2"*")"))))) + #using quotations over numbers allow them to be bold
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/zoox.png", device = "png", width = 5, height = 6)

#analysis for treatment x AFDW graphs

#load data

AFDWdata <- read.csv("Endosymbiont_Coral_Response/Data/AFDW.csv")
View(AFDWdata)

AFDWdata <- left_join(AFDWdata, metadata)

#need to organize data so it recognizes seperate treatments

enriched <- AFDWdata %>%
  filter(treatment == "enriched")


#need to organize data so it recognizes seperate treatments

control <- AFDWdata %>%
  filter(treatment == "control")

#check for normality, use normality plots

qqp(enriched$AFDW, "norm")


#check for normality, use normality plots


qqp(control$AFDW, "norm")

#rename control and enriched data in treatment column in data sheet to have capital letters for plots
AFDWdata <- AFDWdata %>% 
  mutate(treatment = recode(treatment, "enriched='Enriched'; control='Control'"))


#check variance

SummaryByGroup <- AFDWdata %>%
  group_by(treatment) %>%
  summarize(variance=var(AFDW, na.rm=TRUE))
SummaryByGroup

boxplot(AFDW.mg.cm2.~treatment, data=AFDWdata, ylab= "AFDW") 

bartlett.test(AFDW.mg.cm2.~treatment, data=AFDWdata)

#p>0.05 means variances are equal, use students t test

#use two-sample students t-test to test AFDW, students t-test for equal variance t-test

compare_means (AFDW.mg.cm2.~treatment, var.equal=TRUE, data=AFDWdata, method = "t.test")

#mixed effects model with block as random factor 
AFDW.mod <- lmer(AFDW.mg.cm2.~treatment + (1|block), data = AFDWdata)
anova(AFDW.mod)
tab_model(AFDW.mod)

#look at mean difference between treatments and tissue biomass
data.summary<-AFDWdata %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(AFDW.mg.cm2.), se=sd(AFDW.mg.cm2.)/sqrt(n())) #calculates mean 
data.summary

d <- ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=30), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold(paste(atop("Tissue Biomass", "(mg "*cm^"-2"*")"))))) +  #using quotations over numbers allow them to be bold
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/biomass.png", device = "png", width = 5, height = 5)

################################################################################################################################
#load data for tissue N content info

mydata2 <- read.csv("Endosymbiont_Coral_Response/Data/tissueNcontent.Oct.csv")
View(mydata2)

metadata.N <- read.csv("Endosymbiont_Coral_Response/Data/metadata.tissueN.csv")

#left join metadata
mydata2 <- left_join(mydata2, metadata.N)

#rename control and enriched data in treatment column in data sheet to have capital letters for plots
mydata2 <- mydata2 %>% 
  mutate(treatment = recode(treatment, "enriched='Enriched'; control='Control'"))


#need to organize data so it recognizes seperate sample types

AT <- mydata2 %>%
  filter(sample.type == "AT")


#need to organize data so it recognizes seperate sample types

ST <- mydata2 %>%
  filter(sample.type == "ST")

#analysis for sample type x tissue n content
#check for normality, use normality plots

qqPlot(AT$N, "norm")


#check for normality, use normality plots


qqPlot(ST$N, "norm")


#check variance

SummaryByGroup <- mydata2 %>%
  group_by(sample.type) %>%
  summarize(variance=var(N, na.rm=TRUE))
SummaryByGroup

boxplot(N~sample.type, data=mydata2, ylab= "N")

bartlett.test(N~sample.type, data=mydata2)

#p>0.05 means variances are equal, use students t test

#rename AT and ST in sample type column in data sheet 
mydata2 <- mydata2 %>% 
  mutate(sample.type = recode(sample.type, "AT='Coral Tissue'; ST='Algal Endosymbiont'"))

#p>0.05 means variances are equal, use students t test

#use two-sample student t-test to test total zoox, student t-test for pooled variance t-test

compare_means (N~sample.type, var.equal=TRUE, data = mydata2, method = "t.test")

#mixed effects model with block as random factor
tissue.mod <- lmer(N~sample.type + (1|block), data = mydata2)
anova(tissue.mod)


#plot % N values between AT and ST overall
ggplot(mydata2, aes(x=sample.type, y=N, fill=sample.type)) + 
  geom_boxplot(alpha=0.8) +
  scale_color_manual(values = wes_palette("Royal1")) +
  labs(fill = "Sample Type") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Sample Type") + ylab("% Nitrogen Content") +
  theme(axis.text.x=element_text(face="bold", color="black", size=24), axis.text.y=element_text(face="bold", color="black", size=18), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank())+  #adjust themes for chart x and y axis labels and axis tick mark labels 
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.png", device = "png", width = 10, height = 10)

#plot % N values between AT and ST by treatment
ggplot(mydata2, aes(x=treatment, y=N, col = treatment)) + 
  geom_boxplot(alpha=0.8) +
  scale_color_manual(values = wes_palette("Royal1")) +
  labs(fill = "Sample Type") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="none") +
  facet_wrap(~sample.type) +
  xlab("Treatment") + ylab("% Nitrogen Content") +
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=16), axis.title.x = element_text(face="bold", color="black", size=18), axis.title.y = element_text(face="bold", color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size=12, face = "bold"))

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.per treatment.png", device = "png", width = 10, height = 10)

#filter just ST values to compare pgN/cell and % N algal endosymbionts
mydata2.ST <- mydata2 %>%
  filter(sample.type == "Algal Endosymbiont")

#check variance

SummaryByGroup <- mydata2.ST %>%
  group_by(treatment) %>%
  summarize(variance=var(pgN.cell, na.rm=TRUE))
SummaryByGroup

#look at mean difference between ST pgn/cell content
data.summary<-mydata2.ST %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(pgN.cell), se=sd(pgN.cell)/sqrt(n())) #calculates mean 
data.summary

#mixed effects model with block as random factor
pgN.mod <- lmer(pgN.cell~treatment + (1|block), data = mydata2.ST)
anova(pgN.mod)
tab_model(pgN.mod)

#plot pg/cell between treatments
e <- ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Endosymbiont Nitrogen Content", "(pg N"*~ cell^"-1"*")"))))) +
  theme(legend.position = "none") 

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.ST.pgN.png", device = "png", width = 5, height = 6)


#check variance

SummaryByGroup <- mydata2.ST %>%
  group_by(treatment) %>%
  summarize(variance=var(N, na.rm=TRUE))
SummaryByGroup

#mixed effects model with block as random factor
N.ST.mod <- lmer(N~treatment + (1|block), data = mydata2.ST)
anova(N.ST.mod)
tab_model(N.ST.mod)


#look at mean difference between ST % N content
data.summary<-mydata2.ST %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(N), se=sd(N)/sqrt(n())) #calculates mean 
data.summary

#plot pg/cell between treatments
f <- ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=18), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Total Endosymbiont", "% Nitrogen Content"))))) +
  theme(legend.position = "none") 

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.ST.percentN.png", device = "png", width = 10, height = 10)

#filter just AT values to compare % N coral tissue
mydata2.AT <- mydata2 %>%
  filter(sample.type == "Coral Tissue")

#check variance

SummaryByGroup <- mydata2.AT %>%
  group_by(treatment) %>%
  summarize(variance=var(N, na.rm=TRUE))
SummaryByGroup

#mixed effects model with block as random factor
N.AT.mod <- lmer(N~treatment + (1|block), data = mydata2.AT)
anova(N.AT.mod)
tab_model(N.AT.mod)


#look at mean difference between AT % N content
data.summary<-mydata2.AT %>%
  group_by(treatment) %>% #tells to group by treatment
  summarise(mean=mean(N), se=sd(N)/sqrt(n())) #calculates mean 
data.summary

#plot % N AT between treatments
g <- ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=18), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("") + ylab(expression(bold(paste(atop("Coral Tissue", "% Nitrogen Content"))))) +
  theme(legend.position = "none") 

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.AT.percentN.png", device = "png", width = 10, height = 10)

#use patchwork to organize figure

figure <- a + b + c + e + f + g + d  +           #patchwork to combine plots
    plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
    theme(plot.tag = element_text(size = 20, face = "bold"))   #edit the lettered text


figure

ggsave(filename = "Endosymbiont_Coral_Response/Output/physio_graphs.png", device = "png", width = 17, height = 19)


#arrange models to make table of models results 
endo_coral_M2 <- tab_model(chl.mod, chl.cell.mod, zoox.mod,  pgN.mod, N.ST.mod, N.AT.mod, AFDW.mod, pred.labels = c("Intercept", "Treatment"), show.ci = FALSE,
                              dv.labels = c("Chlorophyll Content", "Chlorophyll Content per Cell", "Endosymbiont Densities", "Endosymbiont N Content per Cell", "Endosymbiont % N Content", "Coral % N Content", "Tissue Biomass"), string.ci = "Conf. Int (95%)",
                              string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/M2.endo.coral.doc")




ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold("Endosymbiont Density (" *x*"10"^"6" *~cells *~ cm^"-2"*")"))) + #using quotations over numbers allow them to be bold
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/zoox.png", device = "png", width = 5, height = 7)

ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=28), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold("Total Chlorophyll Content ("*mu*g *~ cm^"-2"*")")))  + #using quotations over numbers allow them to be bold
  theme(legend.position = "none")  

ggsave(filename = "Endosymbiont_Coral_Response/Output/chloro.png", device = "png", width = 5, height = 7)

ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=26), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=30), axis.title.y = element_text(face="bold", color="black", size=24),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold("Tissue Biomass (mg "*cm^"-2"*")")))  +  #using quotations over numbers allow them to be bold
  theme(legend.position = "none")

ggsave(filename = "Endosymbiont_Coral_Response/Output/biomass.png", device = "png", width = 5, height = 7)

ggplot(data.summary, aes(x=treatment, y=mean, col = treatment)) + 
  geom_point(size = 6) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.title = element_blank()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) +
  labs(fill = "Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(face="bold", color="black", size=24), axis.text.y=element_text(face="bold", color="black", size=20), axis.title.x = element_text(face="bold", color="black", size=24), axis.title.y = element_text(face="bold", color="black", size=21),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  xlab("Treatment") + ylab(expression(bold(paste("Endosymbiont Nitrogen Content (pg N"*~ cell^"-1"*")")))) +
  theme(legend.position = "none") 

ggsave(filename = "Endosymbiont_Coral_Response/Output/tissue.n.content.ST.pgN.png", device = "png", width = 5, height = )
