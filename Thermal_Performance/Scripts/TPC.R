##Practicing TPC

rm(list=ls())

##Install packages
# load packages
library(nls.multstart)
library(sjPlot)
library(lsmeans)
library(broom)
library(wesanderson)
library(patchwork)
library(tidyverse)
library(nlstools)
library(dplyr)
library(car)
library(nls2)
library(car)
library(here)
library(lme4)
library(lmerTest)

#set wd
here()


#load data

mydata <- read.csv("../Thermal_Performance/Data/Photo.T.csv")
mydata$X <- NULL
View(mydata)
glimpse(mydata)

mydata <- mydata %>% #filtering out NP, removes it from the list 
  filter(rate.type !="NP") 

mydata$log.rate <- log(mydata$umol.cm2.hr + 1)  #logging and adding 1 because a log of zero does not exist

# convert temp to K
mydata$K<-mydata$Temp.C + 273.15

#change GP and R names
mydata <- mydata %>% 
  mutate(rate.type = recode_factor(rate.type, `GP` = "Gross Photosynthesis", `R` = "Dark Respiration"))

#chnage control and enriched names
mydata <- mydata %>% 
  mutate(treatment = recode_factor(treatment, `control` = "Control", `enriched` = "Enriched"))

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp))) #units are eV/K, electrovolts/Kelvin
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}


# # #run nls_multstart
# fit <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26),
#                       data = PA1,
#                       iter = 500,
#                       start_lower = c(lnc = -2, E = 0.1, Eh = 0.2, Th = 285),
#                       start_upper = c(lnc = 3, E = 2, Eh = 5, Th = 330),
#                       supp_errors = 'Y',
#                       na.action = na.omit,
#                       lower = c(lnc = -2, E = 0, Eh = 0, Th = 0))
# 
#  fit


# fit over each set of groupings
fits <- mydata %>%
  group_by(., rate.type, fragment.ID) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26.8),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                                start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))



#broom, models over and over again and purr for lists
#make dplyr code and export the slope, intercept and R2 values
#get r2, extract predit and pull out and join lists, shows o vs predcited and r2
#PredictedPA1_D <- predict(fits$fit[[1]])
#ObservedPA1_D <- mydata$umol.cm2.hr[mydata$fragment.ID == "PA1_D"]

#po <- lm(PredictedPA1_D ~ ObservedPA1_D)
#summary(po)
#plot(ObservedPA1_D,PredictedPA1_D)
#abline(po)
#legend("topleft", bty="n", legend=paste("r2 =", format(summary(po)$adj.r.squared, digits=4)))



# look at a single fit
summary(fits$fit[[1]])

# look at output object
select(fits, fragment.ID, data, fit)  

# get summary info
info <- fits %>%
  unnest_legacy(fit %>% map(glance))

# get params
params <- fits %>%
  unnest_legacy(fit %>% map(tidy))


#left join params with meta data file to have treatment in data frame
metadata <- read.csv(file="../Thermal_Performance/Data/metadata.csv", header=T) #read in metadata file to add site block and recovery block
params <- left_join(params, metadata)

#left join params with meta data file to have block in file
metadata <- read.csv(file="../Thermal_Performance/Data/metadata.csv", header=T) #read in metadata file to add site block and recovery block
params <- left_join(params, metadata)

# get confidence intervals
CI <- fits %>% 
  unnest_legacy(fit %>% map(~ confint2(.x) %>%
                       data.frame() %>%
                       rename(., conf.A = X2.5.., conf.B = X97.5..))) %>%
  group_by(., fragment.ID) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup()

#colnames(CI)[4:5]<-c("conf.low", "conf.high") # rename columns

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  unnest_legacy(fit %>% map(augment))

#select(info, fragment.ID, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions, do this to set a sequence to make a smooth curve with your prediction points
new_preds <- mydata %>%  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE)) #setting a specific sequence so you can have a smooth curve

# max and min for each curve
max_min <- mydata %>% group_by(fragment.ID) %>%
  summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  unnest_legacy(fit %>% map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = "fragment.ID") %>%
  group_by(., fragment.ID) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  rename(., ln.rate = .fitted) %>%
  ungroup()


#left join preds2 with meta data file to have treatment on dataframe
preds2 <- left_join(preds2, metadata)

#want to do ggplot where we look at individual curves for respiration rates of each fragment
  mydataR <- mydata %>%
    filter(rate.type =="Dark Respiration")
  
  preds2R<- preds2 %>%
    filter(rate.type =="Dark Respiration")
  
  ggplot() +
    geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataR) +
    geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2R) +
    facet_wrap(~ fragment.ID, labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual(values = c('green4', 'blue')) +
    theme_bw(base_size = 12, base_family = 'Helvetica') +
    ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
    xlab('Temperature (ºC)') +
    geom_hline(yintercept=0, color = "red") +
    theme(legend.position = c(0.91, 0.85))+
    labs(color = "Rate Type")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(filename = "../Thermal_Performance/Output/respindiv.curves.pdf", device = "pdf", width = 10, height = 10)
  
#want to do ggplot where we look at individual curves for photosynthesis rates of each fragment
  mydataGP <- mydata %>%
    filter(rate.type =="Gross Photosynthesis")
  
  preds2GP<- preds2 %>%
    filter(rate.type =="Gross Photosynthesis")
  
  ggplot() +
    geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataGP) +
    geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2GP) +
    facet_wrap(~ fragment.ID, labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual(values = c('green4', 'blue')) +
    theme_bw(base_size = 12, base_family = 'Helvetica') +
    ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
    xlab('Temperature (ºC)') +
    geom_hline(yintercept=0, color = "red") +
    theme(legend.position = c(0.91, 0.85))+
    labs(color = "Rate Type")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(filename = "../Thermal_Performance/Output/photoindiv.curves.pdf", device = "pdf", width = 10, height = 10)


  
#chnage control and enriched names
  preds2 <- preds2 %>% 
    mutate(treatment = recode_factor(treatment, `control` = "Control", `enriched` = "Enriched"))
  

# plot all P and R values in TPCs
a <- ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = treatment), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = treatment, group = fragment.ID), alpha = 0.5, preds2) +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression(bold("Rate [Log (" *mu*"mol " *cm^-2 *hr^-1*")]"))) +
  xlab('Temperature (ºC)') +
  facet_grid(~ rate.type) + #rename facet wrap headings
  theme(strip.text.x = element_text(size=18, face = "bold")) +
  theme(legend.position = "none",
        axis.text.x=element_text(face="bold", color="black", size=16), 
        axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=20, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  labs(col = "Treatment") 


ggsave(filename = "../Thermal_Performance/Output/TPCcurves.pdf", device = "pdf")


#want to do ggplot where we look at one single set of curves
mydataindiv <- mydata %>%
  filter(individual.ID =="PV2")

#need to make individual data frames to the join because it isnt pulling from correct data sheet if we dont

preds2_D<- preds2 %>%
  filter(fragment.ID =="PV2_D") 

preds2_L<- preds2 %>%
  filter(fragment.ID =="PV2_L")

#join data frames together to use for plot

preds2indiv <- rbind(preds2_D, preds2_L)

#change names of rate type params in data sheet
preds2indiv <- preds2indiv %>% 
  mutate(rate.type = recode_factor(rate.type, `GP` = "Gross Photosynthesis", `R` = "Respiration"))

mydataindiv <- mydataindiv %>% 
  mutate(rate.type = recode_factor(rate.type, `GP` = "Gross Photosynthesis", `R` = "Respiration"))

ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataindiv) +
  geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2indiv) +
  scale_colour_manual(values = c('green4', 'blue')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  labs(x = "Temperature (ºC)", y = expression(bold("Log Rate (" *mu*"mol " *cm^-2 *hr^-1*")"))) +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(size = 14, color = "black", face = "bold"), axis.text.x=element_text(size = 11, color = "black"), axis.text.y = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 14, color = "black", face = "bold")) +
  geom_point(size = 3) + #increase the point size 
  labs(col = "Rate Type") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggsave(filename = "../Thermal_Performance/Output/singlecurve.CSUN.pdf", device = "pdf", width = 6, height = 5)


# function for calculating Topt

get_topt <- function(E, Th, Eh){
  return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
}

# calc topts for all 
Topt_data <- params %>%
  dplyr::select(fragment.ID,  term, estimate, treatment, rate.type, block) %>%
  spread(term, estimate) %>%
  mutate(Topt = get_topt(E, Th, Eh)) 

write.csv(Topt_data, '../Thermal_Performance/Data/Topt_data.csv') 

#chnage control and enriched names
Topt_data <- Topt_data %>% 
  mutate(treatment = recode_factor(treatment, `control` = "Control", `enriched` = "Enriched"))

#get temerature back in celcius not K
Topt_data$Topt <- Topt_data$Topt - 273.15 

#drop NP
Topt_data$rate.type <- droplevels(Topt_data$rate.type)

#anova function mixed effects model with block as a random factor
#topt model for GP
Topt.GP <- lmer(Topt~treatment + (1|block), data = Topt_data, subset = rate.type=="Gross Photosynthesis")
anova(Topt.GP)
summary(Topt.GP)


#topt model for dark respiration
Topt.R <- lmer(Topt~treatment + (1|block), data = Topt_data, subset = rate.type=="Dark Respiration")
anova(Topt.R)
summary(Topt.R)


#check for normality, use normality plots

qqnorm(resid(Topt.R))
qqline(resid(Topt.R))

qqnorm(resid(Topt.GP))
qqline(resid(Topt.GP))

#check heteroscisity with boxplots

boxplot(resid(Topt.R)~Topt_data$treatment)

#extract model paramaters (means/SE) 
Topt.R.ls <- lsmeans(Topt.R, "treatment")
Topt.GP.ls <- lsmeans(Topt.GP, "treatment")

#make data table into a data frame
Topt.R.dat <- as.data.frame(Topt.R.ls)
Topt.GP.dat<- as.data.frame(Topt.GP.ls)

#make another column in lsmeans data frame to distinguish between R and GP
Topt.R.dat$rate.type  <- c("Topt.R")
Topt.GP.dat$rate.type <- c("Topt.GP")

#combine lsmeans tables for Topt values
data.summary.Topt <- rbind(Topt.GP.dat, Topt.R.dat)

#data.summary$rate.type <- factor(data.summary$rate.type, levels = c("R", "GP")) 

d <- ggplot(data.summary.Topt, aes(x=treatment, y=lsmean, col = treatment)) +
  theme_bw() +
  geom_point(size=6) +
  geom_errorbar(aes(ymax=lsmean+SE, ymin=lsmean-SE), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  labs(x="Treatment", y= expression(bold(paste(atop("Acute Thermal Optimum",  "(°C)"))))) +
  facet_wrap(. ~ rate.type, scales = "free") +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.position="none", 
        axis.text.x=element_text(face="bold", color="black", size=18), 
        axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=20, face="bold"), axis.title.y = element_text(color="black", size=20, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x = element_blank())
 

 
ggsave(filename = "../Thermal_Performance/Output/Topt_graph.pdf", device = "pdf", width = 8, height = 5)



#lnc model for dark respiration
lnc.GP <- lmer(lnc~treatment + (1|block), data = Topt_data, subset = rate.type=="Gross Photosynthesis")
anova(lnc.GP)
summary(lnc.GP)


#lnc model for dark respiration
lnc.R <- lmer(lnc~treatment + (1|block), data = Topt_data, subset = rate.type=="Dark Respiration")
anova(lnc.R)
summary(lnc.R)


#check for normality, use normality plots

qqnorm(resid(lnc.R))
qqline(resid(lnc.R))

qqnorm(resid(lnc.GP))
qqline(resid(lnc.GP))

#check heteroscisity with boxplots

boxplot(resid(lnc.mod)~Topt_data$treatment*Topt_data$rate.type)

#extract model paramaters (means/SE) 
lnc.R.ls <- lsmeans(lnc.R, "treatment")
lnc.GP.ls <- lsmeans(lnc.GP, "treatment")

#make data table into a data frame
lnc.R.dat <- as.data.frame(lnc.R.ls)
lnc.GP.dat<- as.data.frame(lnc.GP.ls)

#make another column in lsmeans data frame to distinguish between R and GP
lnc.R.dat$rate.type  <- c("lnc.R")
lnc.GP.dat$rate.type <- c("lnc.GP")

#combine lsmeans tables for Topt values
data.summary.lnc <- rbind(lnc.GP.dat, lnc.R.dat)


c <- ggplot(data.summary.lnc, aes(x=treatment, y=lsmean, col = treatment)) +
  theme_bw() +
  geom_point(size=6) +
  geom_errorbar(aes(ymax=lsmean+SE, ymin=lsmean-SE), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  labs(x="", y=expression(bold(paste(atop("Rate at a reference temperature",  "(log " *mu*"mol " *cm^-2 *hr^-1*")"))))) +
  facet_wrap(. ~ rate.type, scales = "free") +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.position="none", 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(face="bold", color="black", size=13), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size=20, face="bold"), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  theme(strip.text.x = element_blank())


ggsave(filename = "../Thermal_Performance/Output/lnc_graph.pdf", device = "pdf", width = 8, height = 7)
    

#calculate Pmax values between sites

Pmax_data <- Topt_data %>%
  #select(fragment.ID,  term, estimate, site, rate.type) %>%
  #spread(term, estimate) %>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = Topt + 273.15, Tc = 26)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., rate.type, treatment, fragment.ID)

write.csv(Pmax_data, '../Thermal_Performance/Data/Pmax_data.csv') # export all the uptake rates

#drop NP
Pmax_data$rate.type <- droplevels(Pmax_data$rate.type)

#pmax model for dark respiration
Pmax.GP <- lmer(Pmax~treatment + (1|block), data = Pmax_data, subset = rate.type=="Gross Photosynthesis")
anova(Pmax.GP)
summary(Pmax.GP)


#pmax model for dark respiration
Pmax.R <- lmer(Pmax~treatment + (1|block), data = Pmax_data, subset = rate.type=="Dark Respiration")
anova(Pmax.R)
summary(Pmax.R)


#check for normality, use normality plots

qqnorm(resid(Pmax.R))
qqline(resid(Pmax.R))

#check heteroscisity with boxplots

#boxplot(resid(Pmax.R)~Pmax_data$treatment) 
  
#extract model paramaters (means/SE) 
Pmax.R.ls <- lsmeans(Pmax.R, "treatment")
Pmax.GP.ls <- lsmeans(Pmax.GP, "treatment")

#make data table into a data frame
Pmax.R.dat <- as.data.frame(Pmax.R.ls)
Pmax.GP.dat<- as.data.frame(Pmax.GP.ls)

#make another column in lsmeans data frame to distinguish between R and GP
Pmax.R.dat$rate.type  <- c("Pmax.R")
Pmax.GP.dat$rate.type <- c("Pmax.GP")

#combine lsmeans tables for Topt values
data.summary.Pmax <- rbind(Pmax.GP.dat, Pmax.R.dat)


b <- ggplot(data.summary.Pmax, aes(x=treatment, y=lsmean, col = treatment)) +
  theme_bw() +
  geom_point(size=6) +
  geom_errorbar(aes(ymax=lsmean+SE, ymin=lsmean-SE), position=position_dodge(width=0.9), width=0.1) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  labs(x="", y=expression(bold(paste(atop("Maximal perfromance",  "(log " *mu*"mol " *cm^-2 *hr^-1*")"))))) +
  facet_wrap(. ~ rate.type, scales = "free") +
  scale_color_manual(values = wes_palette("Royal1")) +
  theme(legend.position="none", 
        axis.text.x=element_blank(), 
        axis.text.y=element_text(face="bold", color="black", size=13), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size=20, face="bold"), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  theme(strip.text.x = element_blank()) 


ggsave(filename = "../Thermal_Performance/Output/Pmax_graph.pdf", device = "pdf", width = 8, height = 6)


figure <- a / b / c / d           #patchwork to combine plots
  #plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  #theme(plot.tag = element_text(size = 20, face = "bold"))   #edit the lettered text



figure

ggsave(filename = "../Thermal_Performance/Output/holo_graphs.pdf", device = "pdf", width = 10, height = 20)


#arrange models to make table of models results

TPC_models.M2 <- tab_model(Topt.GP, Topt.R, Pmax.GP, Pmax.R,lnc.GP, lnc.R,  pred.labels = c("Intercept", "Treatment"),  p.val = "satterthwaite",
                              dv.labels = c("Topt_GP", "Topt_R", "umax_GP", "umax_R", "b(Tc)_GP", "b(Tc)_R"), string.ci = "Conf. Int (95%)", show.ci = FALSE,
                              string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/Topt_M2.doc")









