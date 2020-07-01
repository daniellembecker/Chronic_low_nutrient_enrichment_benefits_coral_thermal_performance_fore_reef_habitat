##Processing temperature data for three tanks used in Mo'orea October trip, compare to fore reef temp data

#clear dataframe
rm(list=ls())

##Install packages
# load packages
library("segmented")
library("devtools")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library('tidyverse')
library('lubridate')
library('scales')
library('pastecs')
library(lme4)
library(ggpubr)
library(lmerTest)
library(here)


#set wd
here()


#load temp data

gumptank.1 <- read.csv("October_2019/HOBO/Data/tank.data/GumpTank1.csv")
gumptank.2 <- read.csv("October_2019/HOBO/Data/tank.data/GumpTank2.csv")


#bind all data sets together into one data frame
temp.data <- rbind(gumptank.1, gumptank.2)


#remove any NA's from the data frame
temp.data <- na.omit(temp.data)

view(temp.data)


#modify the date and time structure

temp.data$Date.Time <- mdy_hm(temp.data$Date.Time, quiet=FALSE, tz="UTC", truncated=0) #format date and time 

#save file as csv so transforms to data frame
write.csv(temp.data, 'October_2019/HOBO/Data/temp.data.csv') 

#make tank number a factor
temp.data$tank.number = as.factor(temp.data$tank.number)

#summarize temps for mean, min, max, 90th percentile, etc.
  
temp.data.sum <- temp.data %>%
  group_by(tank.number) %>%
  summarise(mean.T = mean(Temp, na.rm=T), #do this after editing the dates and times
  min.T = min(Temp, na.rm=T),
  max.T = max(Temp, na.rm=T),
  var.T = var(Temp, na.rm=T),
  qnorm.T = quantile(Temp,probs = c(0.90), na.rm=T))

view(temp.data.sum)

# #filter out data from october 27th after 06:00
temp.data <- temp.data%>%
   filter(!(tank.number == "1" & Date.Time > "2019-10-27 05:00:00")) %>%
   filter(!(tank.number == "3" & Date.Time > "2019-10-27 05:00:00")) 

# #filter out data from october 12th BEFORE 18:00
temp.data <- temp.data%>%
  filter(!(tank.number == "1" & Date.Time < "2019-10-12 12:45:00")) %>%
  filter(!(tank.number == "3" & Date.Time < "2019-10-12 06:00:00")) 


#graphing temp data from HOBO loggers for tanks 
ggplot(data = temp.data, aes(x = Date.Time, y = Temp, colour = tank.number)) +
  geom_point() +
  labs(x = "Date.Time", y = "Temperature (°C)") + #label x and y axis
  labs(colour="Tank Number") + #change legend title
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) + #change color scheme for lines
  geom_line(aes(group=tank.number))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_datetime(date_breaks = "12 hour", labels = date_format("%b %d - %H:%M")) +
  labs(x="Date - Time", y="Temperature")

ggsave(filename = "October_2019/HOBO/Output/temp.graph.png", device = "png", width = 20, height = 10)

#variances of temp data between tanks
boxplot(ylab="Temp", xlab= "tank.number", Temp~tank.number, data=temp.data)

#if >0.05 same variance, TRUE students t test
bartlett.test(Temp~tank.number, data=temp.data)

#get data sum with just the mean and st_errors


datasum.temp <- temp.data %>% 
  group_by(tank.number) %>% 
  summarise(mean.T = mean(Temp, na.rm=T), #do this after editing the dates and times
            var.T = var(Temp, na.rm=T), qnorm.T = quantile(Temp,probs = c(0.90), na.rm=T))

datasum.temp

#look at forreef lter site 2 data from 2005-06/19

forereef.temp <- read.csv("October_2019/HOBO/Data/FOR02_10m_temp.csv")

#look at mean temp over 
datasum.temp <- forereef.temp %>% 
  summarise(mean.T = mean(temperature_c, na.rm=T), #do this after editing the dates and times
            var.T = var(temperature_c, na.rm=T), qnorm.T = quantile(temperature_c,probs = c(0.90), na.rm=T))

datasum.temp


