#Photosynthesis Irradiance Curve code October 2019
##created by Danielle Becker 03/02/20
##edited by Danielle Becker 04/02/20

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library('devtools')
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 


#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('tidyverse')


##### PHOTOSYNTHESIS AND RESPIRATION #####

# get the file path
setwd("PI_curve/")
path.p<-"../PI_curve/Data/PI_curve_resp" #the location of all your respirometry files 

# bring in the oxygen files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name (note, these are the same for this example, but I often have lots of subfolders for different Runs)
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate an empty 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names), ncol=5))
colnames(Photo.R) <- c("fragment.ID.full","Intercept", "umol.L.sec","Temp.C","Light_Dark") # name the columns


#Load Sample Info
#Sample.Info <- read.csv(file=paste0(path.p,"/../Panama MetaData/Nubbin_Sample_Info_T0_Panama_QC.csv"), header=T) #read sample.info data
Sample.Info <- read.csv(file=paste0(path.p,"/../resp_data.csv"), header=T) #read sample.info data

# load surface area data
SA <- read.csv(file=paste0(path.p,"/../sample_info_PI.csv"), header=T) #read sample.info data
# add 620 ml to the NAs in volume (the blanks)
#Calculat the volume of water
as.numeric(SA$volume)

#Sample.Info$Volume[which(is.na(Sample.Info$Volume))]<-650
# add 0's for the "not blanks"
View(SA)

# joint the sample info and surface area and volume measurements
Sample.Info<-left_join(Sample.Info, SA)

View(Sample.Info)


# make start and stop times real times
Sample.Info$start.time <- as.POSIXct(Sample.Info$start.time,format="%H:%M:%S", tz = "") #convert time from character to time
Sample.Info$stop.time <- as.POSIXct(Sample.Info$stop.time,format="%H:%M:%S", tz = "") #convert time from character to time


PR<-c('Photo','Resp')

 # for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
for(i in 1:length(file.names.full)) { # for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
  
  #find the lines in sample info that have the same file name that is being brought it
  FRow<-which(Sample.Info$fragment.ID==strsplit(file.names[i],'.csv'))
  
  # read in the O2 data one by one
  Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]), skip = 1, header=T) # skips the first line
  Photo.Data1  <- Photo.Data1[,c("Time","Value","Temp")] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  Photo.Data1 <- na.omit(Photo.Data1)

  
  # clean up some of the data
    n<-dim(Photo.Data1)[1] # length of full data
    Photo.Data1 <-Photo.Data1[(n-120):(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data1)[1] #list length of trimmed data
    Photo.Data1$sec <- (1:n) #set seconds by one from start to finish of run in a new column
  
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub(".csv","", file.names[i]) # remove all the extra stuff in the file name
     
    pdf(paste0("Output/",rename,"thinning.pdf")) # open the graphics device
    
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot (empty plot to fill) data as a function of time
    usr  <-  par('usr') # extract the size of the figure margins
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) # put a grey background on the plot
    whiteGrid() # make a grid
    box() # add a box around the plot
    points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1) # add the x axis
    axis(2, las=1) # add the y-axis
    
    # This the data to make the code run faster
    Photo.Data.orig<-Photo.Data1#save original unthinned data
    Photo.Data1 <-  thinData(Photo.Data1 ,by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data1$sec <- as.numeric(rownames(Photo.Data1 )) #maintain numeric values for time
    Photo.Data1$Temp<-NA # add a new column to fill with the thinned data
    Photo.Data1$Temp <-  thinData(Photo.Data.orig,xy = c(1,3),by=20)$newData1[,2] #thin data by every 20 points for the temp values
   
    # plot the thinned data
    plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    #option to add multiple outputs method= c("z", "eq", "pc")
    Regs  <-  rankLocReg(xall=Photo.Data1$sec, yall=Photo.Data1$Value, alpha=0.5, method="pc", verbose=TRUE)  
   
    # add the regression data
    plot(Regs)
    dev.off()
    
    # fill in all the O2 consumption and rate data
    Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[i,1] <- rename #stores the file name in the Date column
    Photo.R[i,4] <- mean(Photo.Data1$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    #Photo.R[i,5] <- PR[j] #stores whether it is photosynthesis or respiration
 
    
    # rewrite the file everytime... I know this is slow, but it will save the data that is already run
}
write.csv(Photo.R, 'Output/Photo.R.csv')  
View(Photo.R)

# Calculate P and R rate

Photo.R$fragment.ID.full<-Photo.R$fragment.ID
Photo.R$fragment.ID<-NULL
View(Photo.R)

Photo.R<-left_join(Photo.R, Sample.Info)
View(Photo.R)
#Convert sample volume to mL
Photo.R$volume <- Photo.R$volume/1000 #calculate volume

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Photo.R$umol.sec <- Photo.R$umol.L.sec*Photo.R$volume

#Account for blank rate by temperature
#convert character columns to factors
Photo.R <- Photo.R %>%
  mutate_if(sapply(., is.character), as.factor)
View(Photo.R)

#make the blank column a factor
#ask nyssa question here, need to pull the BLANKs from sample.type?
Photo.R$BLANK<-ifelse(Photo.R$treatment=='BLANK', 1,0)
Photo.R$BLANK<-as.factor(Photo.R$BLANK)
View(Photo.R)

photo.blnk <- aggregate(umol.sec ~ species*run*BLANK, data=Photo.R, mean)
# pull out only the blanks
#photo.blnk<-photo.blnk[photo.blnk$Species=='BK',]
photo.blnk<-photo.blnk[photo.blnk$BLANK==1,]

# remove the species column and join with the full data set
#photo.blnk$Species<-NULL
# remove the blank column
photo.blnk$BLANK<-NULL

colnames(photo.blnk)[3]<-'blank.rate' # rename the blank rate 
# join the blank data with the rest of the data
Photo.R<-left_join(Photo.R, photo.blnk)

# subtract the blanks######################
Photo.R$umol.sec.corr<-Photo.R$umol.sec-Photo.R$blank.rate

View(Photo.R)

#### Normalize to organic biomass (ash free dry weight)#####

#Calculate net P and R
Photo.R$umol.cm2.hr <- (Photo.R$umol.sec.corr*3600)/Photo.R$surf.area.cm2 #mmol cm-2 hr-1

#Photo.R<-Photo.R[complete.cases(Photo.R),] # remove NAs and blanks
Photo.R<-Photo.R[Photo.R$BLANK==0,]

#make respiration positive
#Photo.R$umol.cm2.hr[Photo.R$PR=='Respiration']<-abs(Photo.R$umol.cm2.hr[Photo.R$PR=='Respiration'])
Photo.R$umol.cm2.hr<- Photo.R$umol.cm2.hr

# log the rates
Photo.R$Rate.ln<-log(Photo.R$umol.cm2.hr+0.1)
#remove empty rows
#Photo.R<-Photo.R[-which(is.na(Photo.R$ID)),]

#ggplot(Photo.R, aes(x=Temp.C, y=umol.cm2.hr,  group=c(Species), col = Organism.ID))+
#geom_point(aes(shape=Species), position = position_dodge(width = 0.2), size=4)+
#ylim(0,1.5)+  
#facet_wrap(~ Species, labeller = labeller(.multi_line = FALSE))
  

write.csv(Photo.R, 'Output/PI_curvesRates.csv') # export all the uptake rates
View(Photo.R)

PhotoMeans<- Photo.R %>%
  group_by(species, run)%>%
  summarise(rates.mean = mean(umol.cm2.hr), se = sd(umol.cm2.hr)/sqrt(n()))


# plot the raw data with the means on top
ggplot()+
  theme_bw()+  
  #geom_point(data=Photo.R, aes(x=run, y=umol.cm2.hr, alpha = 0.05), position = position_dodge(width = 0.2), size=4)+
  geom_point(data=PhotoMeans, aes(x=run, y=rates.mean),  size=1)+
  geom_line(data = PhotoMeans,  aes(x=run, y=rates.mean), size=1)+
  geom_errorbar(data = PhotoMeans, aes(x = run, ymin=rates.mean-se, ymax=rates.mean+se, width=.2))
  #facet_wrap(~ Species, labeller = labeller(.multi_line = FALSE))+
  ggsave('Output/RespirationRates.pdf')

#Mo'orea PI curve fit
#pulling out numeric for everything, pull put for high and pull out for low and do a curve 
  
PAR <- as.numeric(Photo.R$run)
View(Photo.R)

Pc <- as.numeric(Photo.R$umol.cm2.hr) 
  
plot(PAR,Pc,xlab="",
       ylab="",
       xlim=c(0,max(PAR)),
       ylim=c(-1,1.2),
       cex.lab=0.8,cex.axis=0.8,cex=1, main="",
       adj=0.05) 

  #set plot info
  
  mtext(expression("Photon Flux Density ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1)
  
  #add labels
  
  mtext(expression(Rate* " ("*mu*"mol"*~O[2]*" "*cm^-2 *~h^-1*")"),side=2,line=2,cex=1)
  
  #add labels
  
  #fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

  curve.nlslrc= nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,start=list(Am=(max(Pc)-min(Pc)),AQY=0.001,Rd=-min(Pc),theta=0.8)) 
  
  my.fit <- summary(curve.nlslrc)
  #summary of model fit
  
  
  #draw the curve using the model fit
  #hyperbolic tangent, more common way to fit 
  
  mor.curve.fitting<- curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
  
  
  
  
  #Amax (max gross photosytnthetic rate)
  
  Pmax.gross <- my.fit$parameters[1]
  
  
  #AQY (apparent quantum yield) alpha
  
  AQY <- my.fit$parameters[2]
  
  
  #Rd (dark respiration)
  
  Rd <- my.fit$parameters[3]
  
  
  # Ik light saturation point
  
  Ik <- Pmax.gross/AQY

 # add a line to figure for the Ik or saturation point 
  
  abline(v=Ik, col="red", lty=3, lwd = 3)
  text(x = 410, y = 1.0, label = "Ik", srt = 0)
  
  dev.copy(pdf,'Output/Moorea_Oct19_PI.pdf', width = 5, height = 4)
  dev.off()
  
  # Ic light compensation point
  
  Ic <- Rd/AQY
  
  
  # Net photosynthetic rates
  
  Pmax.net <- Pmax.gross-Rd
  
  
  #output parameters into a table
  
  Mor.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
  
  row.names(Mor.PI.Output)<- c("Pg.max","Pn.max","Rdark","alpha","Ik","Ic")
  
  
  
  pdf("Output/Figures_Values_PICurves.pdf")
  
  getwd()
    
    dev.off()

#randomize my tanks and samples
x3 <- sample (1:10, 10)
x3

