
#author: "HM Putnam"
#date: "6/15/2020"
#edited by: Danielle Becker

# load packages
library(tidyverse)
library(vegan)

#ITS2 Profile level Analysis and Plotting
data <- read.csv("October_2019/bioinformatics/Data/counts.profile.data.csv", sep=",", header=TRUE) #load data
rownames(data) <-data$sample_name #identify sample names
otumat <- as.matrix(data[,c(21:22)]) #counts data
meta <- data[,-c(21:22)] #sample info

#counts descriptive stats
sum(data$D1)+sum(data$C42g.C42a.C42.2.C1.C1b.C42b.C42h) # total number of sequences
range(rowSums(otumat)) #range of sequences per sample
mean(rowSums(otumat)) #mean of sequences per sample

#percent D per sample
which(data$D1>0)
percent.D <- (data[which(data$D1>0),22]/data[which(data$D1>0),21])*100
percent.D

#percent D overall
(sum(data$D1)/sum(data$C42g.C42a.C42.2.C1.C1b.C42b.C42h))*100
sum(data$D1)+sum(data$C42g.C42a.C42.2.C1.C1b.C42b.C42h)

#Check Relative Abundance
rel.profile.mat <- otumat/rowSums(otumat) #calculate relative abundance
N <- rowSums(rel.profile.mat) # calculate row sums of the relative abundnaces, should sum to 1
N #Display Row Sums, should all be 1

#Transform Relative Abundance Matrix using sqrt
Trans.Rel.Sym.Data <- sqrt(rel.profile.mat) #sqrt transform the matrix
Trans.Rel.Sym.Data #view data

#Identify profile names
types <- colnames(Trans.Rel.Sym.Data)
types 

#change format to longer
df <- data %>%
  pivot_longer(cols = C42g.C42a.C42.2.C1.C1b.C42b.C42h:D1,
               names_to = c("Profile"),
               values_to="rel.abun")

#plot relative abundance of profile types
ggplot(data = df, mapping = aes(x = Rep, y = rel.abun, fill = Profile)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(facets =  vars(Treatment), ncol = 2)+ 
  ylab("Relative Abundance") + xlab("Replicate") +
  theme(legend.position="top") + theme(legend.direction = "horizontal") +
  scale_fill_manual("ITS2 Profiles", values = c("C42g.C42a.C42.2.C1.C1b.C42b.C42h" = "grey", "D1" = "black")) +
  theme(axis.text.x = element_text(size=rel(0.3))) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 16), legend.title=element_text(size=20), legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=20), axis.title.y = element_text(color="black", size=20),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave("October_2019/bioinformatics/Output/rel.abund.its2.profiles.pdf", width = 12, height = 10)

#Test for differences by site
dist <- vegdist(Trans.Rel.Sym.Data, method="bray") #calculate distance matrix of sym profiles
permanova <- adonis(dist ~ Treatment, data=data, permutations=9999) #apply permutations test
permanova #view results

