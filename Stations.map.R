###LOAD PACKAGES
library(ggcyto)
library(stringr)
library(car)
library(rgl)
library(ggplot2)
library(geo)
library(devtools)
library(maps)
library(mapdata)
library(lubridate)
library(tidyverse)
library(cowplot)
theme_set(theme_grey())
library(ggmap)


###VARIABLES
	Map.name <- "Spring_2017_Phytoplankton_concentration_around_Iceland.png" #Name of the map
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/figures/'#Path of the folder in which the map will be saved
	
	#INFORMATIONS ABOUT THE STATIONS
	Infos.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/20180704_Summary_after_automatic_gating_QC.csv'#Path of the csv file containing the informations about the stations
	
	Infos1 <- read.csv(Infos.path, sep=';', header=TRUE)
	iceland <- read.csv('/Users/bland/Desktop/Flow-cytometry_data/Summary_stations/iceland.csv', sep=';', header=TRUE) ###Read the csv file containing the coordinate of Iceland
	




head(iceland)
glimpse(iceland)

Infos <- Infos1[c(which(Infos1$Smp.depth=="0"), which(Infos1$Smp.depth=="10"), which(Infos1$Smp.depth=="20"),which(Infos1$Smp.depth=="30"),which(Infos1$Smp.depth=="50"), which(Infos1$Smp.depth=="100"), which(Infos1$Smp.depth=="200"), which(Infos1$Smp.depth=="400"), which(Infos1$Smp.depth=="600")),]
Infos$stn.lon <- gsub(",",".",Infos$stn.lon)
Infos$stn.lat <- gsub(",",".",Infos$stn.lat)
Infos$Smp.depth <- as.numeric(as.vector(Infos$Smp.depth))
Infos$Abundance <- as.numeric(as.vector(Infos$Abundance))
Infos$stn.lat <- as.numeric(as.vector(Infos$stn.lat))
Infos$stn.lon <- as.numeric(as.vector(Infos$stn.lon))


p <- ggplot(iceland, aes(long, lat)) + labs(x = "Longitude [decimal]", y = "Latitude [decimal]") + xlim(-29, -7) + ylim(62.5, 68.5) + theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_line(colour = "grey80"), panel.grid.minor = element_blank())
p1 <- p + geom_point()
p2 <- p + geom_line()
p3 <- p + geom_path() + geom_polygon(fill = "grey60", colour="grey60")



png.file <- png(paste(img.path,Map.name ,sep =""), width=1200, height=750)	
png.file
p3 + geom_point(data=Infos, aes(x=as.numeric(as.character(stn.lon)), y=as.numeric(as.character(stn.lat)), size=Abundance), colour = "yellowgreen", alpha = 0.6) + scale_size_area(max_size = 20) + facet_wrap(~Smp.depth, ncol=2) + labs(title = "Concentration of phytoplankton (events/mL) around Iceland in spring 2017 depending on the depth") + theme(legend.title = element_blank(), legend.key = element_rect(fill = "white", colour = "white"), plot.title = element_text(size = 23, face="bold"), axis.text = element_text(size=15), legend.text = element_text(size=18), axis.title = element_text(size=18), strip.text= element_text(size=18))



dev.off() #close window
