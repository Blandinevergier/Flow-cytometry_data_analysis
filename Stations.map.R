###LOAD PACKAGES
library(flowCore)
library(ggcyto)
library(stringr)
library(car)
library(rgl)

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

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Figures/Plots/3D_plots/' #Path of the folder containing the PDF Files for the results
	csv.name <- "_Abundance_with_all_info_results.csv" #Name of the CSV file containing the results
	pdf.name <- "_Plots_with_gating.pdf" #Name of the pdf containing the plots with the gates
	liste.stations <- c('FX5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	today <- '20180615'
	
	#INFORMATIONS ABOUT THE STATIONS
	Infos.path <- '/Users/bland/Desktop/Flow-cytometry_data/Summary_stations/hydrostn.csv' #Path of the csv file containing the informations about the stations
	
	Infos1 <- read.csv(Infos.path, sep=';', header=TRUE)
	


Infos1 %>% 
 ggplot(aes(stn.lat, stn.lon, size = stn.max.depth)) +
 geom_point(colour = "red")
