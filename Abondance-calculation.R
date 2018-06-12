###LOAD PACKAGES
	library(flowCore)
	library(ggcyto)
	library(stringr)

###VARIABLES

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results
	csv.name <- "_Abundance_with_all_info_results.csv" #Name of the CSV file containing the results
	liste.stations <- c('FX5','FX8','HB1','IH1','LB2','LB8','ST2','KG1','KG2','KG3','KG4','KG5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	today <- '20180606'
	
	#TYPE OF TRANSFORMATION
	Transfo.type <- logTransform(transformationId="LogTransform", logbase=10, r=1, d=1) #Type of transformation
	to.transform <- c('FSC.A', 'SSC.A', 'Chlorophyll.A', 'SybrGreen.A', 'PE.A') #List of the names of the measurement parameters on which the transformations have to be performed
	
	#BACKGROUND NOISE
	NoiseSyb.min <- 0 #Minimal value of the noise in SybrGreen.A
	NoiseSyb.max <- 2 #Maximal value of the noise in SybrGreen.A
	NoisePE.min <- 0 #Minimal value of the noise in PE.A
	NoisePE.max <- 2 #Maximal value of the noise in PE.A
	NoiseChl.min <- 0 #Minimal value of the noise in Chlorophyll.A
	NoiseChl.max <- 3.15 #Maximal value of the noise in Chlorophyll.A
	
	#BEADS GATE
	BeadsSyb.min <- 5 #Minimal value of the beads gate in SybrGreen.A
	BeadsSyb.max <- 7 #Maximal value of the beads gate in SybrGreen.A
	BeadsChl.min <- 4.5 #Minimal value of the beads gate in Chlorophyll.A
	BeadsChl.max <- 5.5 #Maximal value of the beads gate in Chlorophyll.A
	BeadsSSC.min <- 4.5 #Minimal value of the beads gate in SSC.A
	BeadsSSC.max <- 5.5 #Maximal value of the beads gate in SSC.A
	
	#INFORMATIONS ABOUT THE STATIONS
	stn.lane.list <- c("Faxafloi","Faxafloi","Faxafloi","Faxafloi","Faxafloi","Faxafloi","Faxafloi","Faxafloi","Faxafloi","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Latrabjarg","Kogur", "Kogur", "Kogur", "Kogur", "Kogur","Kogur", "Hornbanki","Hornbanki","Hornbanki","Hornbanki","Hornbanki","Hornbanki", "Siglunes","Siglunes","Siglunes","Siglunes","Siglunes","Siglunes","Siglunes","Siglunes", "Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes","Langanes", "Krossanes","Krossanes","Krossanes","Krossanes","Krossanes","Krossanes", "Stokksnes","Stokksnes","Stokksnes","Stokksnes","Stokksnes", "Ingolfshofdi","Ingolfshofdi","Ingolfshofdi", "Selvogsbanki","Selvogsbanki","Selvogsbanki","Selvogsbanki","Selvogsbanki")
	stn.name.list <- c("Faxafloi 1","Faxafloi 2","Faxafloi 3","Faxafloi 4","Faxafloi 5","Faxafloi 6","Faxafloi 7","Faxafloi 8","Faxafloi 9","Latrabjarg 1","Latrabjarg 2","Latrabjarg 3","Latrabjarg 4","Latrabjarg 5","Latrabjarg 6","Latrabjarg 7","Latrabjarg 8","Latrabjarg 9", "Kogur 1", "Kogur 2", "Kogur 3", "Kogur 4", "Kogur 5","Kogur 6", "Hornbanki 1","Hornbanki 2","Hornbanki 3","Hornbanki 4","Hornbanki 5","Hornbanki 6", "Siglunes 1","Siglunes 2","Siglunes 3","Siglunes 4","Siglunes 5","Siglunes 6","Siglunes 7","Siglunes 8", "Langanes 1","Langanes 2","Langanes 3","Langanes 4","Langanes 5","Langanes 6","Langanes 1","Langanes 2","Langanes 3","Langanes 4","Langanes 5","Langanes 6","Langanes 7","Langanes 8", "Krossanes 1","Krossanes 2","Krossanes 3","Krossanes 4","Krossanes 5","Krossanes 6", "Stokksnes 1","Stokksnes 2","Stokksnes 3","Stokksnes 4","Stokksnes 5", "Ingolfshofdi 1","Ingolfshofdi 2","Ingolfshofdi 3", "Selvogsbanki 1","Selvogsbanki 2","Selvogsbanki 3","Selvogsbanki 4","Selvogsbanki 5")
	stn.id.list <- c("FX1", "FX2", "FX3","FX4","FX5","FX6","FX7","FX8","FX9","LB1","LB2","LB3","LB4","LB5","LB6","LB7","LB8","LB9","KG1","KG2","KG3","KG4","KG5","KG6","HB1","HB2","HB3","HB4","HB5","HB6","SI1","SI2","SI3","SI4","SI5","SI6","SI7","SI8","LN1", "LN2","LN3","LN4","LN5","LN6","LA1","LA2","LA3","LA4","LA5","LA6","LA7","LA8","KR1","KR2","KR3","KR4","KR5","KR6","ST1","ST2","ST3","ST4","ST5","IH1","IH2","IH3","SB1","SB2","SB3","SB4","SB5")
	stn.lat.list <- c(64.2,64.2,64.2,64.2,64.2,64.2,64.2,64.2,64.2,65.3,65.35,65.4,65.45,65.5,65.54,66.1,66.5,66.9,66.3,66.41,66.53,67.5,67.2,67.35,66.4,66.5,67,67.1,67.2,67.3,66.16,66.24,66.32,66.44,67,67.2,67.4,68,66.37,67,67.15,67.3,67.45,68,66.22,66.22,66.22,66.22,66.22,66.22,66.22,66.22,65,65,65,65,65,65,64.12,64.2,63.52,63.48,63.38,63.46,63.42,63.34,63.41,63.29,63.19,63.9,63)
	stn.lon.list <- c(-22.25,-22.45,-23.15,-23.45,-24.2,-25,-26,-27,-27.57,-24.34,-24.55,-25.16,-25.39,-26,-26.29,-26.48,-27.3,-27.15,-23,-23.9,-23.19,-23.28,-23.42,-23.56,-21.35,-21.35,-21.35,-21.35,-21.35,-21.35,-18.5,-18.5,-18.5,-18.5,-18.5,-18.5,-18.5,-18.5,-14.16,-13.55,-13.34,-13.16,-12.58,-12.4,-14.22,-14.1,-13.35,-13,-12.5,-11,-10,-9,-13.3,-12.49,-11.4,-11.17,-10.7,-9,-14.5,-14.28,-14.8,-13.58,-13.4,-16.36,-16.3,-16.18,-20.41,-20.54,-21.7,-21.18,-21.28)
	stn.max.depth.list <- c(70,40,94,190,220,235,327,430,1010,34,64,108,256,228,290,410,658,510,45,78,236,240,500,980,100,120,200,230,330,608,80,435,470,700,230,498,408,1045,188,420,1576,1760,1830,1890,60,156,275,160,1094,1400,1300,1260,54,116,230,518,600,1420,84,141,216,546,1192,72,90,108,46,90,155,510,1004)


	
	
###FUNCTIONS =====================================================================


Sort.Files <- function(listouille, max.depth) {  ##This function returns the list of the sorted files

	List.sorted <- c()
	
	Match2 <- grep("BLANK", listouille, value = TRUE) ###Find the file corresponding to the blank (to be sorted at the first position)
	
	if(length(Match2)>0) {
			
		List.sorted <- append(List.sorted, Match2[1])
		
	}
	
	for(dep in 0:max.depth) { ###Sort the file with the depth
		
		dep.match <- paste("_", dep, "m", sep="")
		Match <- grep(dep.match, listouille, value = TRUE)
		
		if(length(Match)>0) {
			
			List.sorted <- append(List.sorted, Match[1])
			
		}
	
	}
	

	return(List.sorted)

}

Find.Depth <- function(listouille, max.depth) {  ###This function returns the list of the depths

	Profondeur <- c()
	
	Match2 <- grep("BLANK", listouille, value = TRUE) ###Find the file corresponding to the blank (to be sorted at the first position)
	
	if(length(Match2)>0) {
			
		Profondeur <- append(Profondeur, "BLANK")
		
	}
	
	for(dep in 0:max.depth) { ###Sort the file with the depth
		
		dep.match <- paste("_", dep, "m", sep="")
		Match <- grep(dep.match, listouille, value = TRUE)
		
		if(length(Match)>0) {
			
			Profondeur <- append(Profondeur, dep)
			
		}
	
	}
	

	return(Profondeur)

}

###===============================================================================
	
	
###INITIALISATION DE LA LISTE DE TRANSFORMATION
	myTrans <- transformList(to.transform, Transfo.type)

###DETERMINATION OF BEADS GATE
	BeadsSyb.Gate <- rectangleGate(filterId="Beads Region","SybrGreen.A"=c(BeadsSyb.min, BeadsSyb.max))
	BeadsChl.Gate <- rectangleGate(filterId="Beads Region","Chlorophyll.A"=c(BeadsChl.min, BeadsChl.max))
	BeadsSSC.Gate <- rectangleGate(filterId="Beads Region","SSC.A"=c(BeadsSSC.min, BeadsSSC.max))
	Beads.Gate <- BeadsSyb.Gate & BeadsChl.Gate & BeadsSSC.Gate


###DETERMINATION OF NOISE GATE
	NoiseSyb.Gate <- rectangleGate(filterId="Noise","SybrGreen.A"=c(NoiseSyb.min, NoiseSyb.max))
	NoisePE.Gate <- rectangleGate(filterId="Noise","PE.A"=c(NoisePE.min, NoisePE.max))
	NoiseChl.Gate <- rectangleGate(filterId="Noise","Chlorophyll.A"=c(NoiseChl.min, NoiseChl.max))
	Noise.Gate <- NoiseSyb.Gate & NoisePE.Gate & NoiseChl.Gate

###CREATION OF FLOWSET (ONE PER STATION)
	Station.frames <- c()
	Beads.frames <- c()
	Noise.frames <- c()
	
	Samp.Name <- c("Name of the FCS file") #List containing the names of all the analysed files
	Smp.depth <- c("Depth (in meters)")
	
	for (station.index in 1:length(liste.stations)) {
		
		setwd(Folder.path)
		list.FCSname <- Sort.Files(list.files(pattern=liste.stations[station.index]), 2000)
		Smp.depth <- append(Smp.depth, Find.Depth(list.files(pattern=liste.stations[station.index]), 2000))
		
		for(truc in 1:length(list.FCSname)) {
		
			Samp.Name <- append(Samp.Name, list.FCSname[truc])
		
		}
		
		###READ FLOWSET
		fs <- read.flowSet(files=list.FCSname, alter.names=TRUE, transformation =FALSE)
		
		###TRANSFORMATION
		fs.trans <- transform(fs, myTrans) #Transformation of the data
		Station.frames <- append(Station.frames, fs.trans)
		
		###SUPPRESSION OF BACKGROUND NOISE
		Noise <- Subset(fs.trans, Noise.Gate)
		Noise.frames <- append(Noise.frames, Noise)
		
		###BEADS GATING
		Beads <- Subset(fs.trans, Beads.Gate)
		Beads.frames <- append(Beads.frames, Beads)
		
		
	}

###INITIALISATION OF THE VECTORS CONTAINING THE DATA (the first value of each vector is the label)
	Nb.Totevent <- c("Total number of events")
	Nb.beads <- c("Number of beads")
	Nb.Noise <- c("Number of events to remove (background noise)")
	Pourc.Noise <- c("Pourcentage of background noise removed (%)")
	Abundance <- c("Concentration of phytoplankton (number of events / mL)")
	stn.lane <- c("Station lane")
	stn.name <- c("Station name")
	stn.id <- c("Station ID")
	stn.lat <- c("Station latitude")
	stn.lon <- c("Station longitude")
	stn.max.depth <- c("Maximal depth (in meters)")
	
###CREATION OF ONE PDF FILE CONTENING PLOTS / STATION
	for (station in 1:length(liste.stations)) {
	
	ind <- match(liste.stations[station],stn.id.list)
	
	
		for (prof in 1:length(Station.frames[[station]])) {
				
				NbEvent <- nrow(Station.frames[[station]][[prof]])
				NbBeads <- nrow(Beads.frames[[station]][[prof]])
				NbNoise <- nrow(Noise.frames[[station]][[prof]])
				PNoise <- (NbNoise*100)/NbEvent
				
				Abond <- ((NbEvent-NbBeads-NbNoise)/(NbBeads/1080000))
				
				Nb.Totevent <- append(Nb.Totevent, NbEvent)
				Nb.beads <- append(Nb.beads, NbBeads)
				Nb.Noise <- append(Nb.Noise, NbNoise)
				Abundance <- append(Abundance, Abond)
				stn.lane <- append(stn.lane, stn.lane.list[ind])
				stn.name <- append(stn.name, stn.name.list[ind])
				stn.id <- append(stn.id, stn.id.list[ind])
				stn.lat <- append(stn.lat, stn.lat.list[ind])
				stn.lon <- append(stn.lon, stn.lon.list[ind])
				stn.max.depth <- append(stn.max.depth, stn.max.depth.list[ind])
				Pourc.Noise <- append(Pourc.Noise, PNoise)
				
				
			}
		
		
	}

###CREATION OF A CSV FILE TO STORE THE DATA	
	csv.name <- paste(csv.path, today, csv.name, sep="")
	results <- cbind(Samp.Name, stn.lane, stn.name, stn.id, stn.lat, stn.lon, stn.max.depth, Smp.depth, Nb.Totevent, Nb.beads, Nb.Noise, Pourc.Noise, Abundance)
	write.csv(results,csv.name)
