###LOAD PACKAGES
	library(flowCore)
	library(ggcyto)
	library(stringr)
	library(car)
	library(rgl)

###VARIABLES

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Figures/Plots/3D_plots/' #Path of the folder containing the PDF Files for the results
	csv.name <- "_Abundance_with_all_info_results.csv" #Name of the CSV file containing the results
	pdf.name <- "_Plots_with_gating.pdf" #Name of the pdf containing the plots with the gates
	liste.stations <- c('FX5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	today <- '20180615'
	
	#MINIMAL NUMBER OF BEADS AND EVENT
	minEvents <- 9999 #minimal number of events
	minBeads <- 999 #minimal number of beads
	
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
	NoiseSSC.min <- 0 #Minimal value of the noise in SSC.A
	NoiseSSC.max <- 2 #Maximal value of the noise in SSC.A
	
	#BEADS GATE
	BeadsSyb.min <- 5 #Minimal value of the beads gate in SybrGreen.A
	BeadsSyb.max <- 7 #Maximal value of the beads gate in SybrGreen.A
	BeadsChl.min <- 4.5 #Minimal value of the beads gate in Chlorophyll.A
	BeadsChl.max <- 5.5 #Maximal value of the beads gate in Chlorophyll.A
	BeadsSSC.min <- 4.5 #Minimal value of the beads gate in SSC.A
	BeadsSSC.max <- 5.5 #Maximal value of the beads gate in SSC.A
	
	#INFORMATIONS ABOUT THE STATIONS
	Infos.path <- '/Users/bland/Desktop/Flow-cytometry_data/Summary_stations/hydrostn.csv' #Path of the csv file containing the informations about the stations
	
	Infos1 <- read.csv(Infos.path, sep=';', header=TRUE)
	Infos <- as.matrix(Infos1)
	stn.lane.list <- Infos[,"stn.lane"]
	stn.name.list <- Infos[,"stn.name"]
	stn.id.list <- Infos[,"stn.id"]
	stn.lat.list <- Infos[,"stn.lat"]
	stn.lon.list <- Infos[,"stn.lon"]
	stn.max.depth.list <- Infos[,"stn.max.depth"]

	
	
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

FFtoDF <-function(FF){

  if(class(FF) == "flowFrame"){
    return(as.data.frame(exprs(FF)))
  }

  if(class(FF) == "list"){
    frameList<-list()
    length(frameList)<-length(FF)
    for(i in 1:length(FF)){
      if(class(FF[[i]]) == "flowFrame"){
        frameList[[i]]<-as.data.frame(flowCore::exprs(FF[[i]]))
        names(frameList)[[i]]<-names(FF)[[i]]
      }
      else{
        warning(paste("Object at index",i,"not of type flowFrame"))
      }
    }
    return(frameList)
  }
  else {
    stop("Object is not of type flowFrame")
  }
}

Make.3Dplot <- function(Station.flowFrame, Beads.flowFrame, Noise.flowFrame, X.index, Y.index, Z.index, xlabel, ylabel, zlabel, titre){ ###This function allows to create a 3D scatter plot from a flowframe

	AllStation <- FFtoDF(Station.flowFrame)
	Allevent <- rep("Entire sample",length(AllStation[,1]))
	AllStation$Gate <- Allevent
	
	AllBeads <- FFtoDF(Beads.flowFrame)
	AllBevent <- rep("Beads",length(AllBeads[,1]))
	AllBeads$Gate <- AllBevent
	
	AllNoise <- FFtoDF(Noise.flowFrame)
	AllNevent <- rep("Noise",length(AllNoise[,1]))
	AllNoise$Gate <- AllNevent
	
	matr <- rbind(AllStation, AllBeads, AllNoise)
	

	list.toRemove <- grep("Inf", matr[,X.index],value = FALSE)
	list.toRemove <- append(list.toRemove, grep("Inf", matr[,Y.index],value = FALSE))
	list.toRemove <- append(list.toRemove, grep("Inf", matr[,Z.index],value = FALSE))
	
	if(length(list.toRemove)>0){
	
		matr <- matr[-list.toRemove,] ### Remove rows containing "Inf" value
		
	}
	
	matr$Gate <- as.factor(matr$Gate)
	
				NbEvent <- nrow(Station.flowFrame) + nrow(Beads.flowFrame) + nrow(Noise.flowFrame)
				NbBeads <- nrow(Beads.flowFrame)
				NbNoise <- nrow(Noise.flowFrame)
				
				if (NbEvent > minEvents && NbBeads > minBeads){
					Abond <- ((NbEvent-NbBeads-NbNoise)/(NbBeads/1080000))
					
					}else{
					Abond <- "ERROR"
					}
				
	
	plt3D <- scatter3d(x = matr[,X.index], y = matr[,Y.index], z = matr[,Z.index], xlab=xlabel, ylab=ylabel, zlab=zlabel, sphere.size=0.1, groups = matr$Gate, surface.col=c("darkorange1","steelblue4","snow3"), axis.col=c("black","black","black"), surface=FALSE)
	+ legend3d("topright", legend = c(paste(titre, ' (', NbEvent,' events)', sep=""), ' ', paste('Beads (',NbBeads, ' events)', sep=""), paste('Microbes communities (', nrow(Station.flowFrame), ' events)', sep=""), paste('Background noise (', NbNoise, ' events)', sep=""), ' ', paste('Cell concentration : ', Abond, ' events/mL', sep="")), pch = 16, col = c("white","white","darkorange1","steelblue4","snow3","white","white"), cex=1, inset=c(0.02))
	
	
	return(plt3D)
 
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
	NoiseSSC.Gate <- rectangleGate(filterId="Noise","SSC.A"=c(NoiseSSC.min, NoiseSSC.max))
	Noise.Gate <- NoiseSyb.Gate & NoisePE.Gate & NoiseChl.Gate & NoiseSSC.Gate 
	
###WITHOUT BEADS AND NOISE

	Sans.Noise.Gate <- !Noise.Gate
	Sans.Beads.Gate <- !Beads.Gate
	
	Filtered.Gate <- Sans.Noise.Gate & Sans.Beads.Gate
	

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
		Without.BeadandNoise <- Subset(fs.trans, Filtered.Gate)
		Station.frames <- append(Station.frames, Without.BeadandNoise)
		
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
	
###FILLING OF THE DIFFERENT LISTS CONTAINING THE INFORMATIONS (station names, depths, localisation, abundance...)
	for (station in 1:length(liste.stations)) {
	
	ind <- match(liste.stations[station],stn.id.list)
	
	
		for (prof in 1:length(Station.frames[[station]])) {
				
				NbEvent <- nrow(Station.frames[[station]][[prof]]) + nrow(Beads.frames[[station]][[prof]]) + nrow(Noise.frames[[station]][[prof]])
				NbBeads <- nrow(Beads.frames[[station]][[prof]])
				NbNoise <- nrow(Noise.frames[[station]][[prof]])
				PNoise <- (NbNoise*100)/NbEvent
				
				if (NbEvent > minEvents && NbBeads > minBeads){
					Abond <- ((NbEvent-NbBeads-NbNoise)/(NbBeads/1080000))
					
					}else{
					Abond <- "ERROR"
					}
				
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
	#write.csv(results,csv.name)
	
	
###CREATION OF 3D PLOTS
index <- 1	
	
	for (station in 1:length(liste.stations)) {

	path3D <- paste(img.path, liste.stations[station], sep="")
	path3D1 <- paste(path3D,"/SSC_vs_SYBRGREEN_vs_PE", sep="")
	path3D2 <- paste(path3D,"/Chlorophyll_vs_SYBRGREEN_vs_PE", sep="")
	
 	print(dir.create(path3D))
	print(dir.create(path3D1))
	print(dir.create(path3D2))
	
		for (prof in 1:length(Station.frames[[station]])) {
		
				index <- index + 1
				
				print(par3d(windowRect = 50 + c(0,0,940,740)))
				
				print(Make.3Dplot(Station.frames[[station]][[prof]], Beads.frames[[station]][[prof]], Noise.frames[[station]][[prof]], "SSC.A", "SybrGreen.A", "PE.A", "log(SSC.A) [arbitratry unit]", "log(SyberGreen.A) [arbitratry unit]", "log(PE.A) [arbitratry unit]", paste("3D plot of the sample ", liste.stations[station], "_", Smp.depth[index], sep="")))
				print(writeWebGL(filename = paste(path3D1, "/",today, "_3Dplot_SSC-SybrGreen-PE_", liste.stations[station], "_", Smp.depth[index], ".html", sep="")))
				
				print(Make.3Dplot(Station.frames[[station]][[prof]], Beads.frames[[station]][[prof]], Noise.frames[[station]][[prof]], "Chlorophyll.A", "SybrGreen.A", "PE.A", "log(Chlorophyll.A) [arbitratry unit]", "log(SyberGreen.A) [arbitratry unit]", "log(PE.A) [arbitratry unit]", paste("3D plot of the sample ", liste.stations[station], "_", Smp.depth[index], sep="")))
				print(writeWebGL(filename = paste(path3D2, "/",today, "_3Dplot_Chlorophyll-SybrGreen-PE_", liste.stations[station], "_", Smp.depth[index], ".html", sep="")))
				
				
			}
		
		
	}

