###LOAD PACKAGES
	library(flowCore)
	library(ggcyto)
	library(stringr)
	library(car)
	library(rgl)
	library(flowPeaks)

###VARIABLES

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Figures/Plots/3D_plots/' #Path of the folder containing the PDF Files for the results
	csv.name <- "_Abundance_with_all_info_results.csv" #Name of the CSV file containing the results
	pdf.name <- "_Plots_with_gating.pdf" #Name of the pdf containing the plots with the gates
	liste.stations <- c('LB2') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	today <- '20180615'
	
	#MINIMAL NUMBER OF BEADS AND EVENT
	minEvents <- 9999 #minimal number of events
	minBeads <- 999 #minimal number of beads
	
	#TYPE OF TRANSFORMATION
	Transfo.type <- logTransform(transformationId="LogTransform", logbase=10, r=1, d=1) #Type of transformation
	to.transform <- c('FSC.A', 'SSC.A', 'Chlorophyll.A', 'SybrGreen.A', 'PE.A') #List of the names of the measurement parameters on which the transformations have to be performed
	
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


Sort.Files <- function(listouille, max.depth) {  ###This function returns the list of the sorted files

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

Cluster.Gating <- function(Flowframe, list.parameters){ ###This function performs automated gating based on a clustering approach

	barcode <- FFtoDF(Flowframe) #The input needs to be converted into a dataframe object

	list.toRemove <- c()

	for(para in 1:length(list.parameters)){
		
		list.toRemove <- append(list.toRemove, grep("Inf", barcode[,list.parameters[para]],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NaN", barcode[,list.parameters[para]],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NA", barcode[,list.parameters[para]],value = FALSE))
		
	}
	
	if(length(list.toRemove)>0){
	
		barcode <- barcode[-list.toRemove,] ### Remove rows containing "Inf", "NaN" and "NA" value
		
	}
	
	fp2 <- flowPeaks(barcode[,list.parameters])
	fpc <- assign.flowPeaks(fp2,fp2$x)
	barcode$gate <- paste("Gate ", fpc, sep="")
	
	for(NG in 1:13){

	barcode[barcode==paste("Gate -",NG, sep="")] <- "A Non-gated"
	
	}
	
	
	return(barcode)

}

Find.BandN.Gate <- function(barcode1){

	barcode <- barcode1
	
	barcode$gate <- as.factor(barcode$gate)
	
	coucou <- aggregate(SybrGreen.A ~ gate, barcode, FUN=mean)
	coucoutruc <- as.vector(coucou$gate)
	levels(barcode$gate)[levels(barcode$gate)==coucoutruc[which.max(coucou[,'SybrGreen.A'])]] <- "B Beads gate" ###determine the beads gate

	coucou <- aggregate(PE.A ~ gate, barcode, FUN=mean)
	coucoutruc <- as.vector(coucou$gate)
	levels(barcode$gate)[levels(barcode$gate)==coucoutruc[which.min(coucou[,'PE.A'])]] <- "C Noise gate" ###determine the noise gate

	return(barcode)
	
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

get.color <- function(Dataframe){ ###This function attribute a color to each gate

color1 <- levels(Dataframe$gate)
color1[color1=="A Non-gated"] <- "black"
color1[color1=="B Beads gate"] <- "orangered"
color1[color1=="C Noise gate"] <- "gray80"

other.gate.colors <- c("royalblue","royalblue1","royalblue2","royalblue3","royalblue4", "dodgerblue4", "dodgerblue3", "slateblue","slateblue1","slateblue2", "slateblue3", "slateblue4", "dodgerblue2") 

for(ind in 1:length(color1)){

color1[color1==paste("Gate ", ind, sep="")] <- other.gate.colors[ind]

}

return(color1)

}

###===============================================================================
	
	
###INITIALISATION DE LA LISTE DE TRANSFORMATION
	myTrans <- transformList(to.transform, Transfo.type)
	

###CREATION OF FLOWSET (ONE PER STATION)
	Station.frames <- c()
	
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
		
	}

gate.dataframe <- Cluster.Gating(Station.frames[[1]][[7]], c("SSC.A","SybrGreen.A", "PE.A","Chlorophyll.A"))

list.toRemove <- c()
		
		list.toRemove <- append(list.toRemove, grep("Inf", gate.dataframe[,"PE.A"],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NaN", gate.dataframe[,"PE.A"],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NA", gate.dataframe[,"PE.A"],value = FALSE))
		
	
	if(length(list.toRemove)>0){
	
		gate.dataframe <- gate.dataframe[-list.toRemove,] ### Remove rows containing "Inf", "NaN" and "NA" value
		
	}

gate.dataframe <- Find.BandN.Gate(gate.dataframe)

scatter3d(x = gate.dataframe[,"SSC.A"], y = gate.dataframe[,"Chlorophyll.A"], z = gate.dataframe[,"SybrGreen.A"], xlab="SSC.A", ylab="Chlorophyll.A", zlab="SybrGreen.A", sphere.size=0.1, groups = gate.dataframe$gate, axis.col=c("black","black","black"), surface.col=get.color(gate.dataframe), surface=FALSE)

