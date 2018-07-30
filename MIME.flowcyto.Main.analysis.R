###LOAD PACKAGES
	library(flowCore)
	library(ggcyto)
	library(stringr)
	library(car)
	library(rgl)
	library(flowPeaks)

###VARIABLES

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input2/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results 
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Figures/Plots/3D_plots/' #Path of the folder containing the PDF Files for the results
	csv.name <- "_Summary_after_automatic_gating.csv" #Name of the CSV file containing the results
	pdf.name <- "_Plots_with_gating.pdf" #Name of the pdf containing the plots with the gates
	liste.stations <- c('FX1', 'FX2', 'FX4', 'FX5', 'FX6', 'FX7','FX8', 'FX9','HB1','HB2','HB3','HB4','HB5','IH1','IH2', 'IH3','KG1','KG2','KG3','KG4','KG5','KG6','KR1','KR2','KR3','KR4','KR5','KR6','LA1','LA2','LA3','LA4','LA5','LA6','LA7','LA8','LB1','LB2','LB3','LB4','LB5','LB6','LB7','LB8','LB9','LB10','LN1','LN2','LN3','LN4','LN5','LN6','MS1','MS2','MS3','MS4','MS5','MS6','SB2','SB3','SB4','SI1','SI2','SI3','SI4','SI5','SI6','SI7','SI8','ST1','ST2','ST3','ST5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	
	today <- '20180730'

	
	#TYPE OF TRANSFORMATION
	Transfo.type <- logTransform(transformationId="LogTransform", logbase=10, r=1, d=1) #Type of transformation
	to.transform <- c('FSC.A', 'SSC.A', 'Chlorophyll.A', 'SybrGreen.A', 'PE.A') #List of the names of the measurement parameters on which the transformations have to be performed
	
	#INFORMATIONS ABOUT THE STATIONS
	Infos.path <- '/Users/bland/Desktop/Flow-cytometry_data/Summary_stations/hydrostndec.csv' #Path of the csv file containing the informations about the stations
	
	Infos1 <- read.csv(Infos.path, sep=';', header=TRUE)
	Infos <- as.matrix(Infos1)
	stn.lane.list <- Infos[,"stn.lane"]
	stn.name.list <- Infos[,"stn.name"]
	stn.id.list <- Infos[,"stn.id"]
	stn.lat.list <- Infos[,"stn.lat"]
	stn.lon.list <- Infos[,"stn.lon"]
	stn.max.depth.list <- Infos[,"stn.max.depth"]
	
	#MANUAL BEADS GATE
	BeadsSyb.min <- 5 #Minimal value of the beads gate in SybrGreen.A
	BeadsSyb.max <- 7 #Maximal value of the beads gate in SybrGreen.A
	BeadsChl.min <- 4.5 #Minimal value of the beads gate in Chlorophyll.A
	BeadsChl.max <- 5.5 #Maximal value of the beads gate in Chlorophyll.A
	BeadsSSC.min <- 4.5 #Minimal value of the beads gate in SSC.A
	BeadsSSC.max <- 5.5 #Maximal value of the beads gate in SSC.A

	
	
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

FFtoDF <-function(FF){ ###The function is able the convert a flowframe into a dataframe

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

	for(para in 1:length(to.transform)){
		
		list.toRemove <- append(list.toRemove, grep("Inf", barcode[,to.transform[para]],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NaN", barcode[,to.transform[para]],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NA", barcode[,to.transform[para]],value = FALSE))
		
	}
	
	if(length(list.toRemove)>0){
	
		barcode <- barcode[-list.toRemove,] ### Remove rows containing "Inf", "NaN" and "NA" value
		
	}
	
	fp2 <- flowPeaks(barcode[,list.parameters])
	fpc <- assign.flowPeaks(fp2,fp2$x)
	barcode$gate <- paste("Gate ", fpc, sep="")
	
	for(NG in 1:13){

	barcode[barcode==paste("Gate -",NG, sep="")] <- "Non-gated"
	
	}
	
	
	return(barcode)

}

Noise.MGating <- function(DF.withNoise){ ###This function is able to gate manually the noise in case of the clustering approach doesn't work

Dataframe.withNoise <- DF.withNoise

levels(Dataframe.withNoise$gate) <- c(levels(Dataframe.withNoise$gate), "Noise gate (manual gating)")

Dataframe.withNoise$gate[which(Dataframe.withNoise$SybrGreen.A < 2 & Dataframe.withNoise$SSC.A < 2 & Dataframe.withNoise$PE.A < 2 & Dataframe.withNoise$Chlorophyll.A < 3.15)] <- "Noise gate (manual gating)"

Dataframe.withNoise$gate <- as.factor(Dataframe.withNoise$gate)

return(Dataframe.withNoise)


}

Find.BandN.Gate <- function(barcode1){ #This function is able to sort the different gates detected by the clustering approach

	barcode <- barcode1
	
	barcode$gate <- as.factor(barcode$gate)
	
	coucou <- aggregate(SybrGreen.A ~ gate, barcode, FUN=mean)
	coucoutruc <- as.vector(coucou$gate)
	levels(barcode$gate)[levels(barcode$gate)==coucoutruc[which.max(coucou[,'SybrGreen.A'])]] <- "Beads gate" ###determine the beads gate
	
	
	coucou <- aggregate(PE.A ~ gate, barcode, FUN=mean)
	coucoutruc <- as.vector(coucou$gate)
	levels(barcode$gate)[levels(barcode$gate)==coucoutruc[which.min(coucou[,'PE.A'])]] <- "Noise gate" ###determine the noise gate

	
	indou <- grep("Noise gate", barcode$gate,value = FALSE)
	subdata <- barcode[indou,]
	propor <- sum(subdata$SSC.A < 2 & subdata$PE.A < 2 & subdata$Chlorophyll.A < 3.15)
	propor2 <- sum(barcode$SSC.A < 2 & barcode$PE.A < 2 & barcode$Chlorophyll.A < 3.15)
	
	if(propor < 0.7*nrow(subdata) | propor < 0.7*propor2){ #Quality control of the noise gate --> The gate will be determined manually if the clustering gating wasn't well performed

	levels(barcode$gate)[levels(barcode$gate)=="Noise gate"] <- "Non-gated"
	
	barcode <- Noise.MGating(barcode)
	
	}
	
	barcode$gate <- as.factor(barcode$gate)

	return(barcode)
	
	}
	
Manual.Gating <- function(DataF,SSCmin, SSCmax, Chlmin, Chlmax, Sybmin, Sybmax){

DataF$ManualGate <- "Gate -1"

levels(DataF$ManualGate) <- c("Gate -1", "Gate 1")

DataF$ManualGate[which(DataF$SybrGreen.A > Sybmin & DataF$SybrGreen.A < Sybmax & DataF$Chlorophyll.A > Chlmin & DataF$Chlorophyll.A < Chlmax & DataF$SSC.A > SSCmin & DataF$SSC.A < SSCmax)] <- "Gate 1"

return(DataF)
}	

Auto.Gating <- function(DataF){

DataF$AutoGate <- "Gate -1"

levels(DataF$AutoGate) <- c("Gate -1", "Gate 1")

DataF$AutoGate[DataF$gate == "Beads gate"] <- "Gate 1"

return(DataF)

}
	
get.color <- function(Dataframe){ ###This function attribute a color to each gate

color1 <- levels(Dataframe$gate)
color1[color1=="Non-gated"] <- "black"
color1[color1=="Beads gate"] <- "orangered"
color1[color1=="Noise gate"] <- "gray80"
color1[color1=="Noise gate (manual gating)"] <- "gray80"

other.gate.colors <- c("darkorchid4","brown4", "royalblue4", "chartreuse4", "hotpink4","mediumpurple4", "deepskyblue4", "slateblue4", "mediumorchid4", "purple4", "steelblue4", "navyblue", "deeppink4", "dodgerblue", "cyan4") 

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
	
	Samp.Name <- c() #List containing the names of all the analysed files
	Smp.depth <- c()
	
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

###INITIALISATION OF THE VECTORS CONTAINING THE DATA 	
beads.count <- c()
noise.count <- c()
total.count <- c()
stn.lane <- c()
stn.name <- c()
stn.id <- c()
stn.lat <- c()
stn.lon <- c()
stn.max.depth <- c()
Noise.gating.type <- c()
Eval.Cluster <- c()	

	
###FILLING OF THE DIFFERENT LISTS CONTAINING THE INFORMATIONS (station names, depths, localisation, abundance...)
	for (station in 1:length(liste.stations)) {
	
	ind <- match(liste.stations[station],stn.id.list)
	
		for (prof in 1:length(Station.frames[[station]])) {
				
				stn.lane <- append(stn.lane, stn.lane.list[ind])
				stn.name <- append(stn.name, stn.name.list[ind])
				stn.id <- append(stn.id, stn.id.list[ind])
				stn.lat <- append(stn.lat, stn.lat.list[ind])
				stn.lon <- append(stn.lon, stn.lon.list[ind])
				stn.max.depth <- append(stn.max.depth, stn.max.depth.list[ind])
				
				
				
				
			}
		
		
	}

Index <- 0

datepath <- paste(img.path, today, sep="")
print(dir.create(datepath))

for (station in 1:length(liste.stations)) {


path3D <- paste(datepath, "/", liste.stations[station], sep="")
path3D1 <- paste(path3D,"/SSC_vs_Chlorophyll_vs_SybrGreen", sep="")
path3D2 <- paste(path3D,"/SSC_vs_SYBRGREEN_vs_PE", sep="")
path3D3 <- paste(path3D,"/Chlorophyll_vs_SYBRGREEN_vs_PE", sep="")
	
	
 	print(dir.create(path3D))
	print(dir.create(path3D1))
	print(dir.create(path3D2))
	print(dir.create(path3D3))
	
	print(par3d(windowRect = 50 + c(0,0,940,740)))
	
	
		for (prof in 1:length(Station.frames[[station]])) {
		
		Index <- Index + 1

			gate.dataframe <- Cluster.Gating(Station.frames[[station]][[prof]], c("PE.A","SSC.A","SybrGreen.A","Chlorophyll.A"))

			gate.dataframe <- Find.BandN.Gate(gate.dataframe)
			
			gate.dataframe <- Manual.Gating(gate.dataframe,BeadsSSC.min, BeadsSSC.max, BeadsChl.min, BeadsChl.max, BeadsSyb.min, BeadsSyb.max)

			gate.dataframe <- Auto.Gating(gate.dataframe)
			
			Eval.Cluster <- append(Eval.Cluster, evalCluster(gate.dataframe$ManualGate,gate.dataframe$AutoGate,method="Vmeasure"))

			tr <- table(gate.dataframe$gate)

			beads.count <- append(beads.count, as.numeric(tr["Beads gate"]))
			
			if(is.element('Noise gate (manual gating)', levels(gate.dataframe$gate))==TRUE){

			noise.count <- append(noise.count, as.numeric(tr['Noise gate (manual gating)']))
			Noise.gating.type <- append(Noise.gating.type, "corrected")

			}else{
			
			noise.count <- append(noise.count, as.numeric(tr["Noise gate"]))
			Noise.gating.type <- append(Noise.gating.type, "Without correction")
			
			}

			
			
			total.count <- append(total.count, nrow(gate.dataframe))
			abond <- (nrow(gate.dataframe) - noise.count[Index] - as.numeric(tr["Beads gate"]))/((as.numeric(tr["Beads gate"]))/1080000)
			
			
				
				plot3dd <- scatter3d(x = gate.dataframe[,"SSC.A"], y = gate.dataframe[,"Chlorophyll.A"], z = gate.dataframe[,"SybrGreen.A"], xlab="log(SSC.A) [arbitratry unit]", ylab="log(Chlorophyll.A) [arbitratry unit]", zlab="log(SybrGreen.A) [arbitratry unit]", sphere.size=0.1, groups = gate.dataframe$gate, axis.col=c("black","black","black"), surface.col=get.color(gate.dataframe), surface=FALSE)
				+ legend3d("topright", legend = c(paste(stn.id[Index],"_", Smp.depth[Index], " (", total.count[Index], " events)", sep=""), " ", levels(gate.dataframe$gate), " ", paste("beads count : ", beads.count[Index], sep=""), paste("Abundance : ", abond, " events/mL", sep="")), pch = 16, col = c("white", "white", get.color(gate.dataframe), "white","white","white"), cex=1, inset=c(0.01))
				print(plot3dd)
				print(writeWebGL(filename = paste(path3D1, "/",today, "_3Dplot_SSC-Chlorophyll-SybrGreen_", liste.stations[station], "_", Smp.depth[Index], ".html", sep="")))
				
				plot3dd <- scatter3d(x = gate.dataframe[,"SSC.A"], y = gate.dataframe[,"PE.A"], z = gate.dataframe[,"SybrGreen.A"], xlab="log(SSC.A) [arbitratry unit]", ylab="log(PE.A) [arbitratry unit]", zlab="log(SybrGreen.A) [arbitratry unit]", sphere.size=0.1, groups = gate.dataframe$gate, axis.col=c("black","black","black"), surface.col=get.color(gate.dataframe), surface=FALSE)
				+ legend3d("topright", legend = c(paste(stn.id[Index],"_", Smp.depth[Index], " (", total.count[Index], " events)", sep=""), " ", levels(gate.dataframe$gate), " ", paste("beads count : ", beads.count[Index], sep=""), paste("Abundance : ", abond, " events/mL", sep="")), pch = 16, col = c("white", "white", get.color(gate.dataframe), "white","white","white"), cex=1, inset=c(0.01))
				print(plot3dd)
				print(writeWebGL(filename = paste(path3D2, "/",today, "_3Dplot_SSC-SYBRGREEN-PE_", liste.stations[station], "_", Smp.depth[Index], ".html", sep="")))
				
				plot3dd <- scatter3d(x = gate.dataframe[,"PE.A"], y = gate.dataframe[,"Chlorophyll.A"], z = gate.dataframe[,"SybrGreen.A"], xlab="log(PE.A) [arbitratry unit]", ylab="log(Chlorophyll.A) [arbitratry unit]", zlab="log(SybrGreen.A) [arbitratry unit]", sphere.size=0.1, groups = gate.dataframe$gate, axis.col=c("black","black","black"), surface.col=get.color(gate.dataframe), surface=FALSE)
				+ legend3d("topright", legend = c(paste(stn.id[Index],"_", Smp.depth[Index], " (", total.count[Index], " events)", sep=""), " ", levels(gate.dataframe$gate), " ", paste("beads count : ", beads.count[Index], sep=""), paste("Abundance : ", abond, " events/mL", sep="")), pch = 16, col = c("white", "white", get.color(gate.dataframe), "white","white","white"), cex=1, inset=c(0.01))
				print(plot3dd)
				print(writeWebGL(filename = paste(path3D3, "/",today, "_3Dplot_Chlorophyll-SYBRGREEN-PE_", liste.stations[station], "_", Smp.depth[Index], ".html", sep="")))
				
			
			
		}
		
		print(rgl.close())

}

Abundance <- (total.count - noise.count - beads.count)/(beads.count/1080000)


###CREATION OF A CSV FILE TO STORE THE DATA	
	csv.name <- paste(csv.path, today, csv.name, sep="")
	results <- cbind(Samp.Name, stn.lane, stn.name, stn.id, stn.lat, stn.lon, stn.max.depth, Smp.depth, total.count, beads.count, noise.count, Noise.gating.type, Eval.Cluster, Abundance)
	write.csv(results,csv.name)	


