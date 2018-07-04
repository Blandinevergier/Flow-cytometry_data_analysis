###LOAD PACKAGES
	library(flowCore)
	library(ggcyto)
	library(stringr)
	library(car)
	library(rgl)
	library(flowPeaks)
	library(gridExtra)

###VARIABLES

	Folder.path <- '/Users/bland/Desktop/Flow-cytometry_data/Input2/' #Path of the folder containing the FCS Files
	csv.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Dataframe/' #Path of the folder containing the csv file for the results 
	img.path <- '/Users/bland/Desktop/Flow-cytometry_data/Output/Figures/Plots/' #Path of the folder containing the PDF Files for the results
	csv.name <- "_Summary_after_automatic_gating.csv" #Name of the CSV file containing the results
	pdf.name <- "_Plots_with_gating.pdf" #Name of the pdf containing the plots with the gates
	#liste.stations <- c('FX1', 'FX2', 'FX4', 'FX5', 'FX6', 'FX7','FX8', 'FX9','HB1','HB2','HB3','HB4','HB5','IH1','IH2', 'IH3','KG1','KG2','KG3','KG4','KG5','KG6','KR1','KR2','KR3','KR4','KR5','KR6','LA1','LA2','LA3','LA4','LA5','LA6','LA7','LA8','LB1','LB2','LB3','LB4','LB5','LB6','LB7','LB8','LB9','LB10','LN1','LN2','LN3','LN4','LN5','LN6','MS1','MS2','MS3','MS4','MS5','MS6','SB2','SB3','SB4','SI1','SI2','SI3','SI4','SI5','SI6','SI7','SI8','ST1','ST2','ST3','ST5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	#liste.stations <- c('SI1','SI2','SI3','SI4','SI5','SI6','SI7','SI8','ST1','ST2','ST3','ST5') #List of the keywords of the stations to analyse ###be sure that all the FCS files corresponding to the stations are in the folder and that the keywords correspond to a unique station
	liste.stations <- c('LB4')
	
	today <- '20180702'
	
	#MINIMAL NUMBER OF BEADS AND EVENT
	minEvents <- 9999 #minimal number of events
	minBeads <- 999 #minimal number of beads
	
	#TYPE OF TRANSFORMATION
	Transfo.quad <- quadraticTransform(transformationId="defaultQuadraticTransform", a = 1, b = 1, c = 0) #Type of transformation
	Transfo.log <- logTransform(transformationId="LogTransform", logbase=10, r=1, d=1) #Type of transformation
	Transfo.lin <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=1, b=1, c=0)
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

Remove.NA <- function(DF.withNA){

	DF.withoutNA <- DF.withNA

		list.toRemove <- c()
		
		for(ii in 1:ncol(DF.withoutNA)){
		list.toRemove <- append(list.toRemove, grep("Inf", DF.withoutNA[,ii],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NaN", DF.withoutNA[,ii],value = FALSE))
		list.toRemove <- append(list.toRemove, grep("NA", DF.withoutNA[,ii],value = FALSE))
		}
	
	
	if(length(list.toRemove)>0){
	
		DF.withoutNA <- DF.withoutNA[-list.toRemove,] ### Remove rows containing "Inf", "NaN" and "NA" value
		
	}
	
	return(DF.withoutNA)
	
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

###===============================================================================
	
	
###INITIALISATION DE LA LISTE DE TRANSFORMATION
	myTranslog <- transformList(to.transform, Transfo.log)
	myTranslin <- transformList(to.transform, Transfo.lin)
	myTransquad <- transformList(to.transform, Transfo.quad)
	

###CREATION OF FLOWSET (ONE PER STATION)
	Station.frames.log <- c()
	Station.frames.quad <- c()
	Station.frames.lin <- c()
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
		Station.frames <- append(Station.frames, fs)
		
		###TRANSFORMATION
		fs.transquad <- transform(fs, myTransquad) #Transformation of the data
		Station.frames.quad <- append(Station.frames.quad, fs.transquad)
		
		fs.translog <- transform(fs, myTranslog) #Transformation of the data
		Station.frames.log <- append(Station.frames.log, fs.translog)
		
		fs.translin <- transform(fs, myTranslin) #Transformation of the data
		Station.frames.lin <- append(Station.frames.lin, fs.translin)
		
	}
	
	

DF <- FFtoDF(Station.frames[[1]][[2]])

DF <- Remove.NA(DF)

DFlog <- FFtoDF(Station.frames.log[[1]][[2]])

DFlog <- Remove.NA(DFlog)

DFquad <- FFtoDF(Station.frames.quad[[1]][[2]])

DFquad <- Remove.NA(DFquad)

DFlin <- FFtoDF(Station.frames.lin[[1]][[2]])

DFlin <- Remove.NA(DFlin)

DF2 <- FFtoDF(Station.frames[[1]][[3]])

DF2 <- Remove.NA(DF2)

DFlog2 <- FFtoDF(Station.frames.log[[1]][[3]])

DFlog2 <- Remove.NA(DFlog2)

DFquad2 <- FFtoDF(Station.frames.quad[[1]][[3]])

DFquad2 <- Remove.NA(DFquad2)

DFlin2 <- FFtoDF(Station.frames.lin[[1]][[3]])

DFlin2 <- Remove.NA(DFlin2)

DFbl <- FFtoDF(Station.frames[[1]][[1]])

DFbl <- Remove.NA(DFbl)

DFlogbl <- FFtoDF(Station.frames.log[[1]][[1]])

DFlogbl <- Remove.NA(DFlogbl)

DFquadbl <- FFtoDF(Station.frames.quad[[1]][[1]])

DFquadbl <- Remove.NA(DFquadbl)

DFlinbl <- FFtoDF(Station.frames.lin[[1]][[1]])

DFlinbl <- Remove.NA(DFlinbl)

	
plot1 <- ggplot(DF, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "blue") + labs(x = 'SSC.A [A.U]', y = 'Chlorophyll.A [A.U]') + ggtitle("0 meter", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))

 

plot2 <- ggplot(DFlog, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "red") + labs(x = 'log(SSC.A) [A.U]', y = 'log(Chlorophyll.A) \n [A.U]')  + ggtitle("0 meter", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot3 <- ggplot(DFquad, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "orange") + labs(x = 'quadraticTransform(SSC.A) [A.U]', y = 'quadraticTransform(Chlorophyll.A) \n [A.U]') + ggtitle("0 meter", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot4 <- ggplot(DFlin, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "green3") + labs(x = 'arcsinhTransform(SSC.A) [A.U]', y = 'arcsinhTransform(Chlorophyll.A) \n [A.U]') + ggtitle("0 meter", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot12 <- ggplot(DF2, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "blue") + labs(x = 'SSC.A [A.U]', y = 'Chlorophyll.A [A.U]') + ggtitle("-10 meters", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))


plot22 <- ggplot(DFlog2, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "red") + labs(x = 'log(SSC.A) [A.U]', y = 'log(Chlorophyll.A) \n [A.U]') + ggtitle("-10 meters", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot32 <- ggplot(DFquad2, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "orange") + labs(x = 'quadraticTransform(SSC.A) [A.U]', y = 'quadraticTransform(Chlorophyll.A) \n [A.U]')  + ggtitle("-10 meters", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot42 <- ggplot(DFlin2, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "green3") + labs(x = 'arcsinhTransform(SSC.A) [A.U]', y = 'arcsinhTransform(Chlorophyll.A) \n [A.U]') + ggtitle("-10 meters", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot1bl <- ggplot(DFbl, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "blue") + labs(x = 'SSC.A [A.U]', y = 'Chlorophyll.A [A.U]') + ggtitle("BLANK", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))


plot2bl <- ggplot(DFlogbl, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "red") + labs(x = 'log(SSC.A) [A.U]', y = 'log(Chlorophyll.A) \n [A.U]') + ggtitle("BLANK", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot3bl <- ggplot(DFquadbl, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "orange") + labs(x = 'quadraticTransform(SSC.A) [A.U]', y = 'quadraticTransform(Chlorophyll.A) \n [A.U]') + ggtitle("BLANK", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



plot4bl <- ggplot(DFlinbl, aes(SSC.A, Chlorophyll.A))	 +
 geom_point(colour = "green3") + labs(x = 'arcsinhTransform(SSC.A) [A.U]', y = 'arcsinhTransform(Chlorophyll.A) \n [A.U]')  + ggtitle("BLANK", subtitle = NULL) + theme(plot.title = element_text(size=30, face="plain", hjust=0.5, vjust=0, colour="grey47"), axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))



png.file <- png(paste(img.path,"LB4_Best_Transformation.png",sep =""), width=900, height=1100)
		
png.file
grid.arrange(plot1bl, plot1,plot12, plot3bl, plot3, plot32, plot4bl, plot4,plot42, plot2bl, plot2,plot22, nrow=4, ncol=3)
dev.off()
