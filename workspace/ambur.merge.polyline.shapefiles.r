###############################################################################


library(rgdal)
library(maptools)
library(rgeos)
library(tcltk)



tkmessageBox(message = "Shapefiles must be in the same directory.")

filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- list(tk_choose.files("","Choose file",multi = TRUE,filetype,1))
#shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata[[1]][1])
setwd(workingdir)



# obtain shapefiles in current directory
files <- getdata

uid<-1

# Get polylines
#-------------------------------------

poly.data<- readOGR(files[[1]][1], gsub(".shp","",basename(files[[1]][1])))
n <- length(slot(poly.data, "lines"))
poly.data <- spChFIDs(poly.data, as.character(uid:(uid+n-1)))
uid <- uid + n

shapename <- rep(gsub(".shp","",basename(files[[1]][1])),n) ###get shapenames to add to append to the column

# ]combine remaining  polylines with first polyline
#-----------------------------------------------------------------
pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(files[[1]])), 50)

for (i in 2:length(files[[1]])) {
     temp.data <- readOGR(files[[1]][i], gsub(".shp","",basename(files[[1]][i])))
     n <- length(slot(temp.data, "lines"))
     temp.data <- spChFIDs(temp.data, as.character(uid:(uid+n-1)))
     uid <- uid + n
     poly.data <- spRbind(poly.data,temp.data)
     temp.name <- rep(gsub(".shp","",basename(files[[1]][i])),n)
     shapename <- c(shapename,temp.name)

Pcnt.Complete <-  round(((i)/ length(files[[1]])) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb, i, sprintf("AMBUR: Merging polyline shapefiles (%s)", info), info)
     
     
}

#names(poly.data)
outputdata <- poly.data

outputdata$shapename <- shapename

   projectionString <- proj4string(poly.data) # contains projection info
  
  proj4string(outputdata) <- projectionString
   
writeOGR(outputdata, ".", "all_polylines_merged", driver="ESRI Shapefile")



