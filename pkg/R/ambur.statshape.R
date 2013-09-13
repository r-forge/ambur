ambur.statshape <-
function(npts=50) {
#require(tcltk)
#require(rgdal)


tkmessageBox(message = "Please select the GIS_stats_table_short.csv file...")


getdata <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
#attach(mydata)

#mydata$ID <- seq(0,length(mydata$Transect)-1,1)

#replace characters in fields that GIS doesn't read
colnames(mydata) <- gsub(".", "_", colnames(mydata),fixed=TRUE)
colnames(mydata) <- gsub(" ", "", colnames(mydata),fixed=TRUE)


dir.path <- dirname(getdata)
setwd(dir.path)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_gisdata", showWarnings=FALSE)
setwd("AMBUR_gisdata")

dir.create(paste(time.stamp2," ","gisdata",sep=""))
setwd(paste(time.stamp2," ","gisdata",sep=""))
################################################################################################################



ddTable <- data.frame(mydata)
coordinates(ddTable) <- data.frame(x=(mydata["Outer_X"]), y=(mydata["Outer_Y"]))
writeOGR(ddTable, ".", "outer_pts", driver="ESRI Shapefile")


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["Inner_X"]), y=(mydata["Inner_Y"]))
writeOGR(ddTable, ".", "inner_pts", driver="ESRI Shapefile")


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["Start_X"]), y=(mydata["Start_Y"]))
writeOGR(ddTable, ".", "start_pts", driver="ESRI Shapefile")


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["End_X"]), y=(mydata["End_Y"]))
writeOGR(ddTable, ".", "end_pts", driver="ESRI Shapefile")


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["Max_DateX"]), y=(mydata["Max_DateY"]))
writeOGR(ddTable, ".", "max_date_pts", driver="ESRI Shapefile")


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["Min_DateX"]), y=(mydata["Min_DateY"]))
writeOGR(ddTable, ".", "min_date_pts", driver="ESRI Shapefile")

################################################################################################################
new_trandata <- data.frame(mydata)

row.names(new_trandata) <- new_trandata$Transect

Transect.Factor <- factor(new_trandata$Transect)    #fixed 20130224 to get proper order of transects to match LineIDs with row.names of new_trandata

shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Start_X[new_trandata$Transect == x], new_trandata$End_X[new_trandata$Transect == x]), y=c(new_trandata$Start_Y[new_trandata$Transect == x],new_trandata$End_Y[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
#edit(data.frame(getSLLinesIDSlots(shape.final2)) )
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
writeOGR(shape.final3, ".", "original_transects", driver="ESRI Shapefile")


shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Outer_X[new_trandata$Transect == x], new_trandata$Inner_X[new_trandata$Transect == x]), y=c(new_trandata$Outer_Y[new_trandata$Transect == x],new_trandata$Inner_Y[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
writeOGR(shape.final3, ".", "envelope_transects", driver="ESRI Shapefile")


shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$Min_DateX[new_trandata$Transect == x], new_trandata$Max_DateX[new_trandata$Transect == x]), y=c(new_trandata$Min_DateY[new_trandata$Transect == x],new_trandata$Max_DateY[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)
writeOGR(shape.final3, ".", "net_min_max__date_transects", driver="ESRI Shapefile")

################################################################################################################


#tidy up and remove all objects

#detach(mydata)
rm(list = ls())

}

