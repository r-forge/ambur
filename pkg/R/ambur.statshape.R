ambur.statshape <-
function(npts=50) {

winDialog("ok","Please select the GIS_stats_table_short.csv file...")

getdata <- choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_gisdata", showWarnings=FALSE)
setwd("AMBUR_gisdata")

dir.create(paste(time.stamp2," ","gisdata",sep=""))
setwd(paste(time.stamp2," ","gisdata",sep=""))

require(shapefiles)

id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field),X=(mydata["Outer_X"]),Y=(mydata["Outer_Y"]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 1)
write.shapefile(ddShapefile, paste("outer_pts",sep=""), arcgis=T)

id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field),X=(mydata["Inner_X"]),Y=(mydata["Inner_Y"]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 1)
write.shapefile(ddShapefile, paste("inner_pts",sep=""), arcgis=T)

id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field),X=(mydata["Start_X"]),Y=(mydata["Start_Y"]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 1)
write.shapefile(ddShapefile, paste("start_pts",sep=""), arcgis=T)

id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field),X=(mydata["End_X"]),Y=(mydata["End_Y"]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 1)
write.shapefile(ddShapefile, paste("end_pts",sep=""), arcgis=T)

id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field),X=(mydata["Max_DateX"]),Y=(mydata["Max_DateY"]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 1)
write.shapefile(ddShapefile, paste("max_date_pts",sep=""), arcgis=T)


id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field[,],id.field[,]),X=c(mydata["Outer_X"][,],mydata["Inner_X"][,]),Y=c(mydata["Outer_Y"][,],mydata["Inner_Y"][,]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable,"Transect", 3)
write.shapefile(ddShapefile, paste("transects_envelope",sep=""), arcgis=T)


id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field[,],id.field[,]),X=c(mydata["Start_X"][,],mydata["End_X"][,]),Y=c(mydata["Start_Y"][,],mydata["End_Y"][,]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable,"Transect", 3)
write.shapefile(ddShapefile, paste("transects_original",sep=""), arcgis=T)


id.field <- mydata["Transect"]
dd <- data.frame(Id=c(id.field[,],id.field[,]),X=c(mydata["Min_DateX"][,],mydata["Max_DateX"][,]),Y=c(mydata["Min_DateY"][,],mydata["Max_DateY"][,]))
ddTable <- data.frame(mydata)
ddShapefile <- convert.to.shapefile(dd, ddTable,"Transect", 3)
write.shapefile(ddShapefile, paste("transects_min_max_date",sep=""), arcgis=T)



#tidy up and remove all objects
detach("package:shapefiles")
detach(mydata)
rm(list = ls())

}

