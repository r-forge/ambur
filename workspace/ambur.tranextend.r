ambur.tranextend <-
function(start.extend=50, end.extend=100) {

require(tcltk)
 require(rgdal)
require(rgeos)


#open ambur plotting file
#choose dbf file to import
tkmessageBox(message = "Please select an ambur generated transects shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)

#mydata <- as(shapedata, "SpatialLinesDataFrame")

mydata <- data.frame(shapedata)

colnames(mydata) <- gsub("Start_X", "StartX", colnames(mydata))
colnames(mydata) <- gsub("End_X", "EndX", colnames(mydata))
colnames(mydata) <- gsub("Start_Y", "StartY", colnames(mydata))
colnames(mydata) <- gsub("End_Y", "EndY", colnames(mydata))
colnames(mydata) <- gsub("T_azimuth", "Azimuth", colnames(mydata))

workingdir <- dirname(getdata)
setwd(workingdir)
path <- getdata



time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)



crdl0 <- coordinates(shapedata)
crdl1 <- sapply(crdl0, function(x) do.call("rbind", x))
crdl2 <- data.frame(t(crdl1))

colnames(crdl2) <- c("x.start","x.end","y.start","y.end")

endx.extend <- sin((mydata[,"Azimuth"] * pi/180)) * (end.extend) + crdl2[,"x.end"]
endy.extend <- cos((mydata[,"Azimuth"] * pi/180)) * (end.extend) + crdl2[,"y.end"]
startx.extend <- sin((mydata[,"Azimuth"] * pi/180)) * (start.extend)*-1 + crdl2[,"x.start"]
starty.extend <- cos((mydata[,"Azimuth"] * pi/180)) * (start.extend)*-1 + crdl2[,"y.start"]

plot(crdl2[,1],crdl2[,3],asp=T,col=0)
segments(startx.extend,starty.extend,endx.extend,endy.extend, col="red")
segments(crdl2[,1],crdl2[,3],crdl2[,2],crdl2[,4], col="green")



################################################################################################################
new_trandata <- data.frame(mydata)

new_trandata$StartX <- startx.extend
new_trandata$StartY <- starty.extend
new_trandata$EndX <- endx.extend
new_trandata$EndY <- endy.extend


row.names(new_trandata) <- new_trandata$Transect

Transect.Factor <- factor(new_trandata$Transect)   

shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$StartX[new_trandata$Transect == x], new_trandata$EndX[new_trandata$Transect == x]), y=c(new_trandata$StartY[new_trandata$Transect == x],new_trandata$EndY[new_trandata$Transect == x])))), ID=(as.numeric(x))))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)

new_trandata$TranDist <- as.vector(SpatialLinesLengths(shape.final2, longlat=FALSE))     #get new line lengths

shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)

  # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info
  
  proj4string(shape.final3) <- projectionString

  outputname <- paste(shapename,"_extended",sep="")

writeOGR(shape.final3, ".", outputname, driver="ESRI Shapefile")




}