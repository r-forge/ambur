require(raster)
require(rgdal) # also loads sp and lattice packages
 require(tcltk)
tkmessageBox(message = "Please select the ESRI grid header file...")

getdata <- tk_choose.files(default = " ",multi = FALSE)
shapename <- gsub(".", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)

 
   shapedata <- raster(workingdir)

spplot(shapedata)

nrow(shapedata) * ncol(shapedata)

zz1 <- getValues(shapedata,1,(nrow(shapedata)/2))
xy1 <- xyFromCell(shapedata,1:length(zz1))

write.csv(cbind(xy1,zz1),row.names=FALSE,"xyz_1.txt") 

rm(zz1)
rm(xy1)

zz1 <- getValues(shapedata,(nrow(shapedata)/2)+1,nrow(shapedata))
xy1 <- xyFromCell(shapedata,(length(zz1)+1):(length(zz1)*2))
write.csv(cbind(xy1,zz1),row.names=FALSE,"xyz_2.txt") 

rm(zz1)
rm(xy1)

zz1 <- getValues(shapedata,1,5)
xy1 <- xyFromCell(shapedata,1:length(zz1))
write.csv(cbind(xy1,zz1),row.names=FALSE,"xyz_small.txt")  
  
  
