ambur.rastertopoly <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1

require(raster)
require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the raster image...")
filters <- matrix(c("Imagine", ".img"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".img", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)

shapedata <- raster(getdata)

#bb <- data.frame(values(shapedata))



#raster.poly <-  rasterToPolygons(shapedata, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
#raster.tran <-  extract(shapedata, shapedata2,fun=mean, na.rm=TRUE, sp=TRUE, progress='text')


