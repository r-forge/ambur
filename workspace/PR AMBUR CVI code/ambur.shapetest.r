ambur.shapetest <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select a shapefile to check its geometry and projection...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
#getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
#attrtable <- data.frame(shapedata)

print(gIsValid(shapedata, reason = TRUE))

print(ifelse(is.projected(shapedata)== TRUE,"Shapefile is projected","Shapefile is not projected. Please update projection info.") )

print(proj4string(shapedata)) # contains projection info





}