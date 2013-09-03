ambur.sovi_inland.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the NOAA SOVI shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, 100, 10)

 ########
time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
setwd("AMBUR_HVA_data_prep")

#dir.create(paste(time.stamp2," ","shorepts",sep=""))
#setwd(paste(time.stamp2," ","shorepts",sep=""))

####create spatialLines only

#shapedataA <- as(shapedata, "SpatialLines")
#shapedataB <- as(shapedata2, "SpatialLines")

 shapedatacull <- shapedata
 shapetablecull <-  attrtable
 
 ######HVA index setup
hva_cat <- numeric(length(shapetablecull[,1] ))

hva_cat[shapetablecull[,"SoVI"] > 5] <- 5
hva_cat[shapetablecull[,"SoVI"] > 1 & shapetablecull[,"SoVI"] < 5] <- 4
hva_cat[shapetablecull[,"SoVI"] > -1 & shapetablecull[,"SoVI"] < 1] <- 3
hva_cat[shapetablecull[,"SoVI"] > -5 & shapetablecull[,"SoVI"] < -1] <- 2
hva_cat[shapetablecull[,"SoVI"] < -5] <- 1

shapetablecull$hva_cat <- hva_cat          



###########



setTkProgressBar(pb,90,"Still working...","AMBUR: Writing new polygons to shapefile")



plys.output <- SpatialPolygonsDataFrame(shapedatacull,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(plys.output) <- projectionString

writeOGR(plys.output, ".", "ambur_hva_sovi_plys", driver="ESRI Shapefile")

setTkProgressBar(pb,100,"Done!", "AMBUR: Writing new polygons to shapefile")































}