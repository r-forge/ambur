ambur.slosh.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the SLOSH MOM-AGL shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)



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
hva_cat <- numeric(length(row.names(shapedatacull) ))


hva_cat[shapetablecull[,"c5_high"] > 0 & shapetablecull[,"c5_high"] != 99.9] <- 1
hva_cat[shapetablecull[,"c4_high"] > 0 & shapetablecull[,"c4_high"] != 99.9] <- 2
hva_cat[shapetablecull[,"c3_high"] > 0 & shapetablecull[,"c3_high"] != 99.9] <- 3
hva_cat[shapetablecull[,"c2_high"] > 0 & shapetablecull[,"c2_high"] != 99.9] <- 4
hva_cat[shapetablecull[,"c1_high"] > 0 & shapetablecull[,"c1_high"] != 99.9] <- 5


shapetablecull$hva_cat <- hva_cat          



###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)



plys.output <- SpatialPolygonsDataFrame(shapedatacull,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(plys.output) <- projectionString

writeOGR(plys.output, ".", "ambur_hva_slosh_plys", driver="ESRI Shapefile")

 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)































}