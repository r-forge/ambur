ambur.q3firm.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the FEMA Q3/FIRM shapefile...")
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


 pb <- tkProgressBar("AMBUR: progress bar", "Preparing data...", 0, 100, 25)

shapetablecull <- attrtable
 
 ######HVA index setup
hva_cat <- numeric(length(shapetablecull[,"FLD_ZONE"] ))

hva_cat[shapetablecull[,"FLD_ZONE"] == "0.2 PCT ANNUAL CHANCE FLOOD HAZARD"] <- 1
hva_cat[shapetablecull[,"FLD_ZONE"] == "X"] <- 1
hva_cat[shapetablecull[,"FLD_ZONE"] == "X_500"] <- 1
hva_cat[shapetablecull[,"FLD_ZONE"] == "B"] <- 1
hva_cat[shapetablecull[,"FLD_ZONE"] == "C"] <- 1
hva_cat[shapetablecull[,"FLD_ZONE"] == "A"] <- 2
hva_cat[shapetablecull[,"FLD_ZONE"] == "AE"] <- 3
hva_cat[shapetablecull[,"FLD_ZONE"] == "AE (Floodway)"] <- 4
hva_cat[shapetablecull[,"FLD_ZONE"] == "AH"] <- 4
hva_cat[shapetablecull[,"FLD_ZONE"] == "VE"] <- 5
hva_cat[shapetablecull[,"FLD_ZONE"] == "V"] <- 5
hva_cat[shapetablecull[,"FLD_ZONE"] == "OPEN WATER"] <- 5

shapetablecull$hva_cat <- hva_cat          




###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new polygons (%s)", info), info)


plys.output <- SpatialPolygonsDataFrame(shapedata,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(plys.output) <- projectionString

writeOGR(plys.output, ".", "ambur_hva_q3_polys", driver="ESRI Shapefile")

 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)































}