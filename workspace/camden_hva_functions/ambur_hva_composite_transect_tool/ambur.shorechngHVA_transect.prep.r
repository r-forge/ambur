ambur.shorechngHVA_transect.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the post-AMBUR HVA ambur_hva_calculations_pts.shp...")
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
###########



 out <- SpatialPoints(cbind(attrtable[,"Max_DtX"],attrtable[,"Max_DtY"])) #calculate point coordinates of youngest shoreline 




Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Preparing data... (%s)", info), info)

######HVA index setup

attrtable$hva_cat <- attrtable$hv_scld         

hva_cat <- attrtable$hva_cat 

esi_hva_table <- cbind(attrtable$Transct, hva_cat)

colnames(esi_hva_table) <- c("Transect","hva_cat")

write.table(esi_hva_table, file = "ambur_hva_shorechng_HVAdata.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


 pts.output <- SpatialPointsDataFrame(out,attrtable)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_hva_shorechngHVA_pts", driver="ESRI Shapefile")

#write Google Earth KML file
#ambur_hva_shorechng_pts.kml <- spTransform(pts.output, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(ambur_hva_shorechng_pts.kml["hva_cat"],"ambur_hva_shorechng_pts.kml","hva_cat","KML")



 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




}