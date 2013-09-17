ambur.shoregeo.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select a post-AMBUR analysis transect shapefile generated from ambur.statshape function...")
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



 out <- SpatialPoints(cbind(attrtable[,"Max_DateX"],attrtable[,"Max_DateY"])) #calculate point coordinates of youngest shoreline 




Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Preparing data... (%s)", info), info)

######HVA index setup
ESI_values <- attrtable[,"MaxClass1"]

hva_cat <- numeric(length(ESI_values ))

ESIcat1 <- c("rocky","cliff")
ESIcat2 <- c("high bluffs")
ESIcat3 <- c("anthropogenic")
ESIcat4 <- c("upland","upland or apparent mean high water")
ESIcat5 <- c("marsh")


hva_cat[attrtable[,"MaxClass1"] %in% ESIcat1] <- 1
hva_cat[attrtable[,"MaxClass1"] %in% ESIcat2] <- 2
hva_cat[attrtable[,"MaxClass1"] %in% ESIcat3] <- 3
hva_cat[attrtable[,"MaxClass1"] %in% ESIcat4] <- 4
hva_cat[attrtable[,"MaxClass1"] %in% ESIcat5] <- 5


attrtable$hva_cat <- hva_cat          

esi_hva_table <- cbind(attrtable$Transect, hva_cat)

colnames(esi_hva_table) <- c("Transect","hva_cat")

write.table(esi_hva_table, file = "ambur_hva_shoregeo_data.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


 pts.output <- SpatialPointsDataFrame(out,attrtable)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_hva_shoregeo_pts", driver="ESRI Shapefile")



 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




}