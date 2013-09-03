ambur.spatvar.prep <-
function(bufferdist=100) {

bufferdist <- 100

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

 outptsdf <- SpatialPointsDataFrame(out,attrtable)



out2 <- gBuffer(out, byid=TRUE, id=NULL, width=bufferdist, quadsegs=5, capStyle="ROUND",joinStyle="ROUND", mitreLimit=1.0)

out2a.df <- data.frame(attrtable$Transect)
colnames(out2a.df) <- "Transect" 
out2a <- SpatialPolygonsDataFrame(out2,out2a.df)
proj4string(out2a) <- projectionString
writeOGR(out2a, ".", "ambur_hva_spatvar_polys", driver="ESRI Shapefile")


#out3 <- over(out2, outptsdf[,"Range_Dst"], returnList = TRUE, fn = mean)

out3 <- over(out2, outptsdf[,"Range_Dst"], returnList = FALSE, fn = sd)


Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Preparing data... (%s)", info), info)


plot(out3[,1])
hist(out3[,1])
quantile(out3[,1],na.rm=TRUE)
mean(out3[,1],na.rm=TRUE)

SpatVar <- out3[,1]

#plot(attrtable[,"Range_Dst"])
#plot(abs(attrtable[,"Net_Chng"])/attrtable[,"Range_Dst"])
#plot(((abs(attrtable[,"Net_Chng"])/attrtable[,"Range_Dst"]) * 100) > 100)
#TempVarPct <- ((abs(attrtable[,"Net_Chng"])/attrtable[,"Range_Dst"]) * 100)

######HVA index setup
ESI_values <- SpatVar

hva_cat <- numeric(length(ESI_values ))

hva_cat[SpatVar > 20] <- 5
hva_cat[SpatVar > 15 & SpatVar <= 20] <- 4
hva_cat[SpatVar > 10 & SpatVar <= 15] <- 3
hva_cat[SpatVar > 5 & SpatVar <= 10] <- 2
hva_cat[SpatVar <= 5] <- 1

attrtable$SpatVar <- SpatVar
attrtable$hva_cat <- hva_cat          

esi_hva_table <- cbind(attrtable$Transect, hva_cat)

colnames(esi_hva_table) <- c("Transect","hva_cat")

write.table(esi_hva_table, file = "ambur_hva_spatvar_data.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


 pts.output <- SpatialPointsDataFrame(out,attrtable)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_hva_spatvar_pts", driver="ESRI Shapefile")



 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




}