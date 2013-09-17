ambur.shorechng.prep <-
function(bufferdist=150) {


# Establish the inputs



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

######HVA index setup for shoreline change rate rank
EPR_values <- attrtable[,"EPR"]

scr_rank <- character(length(EPR_values ))

scr_rank[attrtable[,"EPR"] >= 1] <- "high"
scr_rank[attrtable[,"EPR"] > 0.2 & attrtable[,"EPR"] < 1] <- "medium"
scr_rank[attrtable[,"EPR"] >= -0.2 & attrtable[,"EPR"] <= 0.2] <- "low"
scr_rank[attrtable[,"EPR"] > -1 & attrtable[,"EPR"] < -0.20] <- "medium"
scr_rank[attrtable[,"EPR"] <= -1] <- "high"

attrtable$scr_rank <- scr_rank          


######HVA index setup for absolute value of shoreline change distance variability rank (coef of variation at single transect (uses only the mean at the transect and not the entire population)) NOTE: shoreline change rates could not be used because it is interval data and not ratio data.  Raw distance data is ratio data as long at is is non-negative
tvar_values <- attrtable[,"CVErasDst"]

tscrv_rank <- character(length(tvar_values  ))

tscrv_rank[tvar_values >= 1] <- "high"
tscrv_rank[tvar_values >= 0.5 & tvar_values < 1] <- "medium"
tscrv_rank[tvar_values < 0.5] <- "low"

attrtable$tscrv_rank <- tscrv_rank          

#################get spatial variability at a given buffer distance
 outptsdf <- SpatialPointsDataFrame(out,attrtable)



out2 <- gBuffer(out, byid=TRUE, id=NULL, width=bufferdist, quadsegs=5, capStyle="ROUND",joinStyle="ROUND", mitreLimit=1.0)

out2a.df <- data.frame(attrtable$Transect)
colnames(out2a.df) <- "Transect" 
out2a <- SpatialPolygonsDataFrame(out2,out2a.df)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(out2a) <- projectionString
writeOGR(out2a, ".", "ambur_hva_spatvar_polys", driver="ESRI Shapefile")


out3 <- over(out2, outptsdf[,"CVErasDst"], returnList = FALSE, fn = mean)

svar_values <- out3[,1]


sscrv_rank <- character(length(svar_values  ))

sscrv_rank[svar_values >= 1] <- "high"
sscrv_rank[svar_values >= 0.5 & svar_values < 1] <- "medium"
sscrv_rank[svar_values < 0.5] <- "low"

attrtable$sscrv_rank <- sscrv_rank 

##################HVA Category

hva_cat <- numeric(length(svar_values  ))

hva_cat[scr_rank == "high" & sscrv_rank == "low"] <- 5
hva_cat[(scr_rank == "high" | scr_rank == "medium") & (sscrv_rank == "medium" | sscrv_rank == "high")] <- 4
hva_cat[(scr_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "high")] <- 4
hva_cat[(scr_rank == "medium") & (tscrv_rank == "medium") & (sscrv_rank == "medium")] <- 3
hva_cat[(scr_rank == "low") & (tscrv_rank == "medium" | tscrv_rank == "high") & (sscrv_rank == "medium" | sscrv_rank == "high"  )] <- 2
hva_cat[(scr_rank == "low") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 1

attrtable$hva_cat <- hva_cat 

sc_hva_table <- cbind(attrtable$Transect, hva_cat,scr_rank,tscrv_rank,sscrv_rank)

colnames(sc_hva_table) <- c("Transect","hva_cat","scr_rank","tscrv_rank","sscrv_rank")

write.table(sc_hva_table, file = "ambur_hva_shorechng_data.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


 pts.output <- SpatialPointsDataFrame(out,attrtable)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_hva_shorechng_pts", driver="ESRI Shapefile")

#write Google Earth KML file
#ambur_hva_shorechng_pts.kml <- spTransform(pts.output, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(ambur_hva_shorechng_pts.kml["hva_cat"],"ambur_hva_shorechng_pts.kml","hva_cat","KML")



 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




}