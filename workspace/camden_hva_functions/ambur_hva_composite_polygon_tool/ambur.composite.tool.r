ambur.composite.tool <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)
 
 
 tkmessageBox(message = "Please select the AMBUR prepped NOAA SOVI shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)


 tkmessageBox(message = "Please select the AMBUR HVA Inundation shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata2 <- tk_choose.files(filter = filters,multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
shapedata2@data <- shapedata2@data[,-(1:3)]  #get rid of transID and shoreID
attrtable2 <- data.frame(shapedata2)

     

#test <- gUnion(shapedata, shapedata2, byid=FALSE, id=NULL)



###########

####create spatialLines only

shapedataA <- as(shapedata, "SpatialPolygons")
shapedataB <- as(shapedata2, "SpatialPolygons")

int <- gIntersects(shapedata, shapedata2, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(seq(along=vec))), 50)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   


vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(shapedata[i,], shapedata2[int[,i],], byid=TRUE)  else 0 
setTkProgressBar(pb, i, sprintf("AMBUR: Intersecting data layers (%s)", info), info)

}



Pcnt.Complete <-  40                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Extracting data (%s)", info), info)


#cond <- lapply(vec, function(x) length(x) > 0)

cond <- lapply(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]

out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))


Pcnt.Complete <-  50                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Writing new polygons to shapefile (%s)", info), info)




transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
#POINT_X <-  data.frame(coordinates(out))$x
#POINT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(1,length(transID),1)

inter.data <- data.frame(transID,shoreID,sortID)
tran.data <- data.frame(shapedata@data)
shore.data <- data.frame(shapedata2@data)


tran.data$Id <- as.numeric(row.names(tran.data))

shore.data$mergeID <- as.numeric(row.names(shore.data))


tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)

tet3 <- tet2[ order(tet2[,"sortID"]) , ]
tet3$Id <- tet2[,"sortID"]

row.names(tet3) <- seq(1,length(tet3$Id),1)


#change the shape IDs to match the tet3 IDs
out2 <- spChFIDs(out, as.character(as.character(row.names(tet3))))

hva.calc <- sqrt((tet3$hva_cat * tet3$hva_scaled)/2)

tet3$hva_calc2 <- hva.calc

####scale hva.calc values back to 1 to 5
n.parameters <- 2
cat1 <- sqrt((1^n.parameters)/n.parameters)
cat2 <- sqrt((2^n.parameters)/n.parameters)
cat3 <- sqrt((3^n.parameters)/n.parameters)
cat4 <- sqrt((4^n.parameters)/n.parameters)
cat5 <- sqrt((5^n.parameters)/n.parameters)

break.pts <- c(cat1,cat2,cat3,cat4,cat5)

break.pts2 <- ((break.pts[-1] - break.pts[-5])/2) + break.pts[-5]

break.pts3 <- c(break.pts[1],break.pts2,break.pts[5])

#plot(break.pts3,asp=1)
#points(break.pts,col="blue")


hva_scaled <- numeric(length(hva.calc ))

hva_scaled[hva.calc <= break.pts3[2]] <- 1
hva_scaled[hva.calc > break.pts3[2] & hva.calc <= break.pts3[3]] <- 2
hva_scaled[hva.calc > break.pts3[3] & hva.calc <= break.pts3[4]] <- 3
hva_scaled[hva.calc > break.pts3[4] & hva.calc <= break.pts3[5]] <- 4
hva_scaled[hva.calc > break.pts3[5]] <- 5
#############################################

 tet3$hva_scaled2 <- hva_scaled




#tet3$Distance <- (((tet3$POINT_X - tet3$StartX)^2 +  (tet3$StartY - tet3$POINT_Y)^2)^(1/2))

#tet3$BASE_LOC <- ifelse(tet3$BASE_LOC == "NA", tet3$BaseOrder, tet3$BASE_LOC)

outputdata <- SpatialPolygonsDataFrame(out2,tet3)



Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Writing new polygons to shapefile (%s)", info), info)


 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(outputdata) <- projectionString

writeOGR(outputdata, ".", "composite_composite_results", driver="ESRI Shapefile")


Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new polygons to shapefile (%s)", info), info)


 #getSpPPolygonsIDSlots(out) 
 
 }