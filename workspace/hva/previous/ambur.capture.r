ambur.capture <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the shoreline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)


tkmessageBox(message = "Please select the transect shapefile...")
getdata2 <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))

shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)

 ########
time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_capture", showWarnings=FALSE)
setwd("AMBUR_capture")

dir.create(paste(time.stamp2," ","shorepts",sep=""))
setwd(paste(time.stamp2," ","shorepts",sep=""))


###########

####create spatialLines only

shapedataA <- as(shapedata, "SpatialLines")
shapedataB <- as(shapedata2, "SpatialLines")

int <- gIntersects(shapedata2, shapedata, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(seq(along=vec))), 50)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb, i, sprintf("AMBUR: Capture shoreline positions (%s)", info), info)

vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(shapedata2[i,], shapedata[int[,i],], byid=TRUE)  else 0 }






#cond <- lapply(vec, function(x) length(x) > 0)

cond <- lapply(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]

out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))


Pcnt.Complete <-  50                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
POINT_X <-  data.frame(coordinates(out))$x
POINT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(0,length(POINT_X)-1,1)

inter.data <- data.frame(POINT_X,POINT_Y,transID,shoreID,sortID)
tran.data <- data.frame(shapedata2@data)
shore.data <- data.frame(shapedata@data)


tran.data$Id <- as.numeric(row.names(tran.data))

shore.data$mergeID <- as.numeric(row.names(shore.data))


tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)

tet3 <- tet2[ order(tet2[,"sortID"]) , ]

tet3$Id <- tet2[,"sortID"]

#tet3$Distance <- (((tet3$POINT_X - tet3$StartX)^2 +  (tet3$StartY - tet3$POINT_Y)^2)^(1/2))

#tet3$BASE_LOC <- ifelse(tet3$BASE_LOC == "NA", tet3$BaseOrder, tet3$BASE_LOC)

outputdata <- SpatialPointsDataFrame(out,tet3)



Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(outputdata) <- projectionString

writeOGR(outputdata, ".", "shore_pts6", driver="ESRI Shapefile")


Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)


}