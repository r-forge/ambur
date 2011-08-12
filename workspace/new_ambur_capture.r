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



int <- gIntersects(shapedata2, shapedata, byid=TRUE) 
vec <- vector(mode="list", length=dim(int)[2]) 
for (i in seq(along=vec)) vec[[i]] <- gIntersection(shapedata2[i,], shapedata[int[,i],], byid=TRUE) 
out <- do.call("rbind", vec) 
rn <- row.names(out) 
nrn <- do.call("rbind", strsplit(rn, " ")) 




transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
POINT_X <-  data.frame(coordinates(out))$x
POINT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(1,length(POINT_X),1)

inter.data <- data.frame(POINT_X,POINT_Y,transID,shoreID,sortID)
tran.data <- data.frame(shapedata2@data)
shore.data <- data.frame(shapedata@data)


tran.data$Id <- as.numeric(row.names(tran.data))

shore.data$mergeID <- as.numeric(row.names(shore.data))


tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)

tet3 <- tet2[ order(tet2[,"sortID"]) , ]

tet3$Id <- tet2[,"sortID"]   

tet3$Distance <- (((tet3$POINT_X - tet3$StartX)^2 +  (tet3$StartY - tet3$POINT_Y)^2)^(1/2))

outputdata <- SpatialPointsDataFrame(out,tet3)
      

   
writeOGR(outputdata, ".", "shore_pts", driver="ESRI Shapefile")