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
 
ttt <- gIntersection(shapedata2, shapedata, byid=TRUE)
POINT_X <-  coordinates(ttt)[,1]
POINT_Y <-  coordinates(ttt)[,2]

rn1 <- row.names(ttt) 
nrn1 <- data.frame(do.call("rbind", strsplit(rn1, " ")),POINT_X,POINT_Y)



dataf1 <- merge(nrn1,attrtable2,by.x="X1",by.y="row.names",all=FALSE) 
dataf2 <- merge(dataf1,attrtable,by.x="X2",by.y="row.names",all=FALSE) 




pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(seq(along=vec))), 50)


Pcnt.Complete <-  round(((i)/ length(seq(along=vec))) * 100, 0)/2
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="")
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i/2, sprintf("AMBUR: Capture shoreline positions (%s)", info), info)

tet3 <- dataf2

tet3$Id <- as.numeric(row.names(dataf2))

tet3$Distance <- (((tet3$POINT_X - tet3$StartX)^2 +  (tet3$StartY - tet3$POINT_Y)^2)^(1/2))

tet3$BASE_LOC <- ifelse(tet3$BASE_LOC == "NA", tet3$BaseOrder, tet3$BASE_LOC)

interpts <- SpatialPoints(data.frame(x=tet3$POINT_X,y=tet3$POINT_Y))


outputdata <- SpatialPointsDataFrame(interpts,tet3)

Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Capture shoreline positions (%s)", info), info) 


 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info
  
  proj4string(outputdata) <- projectionString
   
writeOGR(outputdata, ".", "shore_pts", driver="ESRI Shapefile")


Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Capture shoreline positions (%s)", info), info)


} 