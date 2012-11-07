ambur.groupcull <-
function(indv=1) {

require(tcltk)
require(rgdal)
require(rgeos)

#choose file to import
tkmessageBox(message = "Please select the capture points shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
mydata <- data.frame(shapedata)


workingdir <- dirname(getdata)
setwd(workingdir)

#Intersect.Position <- "max"

#make sure all column names are upper case
colnames(mydata) <- toupper(colnames(mydata))

mydata[,"ID"] <- seq(1:length(mydata[,"ID"]))

 ##########add a distance field if one doesn't exist (for use with ARCGIS intersection)

crdl0 <- coordinates(shapedata)

basexa <- crdl0[,1]
baseya <- crdl0[,2]


mydata$COORDS.X1 - crdl0[,1]




Cx <- mydata$STARTX
Cy <- mydata$STARTY

Cx2 <- mydata$COORDS.X1
Cy2 <- mydata$COORDS.X2

ptdisttoorgin <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

mydata$DISTANCE <- ptdisttoorgin


colnames(mydata) <- gsub("COORDS.X1", "X_COORD", colnames(mydata))
colnames(mydata) <- gsub("COORDS.X2", "Y_COORD", colnames(mydata))
##########################################################################################

test <- data.frame(mydata[ ,"TRANSECT"], mydata[ ,"GROUP"],mydata[ ,"DISTANCE"],mydata[ ,"ID"])[ order(mydata[ ,"TRANSECT"], mydata[ ,"GROUP"],mydata[ ,"DISTANCE"]),]   #sort to get proper order

colnames(test) <- c("TRANSECT",  "GROUP", "DISTANCE" ,"ID")


testblah <- data.frame(test)

testblah2 <- aggregate(. ~ testblah$TRANSECT+testblah$GROUP,data = testblah[,3:4],FUN=function(x) c(mn =min(x), n=x[1] ) )  #get first occurence of the min distance ID
colnames(testblah2) <- c("TRANSECT",  "GROUP", "DISTANCE" ,"DISTANCE2")

testblah3 <- testblah2[,4]

#aa <- as.character(testblah3)
#unlist(strsplit(aa, " .*"))
######################################################################################################




 valid.pts <-  as.numeric(testblah3[,2])

 shapedatacull <- shapedata[valid.pts,]
 shapetablecull <-  data.frame(shapedatacull)

 #mydata[,"ID"] %in% valid.pts    mydata[valid.pts+1,]

 plot(shapedatacull)

 pts.output <- SpatialPointsDataFrame(shapedatacull,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_group_pts_prep0", driver="ESRI Shapefile")

 plot(mydata$STARTX, mydata$STARTY,asp=1)
 points(shapedatacull)
 
 }


