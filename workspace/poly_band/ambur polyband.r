
ambur.polyband <- function(pwidth=50) {

#pwidth=500
require(shapefiles)
winDialog("ok","Please select the transects shapefile...")
path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)


mydata_trans <- data.frame(read.dbf(paste(substr(path1,1,max(del.ext)-4),".dbf",sep="")))

colnames(mydata_trans) <- gsub("ID", "Id", colnames(mydata_trans))

mydata_trans$dbf.Id <- seq(1,length(mydata_trans$dbf.Id),by=1)



winDialog("ok","Please select the GIS_stats_table_short.csv file...")

getdata <- choose.files(default = "*.csv",multi = FALSE)

mydata <- read.table(getdata, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(getdata)
setwd(dir.path)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_gisdata", showWarnings=FALSE)
setwd("AMBUR_gisdata")

dir.create(paste(time.stamp2," ","gisdata",sep=""))
setwd(paste(time.stamp2," ","gisdata",sep=""))

#set the variables
basecastdirx <- ifelse(min(mydata[,"Base_Off"] == 1), -1, 1)
tranlength <- pwidth


#### start code
t.startx <- mydata_trans$dbf.StartX
t.starty <- mydata_trans$dbf.StartY
t.azimuth <- mydata_trans$dbf.Azimuth
t.length <- tranlength * basecastdirx

t.endx <- sin((t.azimuth * pi/180)) * t.length + t.startx
t.endy <- cos((t.azimuth * pi/180)) * t.length + t.starty


Cx <- t.startx
Cy <- t.starty

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

t.startx2 <- c(t.startx[1],(sin((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.startx)[-length(t.startx)])
t.starty2 <- c(t.starty[1],(cos((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.starty)[-length(t.starty)])

t.startx2a <- c((sin((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.startx)[-length(t.startx)],t.startx[length(t.startx)])
t.starty2a <- c((cos((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.starty)[-length(t.starty)],t.starty[length(t.starty)])



Cx <- t.endx
Cy <- t.endy

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

t.endx2 <- c(t.endx[1],(sin((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.endx)[-length(t.endx)])
t.endy2 <- c(t.endy[1],(cos((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.endy)[-length(t.endy)])

t.endx2a <- c((sin((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.endx)[-length(t.endx)],t.endx[length(t.endx)])
t.endy2a <- c((cos((Segment.Azimuth * pi/180)) * (Segment.Length/2) + t.endy)[-length(t.endy)],t.endy[length(t.endy)])




#build points table for polygons

polypt1x <-  t.startx2
polypt1y <-  t.starty2

polypt2x <-  t.startx2a
polypt2y <-  t.starty2a

polypt3x <-  t.endx2a
polypt3y <-  t.endy2a

polypt4x <-  t.endx2
polypt4y <-  t.endy2

testxx <- c(polypt1x,polypt2x,polypt3x,polypt4x)
testyy <- c(polypt1y,polypt2y,polypt3y,polypt4y)
testid <- rep(mydata_trans$dbf.Transect,4)
testid2 <- c(rep(1,length(mydata_trans$dbf.Transect)),rep(2,length(mydata_trans$dbf.Transect)),rep(3,length(mydata_trans$dbf.Transect)),rep(4,length(mydata_trans$dbf.Transect)))



polycoords <- data.frame(testid,testid2,testxx,testyy)[order(testid,testid2), ]





#build gis table
gistable <- merge(mydata,mydata_trans,by.x = "Transect", by.y = "dbf.Transect")
colnames(gistable) <- gsub("dbf.", "t", colnames(gistable))

id.field <- gistable$Transect
dd <- data.frame(Id=polycoords$testid,X=polycoords$testxx,Y=polycoords$testyy)
ddTable <- data.frame(gistable)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Transect", 5)
write.shapefile(ddShapefile, paste("polyband",sep=""), arcgis=T)




detach("package:shapefiles")

#test plot
plot(c(t.startx,t.endx),c(t.starty,t.endy), type="p",asp=1,col="white",main="Perpendicular Transects",xlab="X",ylab="Y")
lines(mydata[,"Start_X"],mydata[,"Start_Y"], type="l",col="gray")
segments(t.startx,t.starty,t.endx,t.endy,col="blue")
points(t.startx,t.starty,col="green")
points(mydata[,"Start_X"],mydata[,"Start_Y"],col="red")
points(polypt1x,polypt1y,col="purple")
points(polypt2x,polypt2y,col="orange")
points(polypt3x,polypt3y,col="purple")
points(polypt4x,polypt4y,col="orange")
polygon(c(polypt1x,polypt2x,polypt3x,polypt4x), c(polypt1y,polypt2y,polypt3y,polypt4y), col="red")



#tidy up and remove all objects
detach(mydata)
rm(list = ls())

}