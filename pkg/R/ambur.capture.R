ambur.capture <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(shapefiles)
require(locfit)
require(spatstat)
require(tcltk)


winDialog("ok","Please select the shorelines shapefile...")

path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

shore.dbf <- data.frame(read.dbf(paste(substr(path1,1,max(del.ext)-4),".dbf",sep="")))

colnames(shore.dbf) <- gsub("ID", "Id", colnames(shore.dbf))

shore.dbf$dbf.Id <- seq(1,length(shore.dbf$dbf.Id),by=1)

my.shapefile <- read.shp(path1a)

mydata_inner <- convert.to.simple(my.shapefile)

path2 <- dirname(path1)
setwd(path2)

winDialog("ok","Please select the transects shapefile...")

path3 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path3)

path3a <- paste(substr(path3,1,max(del.ext)-4),".shp",sep="")

transect.dbf <- data.frame(read.dbf(paste(substr(path3,1,max(del.ext)-4),".dbf",sep="")))

my.shapefile2 <- read.shp(path3a)

mydata_outer <- convert.to.simple(my.shapefile2)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_capture", showWarnings=FALSE)
setwd("AMBUR_capture")

dir.create(paste(time.stamp2," ","shorepts",sep=""))
setwd(paste(time.stamp2," ","shorepts",sep="")) 

##status plot window
plot(mydata_outer$X,mydata_outer$Y,col="white", asp=1,xlab="X",ylab="Y")
segments(transect.dbf$dbf.StartX,transect.dbf$dbf.StartY,transect.dbf$dbf.EndX,transect.dbf$dbf.EndY,col="green")



####start analyses
nsegments <- sort(unique(c(mydata_inner$Id)))

pb <- tkProgressBar("AMBUR: progress bar", "Some information in %", 0, max(length(unique(nsegments))), 50)

Cx <- 0
Cy <- 0
Cid <- 0

Cx2 <- 0
Cy2 <- 0
 
 pts.id <- 0 
  int.ptsX <- 0
    int.ptsY <- 0
      trans.id <- 0
      
      
for (b in 1:length(nsegments)) {

Cx <- mydata_inner$X[mydata_inner$Id == b]
Cy <- mydata_inner$Y[mydata_inner$Id == b]
Cid <- mydata_inner$Id[mydata_inner$Id == b]

Cx2 <- c(mydata_inner$X[mydata_inner$Id == b][-1],mydata_inner$X[mydata_inner$Id == b][length(mydata_inner$X[mydata_inner$Id == b])])
Cy2 <- c(mydata_inner$Y[mydata_inner$Id == b][-1],mydata_inner$Y[mydata_inner$Id == b][length(mydata_inner$Y[mydata_inner$Id == b])])




Cx <- Cx
Cy <- Cy
Cid <- Cid

Cx2 <- Cx2
Cy2 <- Cy2



test.wx <- c(mydata_inner$X,mydata_outer$X)
test.wy <- c(mydata_inner$Y,mydata_outer$Y)


Cxo <- transect.dbf$dbf.StartX
Cyo <- transect.dbf$dbf.StartY
Cido <- transect.dbf$dbf.Transect

Cx2o <- transect.dbf$dbf.EndX
Cy2o <- transect.dbf$dbf.EndY
Cido2 <- transect.dbf$dbf.Transect





Test.w <- owin()
Test.w <- owin(c(min(test.wx-100000),max(test.wx+100000)), c(min(test.wy-100000),max(test.wy+100000)))

TY.w <- owin()
TY.w <- owin(c(min(test.wx-100000),max(test.wx+100000)), c(min(test.wy-100000),max(test.wy+100000)))
TY <- psp(Cx,Cy,Cx2,Cy2,mark=Cid,window=Test.w)




 all.ptsX <- 0
 all.ptsY <- 0
 pts.id2 <- 0
 trans.id2 <- 0

for (i in 1:length(Cxo)) {
 
b.x <- psp(Cxo[i], Cyo[i], Cx2o[i], Cy2o[i], window=Test.w)

 inner.trim <- crossing.psp(b.x,TY)

 pts.dists <- ((inner.trim$x - Cxo[i])^2 +  (inner.trim$y - Cyo[i])^2)^(1/2)



 all.ptsX <- c(all.ptsX,inner.trim$x)
 all.ptsY <- c(all.ptsY,inner.trim$y)
 pts.id2 <- c(pts.id2,numeric(length(inner.trim$x))+b)
 trans.id2 <- c(trans.id2,numeric(length(inner.trim$x))+ Cido[i])

 } 


int.ptsX2 <- all.ptsX[-1]
int.ptsY2 <- all.ptsY[-1]
pts.id2 <- pts.id2[-1]
trans.id2 <- trans.id2[-1]


pts.id <- c(pts.id,pts.id2)
int.ptsX <- c(int.ptsX,int.ptsX2)
int.ptsY <- c(int.ptsY,int.ptsY2)
trans.id <- c(trans.id,trans.id2)


#update plot
segments(Cx,Cy,Cx2,Cy2,col="gray")
points(int.ptsX[-1],int.ptsY[-1],col="blue")


#status update: add progress bar, estimate percent completion and map

Pcnt.Complete <-  round(((b)/ length(nsegments)) * 100, 0)

Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="")



    info <- sprintf("%d%% done", Pcnt.Complete)
    setTkProgressBar(pb, b, sprintf("AMBUR: Capture shoreline positions (%s)", info), info)



}

pts.id <- pts.id[-1]
int.ptsX <- int.ptsX[-1]
int.ptsY <- int.ptsY[-1]
trans.id <- trans.id[-1]

#####################################################################

#plot(int.ptsX,int.ptsY)
#segments(Cxo, Cyo, Cx2o, Cy2o,col="green")
#segments(Cx,Cy,Cx2,Cy2,col="gray")

int.tab <- cbind(trans.id,int.ptsX,int.ptsY,pts.id)
colnames(int.tab) <- c("Transect","Point_X","Point_Y","pointID")

tet <- merge(int.tab,transect.dbf, by.x = "Transect", by.y = "dbf.Transect")
tet2 <- merge(tet,shore.dbf, by.x = "pointID", by.y = "dbf.Id")
colnames(tet2) <- gsub("dbf.", "", colnames(tet2))
colnames(tet2) <- gsub("header.y", "header_y", colnames(tet2))

trandist <- (((tet2$Point_X - tet2$StartX)^2 +  (tet2$StartY - tet2$Point_Y)^2)^(1/2))


#write shapefiles of the transects
#will produce warning messages if CLASS_1,2, or 3 is empty: (In max(nchar(x[!is.na(x)], "b"))no non-missing arguments to max; returning -Inf ....the warning message doesn't hurt anything. 

library(shapefiles)
id.field <- seq(1,length(tet2[,1]),by=1)
dd <- data.frame(Id=c(id.field),X=c(tet2$Point_X),Y=c(tet2$Point_Y))
ddTable <- data.frame(Id=id.field,Transect = tet2$Transect,TranSpace = tet2$TranSpace,TranDist = tet2$TranDist,Location = tet2$Location,MaxBNum = tet2$MaxBNum,BaseOrder = tet2$BaseOrder,OFFshore = tet2$OFFshore,CastDir = tet2$CastDir,BASE_LOC = tet2$BASE_LOC,StartX = tet2$StartX,StartY = tet2$StartY,EndX = tet2$EndX,EndY = tet2$EndY,Azimuth = tet2$Azimuth,ACCURACY = tet2$ACCURACY,SHORE_LOC = tet2$SHORE_LOC,CLASS_1 = tet2$CLASS_1,CLASS_2 = tet2$CLASS_2,CLASS_3 = tet2$CLASS_3,DATE_= tet2$DATE_,POINT_X=tet2$Point_X,POINT_Y=tet2$Point_Y,Distance = trandist,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
write.shapefile(ddShapefile, paste("shore_pts",sep=""), arcgis=T)


detach("package:shapefiles")

}

