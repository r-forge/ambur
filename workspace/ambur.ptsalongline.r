ambur.ptsalongline <-
function(ptspace=50,offsetdist=-25) {

###enter a negative offsetdist number to offset left of the polyline,  positive for right of line

require(tcltk)
require(rgdal)
require(rgeos)


#ptspace <- 50 # for testing
#offsetdist <- -10 # for testing


tkmessageBox(message = "Please select the polyline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
shapedata <- as(shapedata, "SpatialLinesDataFrame")
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)
time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


#dir.create("AMBUR_points_along", showWarnings=FALSE)
#setwd("AMBUR_points_along")

#dir.create(paste(time.stamp2," ","points_along",sep=""))
#setwd(paste(time.stamp2," ","points_along",sep=""))

############################## build functions
#build function for points along a line (modified to add the start and end points)
###################################################################################
pts.along <- function(x,y,pspace,baseshapeID) {

#VARIABLES: x = x coord, y = y coord, pspace = spacing between points


baselineID <- unique(baseshapeID)[1]

Cx <- x
Cy <- y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

Segment.Azimuth[length(Segment.Azimuth)] <- Segment.Azimuth[length(Segment.Azimuth)-1]

#changed to cumsum function 2012
Cumltiv.Sum <- cumsum(Segment.Length)

pts.int <- pspace  ##defines the point spacing along the line

# added code to catch point spacing that is greater than the length of the polyline
pts.int <- ifelse(max(Cumltiv.Sum) > pspace, pspace, max(Cumltiv.Sum)/3)



pts.bin <- seq(from=pts.int,to=max(Cumltiv.Sum),by=pts.int)

Cumltiv.Sum2 <- c(0,Cumltiv.Sum[-max(length(Cumltiv.Sum))])


pts.bin.up <- length(pts.bin)
pts.bin.upx <- length(pts.bin)
pts.bin.upy <- length(pts.bin)
pts.bin.upx2 <- length(pts.bin)
pts.bin.upy2 <- length(pts.bin)
pts.bin.diff <- length(pts.bin)

pts.bin.az <- length(pts.bin)
Dx.bin  <- length(pts.bin)
Dy.bin <- length(pts.bin)

for (i in 1:length(pts.bin)) {

pts.bin.up[i] <-  which.max(Cumltiv.Sum2[Cumltiv.Sum2  <=  pts.bin[i]])
pts.bin.upx[i] <- x[pts.bin.up[i]]
pts.bin.upy[i] <- y[pts.bin.up[i]]

pts.bin.diff[i] <- pts.bin[i] - Cumltiv.Sum2[pts.bin.up[i]]

pts.bin.az[i] <- Segment.Azimuth[pts.bin.up[i]]
Dx.bin[i]  <- Dx[pts.bin.up[i]]
Dy.bin[i] <- Dy[pts.bin.up[i]]

}


t.azimuth <- pts.bin.az

t.startx <- c(sin((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upx)
t.starty <- c(cos((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upy)


baselineid <- rep(baselineID, length(t.azimuth))

cbind(baselineid, t.startx,t.starty,t.azimuth)



}


 ###### break down outer baseline into simply points with IDs
crdl0 <- coordinates(shapedata)
crd.1 <- sapply(crdl0, function(x) do.call("rbind", x),simplify = FALSE)
crd.2    <- do.call("rbind", crd.1)
crd.3 <- as.numeric(sapply(crd.1, function(x) max(row(x)) ,simplify = TRUE))
crd.len.test <- as.numeric(length(crd.3))
if(crd.len.test <= 1) crd.rep <-  1 else crd.rep <- seq(1, length(crd.3),1)
basepointIDs <- rep(crd.rep,crd.3)
baseshapeIDs <- basepointIDs - 1
sortshapeIDs <- seq(1,length(basepointIDs),1)
basex <- crd.2[,1]
basey <- crd.2[,2]

outerbase.tab <- data.frame(sortshapeIDs,baseshapeIDs,basepointIDs,basex,basey)
colnames(outerbase.tab) <- c("sortshapeID","shapeID","baseID","baseX", "baseY")



blah <- data.frame(sapply(1:4, function(x) data.frame(x)))

colnames(blah) <- c("BaseID","StartX","StartY","Azimuth")

 Baseline.Factor <- unique(outerbase.tab$baseID)

for (i in Baseline.Factor) {

trandata <- pts.along(outerbase.tab$baseX[outerbase.tab$baseID == Baseline.Factor[i]],outerbase.tab$baseY[outerbase.tab$baseID == Baseline.Factor[i]],ptspace,outerbase.tab$baseID[outerbase.tab$baseID == Baseline.Factor[i]])

colnames(trandata) <- c("BaseID","StartX","StartY","Azimuth")   ### added to fix mismatch names



blah <- rbind(blah, trandata)

}


pts.indv <- blah[-1,]


#########################################



sortID <- seq(1,length(pts.indv[,1]),1)

inter.data <- data.frame(pts.indv, sortID)
tran.data <- attrtable



tran.data$Id <- as.numeric(row.names(tran.data)) +1




tet <- merge(inter.data,tran.data , by.x = "BaseID", by.y = "Id", sort=FALSE)
tet2 <- tet

tet3 <- tet2[ order(tet2[,"sortID"]) , ]

tet3$Id <- tet2[,"sortID"]

tet3$TranSpace <- ptspace

####offset the points perpendicular

fsamp <- 90

aztable <- matrix(data = fsamp, nrow = length(tet3$StartX), ncol = length(fsamp), byrow = TRUE,dimnames = NULL)

tstarttable <- matrix(data = tet3$StartX, nrow = length(tet3$StartX), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

tendtable  <- matrix(data = tet3$StartY, nrow = length(tet3$StartY), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

t.startx <- tstarttable
t.starty <- tendtable
t.azimuth <- ifelse(tet3$Azimuth + aztable >= 360, tet3$Azimuth + aztable - 360, tet3$Azimuth + aztable)

t.endx <- sin((t.azimuth * pi/180)) * offsetdist + t.startx
t.endy <- cos((t.azimuth * pi/180)) * offsetdist + t.starty

tet3$StartX <- as.vector(t.endx)
tet3$StartY <- as.vector(t.endy)

tet3$CastDir <- offsetdist / abs(offsetdist)

#########################get curvature?
###established azimuth filter functions

winsize <- 10 #catch curves using window of 10

move.avg <- function(x,n=winsize){filter(x,rep(1/n,n), sides=2)}
move.sd <- function(x,n=winsize){filter(x,rep(1/n,n), sides=2)}


filter.azimuths <- function(az.data){
cos.az <- cos(az.data * pi/180)
sin.az <- sin(az.data * pi/180)
avg.cosx <- cos.az
avg.siny <- sin.az
filter.az.rads <- atan2(avg.siny, avg.cosx)

test.pos.neg <- filter.az.rads - c(filter.az.rads[-1],0)
avg.rad <- move.avg(filter.az.rads)
stdv.rad <- move.sd(filter.az.rads)

#filter.az.deg <- ifelse( (filter.az.rads *(180/pi)) < 0, (filter.az.rads *(180/pi)) + 360, (filter.az.rads *(180/pi)))
#filter.az.final <- ifelse(is.na(filter.az.deg) == TRUE, trandata$Azimuth, filter.az.deg)
return(cbind(filter.az.rads,test.pos.neg,avg.rad,test.pos.neg/abs(test.pos.neg)))}

curve.test <- filter.azimuths(tet3$Azimuth)


#####get sinuosity

dist.from.origin <-  ((tet3$StartX - tet3$StartX[1])^2 +  (tet3$StartY - tet3$StartY[1])^2)^(1/2)
dist.from.end <-  ((tet3$StartX - tet3$StartX[length(tet3$StartX)])^2 +  (tet3$StartY - tet3$StartY[length(tet3$StartY)])^2)^(1/2)

tet3$CmlDist <- cumsum(tet3$TranSpace)
tet3$Sinuos  <-   dist.from.origin/tet3$CmlDist
tet3$SinuosRev  <-   dist.from.end/rev(tet3$CmlDist)
tet3$SinuosDiff  <-  abs(tet3$Sinuos - tet3$SinuosRev)
tet3$Meander<-  tet3$Sinuos * 0
tet3$Meander[tet3$SinuosRev <0.93]  <-  1

plot(tet3$SinuosRev,type="l")
points(tet3$Sinuos[tet3$Sinuos <0.93], col="red")

####################################





pts.output <- SpatialPointsDataFrame(cbind(x=tet3$StartX,y=tet3$StartY),tet3)

###############################
#create opposite side
fsamp <- 90

aztable <- matrix(data = fsamp, nrow = length(tet3$StartX), ncol = length(fsamp), byrow = TRUE,dimnames = NULL)

tstarttable <- matrix(data = tet3$StartX, nrow = length(tet3$StartX), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

tendtable  <- matrix(data = tet3$StartY, nrow = length(tet3$StartY), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

t.startx <- tstarttable
t.starty <- tendtable
t.azimuth <- ifelse(tet3$Azimuth + aztable >= 360, tet3$Azimuth + aztable - 360, tet3$Azimuth + aztable)

t.endx <- sin((t.azimuth * pi/180)) * (offsetdist *-2) + t.startx
t.endy <- cos((t.azimuth * pi/180)) * (offsetdist *-2) + t.starty


Cx <- tet3$StartX
Cy <- tet3$StartY

Cx2 <- t.endx
Cy2 <- t.endy

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy



Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

Segment.Azimuth[length(Segment.Azimuth)] <- Segment.Azimuth[length(Segment.Azimuth)-1]

tet3b <- tet3

tet3b$StartX <- as.vector(t.endx)
tet3b$StartY <- as.vector(t.endy)

tet3b$CastDir <- offsetdist*-1 / abs(offsetdist)


pts.output2 <- SpatialPointsDataFrame(cbind(x=tet3b$StartX,y=tet3b$StartY),tet3b)


#final.output <- rbind(pts.output,pts.output2)


##########################################

 locname1 <- tet3b$Location[1]
  locname <- gsub(" ", "_", locname1)


 # Note that readOGR method reads the .prj file when it exists
  outputname <- paste(locname,"_basepts_left",sep="")

   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", outputname, driver="ESRI Shapefile")




 locname1 <- tet3b$Location[1]
  locname <- gsub(" ", "_", locname1)



 # Note that readOGR method reads the .prj file when it exists
  outputname <- paste(locname,"_basepts_right",sep="")

   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output2) <- projectionString

writeOGR(pts.output2, ".", outputname, driver="ESRI Shapefile")

plot(pts.output2)



}