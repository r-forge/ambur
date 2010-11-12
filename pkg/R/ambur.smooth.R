ambur.smooth <-
function(alpha=0.7, degree=2, sampledist=5) {


# Establish the inputs

outer.alpha.rec <- alpha
outer.deg.rec <- degree
outersample <- sampledist





require(shapefiles)
require(locfit)
require(spatstat)


winDialog("ok","Please select a baseline to smooth...")

path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)

mydata_outer <- convert.to.simple(my.shapefile)

path2 <- dirname(path1)
setwd(path2)



time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_smoothing", showWarnings=FALSE)
setwd("AMBUR_smoothing")

dir.create(paste(time.stamp2," ","smooth",sep=""))
setwd(paste(time.stamp2," ","smooth",sep="")) 

baseline.dbf <- data.frame(read.dbf(paste(substr(path1,1,max(del.ext)-4),".dbf",sep="")))

colnames(baseline.dbf) <- gsub("ID", "Id", colnames(baseline.dbf))

baseline.dbf$dbf.Id <- seq(1,length(baseline.dbf$dbf.Id),by=1)


#build function for points along a line (modified to add the start and end points)

pts.along <- function(x,y,pspace) {

#VARIABLES: x = x coord, y = y coord, pspace = spacing between points

Cx <- x
Cy <- y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))

Cumltiv.Sum <- ave(Segment.Length, FUN=cumsum)

pts.int <- pspace  ##defines the point spacing along the line

pts.bin <- seq(from=pts.int,to=max(Cumltiv.Sum),by=pts.int)

Cumltiv.Sum2 <- c(0,Cumltiv.Sum[-max(length(Cumltiv.Sum))])


pts.bin.up <- length(pts.bin)
pts.bin.upx <- length(pts.bin)
pts.bin.upy <- length(pts.bin)
pts.bin.upx2 <- length(pts.bin)
pts.bin.upy2 <- length(pts.bin)
pts.bin.diff <- length(pts.bin)

for (i in 1:length(pts.bin)) {

pts.bin.up[i] <-  which.max(Cumltiv.Sum2[Cumltiv.Sum2  <=  pts.bin[i]])
pts.bin.upx[i] <- x[pts.bin.up[i]]
pts.bin.upy2[i] <- y[pts.bin.up[i]+1]
pts.bin.upx2[i] <- x[pts.bin.up[i]+1]
pts.bin.upy[i] <- y[pts.bin.up[i]]
pts.bin.diff[i] <- pts.bin[i] - Cumltiv.Sum2[pts.bin.up[i]]



}

Dx2 <- pts.bin.upx2 - pts.bin.upx
Dy2 <- pts.bin.upy2 - pts.bin.upy


pts.bin.az  <- ifelse(Dx2 >= 0, 90 -(180/pi) * atan(Dy2/Dx2),270 -(180/pi) * atan(Dy2/Dx2))

t.azimuth <- c(pts.bin.az[1],pts.bin.az,pts.bin.az[max(length(pts.bin.az))])

t.startx <- c(Cx[1],sin((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upx,Cx[max(length(Cx))])
t.starty <- c(Cy[1],cos((pts.bin.az * pi/180)) * pts.bin.diff + pts.bin.upy,Cy[max(length(Cy))])


cbind(t.startx,t.starty,t.azimuth)



}



  #build the locfit smoothing function
smlocfit.line <- function(x,y,n=20,alpha.z =0.7,kern.z="tricube",deg.z = 2,span.z = 2/3, fam.z = "gaussian") {
#VARIABLES: x = x coord, y = y coord, z = degrees of freedom, n  = densify number of points

Ox <- x
Oy <- y
n <- n #densification threshold
span.z <- span.z
deg.z <- deg.z
fam.z <- fam.z
alpha.z <- alpha.z
kern.z <- kern.z

   xy <- cbind(Ox,Oy)

   nP <- length(xy[,1])
   if((nP < 3) || (is.na(any(as.logical(xy[,1])))))
       return(print("need more than 3 points"))
       #return(list(x=numeric(0), y=numeric(0))) #old return call

  ## else :

   nP <- length(xy[,1])
   z <- n*(nP-1)

   i <- 1:nP
   Sx <- function(x) {fitted.values(locfit.raw(i, xy[,1],alpha=alpha.z,deg=deg.z,kern=kern.z,kt="sph"))}
   Sy <- function(x) {fitted.values(locfit.raw(i, xy[,2],alpha=alpha.z,deg=deg.z,kern=kern.z,kt="sph"))}
   ti <- seq(1, nP, length = z)
   opspl <- cbind(x = Sx(ti), y = Sy(ti))
   opspl
   }



##############################################################################################################
#set up the variables for the analyses

nbaselines <- sort(unique(c(mydata_outer$Id)))


finaltable2 <- matrix(ncol=3)
colnames(finaltable2) <- c("smx","smy","baseid")



for (b in 1:length(nbaselines)) {


Bx <- mydata_outer$X[mydata_outer$Id == nbaselines[b]]
By <- mydata_outer$Y[mydata_outer$Id == nbaselines[b]]

##############################################################################################################



Cx <- Bx
Cy <- By

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))




outer.basepts4sm <- pts.along(Bx, By,outersample) #to even out the points for smoothing



outer.baseptsSmooth <- smlocfit.line(outer.basepts4sm[,1], outer.basepts4sm[,2],alpha=outer.alpha.rec,kern="epan",n=1,deg=outer.deg.rec)




Cxo <- outer.baseptsSmooth[,1]
Cyo <- outer.baseptsSmooth[,2]



finaltable1 <- cbind(Cxo,Cyo,b)
colnames(finaltable1) <- c("smx","smy","baseid")

finaltable2 <- rbind(finaltable2,finaltable1)


}

finaltable3 <- finaltable2[-1,]



gistable <- data.frame(baseline.dbf)
colnames(gistable) <- gsub("dbf.", "", colnames(gistable))


gistable2 <- data.frame(cbind(gistable,ALPHA=outer.alpha.rec,DEGREE=outer.deg.rec))

#write shapefiles of the transects
library(shapefiles)


id.field <- finaltable3[,3]
dd <- data.frame(Id=c(id.field),X=c(finaltable3[,1]),Y=c(finaltable3[,2]))
ddTable <- data.frame(gistable2)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"baseline_smooth",sep=""), arcgis=T)


plot(mydata_outer$X,mydata_outer$Y,asp=1,col="gray",xlab="X",ylab="Y",type="p",cex=0.25,main="AMBUR-Smooth")
points(finaltable3,cex=0.25,col="blue")


detach("package:shapefiles")



}

