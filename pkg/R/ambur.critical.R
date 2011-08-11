ambur.critical <-
function(ncritpts=50,sampledist=5) {

require(tcltk)

sample.distance <- sampledist
n.points <- ncritpts

#open baseline file

tkmessageBox(message = "Please select a baseline shapefile...")

library(shapefiles)

shape.path <- tk_choose.files(default = "*.shp",multi = FALSE)

dir.path <- dirname(shape.path )
setwd(dir.path)

del.ext <- nchar(shape.path)

shape.path2 <- paste(substr(shape.path,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(shape.path2)

simpleShpFormat <- convert.to.simple(my.shapefile)


time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_critical", showWarnings=FALSE)
setwd("AMBUR_critical")

dir.create(paste(time.stamp2," ","critical_points",sep=""))
setwd(paste(time.stamp2," ","critical_points",sep=""))




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


outer.basepts4sm <- pts.along(simpleShpFormat$X, simpleShpFormat$Y,sample.distance) #to even out the points for smoothing


#find a linear spline's critical points
Ox <- outer.basepts4sm[,1]
Oy <- outer.basepts4sm[,2]

   n <- 1 #densification threshold
   xy <- cbind(Ox,Oy)

   nP <- length(xy[,1])
   if((nP < 3) || (is.na(any(as.logical(xy[,1])))))
       return(list(x=numeric(0), y=numeric(0)))

  ## else :

   nP <- length(xy[,1])
   z <- n*(nP-1)

   i <- 1:nP
   Sx <- function(x) {approx(i, xy[,1],method="linear",n=n.points)}
   Sy <- function(x) {approx(i, xy[,2],method="linear",n=n.points)}
   ti <- seq(1, nP, length = z)
   opspl <- cbind(x = Sx(ti)$y, y = Sy(ti)$y)


   #plot the results and write them to a file
pdf(paste("plot_critical_points.pdf",sep=""),width = 6.5, height = 6.5, bg="white", paper= "letter" )
 plot(Ox,Oy, col="light gray",type="l",asp=1,main="Linear Spline Critical Points",xlab="X",ylab="Y")
 lines(opspl, col="green")
 points(opspl, col="orange")
#turn the device off
dev.off()

 plot(Ox,Oy, col="light gray",type="l",asp=1,main="Linear Spline Critical Points",xlab="X",ylab="Y")
 lines(opspl, col="green")
 points(opspl, col="orange")



id.field <- seq(1,length(opspl[,1]),by=1)
dd <- data.frame(Id=c(id.field),X=c(opspl[,1]),Y=c(opspl[,2]))
ddTable <- data.frame(Id=id.field,PointX=opspl[,1],PointY=opspl[,2],Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
write.shapefile(ddShapefile, paste("critical_pts",sep=""), arcgis=T)

id.field <- 1
dd <- data.frame(Id=c(id.field),X=c(opspl[,1]),Y=c(opspl[,2]))
ddTable <- data.frame(Id=id.field,Source=shape.path,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("critical_pts_line",sep=""), arcgis=T)

#convert the original shapefile to points
#id.field <- seq(1,length(simpleShpFormat$X),by=1)
#dd <- data.frame(Id=c(id.field),X=c(simpleShpFormat$X),Y=c(simpleShpFormat$Y))
#ddTable <- data.frame(Id=id.field,PointX=simpleShpFormat$X,PointY=simpleShpFormat$Y,Creator="R - AMBUR")
#ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
#write.shapefile(ddShapefile, paste("original_pts",sep=""), arcgis=T)

}

