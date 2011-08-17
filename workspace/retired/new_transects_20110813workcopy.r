require(tcltk)
require(rgdal)
require(rgeos)


transpace <- 50
tranlength <- 2000
innersample <- 5
outersample <- 5
fsamp <- 90
#fsamp <- c(seq(2,178,by=5),178) # radiating transects

tkmessageBox(message = "Please select the outer baseline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the inner baseline shapefile...")
getdata2 <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create("AMBUR_transects", showWarnings=FALSE)
setwd("AMBUR_transects")

dir.create(paste(time.stamp2," ","transects",sep=""))
setwd(paste(time.stamp2," ","transects",sep=""))

############################## build functions
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


#######################################################################
###########build perpendicular transect function

perp.trans <- function(OutX,OutY,transpace,tranlength,fsamp,castdir) {

castdir <- unique(castdir)[1]

basecastdirx <- castdir

Bx <- OutX
By <- OutY


Cx <- Bx
Cy <- By

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))


outer.basepts <- pts.along(Bx,By,transpace)  #for transects



aztable <- matrix(data = fsamp, nrow = length(outer.basepts[,1]), ncol = length(fsamp), byrow = TRUE,dimnames = NULL)

tstarttable <- matrix(data = outer.basepts[,1], nrow = length(outer.basepts[,1]), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

tendtable  <- matrix(data = outer.basepts[,2], nrow = length(outer.basepts[,1]), ncol = length(fsamp), byrow = FALSE,dimnames = NULL)

t.startx <- tstarttable
t.starty <- tendtable
t.azimuth <- ifelse(outer.basepts[,3] + aztable >= 360, outer.basepts[,3] + aztable - 360, outer.basepts[,3] + aztable)
t.length <- tranlength * basecastdirx

t.endx <- sin((t.azimuth * pi/180)) * t.length + t.startx
t.endy <- cos((t.azimuth * pi/180)) * t.length + t.starty

transectid <- seq(1,length(outer.basepts[,1]),by=1)

 perp.table <- cbind(transectid,t.startx,t.starty,t.endx,t.endy,t.azimuth,t.length)
 colnames(perp.table) <- c("Transect","StartX","StartY","EndX","EndY","Azimuth","TranLength")

 return(perp.table)
 
 }

###### break down outer baseline into simply points with IDs
crdl0 <- coordinates(shapedata)
crd.1 <- sapply(crdl0, function(x) do.call("rbind", x),simplify = FALSE)
crd.2    <- do.call("rbind", crd.1)
crd.3 <- as.numeric(sapply(crd.1, function(x) max(row(x)) ,simplify = TRUE))
crd.len.test <- as.numeric(length(crd.3))
if(crd.len.test <= 1) crd.rep <-  1 else crd.rep <- c(1, length(crd.3))
basepointIDs <- rep(crd.rep,crd.3)
baseshapeIDs <- basepointIDs - 1
sortshapeIDs <- seq(1,length(basepointIDs),1)
basex <- crd.2[,1]
basey <- crd.2[,2]

outerbase.tab <- data.frame(sortshapeIDs,baseshapeIDs,basepointIDs,basex,basey)
colnames(outerbase.tab) <- c("sortshapeID","shapeID","baseID","baseX", "baseY")

attrtable$Id <- seq(0,max(length(attrtable[,1]))-1,1)
bbb <- data.frame(merge(outerbase.tab,attrtable,by.x= "shapeID" ,by.y = "Id", all.x = TRUE,sort=FALSE))


### cast transects for individual polylines based on unique IDs
Baseline.Factor <- factor(bbb$BaseOrder)

trans.indv <- sapply(levels(Baseline.Factor), function(x) data.frame(perp.trans(bbb$baseX[bbb$BaseOrder == x],bbb$baseY[bbb$BaseOrder == x],transpace=transpace,tranlength=tranlength,fsamp=fsamp,castdir=bbb$CastDir[bbb$BaseOrder == x])) ,simplify = FALSE)
trans.indv2 <- data.frame(do.call("rbind", trans.indv))

trans.indv2$Transect <- seq(1,length(trans.indv2$Transect),1)

trans.indv2 <- data.frame(do.call("rbind", trans.indv))

trans.all <- data.frame(sapply(levels(Baseline.Factor), function(x) data.frame(perp.trans(bbb$baseX,bbb$baseY,transpace=transpace,tranlength=tranlength,fsamp=fsamp,castdir=bbb$CastDir)) ,simplify = FALSE) )
colnames(trans.all) <- c("Transect","StartX","StartY","EndX","EndY","Azimuth","TranLength")


if(max(unique(bbb$baseID)) > 1)  test.tran <- trans.indv2 else test.tran <- trans.all


 ##### plot the results of transect casting   
    plot(c(test.tran$StartX,test.tran$EndX),c(test.tran$StartY,test.tran$EndY),col="white",asp=1,xlab="X",ylab="Y")
    segments(test.tran$StartX,test.tran$StartY,test.tran$EndX,test.tran$EndY)
    lines(bbb$baseX,bbb$baseY,col="red")
    
       plot(c(trans.all$StartX,trans.all$EndX),c(trans.all$StartY,trans.all$EndY),col="white",asp=1,xlab="X",ylab="Y")
    segments(trans.all$StartX,trans.all$StartY,trans.all$EndX,trans.all$EndY)
    lines(bbb$baseX,bbb$baseY,col="red")
    
       plot(c(trans.indv2$StartX,trans.indv2$EndX),c(trans.indv2$StartY,trans.indv2$EndY),col="white",asp=1,xlab="X",ylab="Y")
    segments(trans.indv2$StartX,trans.indv2$StartY,trans.indv2$EndX,trans.indv2$EndY)
    lines(bbb$baseX,bbb$baseY,col="red")