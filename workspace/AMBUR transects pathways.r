#must have all fields and a CastDir field with a "1" to cast to the right flow & "-1" to cast to the left

#Start the function 
ambur.path <- function(userinput1=50, userinput2=500, userinput3=5, userinput4=5) {


# Establish the inputs
transpace <- userinput1
tranlength <- userinput2
innersample <- userinput3
outersample <- userinput4

#transpace <- 50         
#tranlength <- 500
#innersample <- 5
#outersample <- 5




require(shapefiles)
require(locfit)
require(spatstat)


winDialog("ok","Please select the inner baseline...")

path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)

mydata_inner <- convert.to.simple(my.shapefile)

path2 <- dirname(path1)
setwd(path2)

winDialog("ok","Please select the outer baseline...")
path3 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path3)

path3a <- paste(substr(path3,1,max(del.ext)-4),".shp",sep="")

my.shapefile2 <- read.shp(path3a)

mydata_outer <- convert.to.simple(my.shapefile2)

baseline.dbf <- data.frame(read.dbf(paste(substr(path3,1,max(del.ext)-4),".dbf",sep="")))

colnames(baseline.dbf) <- gsub("ID", "Id", colnames(baseline.dbf))

baseline.dbf$dbf.Id <- seq(1,length(baseline.dbf$dbf.Id),by=1)


time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_transects", showWarnings=FALSE)
setwd("AMBUR_transects")

dir.create(paste(time.stamp2," ","transects",sep=""))
setwd(paste(time.stamp2," ","transects",sep="")) 


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



 
##############################################################################################################
#set up the variables for the analyses


for (m in 1:length(mydata_outer$Id)) {

test1 <- which(mydata_outer$Id[m] == baseline.dbf$dbf.Id)

mydata_outer$Id[m] <- baseline.dbf$dbf.BaseOrder[test1]


}


#setup table to hold the results
finaltable2 <- matrix(ncol=16)
colnames(finaltable2) <- c("startx","starty","perpx","perpy","perpaz","trimx","trimy","trimdist","nearx","neary","nearaz","neardist","basemaxbnum","baseorder","baseoffshore","basecastdir")

texttable2 <- matrix(ncol=3)
colnames(texttable2) <- c("startx","baseloc1","baseloc2")


nbaselines <- sort(unique(c(mydata_outer$Id)))

transtart <- 0

Cx <- 0
Cy <- 0
Cid <- 0

Cx2 <- 0
Cy2 <- 0
 
 pts.id <- 0 
  int.ptsX <- 0
    int.ptsY <- 0
      trans.id <- 0
      


for (b in 1:length(nbaselines)) {


Ax <- mydata_inner$X
Ay <- mydata_inner$Y

Bx <- mydata_outer$X[mydata_outer$Id == nbaselines[b]]
By <- mydata_outer$Y[mydata_outer$Id == nbaselines[b]]

baseloc1x <- as.character(baseline.dbf$dbf.Location[baseline.dbf$dbf.BaseOrder == nbaselines[b]])
basemaxbnumx <-  baseline.dbf$dbf.MaxBNum[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseorderx <-  baseline.dbf$dbf.BaseOrder[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseoffshorex <- baseline.dbf$dbf.OFFshore[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
basecastdirx <-   baseline.dbf$dbf.CastDir[baseline.dbf$dbf.BaseOrder == nbaselines[b]]
baseloc2x <- as.character(baseline.dbf$dbf.BASE_LOC[baseline.dbf$dbf.BaseOrder == nbaselines[b]])

##############################################################################################################



Cx <- Bx
Cy <- By

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Segment.Length <-  ((Cx2- Cx)^2 +  (Cy2 - Cy)^2)^(1/2)

Dx <- Cx2 - Cx
Dy <- Cy2 - Cy

Segment.Azimuth <- ifelse(Dx >= 0, 90 -(180/pi) * atan(Dy/Dx),270 -(180/pi) * atan(Dy/Dx))



####cast transects

outer.basepts <- pts.along(Bx,By,transpace)  #for transects



#cast near transects
Cx <- mydata_inner$X
Cy <- mydata_inner$Y

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Cx3 <- c(Cx[1],Cx[-length(Cx)])
Cy3 <- c(Cy[1],Cy[-length(Cy)])

Cxo <- outer.basepts[,1]
Cyo <- outer.basepts[,2]

Cx2o <- c(Cxo[-1],Cxo[length(Cxo)])
Cy2o <- c(Cyo[-1],Cyo[length(Cyo)])

Cx3o <- c(Cxo[1],Cxo[-length(Cxo)])
Cy3o <- c(Cyo[1],Cyo[-length(Cyo)])


TY.w <- owin()
TY.w <- owin(c(min(Cx),max(Cx)), c(min(Cy),max(Cy)))
TY <- psp(Cx,Cy,Cx2,Cy2,window=TY.w)

TX.w <- owin()
TX.w <- owin(c(min(outer.basepts[,1]),max(outer.basepts[,1])), c(min(outer.basepts[,2]),max(outer.basepts[,2])))
TX <- ppp(outer.basepts[,1],outer.basepts[,2],window=TX.w)

v <- project2segment(TX,TY)
  Xproj <- v$Xproj

Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

#transect filter using average window of azimuths
tfilter.tab <- cbind(TX$x,TX$y,Xproj$x,Xproj$y,v$d,v$mapXY,v$tp)
colnames(tfilter.tab) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos")


Dx2_in <- tfilter.tab[,3] - tfilter.tab[,1]
Dy2_in <- tfilter.tab[,4] - tfilter.tab[,2]



in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tfilter.tab[,3]- tfilter.tab[,1])^2 +  (tfilter.tab[,4] - tfilter.tab[,2])^2)^(1/2)


Inx2 <- c(tfilter.tab[-1,3],tfilter.tab[length(tfilter.tab[,3]),3])
Iny2 <- c(tfilter.tab[-1,4],tfilter.tab[length(tfilter.tab[,4]),4])

Inx3 <- c(tfilter.tab[1,3],tfilter.tab[-length(tfilter.tab[,3]),3])
Iny3 <- c(tfilter.tab[1,4],tfilter.tab[-length(tfilter.tab[,4]),4])

end.tspace <-  ((Inx2- tfilter.tab[,3])^2 +  (Iny2 - tfilter.tab[,4])^2)^(1/2)
end.tspace2 <-  ((Inx3- tfilter.tab[,3])^2 +  (Iny3 - tfilter.tab[,4])^2)^(1/2)



#cast perpendicular transects along outer baseline
t.startx <- outer.basepts[,1]
t.starty <- outer.basepts[,2]
t.azimuth <- ifelse(outer.basepts[,3] + 90 >= 360, outer.basepts[,3] + 90 - 360, outer.basepts[,3] + 90) 
t.length <- tranlength * basecastdirx

t.endx <- sin((t.azimuth * pi/180)) * t.length + t.startx
t.endy <- cos((t.azimuth * pi/180)) * t.length + t.starty


####################################test to trim transects

test.wx <- c(t.startx,t.endx,Cx,Cx2,mydata_outer$X)
test.wy <- c(t.starty,t.endy,Cy,Cy2,mydata_outer$Y)

Test.w <- owin()
Test.w <- owin(c(min(test.wx-100),max(test.wx+100)), c(min(test.wy-100),max(test.wy+100)))

TY.w <- owin()
TY.w <- owin(c(min(Cx),max(Cx)), c(min(Cy),max(Cy)))
TY <- psp(Cx,Cy,Cx2,Cy2,window=Test.w)

trim.x <- numeric(length(tfilter.tab[,1]))
trim.y  <- numeric(length(tfilter.tab[,1]))
trim.length  <- numeric(length(tfilter.tab[,1]))

for (i in 1:length(tfilter.tab[,1])) {
 
b.x <- psp(t.startx[i], t.starty[i], t.endx[i], t.endy[i], window=Test.w)

 inner.trim <- crossing.psp(b.x,TY)

pts.dists <- ((inner.trim$x - tfilter.tab[,1][i])^2 +  (inner.trim$y - tfilter.tab[,2][i])^2)^(1/2)

 trim.x[i] <- ifelse(inner.trim$n == 0, t.endx[i],inner.trim$x[which.min(pts.dists)] )
 trim.y[i] <- ifelse(inner.trim$n == 0, t.endy[i],inner.trim$y[which.min(pts.dists)] )
 trim.length[i] <- ifelse(inner.trim$n == 0,t.length,min(pts.dists))
 } 

filter.length6 <- trim.length
filter.x6 <- trim.x
filter.y6 <- trim.y
#####################################################################

#construct master data table

finaltable1 <- cbind(t.startx,t.starty,t.endx,t.endy,t.azimuth,trim.x,trim.y,trim.length,Xproj$x,Xproj$y,in.az,in.length,basemaxbnumx,baseorderx,baseoffshorex,basecastdirx)
colnames(finaltable1) <- c("startx","starty","perpx","perpy","perpaz","trimx","trimy","trimdist","nearx","neary","nearaz","neardist","basemaxbnum","baseorder","baseoffshore","basecastdir")

finaltable2 <- rbind(finaltable2,finaltable1)


texttable1 <- cbind(t.startx,baseloc1x,baseloc2x)

colnames(texttable1) <- c("ttx","baseloc1","baseloc2")

texttable2 <- rbind(texttable2,texttable1)

transtart <- max(length(t.startx))

}
finaltable3 <- finaltable2[-1,]

texttable3 <- texttable2[-1,]

attach(data.frame(finaltable3))



#plots for fun
par(mfrow=(c(3,1)))
par(pty= "m")
plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Perpendicular Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,perpx,perpy,col="blue")

plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Perpendicular Trimmed Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,trimx,trimy,col="blue")

plot(c(startx,perpx),c(starty,perpy), type="l",asp=1,col="white",main="Near Transects",xlab="X",ylab="Y")
lines(mydata_inner$X,mydata_inner$Y, type="l",col="gray")
segments(startx,starty,nearx,neary,col="blue")

return(finaltable3)

}


outer.init <- ambur.path(50,100,5,5)



 winDialog("ok","Please select a inner baselines shapefile...")
shape.path <- choose.files(default = "*.shp",multi = FALSE)

dir.path <- dirname(shape.path )
#setwd(dir.path)

del.ext <- nchar(shape.path)

shape.path2 <- paste(substr(shape.path,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(shape.path2)

simpleShpFormat <- convert.to.simple(my.shapefile)


mydata1 <- data.frame(read.dbf(paste(substr(shape.path,1,max(del.ext)-4),".dbf",sep="")))



mydata1$Id <- seq(1,length(unique(simpleShpFormat$Id)))




nbaselines <- mydata1$Id

Bx <- outer.init[,"nearx"]

By <- outer.init[,"neary"]

transpace <- 50

nearx <- 0

neary <- 0      

mydata1$dbf.level <- as.numeric(as.character(mydata1$dbf.level))
mydata1$seqn <- seq(1,length(mydata1$dbf.level),1)

repair.order <- numeric(length(simpleShpFormat$Id))

for (a in 1:length(simpleShpFormat$Id)) {

repair.order[a] <- mydata1$dbf.level[mydata1$seqn == simpleShpFormat$Id[a]] 
}


mod.lines <- cbind(simpleShpFormat$Id,simpleShpFormat$X,simpleShpFormat$Y,repair.order)[ order(repair.order) ,]
colnames(mod.lines) <- c("modid","modx","mody","level")

#mod.lines <- mod.lines[mod.lines[,"level"] > 1,]
#mod.adj <- length(mydata1$seqn) - (mod.lines[,"level"] -1) 

mod.lines <- mod.lines[mod.lines[,"level"] < length(mydata1$dbf.level),]
mod.adj <- length(mydata1$seqn) - (mod.lines[,"level"]) 


for (b in 1:(length(nbaselines)-1)) {



#cast near transects
Cx <- mod.lines[,"modx"][mod.adj == b]
Cy <- mod.lines[,"mody"][mod.adj == b]

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Cx3 <- c(Cx[1],Cx[-length(Cx)])
Cy3 <- c(Cy[1],Cy[-length(Cy)])

Cxo <- Bx
Cyo <- By

Cx2o <- c(Cxo[-1],Cxo[length(Cxo)])
Cy2o <- c(Cyo[-1],Cyo[length(Cyo)])

Cx3o <- c(Cxo[1],Cxo[-length(Cxo)])
Cy3o <- c(Cyo[1],Cyo[-length(Cyo)])


TY.w <- owin()
TY.w <- owin(c(min(Cx),max(Cx)), c(min(Cy),max(Cy)))
TY <- psp(Cx,Cy,Cx2,Cy2,window=TY.w)

TX.w <- owin()
TX.w <- owin(c(min(Bx),max(Bx)), c(min(By),max(By)))
TX <- ppp(Bx,By,window=TX.w)

v <- project2segment(TX,TY)
  Xproj <- v$Xproj

Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

nearx <- c(nearx,Xproj$x)
neary <- c(neary,Xproj$y)

Bx <- Xproj$x
By <- Xproj$y

}

trnx <- c(outer.init[,"startx"],nearx[-1])
trny <- c(outer.init[,"starty"],neary[-1])

plot(trnx,trny,asp=T)
segments(outer.init[,"startx"],outer.init[,"starty"],outer.init[,"nearx"],outer.init[,"neary"],asp=T)


transect.num <- rep(seq(1,length(outer.init[,"startx"]),1),length(mydata1$seqn))


#write shapefiles of the transects
library(shapefiles)

dd <- data.frame(Id=transect.num,X=trnx,Y=trny)
ddTable <- data.frame(Id=unique(transect.num),Transect=unique(transect.num))
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("pathway_transects",sep=""), arcgis=T)




detach("package:shapefiles")
detach(data.frame(finaltable3))