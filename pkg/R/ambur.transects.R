ambur.transects <-
function(userinput1=50, userinput2=500, userinput3=5, userinput4=5) {


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
require(tcltk)


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

###########build the status bar
pb <- tkProgressBar("AMBUR: progress bar", "Some information in %", 0, 100, 50)

    info <- sprintf("%d%% done", 0)
    setTkProgressBar(pb, 0, sprintf("AMBUR: Transect casting (%s)", info), info)

#####################



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

#######status checkpoint (5%)
    info <- sprintf("%d%% done", 5)
    setTkProgressBar(pb, 5, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################    




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

#######status checkpoint (10%)
    info <- sprintf("%d%% done", 10)
    setTkProgressBar(pb, 10, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################      


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
TY.w <- owin(c(min(Cx-100000),max(Cx+100000)), c(min(Cy-100000),max(Cy+100000)))
TY <- psp(Cx,Cy,Cx2,Cy2,window=TY.w)

TX.w <- owin()
TX.w <- owin(c(min(outer.basepts[,1]-100000),max(outer.basepts[,1]+100000)), c(min(outer.basepts[,2]-100000),max(outer.basepts[,2]+100000)))
TX <- ppp(outer.basepts[,1],outer.basepts[,2],window=TX.w)

################################### adjusted project2segment function to handle NA values
adj.project2segment <-       function (X, Y, action = "project", check = FALSE) 
{
    stopifnot(is.ppp(X))
    stopifnot(is.psp(Y))
    stopifnot(action %in% c("distance", "identify", "project"))
    if (Y$n == 0) 
        stop("Segment pattern Y contains 0 segments; projection undefined")
    if (X$n == 0) {
        nowt <- numeric(0)
        none <- integer(0)
        switch(action, identify = return(none), distance = return(list(dist = nowt, 
            which = none)), project = return(list(Xproj = X, 
            mapXY = none, d = nowt, tp = nowt)))
    }
    XX <- as.matrix(as.data.frame(unmark(X)))
    YY <- as.matrix(as.data.frame(unmark(Y)))
    d <- distppllmin(XX, YY)
    mapXY <- d$min.which
    if (action == "identify") 
        return(mapXY)
    else if (action == "distance") 
        return(data.frame(dist = d$min.d, which = mapXY))
    alldata <- as.data.frame(cbind(XX, YY[mapXY, , drop = FALSE]))
    colnames(alldata) <- c("x", "y", "x0", "y0", "x1", "y1")
    dx <- with(alldata, x1 - x0)
    dy <- with(alldata, y1 - y0)
    leng <- sqrt(dx^2 + dy^2)
    
    co <- dx/leng
    co[is.na(co)==TRUE] <- 0.000001  #added 6/27/2011 to remove NA values from crashing the function  
    
    si <- dy/leng
    si[is.na(si)==TRUE] <- 0.000001   #added 6/27/2011 to remove NA values from crashing the function  
    
    xv <- with(alldata, x - x0)
    yv <- with(alldata, y - y0)
    xpr <- xv * co + yv * si
    ypr <- -xv * si + yv * co
    left <- (xpr <= 0)
    right <- (xpr >= leng)
    xr <- with(alldata, ifelse(left, 0, ifelse(right, leng, xpr)))
    xproj <- with(alldata, x0 + xr * co)
    yproj <- with(alldata, y0 + xr * si)
    Xproj <- ppp(xproj, yproj, window = X$window, marks = X$marks, 
        check = check)
    tp <- xr/leng
    tp[!is.finite(tp)] <- 0
    return(list(Xproj = Xproj, mapXY = mapXY, d = d$min.d, tp = tp))
}


####################################end adj.project2segment function
v <- adj.project2segment(TX,TY)
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

#######status checkpoint (25%)
    info <- sprintf("%d%% done", 25)
    setTkProgressBar(pb, 25, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################

#cast perpendicular transects along outer baseline
t.startx <- outer.basepts[,1]
t.starty <- outer.basepts[,2]
t.azimuth <- ifelse(outer.basepts[,3] + 90 >= 360, outer.basepts[,3] + 90 - 360, outer.basepts[,3] + 90) 
t.length <- tranlength * basecastdirx

t.endx <- sin((t.azimuth * pi/180)) * t.length + t.startx
t.endy <- cos((t.azimuth * pi/180)) * t.length + t.starty


#######status checkpoint (50%)
    info <- sprintf("%d%% done", 50)
    setTkProgressBar(pb, 50, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################


####################################test to trim transects

test.wx <- c(t.startx,t.endx,Cx,Cx2,mydata_outer$X)
test.wy <- c(t.starty,t.endy,Cy,Cy2,mydata_outer$Y)

Test.w <- owin()
Test.w <- owin(c(min(test.wx-100000),max(test.wx+100000)), c(min(test.wy-100000),max(test.wy+100000)))

TY.w <- owin()
TY.w <- owin(c(min(Cx-100000),max(Cx+100000)), c(min(Cy-100000),max(Cy+100000)))
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

#######status checkpoint (75%)
    info <- sprintf("%d%% done", 75)
    setTkProgressBar(pb, 75, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################

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


#write shapefiles of the transects
library(shapefiles)
id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,perpx),Y=c(starty,perpy))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=tranlength,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=perpx,EndY=perpy,Azimuth=perpaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_perp",sep=""), arcgis=T)

id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,trimx),Y=c(starty,trimy))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=trimdist,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=trimx,EndY=trimy,Azimuth=perpaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_perp_trim",sep=""), arcgis=T)


id.field <- seq(1,length(startx),by=1)
dd <- data.frame(Id=c(id.field,id.field),X=c(startx,nearx),Y=c(starty,neary))
ddTable <- data.frame(Id=id.field,Transect=id.field,TranSpace=transpace,TranDist=neardist,Location=texttable3[,"baseloc1"],MaxBNum=basemaxbnum,BaseOrder=baseorder,OFFshore=baseoffshore,CastDir=basecastdir,BASE_LOC=texttable3[,"baseloc2"],StartX=startx,StartY=starty,EndX=nearx,EndY=neary,Azimuth=nearaz,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects_near_inner",sep=""), arcgis=T)




#detach("package:shapefiles")
detach(data.frame(finaltable3))

#######status checkpoint (100%)
    info <- sprintf("%d%% done", 100)
    setTkProgressBar(pb, 100, sprintf("AMBUR: Trnasect casting (%s)", info), info)

###########################


}

