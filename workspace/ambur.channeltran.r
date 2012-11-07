ambur.channeltran <-
function(ptspace=50,offsetdist=-5,buffdist=0,buffnum=20,winsize=5, indv=1) {

###enter a negative offsetdist number to offset left of the polyline,  positive for right of line

require(tcltk)
require(rgdal)
require(rgeos)


#ptspace <- 50 # for testing
#offsetdist <- -25 # for testing
#buffdist <- 100 #for testing
#buffnum <- 20 #for testing
#winsize <- 5 #for testing
#indv <- 1 #for testing

tkmessageBox(message = "Please select the polyline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
shapedata <- as(shapedata, "SpatialLinesDataFrame")
attrtable <- data.frame(shapedata)

####get the MaxBNum value from baseline to get the max distance of the outer contour
buffdist <-  ifelse(buffdist == 0, attrtable$MaxBNum/buffnum, buffdist)

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


pts.output3 <- rbind(pts.output,pts.output2)

  outputname <- paste(locname,"_basepts_both",sep="")

   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output3) <- projectionString

writeOGR(pts.output3, ".", outputname, driver="ESRI Shapefile")

plot(pts.output3)

###################################################################################################################################################

###########################   start buffers function

###################################################################################################################################################



totbuffers <- seq(buffdist,buffdist*buffnum,buffdist)

#tkmessageBox(message = "Please select the baseline polyline shapefile...")
#getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
#shapename <- gsub(".shp", "", basename(getdata))
#shapedata <- readOGR(getdata,layer=shapename)
#shapedata <- as(shapedata, "SpatialLinesDataFrame")
#attrtable <- data.frame(shapedata)

#workingdir <- dirname(getdata)
#setwd(workingdir)
#time.stamp1 <- as.character(Sys.time())
#time.stamp2 <- gsub("[:]", "_", time.stamp1)


#dir.create("AMBUR_buffer", showWarnings=FALSE)
#setwd("AMBUR_buffer")

#dir.create(paste(time.stamp2," ","buffers",sep=""))
#setwd(paste(time.stamp2," ","buffers",sep=""))



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


 Baseline.Factor <- unique(outerbase.tab$baseID)
 Buffer.Factor <- totbuffers


width.factor <- rep(Buffer.Factor, nrow(shapedata))
line.id <- sort(rep(1:nrow(shapedata),buffnum))

l <- vector("list", nrow(shapedata)*buffnum)

for (i in seq_len(nrow(shapedata)*buffnum)) {

        l[[i]] <- gBuffer(shapedata[line.id[i],], width = width.factor[i], byid=TRUE, capStyle="SQUARE")
        l[[i]] <- spChFIDs(l[[i]], as.character(i))

}

vec2 <- unlist(l)


out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))

final_buffers <- as(out, "SpatialLines")

plot(final_buffers)


buff.tab <- data.frame(line.id,width.factor)
colnames(buff.tab) <- c("baseID","distance")

out.buffers <- SpatialLinesDataFrame(final_buffers, buff.tab)

 locname1 <- attrtable$Location[1]
  locname <- gsub(" ", "_", locname1)
      outputname <- paste(locname,"_buffers_",buffdist,sep="")

writeOGR(out.buffers, ".", outputname, driver="ESRI Shapefile")




#get.exp <- as.character(round(length(buff.tab$distance)/c(length(buff.tab$distance),8,4,2,1),0))  # try to get exponential values
#buff.tab2 <- buff.tab[rownames(buff.tab) %in% get.exp,]


#out.buffers2 <- SpatialLinesDataFrame(final_buffers, buff.tab)

#out.buffers3 <- out.buffers2[as.numeric(get.exp),]

#writeOGR(out.buffers3, ".", "polyline_buffers2", driver="ESRI Shapefile")


###################################################################################################################################################

###########################   ambur pathway transect function

###################################################################################################################################################

require(spatstat)
require(maptools)

#tkmessageBox(message = "Please select the baseline points shapefile...")
#getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
#shapename <- gsub(".shp", "", basename(getdata))
shapedata <- pts.output3
attrtable <- data.frame(shapedata)

#workingdir <- dirname(getdata)
#setwd(workingdir)


#tkmessageBox(message = "Please select the buffers polyline shapefile...")
#getdata2 <- tk_choose.files(default = "*.shp",multi = FALSE)
#shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- out.buffers
attrtable2 <- data.frame(shapedata2)


#time.stamp1 <- as.character(Sys.time())
#time.stamp2 <- gsub("[:]", "_", time.stamp1)


#dir.create("AMBUR_near", showWarnings=FALSE)
#setwd("AMBUR_near")

#dir.create(paste(time.stamp2," ","near",sep=""))
#setwd(paste(time.stamp2," ","near",sep=""))



###### break down outer baseline into simply points with IDs
crdl0 <- coordinates(shapedata)



###### break down polyline into segments with IDs
crdl0a <- coordinates(shapedata2)
crd.1a <- sapply(crdl0a, function(x) do.call("rbind", x),simplify = FALSE)
crd.2a    <- do.call("rbind", crd.1a)
crd.3a <- as.numeric(sapply(crd.1a, function(x) max(row(x)) ,simplify = TRUE))
crd.len.testa <- as.numeric(length(crd.3a))
if(crd.len.testa <= 1) crd.repa <-  1 else crd.repa <- seq(1, length(crd.3a),1)
basepointIDsa <- rep(crd.repa,crd.3a)
baseshapeIDsa <- basepointIDsa - 1
sortshapeIDsa <- seq(1,length(basepointIDsa),1)
basexa <- crd.2a[,1]
baseya <- crd.2a[,2]

innerbase.tab <- data.frame(sortshapeIDsa,baseshapeIDsa,basepointIDsa,basexa,baseya)
colnames(innerbase.tab) <- c("sortshapeID","shapeID","baseID","baseX", "baseY")
innerbase.tab$shapeID <-  innerbase.tab$shapeID + 1



###build function to get segment coordinates
conv.segs <- function(InX,InY) {

Cx <- InX
Cy <- InY

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

cbind(Cx,Cy,Cx2,Cy2)

}




#####build spatstat functions and objects
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




#set up the points in spatstat format
TX.w <- owin()
TX.w <- owin(c(min(crdl0[,1]-100000),max(crdl0[,1]+100000)), c(min(crdl0[,2]-100000),max(crdl0[,2]+100000)))
TX <- ppp(crdl0[,1],crdl0[,2], window=TX.w)




blah <- data.frame(TX$x,TX$y,TX$x,TX$y,sapply(1:3, function(x) data.frame(x)),attrtable$BaseID,sapply(1:2, function(x) data.frame(x)))

colnames(blah) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos","baseID","contourID","tranID")

row.removal <- length(blah[,1])

contourdiff <- attrtable2$distance[2] -attrtable2$distance[1]
minicontour <- attrtable2$distance[1]

blah$contourID <- minicontour
blah$tranID <- as.numeric(rownames(blah))

 ##################function

 Baseline.Factor <- unique(innerbase.tab$baseID)

for (i in Baseline.Factor) {

TY <- as.psp.SpatialLines(shapedata2[i,])




#calculate nearest points
near.analysis <- adj.project2segment(ppp(blah$In_x[blah$contourID == attrtable2$distance[i]],blah$In_y[blah$contourID == attrtable2$distance[i]], window=TX.w),TY)

Xproj <- near.analysis$Xproj
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

#create the data table
tnear.tab <- data.frame(blah$In_x[blah$contourID == attrtable2$distance[i]],blah$In_y[blah$contourID == attrtable2$distance[i]],Xproj$x,Xproj$y,near.analysis$d,near.analysis$mapXY,near.analysis$tp,sapply(1:3, function(x) data.frame(x)))
colnames(tnear.tab) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos","baseID","contourID","tranID")

tnear.tab$baseID <- i
tnear.tab$contourID <- ifelse(attrtable2$distance[i] == max(attrtable2$distance), max(attrtable2$distance), attrtable2$distance[i] + minicontour)
tnear.tab$tranID <- blah$tranID[blah$contourID == attrtable2$distance[i]]
#tnear.tab$T_x <- Xproj$x
#tnear.tab$T_y <- Xproj$y


blah <- rbind(blah, tnear.tab)


}

blah$contourID[is.na(blah$contourID)] <- attrtable2$distance[length(attrtable2$distance)]




 #########################################

 tnear.tabA <- blah[-(1:row.removal),]

  #tnear.tab <- blah[(length(tnear.tabA[,1])-max(row.removal)+1):length(tnear.tabA[,1]),]      ####edited this line of code

 tnear.tab[,1] <- crdl0[,1]
   tnear.tab[,2]<- crdl0[,2]

plot(blah$In_x,blah$In_y)
points(tnear.tabA$In_x,tnear.tabA$In_y,col="green")
points(tnear.tab$In_x,tnear.tab$In_y,col="red")

 #set up for distance an azimuth calcs
Dx2_in <- tnear.tab[,3] - tnear.tab[,1]
Dy2_in <- tnear.tab[,4] - tnear.tab[,2]

in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tnear.tab[,3]- tnear.tab[,1])^2 +  (tnear.tab[,4] - tnear.tab[,2])^2)^(1/2)


####setup the attribute table and adjust values
near.transects <- attrtable
near.transects$EndX <- tnear.tab$In_x
near.transects$EndY <- tnear.tab$In_y
near.transects$TranDist <- in.length
near.transects$Azimuth <- in.az
row.names(near.transects) <- as.character(row.names(near.transects))
near.transects$ID <- as.numeric(row.names(near.transects))
near.transects$Transect <- near.transects$ID

#reqfields <- c("Id","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth")

### build spatial lines for final near transects shapefile
#Transect.Factor <- factor(near.transects$ID)
#shape.near <- sapply(levels(Transect.Factor), function(x)
#list(Lines(list(Line(list(x=c(near.transects$T_x[near.transects$ID == x], near.transects$In_x[near.transects$ID == x]), y=c(near.transects$T_y[near.transects$ID == x],near.transects$In_y[near.transects$ID == x])))), ID=x))
#,simplify = TRUE)
#shape.near2 <- SpatialLines(shape.near)
#shape.near3 <- SpatialLinesDataFrame(shape.near2, near.transects)

 # Note that readOGR method reads the .prj file when it exists
  # projectionString <- proj4string(shapedata2) # contains projection info

 # proj4string(shape.near3) <- projectionString


#create shapefile and write it to the working directory
#writeOGR(shape.near3, ".", "near_transects9", driver="ESRI Shapefile")

#create csv file and write it to the working directory
#write.table(shape.near3, file = "near_results.csv", quote = FALSE, sep = ",", row.names = FALSE)

near.transects$Baseline <- near.transects$BaseID

reqfields <- c("ID","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth")

#remove.fields <- c("BaseID","sortID","OBJECTID","SHAPE_Leng","coords.x1","coords.x2","Id")

near.transects <- near.transects[,colnames(near.transects) %in% reqfields]


#try removing duplicate
rrr <- blah[!duplicated(blah$In_x), ]

rrr.debug <- blah[duplicated(blah$T_x), ]

plot(blah$In_x,blah$In_y)
points(rrr$T_x,rrr$T_y,col="green")
points(rrr$In_x,rrr$In_y,col="red")

transect.num2 <- rrr$tranID
transect.num <- blah$tranID


#write shapefiles of the transects
library(shapefiles)

#dd <- data.frame(Id=transect.num,X=blah$T_x,Y=blah$T_y)
#ddTable <- data.frame(near.transects)
#ddShapefile <- convert.to.shapefile(dd, ddTable, "ID", 3)
#write.shapefile(ddShapefile, paste("pathway_transects",sep=""), arcgis=T)



near.transects.rep <-  near.transects[unique(near.transects$ID) %in%  unique(transect.num2),]   ### added to remove missing transects from original data


 locname1 <- attrtable$Location[1]
  locname <- gsub(" ", "_", locname1)
      outputname <- paste(locname,"_transects_trim",sep="")


dd <- data.frame(Id=transect.num2,X=rrr$In_x,Y=rrr$In_y)
ddTable <- data.frame(near.transects.rep)
ddShapefile <- convert.to.shapefile(dd, ddTable, "ID", 3)
write.shapefile(ddShapefile, outputname, arcgis=T)


outputname2 <- paste(locname,"_transects_path",sep="")

transect.num2 <- c(near.transects$Transect, near.transects$Transect)
dd <- data.frame(Id=near.transects$Transect,X=c(near.transects$StartX,near.transects$EndX),Y=c(near.transects$StartY,near.transects$EndY))
ddTable <- data.frame(near.transects)
ddShapefile <- convert.to.shapefile(dd, ddTable, "ID", 3)
write.shapefile(ddShapefile, outputname2, arcgis=T)

detach("package:shapefiles")


#tidy up and remove objects
#rm(list = ls())


  ####################################################################################################################################################################################################################

  # Filter.transects

 ####################################################################################################################################################################################################################





#set up master data table
trandata <- data.frame(near.transects)

trandata$BaseOrder <- trandata$BaseOrder *  trandata$CastDir  ##for doublesided casted transects on one baseline

##establish moving average function
move.avg <- function(x,n=winsize){filter(x,rep(1/n,n), sides=2)}


###established azimuth filter function
filter.azimuths <- function(az.data){
cos.az <- cos(az.data * pi/180)
sin.az <- sin(az.data * pi/180)
avg.cosx <- move.avg(cos.az)
avg.siny <- move.avg(sin.az)
filter.az.rads <- atan2(avg.siny, avg.cosx)
filter.az.deg <- ifelse( (filter.az.rads *(180/pi)) < 0, (filter.az.rads *(180/pi)) + 360, (filter.az.rads *(180/pi)))
filter.az.final <- ifelse(is.na(filter.az.deg) == TRUE, trandata$Azimuth, filter.az.deg)
return(filter.az.final)}

### filter across all baselines
filter.az.all <- filter.azimuths(trandata$Azimuth)

### filter across individual baselines
Baseline.Factor <- factor(trandata$BaseOrder)
filter.az.indv <- sapply(levels(Baseline.Factor), function(x) filter.azimuths(trandata$Azimuth[trandata$BaseOrder == x]) ,simplify = TRUE)
filter.az.indv <- unlist(filter.az.indv, use.names = FALSE)

### get the appropriate filter azimuth based on user option
filter.az.sel <- if(indv == 1) filter.az.indv else filter.az.all

filter.az.sel[1:floor(winsize/2),] <- 9999  #####identify azimuths that are not filtered base on start and end windows
filter.az.sel[(length(filter.az.sel[,1]) - floor(winsize/2)+1):(length(filter.az.sel[,1])),] <- 9999

filter.az.sel <- ifelse(filter.az.sel == 9999, trandata$Azimuth,filter.az.sel)

### get filter transects ending XY coordinates
filter2.x <- sin((filter.az.sel * pi/180)) * (trandata$TranDist) + trandata$StartX
filter2.y <- cos((filter.az.sel * pi/180)) * (trandata$TranDist) + trandata$StartY


######repair the datafile
### make new attribute table with filtered values
new_trandata <-  trandata
new_trandata$Azimuth <-  as.vector(filter.az.sel)
new_trandata$EndX  <- as.vector(filter2.x)
new_trandata$EndY <-  as.vector(filter2.y)




row.names(new_trandata) <- seq(0,length(trandata[,1])-1,1)

### build spatial lines for final filtered transects shapefile
Transect.Factor <- factor(new_trandata$Transect)
shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$StartX[new_trandata$Transect == x], new_trandata$EndX[new_trandata$Transect == x]), y=c(new_trandata$StartY[new_trandata$Transect == x],new_trandata$EndY[new_trandata$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)

   # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info

  proj4string(shape.final3) <- projectionString

  outputname3 <- paste(locname,"_filtered_transects",sep="")

#create shapefile and write it to the working directory
writeOGR(shape.final3, ".", outputname3 , driver="ESRI Shapefile")


}