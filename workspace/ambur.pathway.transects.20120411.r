# Get required packages
require(tcltk)
require(rgdal)
require(rgeos)
require(spatstat)
require(maptools)

tkmessageBox(message = "Please select the points shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the polyline shapefile...")
getdata2 <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create("AMBUR_near", showWarnings=FALSE)
setwd("AMBUR_near")

dir.create(paste(time.stamp2," ","near",sep=""))
setwd(paste(time.stamp2," ","near",sep=""))



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
tnear.tab$contourID <- attrtable2$distance[i] + contourdiff
tnear.tab$tranID <- blah$tranID[blah$contourID == attrtable2$distance[i]]
#tnear.tab$T_x <- Xproj$x
#tnear.tab$T_y <- Xproj$y


blah <- rbind(blah, tnear.tab)


}

 #########################################

 tnear.tabA <- blah[-(1:row.removal),]

  tnear.tab <- blah[(length(tnear.tabA[,1])-max(row.removal)+1):length(tnear.tabA[,1]),]

 tnear.tab[,1] <- crdl0[,1]
   tnear.tab[,2]<- crdl0[,2]
  

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


transect.num <- blah$tranID


#write shapefiles of the transects
library(shapefiles)

dd <- data.frame(Id=transect.num,X=blah$T_x,Y=blah$T_y)
ddTable <- data.frame(near.transects)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("pathway_transects",sep=""), arcgis=T)


transect.num2 <- c(near.transects$Transect, near.transects$Transect)
dd <- data.frame(Id=near.transects$Transect,X=c(near.transects$StartX,near.transects$EndX),Y=c(near.transects$StartY,near.transects$EndY))
ddTable <- data.frame(near.transects)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("pathway_transects_endpts",sep=""), arcgis=T)



