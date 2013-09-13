ambur.neartran <-      

function(sampledist=50) {

#sampledist=50

# Get required packages
#require(tcltk)
#require(rgdal)
#require(rgeos)
#require(spatstat)

#require(maptools)  #don't need as of yet

tkmessageBox(message = "Please select the outer polyline shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the inner polyline shapefile...")
getdata2 <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create("AMBUR_near", showWarnings=FALSE)
setwd("AMBUR_near")

dir.create(paste(time.stamp2," ","near",sep=""))
setwd(paste(time.stamp2," ","near",sep=""))

####################set up a progress bar for the sapply function
############build progress bar for building transect lines
sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}
########################end function



 pb <- tkProgressBar("AMBUR: progress bar", "Creating near transects...", 0, 100, 1)


##########################################
Baseline.Factor <- factor(shapedata$BaseOrder)

setTkProgressBar(pb, 10 , "AMBUR: progress bar", "Calculating transect locations...")

trans.indv <- sapply_pb(levels(Baseline.Factor), function(x) data.frame("BaseOrder"=unique(shapedata$BaseOrder[shapedata$BaseOrder == x]),coordinates(spsample(shapedata[shapedata$BaseOrder == x,], round(sum(SpatialLinesLengths(shapedata[shapedata$BaseOrder == x,]))/sampledist,0), offset=0.000000, "regular") )) ,simplify = FALSE)


trans.indv2 <- data.frame(do.call("rbind", trans.indv))

crdl0 <- trans.indv2
crdl0$BaseOrder <- as.numeric(trans.indv2$BaseOrder)

new_attributes <- merge(crdl0, attrtable, by="BaseOrder",all.x=TRUE,sort=TRUE)

#shapedata[shapedata$BaseOrder == x,]


###### break down polyline into segments with IDs to convert to segments for near analysis
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

setTkProgressBar(pb, 30 , "AMBUR: progress bar", "Preparing inner baseline...")

Baseline.Factor2 <- factor(innerbase.tab$shapeID)
inner.segs <- (sapply_pb(levels(Baseline.Factor2), function(x) conv.segs(innerbase.tab$baseX[innerbase.tab$shapeID == x], innerbase.tab$baseY[innerbase.tab$shapeID == x])  ,simplify = FALSE))

crd.segs2   <- do.call("rbind", inner.segs)   #repair to handle multiple inner baselines better

inner.segs <- data.frame(crd.segs2)
 inner.segs[,1] <- as.numeric(inner.segs[,1])
    inner.segs[,2] <- as.numeric(inner.segs[,2])
       inner.segs[,3] <- as.numeric(inner.segs[,3])
           inner.segs[,4] <- as.numeric(inner.segs[,4])

colnames(inner.segs) <- c("Cx","Cy","Cx2","Cy2")



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

setTkProgressBar(pb, 50, "AMBUR: progress bar", "Performing near analysis...")

#set up the points in spatstat format
TX.w <- owin()
TX.w <- owin(c(min(crdl0[,2]-100000),max(crdl0[,2]+100000)), c(min(crdl0[,3]-100000),max(crdl0[,3]+100000)))
TX <- ppp(crdl0[,2],crdl0[,3], window=TX.w)

#set up line segments in spatstat format
#convert polylines into individual line segments for analysis
TY.w <- owin()
TY.w <- owin(c(min(inner.segs$Cx-100000),max(inner.segs$Cx+100000)), c(min(inner.segs$Cy-100000),max(inner.segs$Cy+100000)))
TY <- psp(inner.segs$Cx,inner.segs$Cy,inner.segs$Cx2,inner.segs$Cy2,window=TY.w)

#calculate nearest points
near.analysis <- adj.project2segment(TX,TY)
  
Xproj <- near.analysis$Xproj 
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

#create the data table
tnear.tab <- data.frame(cbind(TX$x,TX$y,Xproj$x,Xproj$y,near.analysis$d,near.analysis$mapXY,near.analysis$tp))
colnames(tnear.tab) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos")

 #set up for distance an azimuth calcs
Dx2_in <- tnear.tab[,3] - tnear.tab[,1]
Dy2_in <- tnear.tab[,4] - tnear.tab[,2]

in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tnear.tab[,3]- tnear.tab[,1])^2 +  (tnear.tab[,4] - tnear.tab[,2])^2)^(1/2)


#reqfields <- c("Id","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth")

setTkProgressBar(pb, 70 , "AMBUR: progress bar", "Building final data tables")

####setup the attribute table and adjust values
near.transects <- data.frame("Id" = seq(1,length(tnear.tab[,1]),1),"Transect" = seq(1,length(tnear.tab[,1]),1),"TranSpace"=sampledist)
near.transects$TranDist <- tnear.tab$Dist
near.transects$Location <- new_attributes$Location
near.transects$MaxBNum <- new_attributes$MaxBNum
near.transects$BaseOrder <- new_attributes$BaseOrder
near.transects$OFFshore <- new_attributes$OFFshore
near.transects$CastDir <- new_attributes$CastDir
near.transects$BASE_LOC <- new_attributes$BASE_LOC
near.transects$StartX <- tnear.tab$T_x
near.transects$StartY <- tnear.tab$T_y
near.transects$EndX <- tnear.tab$In_x
near.transects$EndY <- tnear.tab$In_y
near.transects$Azimuth <- in.az
row.names(near.transects) <- as.character(row.names(near.transects))
#near.transects$ID <- as.numeric(row.names(near.transects))

  #test.validate <- data.frame(crdl0,near.transects)

setTkProgressBar(pb, 90 ,  "AMBUR: progress bar", "Building near transects shapefile...")

### build spatial lines for final near transects shapefile
Transect.Factor <- factor(near.transects$Transect)
shape.near <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(near.transects$StartX[near.transects$Transect == x], near.transects$EndX[near.transects$Transect == x]), y=c(near.transects$StartY[near.transects$Transect == x],near.transects$EndY[near.transects$Transect == x])))), ID=x))
,simplify = TRUE)
shape.near2 <- SpatialLines(shape.near)
shape.near3 <- SpatialLinesDataFrame(shape.near2, near.transects)

 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info
  
  proj4string(shape.near3) <- projectionString

setTkProgressBar(pb, 95 , "AMBUR: progress bar", "Writing final shapefile...")
#create shapefile and write it to the working directory
writeOGR(shape.near3, ".", "near_transects", driver="ESRI Shapefile")

#create csv file and write it to the working directory
#write.table(shape.near3, file = "near_results.csv", quote = FALSE, sep = ",", row.names = FALSE)


  setTkProgressBar(pb, 100 , "AMBUR: progress bar", "Done!")


 }