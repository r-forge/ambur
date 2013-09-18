ambur.transects <-
function(userinput1=50,userinput2=500,userinput3=5,userinput4=5,userinput5=90) {

#require(tcltk)
#require(rgdal)
#require(rgeos)
#require(spatstat)
#require(sp)

transpace <- userinput1
tranlength <- userinput2
innersample <- userinput3
outersample <- userinput4
fsamp <- userinput5


#transpace <- 50
#tranlength <- 500
#innersample <- 5
#outersample <- 5
#fsamp <- 90

#fsamp <- c(seq(2,178,by=5),178) # radiating transects

tkmessageBox(message = "Please select the outer baseline shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the inner baseline shapefile...")
getdata2 <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create("AMBUR_transects", showWarnings=FALSE)
setwd("AMBUR_transects")

dir.create(paste(time.stamp2," ","transects",sep=""))
setwd(paste(time.stamp2," ","transects",sep=""))


pb <- tkProgressBar("AMBUR: progress bar", "Reading GIS data...", 0, 100, 10)

############################## build functions
#build function for points along a line (modified to add the start and end points)
###################################################################################
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
########################################################################
perp.trans <- function(OutX,OutY,transpace,tranlength,fsamp,castdir,baseshapeID) {

castdir <- unique(castdir)[1]

basecastdirx <- castdir

baselineID <- unique(baseshapeID)[1]

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
baselineid <- rep(baselineID, length(outer.basepts[,1]))


 perp.table <- cbind(baselineid,transectid,t.startx,t.starty,t.endx,t.endy,t.azimuth,t.length)
 colnames(perp.table) <- c("Base_ID","Transect","StartX","StartY","EndX","EndY","Azimuth","TranDist")

 return(perp.table)
 
 }

 Pcnt.Complete <-  25
info <- sprintf("%d%% Casting perpendicular transects...", Pcnt.Complete)
setTkProgressBar(pb, 25 , sprintf("AMBUR: Transects (%s)", info), info)

#################################################################################################################################
### cast perpendicular transects
###################################################################



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

attrtable$Id <- seq(0,max(length(attrtable[,1]))-1,1)  #repair baseline attr table for sequential IDs to match with bbb$shapeID
bbb <- data.frame(merge(outerbase.tab,attrtable,by.x= "shapeID" ,by.y = "Id", all.x = TRUE,sort=FALSE))


### cast transects for individual polylines based on unique IDs
Baseline.Factor <- factor(bbb$BaseOrder)

trans.indv <- sapply(levels(Baseline.Factor), function(x) data.frame(perp.trans(bbb$baseX[bbb$BaseOrder == x],bbb$baseY[bbb$BaseOrder == x],transpace=transpace,tranlength=tranlength,fsamp=fsamp,castdir=bbb$CastDir[bbb$BaseOrder == x],baseshapeID=bbb$shapeID[bbb$BaseOrder == x])) ,simplify = FALSE)
trans.indv2 <- data.frame(do.call("rbind", trans.indv))

trans.indv2 <- data.frame(do.call("rbind", trans.indv))
trans.indv2$Transect <- seq(1,length(trans.indv2$Transect),1)
trans.indv2$TranDist <- abs(trans.indv2$TranDist)


#trans.all <- data.frame(sapply(levels(Baseline.Factor), function(x) data.frame(perp.trans(bbb$baseX,bbb$baseY,transpace=transpace,tranlength=tranlength,fsamp=fsamp,castdir=bbb$CastDir,baseshapeID=bbb$shapeID)) ,simplify = FALSE) )
#colnames(trans.all) <- c("Base_ID","Transect","StartX","StartY","EndX","EndY","Azimuth","TranDist")
#if(max(unique(bbb$baseID)) > 1)  test.tran <- trans.indv2 else test.tran <- trans.all


test.tran <- trans.indv2
test.tran$TranSpace <- transpace


perp.trans.tab <- data.frame(merge(test.tran,attrtable,by.x= "Base_ID" ,by.y = "Id", all.x = TRUE,sort=FALSE))
perp.trans.tab$Id <- perp.trans.tab$Transect

#establish the required fields
reqfields <- c("Id","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth")

#cull extraneous fields from merging attribute tables and order the fields by reqfields and fix the row names for SpatialLinesDataFrame
perp.transects <- data.frame(perp.trans.tab[match(toupper(reqfields), toupper(names(perp.trans.tab)))])
perp.transects$Creator <- "R AMBUR"
row.names(perp.transects) <- as.character(seq(0,length(perp.transects$Transect)-1,1))

### build spatial lines for final perpendicular transects shapefile
Transect.Factor <- factor(perp.transects$Transect)
shape.final <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(perp.transects$StartX[perp.transects$Transect == x], perp.transects$EndX[perp.transects$Transect == x]), y=c(perp.transects$StartY[perp.transects$Transect == x],perp.transects$EndY[perp.transects$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, perp.transects)

 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info
  
  proj4string(shape.final3) <- projectionString


#create shapefile and write it to the working directory
writeOGR(shape.final3, ".", "perp_transects", driver="ESRI Shapefile")



Pcnt.Complete <-  50
info <- sprintf("%d%% Casting near transects...", Pcnt.Complete)
setTkProgressBar(pb, 50 , sprintf("AMBUR: Transects (%s)", info), info)

#################################################################################################################################
### cast near transects
###################################################################

 ###### break down inner baseline into simply points with IDs
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
innerbase.tab$dummyID <-  1  #initially treat multiple polylines as one line and break up later


###build function to get segment coordinates
conv.segs <- function(InX,InY) {

Cx <- InX
Cy <- InY

Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

cbind(Cx,Cy,Cx2,Cy2)

}


Baseline.Factor2 <- factor(innerbase.tab$shapeID)
inner.segs <- (sapply(levels(Baseline.Factor2), function(x) conv.segs(innerbase.tab$baseX[innerbase.tab$shapeID == x], innerbase.tab$baseY[innerbase.tab$shapeID == x])  ,simplify = FALSE))

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


####################################end adj.project2segment function


Pcnt.Complete <-  75
info <- sprintf("%d%% Calculating near distances...", Pcnt.Complete)
setTkProgressBar(pb, 75 , sprintf("AMBUR: Transects (%s)", info), info)

###build objects
TX.w <- owin()
TX.w <- owin(c(min(perp.transects$StartX-100000),max(perp.transects$StartX+100000)), c(min(perp.transects$StartY-100000),max(perp.transects$StartY+100000)))
TX <- ppp(perp.transects$StartX,perp.transects$StartY, window=TX.w)

TY.w <- owin()
TY.w <- owin(c(min(inner.segs$Cx-100000),max(inner.segs$Cx+100000)), c(min(inner.segs$Cy-100000),max(inner.segs$Cy+100000)))
TY <- psp(inner.segs$Cx,inner.segs$Cy,inner.segs$Cx2,inner.segs$Cy2,window=TY.w)


near.analysis <- adj.project2segment(TX,TY)
  
Xproj <- near.analysis$Xproj 
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #2nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)
Xproj$x <- ifelse(is.na(Xproj$x),Xproj$x[-1],Xproj$x) #3nd iteration:this code finds missing NA values and forces line to next pt 
Xproj$y <- ifelse(is.na(Xproj$y),Xproj$y[-1],Xproj$y)

#transect filter using average window of azimuths
tnear.tab <- data.frame(cbind(TX$x,TX$y,Xproj$x,Xproj$y,near.analysis$d,near.analysis$mapXY,near.analysis$tp))
colnames(tnear.tab) <- c("T_x","T_y","In_x","In_y","Dist","Near_Seg","Pos")


Dx2_in <- tnear.tab[,3] - tnear.tab[,1]
Dy2_in <- tnear.tab[,4] - tnear.tab[,2]

in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tnear.tab[,3]- tnear.tab[,1])^2 +  (tnear.tab[,4] - tnear.tab[,2])^2)^(1/2)


####setup the attribute table and adjust values
near.transects <- perp.transects
near.transects$EndX <- tnear.tab$In_x
near.transects$EndY <- tnear.tab$In_y
near.transects$TranDist <- tnear.tab$Dist
near.transects$Azimuth <- in.az
row.names(near.transects) <- as.character(seq(0,length(near.transects$Transect)-1,1))

### build spatial lines for final near transects shapefile
Transect.Factor <- factor(near.transects$Transect)
shape.near <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(near.transects$StartX[near.transects$Transect == x], near.transects$EndX[near.transects$Transect == x]), y=c(near.transects$StartY[near.transects$Transect == x],near.transects$EndY[near.transects$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.near2 <- SpatialLines(shape.near)
shape.near3 <- SpatialLinesDataFrame(shape.near2, near.transects)

 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info
  
  proj4string(shape.near3) <- projectionString


#create shapefile and write it to the working directory
writeOGR(shape.near3, ".", "near_transects", driver="ESRI Shapefile")





Pcnt.Complete <-  90
info <- sprintf("%d%% Trimming perpendicular transects...", Pcnt.Complete)
setTkProgressBar(pb, 90 , sprintf("AMBUR: Transects (%s)", info), info)
#################################################################################################################################
##### create trimmed transects 
##################################################################


shape.prep3 <- shape.final3

int <- gIntersects(shapedata2, shape.prep3, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

for (i in seq(along=vec)) vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(shapedata2[i,], shape.prep3[int[,i],], byid=TRUE) else 0

cond <- lapply(vec, function(x) class(x) == "SpatialPoints")   # fixed to get single intersections and remove "zero" values from list
vec2 <- vec[unlist(cond)]

out <- do.call("rbind", vec2) 
rn <- row.names(out) 
nrn <- do.call("rbind", strsplit(rn, " ")) 

transID <- data.frame(nrn)[,2]
baseID <- data.frame(nrn)[,1]
INT_X <-  data.frame(coordinates(out))$x
INT_Y <-  data.frame(coordinates(out))$y




sortID <- seq(1,length(INT_X),1)
inter.data <- data.frame(INT_X,INT_Y,transID,baseID,sortID)




tran.data <- data.frame(perp.transects)

## fixed that transects are off by 1 in the intersection matrix because column id start with 0  (8-20-2011) 
tran.data$transdataID <-  tran.data$Transect - 1


tet <- merge(tran.data,inter.data , by.x = "transdataID", by.y = "transID", sort=FALSE, all.x=TRUE)
tet2 <- data.frame(tet[ order(tet[,"Transect"]) , ])



###added to correct for multiple interestions with the baseline  (8-20-2011) start:
Transect.Factor <- factor(tran.data$Transect)
tet2dist <- (((tet2[,"INT_X"]- tet2[,"StartX"])^2 +  (tet2[,"INT_Y"] - tet2[,"StartY"])^2)^(1/2))
tet2disttab <- data.frame(sapply(levels(Transect.Factor), function(x) min(tet2dist[tet2$Transect == x],na.rm=FALSE) ,simplify = TRUE))
tet2disttab2 <- data.frame(sapply(levels(Transect.Factor), function(x) tet2dist[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))
tet2intx <- data.frame(sapply(levels(Transect.Factor), function(x) tet2$INT_X[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))
tet2inty <- data.frame(sapply(levels(Transect.Factor), function(x) tet2$INT_Y[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))
tet3 <- data.frame(tran.data$Transect,tet2disttab,tet2disttab,tet2intx,tet2inty)
colnames(tet3) <- c("Transect","MinDist1","MinDist_check","INT_X","INT_Y")
##############end

### make new attribute table with filtered values
new_trandata <-  tran.data
new_trandata[,"EndX"]<- ifelse(is.na(tet3$INT_X) == TRUE, as.numeric(tran.data$EndX), as.numeric(tet3$INT_X))
new_trandata[,"EndY"] <- ifelse(is.na(tet3$INT_Y) == TRUE, as.numeric(tran.data$EndY), as.numeric(tet3$INT_Y))
new_trandata[,"TranDist"] <- (((new_trandata[,"EndX"]- new_trandata[,"StartX"])^2 +  (new_trandata[,"EndY"] - new_trandata[,"StartY"])^2)^(1/2))




if (sum(as.numeric(int)) == 0) new_trandata <- perp.transects






### build spatial lines for final trim transects shapefile
Transect.Factor3 <- factor(new_trandata$Transect)
shape.trim <- sapply(levels(Transect.Factor3), function(x)
list(Lines(list(Line(list(x=c(new_trandata$StartX[new_trandata$Transect == x], new_trandata$EndX[new_trandata$Transect == x]), y=c(new_trandata$StartY[new_trandata$Transect == x],new_trandata$EndY[new_trandata$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.trim2 <- SpatialLines(shape.trim)
shape.trim3 <- SpatialLinesDataFrame(shape.trim2, new_trandata)

   # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info
  
  proj4string(shape.trim3) <- projectionString


#create shapefile and write it to the working directory
writeOGR(shape.trim3, ".", "trim_perp_transects", driver="ESRI Shapefile")


Pcnt.Complete <-  100
info <- sprintf("%d%% Done!", Pcnt.Complete)
setTkProgressBar(pb, 100 , sprintf("AMBUR: Transects (%s)", info), info)



#################################################################################################################################
##### plot the results of transect casting   
##################################################################
    #plot(c(test.tran$StartX,test.tran$EndX),c(test.tran$StartY,test.tran$EndY),col="white",asp=1,xlab="X",ylab="Y")
    #segments(test.tran$StartX,test.tran$StartY,test.tran$EndX,test.tran$EndY)
    #lines(bbb$baseX,bbb$baseY,col="red")
    
       #plot(c(trans.all$StartX,trans.all$EndX),c(trans.all$StartY,trans.all$EndY),col="white",asp=1,xlab="X",ylab="Y")
    #segments(trans.all$StartX,trans.all$StartY,trans.all$EndX,trans.all$EndY)
    #lines(bbb$baseX,bbb$baseY,col="red")
    
       #plot(c(trans.indv2$StartX,trans.indv2$EndX),c(trans.indv2$StartY,trans.indv2$EndY),col="white",asp=1,xlab="X",ylab="Y")
    #segments(trans.indv2$StartX,trans.indv2$StartY,trans.indv2$EndX,trans.indv2$EndY)
    #lines(bbb$baseX,bbb$baseY,col="red")
    
   ####intersect




 }