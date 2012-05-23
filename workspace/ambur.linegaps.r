ambur.linegaps <-
function(gapdist=50) {
require(rgdal)
require(rgeos)
require(tcltk)
require(spatstat)
require(maptools)

gapsize <- gapdist   

#gapsize <- 50

#open baseline file

tkmessageBox(message = "Please select a polyline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_Line_to_CSV", showWarnings=FALSE)
setwd("AMBUR_Line_to_CSV")

dir.create(paste(time.stamp2," ","vertices",sep=""))
setwd(paste(time.stamp2," ","vertices",sep=""))

#merge polyines with common end vertices
#shapedata <- gLineMerge(shapedata)




#################################################################
###### break down outer polyline into simply points with IDs
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

Baseline.Factor <- factor(bbb$shapeID)






#to even out the points for critical point analysis
startvertex <- sapply(levels(Baseline.Factor), function(x) data.frame(bbb$baseX[bbb$shapeID == x][1],bbb$baseY[bbb$shapeID == x][1]) ,simplify = FALSE)
startvertex.tab <- data.frame(do.call("rbind", startvertex))
colnames(startvertex.tab) <- c("x","y")

endvertex <- endvertex <- sapply(levels(Baseline.Factor), function(x) data.frame(bbb$baseX[bbb$shapeID == x][length(bbb$baseX[bbb$shapeID == x])],bbb$baseY[bbb$shapeID == x][length(bbb$baseY[bbb$shapeID == x])]) ,simplify = FALSE)
endvertex.tab <- data.frame(do.call("rbind", endvertex))
colnames(endvertex.tab) <- c("x","y")


#round down to get ID numbers
#outer.basepts4sm$t.ID <- floor(as.numeric(row.names(outer.basepts4sm)))

plot(startvertex.tab[,1],startvertex.tab[,2],col="blue")
points(endvertex.tab[,1],endvertex.tab[,2],col="red",pch=16)



test <- data.frame(rbind(startvertex.tab,endvertex.tab))

#remove duplicate vertices
#points.test <-  subset(test, !duplicated(x,y))   #works good but leaves one duplicate

points.test <- test[!duplicated(test$x,test$y) == TRUE,]



#set up the points in spatstat format
TX.w <- owin()
TX.w <- owin(c(min(points.test$x-1),max(points.test$x+1)), c(min(points.test$y-1),max(points.test$y+1)))
TX <- ppp(points.test$x,points.test$y, window=TX.w)


neardist <- nndist(TX)
 nearwhich <- nnwhich(TX)
 
   self <- (nearwhich[nearwhich] == seq(nearwhich))
   # plot them
   A <- TX[self]
   B <- TX[nearwhich[self]]
   plot(TX)
   segments(A$x, A$y, B$x, B$y)

 
 #plot(TX %mark% (nndist(TX)), markscale=1)




 #create the data table
tnear.tab <- data.frame(cbind(A$x,A$y,B$x,B$y))
colnames(tnear.tab) <- c("T_x","T_y","In_x","In_y")

 #set up for distance an azimuth calcs
Dx2_in <- tnear.tab[,3] - tnear.tab[,1]
Dy2_in <- tnear.tab[,4] - tnear.tab[,2]

in.az  <- ifelse(Dx2_in >= 0, 90 -(180/pi) * atan(Dy2_in/Dx2_in),270 -(180/pi) * atan(Dy2_in/Dx2_in))

in.length <- ((tnear.tab[,3]- tnear.tab[,1])^2 +  (tnear.tab[,4] - tnear.tab[,2])^2)^(1/2)

tnear.tab$in.length <- in.length

tnear.tab$in.az <- in.az

#filter out lines greater than gap size
tnear.tab <- tnear.tab[tnear.tab$in.length <= gapsize,]


#filter out duplicate lines
tnear.tab <- tnear.tab[!duplicated(tnear.tab$in.length) == TRUE,]

####setup the attribute table and adjust values
near.transects <- data.frame(seq(0,length(tnear.tab$T_x)-1,1))
colnames(near.transects) <- "ID2"
near.transects$StartX <- tnear.tab$T_x
near.transects$StartY <- tnear.tab$T_y
near.transects$EndX <- tnear.tab$In_x
near.transects$EndY <- tnear.tab$In_y
near.transects$TranDist <- tnear.tab$in.length
near.transects$Azimuth <- tnear.tab$in.az
row.names(near.transects) <- as.character(near.transects$ID2)


 Transect.Factor <- factor(near.transects$ID)
shape.near <- sapply(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(near.transects$StartX[near.transects$ID == x], near.transects$EndX[near.transects$ID == x]), y=c(near.transects$StartY[near.transects$ID == x],near.transects$EndY[near.transects$ID == x])))), ID=x))
,simplify = TRUE)
shape.near2 <- SpatialLines(shape.near)
shape.near3 <- SpatialLinesDataFrame(shape.near2, near.transects)

 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info
  
  proj4string(shape.near3) <- projectionString


#create shapefile and write it to the working directory
writeOGR(shape.near3, ".", "polyline_gaps", driver="ESRI Shapefile")







}

