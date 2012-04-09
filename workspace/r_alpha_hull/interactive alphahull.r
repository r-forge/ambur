require(rgdal)
require(tcltk)
library(alphahull)
require(rgeos)
require(rpanel)

alpha.h <- 10
bwidth <- 1
sdist <- 1
presample = 1
sample.distance <- sdist
buffer.width <- bwidth
alpha.units <- alpha.h


tkmessageBox(message = "Please select a polyline shapefile...")
getdata <- tk_choose.files(default = "*.shp",multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

mydata_outer <- attrtable 

workingdir <- dirname(getdata)
setwd(workingdir)



time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_envelope", showWarnings=FALSE)
setwd("AMBUR_envelope")

dir.create(paste(time.stamp2," ","envelope",sep=""))
setwd(paste(time.stamp2," ","envelope",sep=""))



#############################################################
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

####added if statements to keep function from bailing out on a polyline less than point spacing
if (sum(Segment.Length,na.rm = TRUE) < pspace)  {

t.startx <- Cx
t.starty <- Cy
t.azimuth <- Segment.Azimuth

}


if (sum(Segment.Length,na.rm = TRUE) > pspace)  {

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

}

cbind(t.startx,t.starty,t.azimuth)



}
###########################################################
ashape_to_SPLDF <- function(x)
	{
	if(class(x) != 'ashape')
		stop('this function only works with `ashape` class objects')


	
	# convert ashape edges to DF
	x.as.df <- as.data.frame(x$edges)
	
	# convert each edge to a line segment
	l.list <- list()
	l.list2 <- list()
	for(i in 1:nrow(x.as.df))
		{
		# extract line start and end points as 1x2 matrices
		p1 <- cbind(x.as.df$x1[i], x.as.df$y1[i])
		p2 <- cbind(x.as.df$x2[i], x.as.df$y2[i])
		# row-bind into 2x3 matrix
		l.list[[i]] <- Line(rbind(p1, p2))
		l.list2[[i]] <- Lines(l.list[[i]], ID=i)
		
    }
		
	# promote to Lines class, then to SpatialLines class

	
	# copy over CRS data from original point data
	l.spl <- SpatialLines(l.list2)
	
	# promote to SpatialLinesDataFrame, required for export to GRASS / OGR
	l.spldf <- SpatialLinesDataFrame(l.spl, data=data.frame(ID=1:nrow(x.as.df)), match.ID=TRUE)
	
	return(l.spldf)
	}

############################################################

###### break down outer polyline into simply points with IDs

if (presample == 1)  {

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
outer.sample <- sapply(levels(Baseline.Factor), function(x) data.frame(pts.along(bbb$baseX[bbb$baseID == x],bbb$baseY[bbb$baseID == x],sample.distance)) ,simplify = FALSE)

outer.sample <- sapply(levels(Baseline.Factor), function(x) data.frame(pts.along(bbb$baseX[bbb$baseID == 17],bbb$baseY[bbb$baseID == 17],sample.distance)) ,simplify = FALSE)

outer.basepts4smA <- data.frame(do.call("rbind", outer.sample))

#round down to get ID numbers
outer.basepts4smA$t.ID <- floor(as.numeric(row.names(outer.basepts4smA)))
outer.basepts4sm <- data.frame(c(basex,outer.basepts4smA$t.startx),c(basey,outer.basepts4smA$t.starty))
colnames(outer.basepts4sm) <- c("t.startx","t.starty")

}



if (presample != 1)  {
crdl0 <- coordinates(shapedata)
crd.1 <- sapply(crdl0, function(x) do.call("rbind", x),simplify = FALSE)
crd.2    <- do.call("rbind", crd.1)
basex <- crd.2[,1]
basey <- crd.2[,2]

outer.basepts4sm <- data.frame(basex,basey)
colnames(outer.basepts4sm) <- c("t.startx","t.starty")
}

#establish xy values
Cx <- outer.basepts4sm$t.startx
Cy <- outer.basepts4sm$t.starty



####alpha shape
test.spatpts1 <-  data.frame(Cx,Cy)
coordinates(test.spatpts1) <- data.frame(x=Cx, y=Cy)
test.spatpts2 <- remove.duplicates(test.spatpts1, zero = 0.0, remove.second = TRUE)



#############   interactive plot

if (interactive()) {

qq.draw <- function(panel) {
with(panel, {

   alpha   <- as.numeric(alpha)
   z <- ashape(panel$x, panel$y,alpha[1])
   x.as.spldf <- ashape_to_SPLDF(z)
   plot(x.as.spldf, main = paste("alpha =", round(alpha[1], 0)),asp=T)
     points(panel$x, panel$y, cex=0.1, pch=1, col='grey')
   lines(x.as.spldf, col='blue')
legend('bottomright', legend=(paste(alpha[1], "map units", " alpha-shape")), lty=c(1), pch=c(NA), col=c('blue'), bty='n')
 
   }) 
   panel
 
    }
   
panel <- rp.control(x = test.spatpts2$Cx, y = test.spatpts2$Cy, alpha = 500)
#rp.slider(panel, alpha, 0, 2000, qq.draw)

rp.textentry(panel, alpha, qq.draw, labels = c("alpha"), initval = c(500))

rp.do(panel, qq.draw)


} 



####calculate alpha shape
test.ashape <- ashape(test.spatpts2$Cx,test.spatpts2$Cy,alpha.units)   #alpha is in number of map units
x.as.spldf <- ashape_to_SPLDF(test.ashape)
writeOGR(x.as.spldf, ".", "alpha_shape", driver="ESRI Shapefile")


plot(0,0, type="n", axes=F, xlab="", ylab="")
text(0,0,"finished writing alpha_shape")


 
   
   
   


###########

####try to make a polygon
#################################################################
hullbuffer <- gBuffer(x.as.spldf, byid=FALSE, id=NULL, width=buffer.width)

buff.tab <- data.frame(ID="1")
row.names(buff.tab) <- "buffer"
hullbuffer2 <- SpatialPolygonsDataFrame(hullbuffer, buff.tab)


writeOGR(hullbuffer2, ".", "alpha_shape_buffer", driver="ESRI Shapefile")

plot(0,0, type="n", axes=F, xlab="", ylab="")
text(0,0,"finished writing alpha_shape_poly")


polyareas <- numeric(length(hullbuffer@polygons[[1]]@plotOrder))
polyareasid <- numeric(length(hullbuffer@polygons[[1]]@plotOrder))
	pl.list <- list()
	pl.list2 <- list()


for(i in 1:length(hullbuffer@polygons[[1]]@plotOrder)) {

polyareas[i] <-  hullbuffer@polygons[[1]]@Polygons[[i]]@area

polyareasid[i] <- i

		pl.list[[i]] <- Line(hullbuffer@polygons[[1]]@Polygons[[i]]@coords)
		pl.list2[[i]] <- Lines(pl.list[[i]], ID=i)

 }

areatab <- data.frame(polyareas,polyareasid)
areatabsort <- areatab[order(areatab$polyareas,areatab$polyareasid,decreasing = TRUE),]

pl.list3 <- SpatialLines(pl.list2)
pl.list4 <- SpatialLinesDataFrame(pl.list3, areatab)
writeOGR(pl.list4, ".", "alpha_shape_poly_all", driver="ESRI Shapefile")

plot(0,0, type="n", axes=F, xlab="", ylab="")
text(0,0,"finished writing alpha_shape_poly_all polyline")

buff1xy <- hullbuffer@polygons[[1]]@Polygons[[areatabsort[1,2]]]@coords
buff1xy2 <- Lines(Line(buff1xy),ID="1")
buff1xy3 <- SpatialLines(list(buff1xy2))
buff1xy3.tab <- data.frame(ID <- seq(1,length(buff1xy)/2,1))
hullbuffer3 <- SpatialLinesDataFrame(buff1xy3, buff1xy3.tab)
writeOGR(hullbuffer3, ".", "alpha_shape_poly_max1", driver="ESRI Shapefile")


buff1xy <- hullbuffer@polygons[[1]]@Polygons[[areatabsort[2,2]]]@coords
buff1xy2 <- Lines(Line(buff1xy),ID="1")
buff1xy3 <- SpatialLines(list(buff1xy2))
buff1xy3.tab <- data.frame(ID <- seq(1,length(buff1xy)/2,1))
hullbuffer3 <- SpatialLinesDataFrame(buff1xy3, buff1xy3.tab)
writeOGR(hullbuffer3, ".", "alpha_shape_poly_max2", driver="ESRI Shapefile")


buff1xy <- hullbuffer@polygons[[1]]@Polygons[[areatabsort[3,2]]]@coords
buff1xy2 <- Lines(Line(buff1xy),ID="1")
buff1xy3 <- SpatialLines(list(buff1xy2))
buff1xy3.tab <- data.frame(ID <- seq(1,length(buff1xy)/2,1))
hullbuffer3 <- SpatialLinesDataFrame(buff1xy3, buff1xy3.tab)
writeOGR(hullbuffer3, ".", "alpha_shape_poly_max3", driver="ESRI Shapefile")

buff1xy <- hullbuffer@polygons[[1]]@Polygons[[areatabsort[4,2]]]@coords
buff1xy2 <- Lines(Line(buff1xy),ID="1")
buff1xy3 <- SpatialLines(list(buff1xy2))
buff1xy3.tab <- data.frame(ID <- seq(1,length(buff1xy)/2,1))
hullbuffer3 <- SpatialLinesDataFrame(buff1xy3, buff1xy3.tab)
writeOGR(hullbuffer3, ".", "alpha_shape_poly_max4", driver="ESRI Shapefile")

plot(0,0, type="n", axes=F, xlab="", ylab="")
text(0,0,"finished!")
