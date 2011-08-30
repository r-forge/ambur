#use this for points and polygons
require(shapefiles)

winDialog("ok","Please select the transect points...")

path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)

mydata <- data.frame(convert.to.simple(my.shapefile))

path1aa <- dirname(path1)
setwd(path1aa)

winDialog("ok","Please select the polygon area...")

path2 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path2)

path2a <- paste(substr(path2,1,max(del.ext)-4),".shp",sep="")

my.shapefile2 <- read.shp(path2a)

mydata2 <- convert.to.simple(my.shapefile2)

polyboundary <- as.matrix(cbind(mydata2[,2],mydata2[,3]))






arealSubPolygons <- function (x, y = NULL, IDs = row.names(x), boundary, min.area.pct = 0.0001) 
{

#x = mydata[,"X"]
#y = mydata[,"Y"]
#IDs = row.names(mydata)
#boundary = polyboundary
#min.area.pct = 0.0001

    stopifnot(require(tripack))
    stopifnot(require(gpclib))
    stopifnot(require(sp))
    stopifnot(require(maptools))

as.Polygons.gpc.poly <- function(x, ID) {
	thisPolys <- lapply(get.pts(x), function(p) {
		Polygon(rbind(cbind(p$x,p$y),c(p$x[1],p$y[1])), hole=p$hole)
	})
	Polygons(thisPolys, ID)
}    
    
    
    xy <- xy.coords(x, y)
    stopifnot(length(IDs) == length(xy$x))
    xy <- data.frame(x = xy$x, y = xy$y)
    boundary <- as(boundary, "gpc.poly")
    dummies <- data.frame(x = c(-1, -1, 1, 1), y = c(-1, 1, -1, 
        1)) * 10 * max(abs(xy))
    xy <- rbind(xy, dummies)
    vpolys <- voronoi.polygons(voronoi.mosaic(xy))
    vpolys <- lapply(vpolys, as, "gpc.poly")
    subpolys <- lapply(vpolys, intersect, boundary)
    if (min.area.pct > 0) {
        min.area.frac <- min.area.pct/100
        totalArea <- area.poly(boundary)
        ok <- rep(TRUE, length(IDs))
        for (i in seq(along = IDs)) {
            if (area.poly(subpolys[[i]])/totalArea < min.area.frac) 
                ok[i] <- FALSE
        }
        if (any(!ok)) 
            return(arealSubPolygons(xy$x[ok], xy$y[ok], IDs[ok], 
                boundary = boundary, min.area.pct = min.area.pct))
    }
    thisSPs <- list()
    for (i in seq(along = IDs)) {
        if (length(get.pts(subpolys[[i]])) == 0) 
            next
        ii <- length(thisSPs) + 1
        thisSPs[[ii]] <- as.Polygons.gpc.poly(subpolys[[i]], 
            IDs[i])
    }
    SpatialPolygons(thisSPs)
}

#remove duplicate points
my.df2 <- mydata[order(mydata[,"Id"], decreasing=FALSE),] 
my.df3 <- my.df2[ !duplicated( my.df2[,c('X','Y')] ), ] 


test <- arealSubPolygons(my.df3[,"X"], my.df3[,"Y"], IDs = row.names(my.df3), polyboundary, min.area.pct = 0.0001)

centroids <- getSpPPolygonsLabptSlots(test)
x <- centroids[,1]
y <- centroids[,2]
testshape <- SpatialPolygonsDataFrame(test, data=data.frame(x=x, y=y, row.names=getSpPPolygonsIDSlots(test)))


writePolyShape(testshape, "theissen_polys")