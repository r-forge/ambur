require(tcltk)
require(akima)
require(rgdal)
require(maptools)

cellsize <- 0.2


tkmessageBox(message = "Please select the '*.csv' file...")


data.path <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.csv(data.path, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(data.path)
setwd(dir.path)


#establish xyz values from inner to outer
Cx <- mydata$adjE
Cy <- mydata$adjN
Cz <- mydata$Mag




akima.li <- interp(Cx, Cy, Cz,linear=TRUE,extrap=FALSE,xo=seq(from=min(Cx),to=max(Cx),by=cellsize),yo=seq(from=min(Cy),to=max(Cy),by=cellsize),duplicate="mean")



###plot the results of gridding and contouring
rgb.palette <- colorRampPalette(c("blue","green","yellow","orange","red"), space = "rgb")

image  (akima.li, main = "Interpolation Grid and Contours", useRaster = TRUE,col= rgb.palette(30))
contour(akima.li, add = TRUE, col = "white",nlevels = 10)
points(Cx,Cy, pch = 1, cex = 1, col = "white")

plot(Cx,Cy, pch = 1, cex = 1, col = "grey",xlab="",ylab="", asp=F)

Contours.GIS <-  ContourLines2SLDF(contourLines(akima.li,nlevels=10))
writeLinesShape(Contours.GIS, "contours")


###convert the akima image to grid  SpatialPixelsDataFrame

image2Grid <-
function (x, p4 = NA)
{
   if (!all(names(x) %in% c("x", "y", "z")))
        stop("image must have components x, y, and z")
   cells.dim <- dim(x$z)
   xx <- x$x
    yy <- x$y
    lx <- length(xx)
    ly <- length(yy)
    if (all(c(lx, ly) == (cells.dim + 1))) {
        print("corners")
        xx <- xx[-1] - diff(xx[1:2])/2
        yy <- yy[-1] - diff(yy[1:2])/2
    }
    SpatialGridDataFrame(GridTopology(c(xx[1], yy[1]), c(diff(xx[1:2]),
        diff(yy[1:2])), cells.dim), data.frame(z = as.vector(x$z[,
        ncol(x$z):1])), proj4string = CRS(p4))
 }
 
 akima.li.grid <- image2Grid(akima.li)



#writeOGR(akima.li.grid, ".", "ambur_grid", driver="ESRI Shapefile")
writeGDAL(akima.li.grid, "ambur_grid_int.tif", drivername="GTiff", type="Int32",mvFlag=-9999) 
writeGDAL(akima.li.grid, "ambur_grid_flt.tif", drivername="GTiff", type="Float32",mvFlag=-9999)