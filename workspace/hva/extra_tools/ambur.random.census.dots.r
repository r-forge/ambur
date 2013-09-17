ambur.random.census.dots <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)
 require(maptools)

tkmessageBox(message = "Please select the Census Block shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

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



 ########
time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

#dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
#setwd("AMBUR_HVA_data_prep")

 attrtable$POP10[attrtable$POP10 == 0] <- 1

plot(shapedata, col="white") 
plot(spsample(shapedata, n=attrtable$POP10, type='random'), col='red', pch=3, cex=0.5)


#sa <- sapply_pb(slot(shapedata, 'polygons'), function(i) spsample(i, n=10,type='random'))

#sa <- sapply_pb(attrtable$POP10, function(i) spsample(shapedata[i], n=i,type='random'))

 sa <- dotsInPolys(shapedata, attrtable$POP10, f = "random",compatible = FALSE)

# stack into a single SpatialPoints object
#s.merged <- do.call('rbind', sa)

#ids <- sapply(slot(shapedata, 'polygons'), function(i) slot(i, 'ID'))

# determine the number of ACTUAL sample points generated for each poly
#npts <- sapply(sa, function(i) nrow(i@coords))

# generate a reconstituted vector of point IDs
#pt_id <- rep(ids, npts)

# promote to SpatialPointsDataFrame
#s.final <- SpatialPointsDataFrame(s.merged, data=data.frame(poly_id=pt_id))

# check:
#plot(shapedata) ; points(s.final, col=s.final$poly_id, pch=3, cex=0.5)

# copy source data spatial reference system to new object
#s.final@proj4string <- shapedata@proj4string

# write out to new file

writeOGR(sa, ".", "random_pts", driver="ESRI Shapefile")










}