merge.noaa.lidar.xyz.txt.to.shp <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1



require(rgdal)
require(rgeos)
require(tcltk)


tkmessageBox(message = "Please select the directory with the ascii txt (comma separated) that files...")
getdata <- tk_choose.dir(default = "", caption = "Select directory")
workingdir <- getdata
setwd(workingdir)
############################################################### setup a progress bar

lapply_pb <- function(X, FUN, ...)
{
 env <- environment()
 pb_Total <- length(X)
 counter <- 0
 pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   

 # wrapper around FUN
 wrapper <- function(...){
   curVal <- get("counter", envir = env)
   assign("counter", curVal +1 ,envir=env)
   setTxtProgressBar(get("pb", envir=env), curVal +1)
   FUN(...)
 }
 res <- lapply(X, wrapper, ...)
 close(pb)
 res
}


###############################################################





filenames <- list.files(path = getdata)
merge.tab <- do.call("rbind", lapply_pb(filenames, read.csv, header = TRUE))

colnames(merge.tab) <- c("Easting","Northing","Elevation")


coordinates(merge.tab) <- data.frame(x=(merge.tab[,1]), y=(merge.tab[,2]))
 #projectionString <- proj4string(shapedata) # contains projection info
  #proj4string(merge.tab) <- projectionString
writeOGR(merge.tab, ".", "xyz_pts", driver="ESRI Shapefile")


#write.table(merge.tab, file = "merged_csv_data.csv", sep = ",", row.names = FALSE)


}