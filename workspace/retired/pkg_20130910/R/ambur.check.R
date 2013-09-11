ambur.check <-
function(checktype=1) {
require(tcltk)
library(foreign)

path <- tk_choose.files(default = "*.dbf",multi = FALSE)
workingdir <- dirname(path)
setwd(workingdir)


mydata <- foreign::read.dbf(path)

checktype <- checktype

if (checktype == 1) shapetype <- "shorelines"
if (checktype == 2) shapetype <- "transects"
if (checktype == 3) shapetype <- "capture points"

cat("Checking",shapetype,"shapefile table...","\n",sep=" ")


if (checktype == 1) reqfields <- c("Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP")

if (checktype == 2) reqfields <- c("Id","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth")
 
if (checktype == 3) reqfields <- c("Id","DATE_","ACCURACY","SHORE_LOC","CLASS_1","CLASS_2","CLASS_3","GROUP","Transect","TranSpace","TranDist","Location","MaxBNum","BaseOrder","OFFshore","CastDir","BASE_LOC","StartX","StartY","EndX","EndY","Azimuth","Point_X","Point_Y","Distance")

#check for missing fields
fieldcheck <- toupper(colnames(mydata))

 presentfields <-  which(toupper(reqfields) %in% fieldcheck)

   missingfields <- as.character(reqfields[-presentfields])



cat("The following fields are missing:","\n")
cat(missingfields, "\n",sep = ", ")


#detach("package:foreign")
}

