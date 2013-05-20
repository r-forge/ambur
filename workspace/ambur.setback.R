ambur.setback <-
function(multiplier=30, rate.field="EPR") {

require(tcltk)
 require(rgdal)
require(rgeos)


#open ambur plotting file
#choose dbf file to import
tkmessageBox(message = "Please select an points shapefile with ambur.analysis results in the attribute table, or, an ambur.statshape generated points shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
mydata <- data.frame(shapedata)


workingdir <- dirname(getdata)
setwd(workingdir)
path <- getdata



time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

#dir.create("AMBUR_setback", showWarnings=FALSE)
#setwd("AMBUR_setback")

if(is.character(rate.field) == TRUE) field <- mydata[,rate.field] else field <- rate.field



multiplier <- multiplier
offshore.cor <- ifelse(mydata[,"Base_Off"] == 1, -1,1)



field.predx <- sin((mydata[,"Azimuth"] * pi/180)) * (field * multiplier * offshore.cor) + mydata[,"POINT_X"]

field.predy <- cos((mydata[,"Azimuth"]  * pi/180)) * (field * multiplier * offshore.cor) + mydata[,"POINT_Y"]





##############################################################################
colnames(mydata) <- gsub("Base_Loc", "BASE_LOC", colnames(mydata))
mydata[,"BASE_LOC"] <- ifelse(is.na(mydata[,"BASE_LOC"]) ==TRUE, 1, mydata[,"BASE_LOC"])

mydata$SetbackX <- field.predx
mydata$SetbackY <- field.predy


ddTable <- data.frame(mydata) 
coordinates(ddTable) <- data.frame(x=(mydata["SetbackX"]), y=(mydata["SetbackY"]))
outputname_a <- paste(shapename,"_setback_points",sep="")
writeOGR(ddTable, ".", outputname_a, driver="ESRI Shapefile")


##field
Baseline.Factor <- factor(mydata[,"BASE_LOC"])
field.final <- sapply(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(field.predx[mydata[,"BASE_LOC"] == x]), y=c(field.predy[mydata[,"BASE_LOC"] == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
field.final2 <- SpatialLines(field.final)

field.tab <- data.frame(Multiplier=multiplier,RateMeth="field",Source=getdata,Creator="R - AMBUR")

field.tab2 <-  field.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(field.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

field.final3 <- SpatialLinesDataFrame(field.final2, field.tab2)
#create shapefile and write it to the working directory

 # Note that readOGR method reads the .prj file when it exists
  projectionString <- proj4string(shapedata) # contains projection info
  
 proj4string(field.final3) <- projectionString

 outputname <- paste(shapename,"_setback_lines",sep="")

writeOGR(field.final3, ".", outputname, driver="ESRI Shapefile")







}

