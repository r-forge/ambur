ambur.forecast <-
function(years=50) {

#require(tcltk)
 #require(rgdal)



tkmessageBox(message = "Please select a post-AMBUR analysis transect shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

data.path <- getdata 

mydata <- attrtable



time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_forecast", showWarnings=FALSE)
setwd("AMBUR_forecast")

dir.create(paste(time.stamp2," ","forecast",sep=""))
setwd(paste(time.stamp2," ","forecast",sep=""))
mydata$Baseline.Location <- mydata$Base_Loc
mydata$Baseline.Offshore <- mydata$Base_Off
mydata$Baseline.Order <- as.character(mydata$BaseOrder)
mydata$Max.Date.Xcoord <-  mydata$Max_DateX
mydata$Max.Date.Ycoord <-  mydata$Max_DateY 
mydata$Transect.Azimuth <- mydata$T_azimuth

forecast.years <- years
offshore.cor <- ifelse(mydata$Baseline.Offshore == 1, -1,1)



EPR.predx <- ifelse(is.na(mydata$EPR) == FALSE,sin((mydata$Transect.Azimuth * pi/180)) * (mydata$EPR * forecast.years * offshore.cor) + mydata$Max.Date.Xcoord, mydata$Max.Date.Xcoord)
EPR.predy <- ifelse(is.na(mydata$EPR) == FALSE,cos((mydata$Transect.Azimuth * pi/180)) * (mydata$EPR * forecast.years * offshore.cor) + mydata$Max.Date.Ycoord, mydata$Max.Date.Ycoord)

LRR.predx <- ifelse(is.na(mydata$LRR) == FALSE,sin((mydata$Transect.Azimuth * pi/180)) * (mydata$LRR * forecast.years * offshore.cor) + mydata$Max.Date.Xcoord, mydata$Max.Date.Xcoord)
LRR.predy <- ifelse(is.na(mydata$LRR) == FALSE,cos((mydata$Transect.Azimuth * pi/180)) * (mydata$LRR * forecast.years * offshore.cor) + mydata$Max.Date.Ycoord, mydata$Max.Date.Ycoord)

WLR.predx <- ifelse(is.na(mydata$WLR) == FALSE,sin((mydata$Transect.Azimuth * pi/180)) * (mydata$WLR * forecast.years * offshore.cor) + mydata$Max.Date.Xcoord, mydata$Max.Date.Xcoord)
WLR.predy <- ifelse(is.na(mydata$WLR) == FALSE,cos((mydata$Transect.Azimuth * pi/180)) * (mydata$WLR * forecast.years * offshore.cor) + mydata$Max.Date.Ycoord, mydata$Max.Date.Ycoord)

#RLR.slope
#LMS.slope
#JK.avg
#JK.min
#JK.max

#plot(EPR.predx,EPR.predy,type="l",asp="1",col="white",main="Shoreline Forecast",xlab="X",ylab="Y")
#lines(mydata$Max.Date.Xcoord,mydata$Max.Date.Ycoord,col="gray")
#lines(EPR.predx,EPR.predy,col="green")
#lines(LRR.predx,LRR.predy,col="blue")
#lines(WLR.predx,WLR.predy,col="red")



##############################################################################

mydata$Baseline.Location <- ifelse(is.na(mydata$Baseline.Location) ==TRUE, mydata$Baseline.Order, mydata$Baseline.Location)


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






cat("writing EPR prediction shapefile")

##EPR
Baseline.Factor <- factor(mydata$Baseline.Location)
EPR.final <- sapply_pb(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(EPR.predx[mydata$Baseline.Location == x]), y=c(EPR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
EPR.final2 <- SpatialLines(EPR.final)

EPR.tab <- data.frame(nYears=forecast.years,RateMeth="EPR",Source=data.path,Creator="R - AMBUR")

EPR.tab2 <-  EPR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(EPR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

EPR.final3 <- SpatialLinesDataFrame(EPR.final2, EPR.tab2)
#create shapefile and write it to the working directory
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(EPR.final3) <- projectionString

writeOGR(EPR.final3, ".", "EPR_forecast", driver="ESRI Shapefile")


cat("writing LRR prediction shapefile")
##LRR
Baseline.Factor <- factor(mydata$Baseline.Location)
LRR.final <- sapply_pb(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(LRR.predx[mydata$Baseline.Location == x]), y=c(LRR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
LRR.final2 <- SpatialLines(LRR.final)

LRR.tab <- data.frame(nYears=forecast.years,RateMeth="LRR",Source=data.path,Creator="R - AMBUR")

LRR.tab2 <-  LRR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(LRR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

LRR.final3 <- SpatialLinesDataFrame(LRR.final2, LRR.tab2)
#create shapefile and write it to the working directory
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(LRR.final3) <- projectionString

writeOGR(LRR.final3, ".", "LRR_forecast", driver="ESRI Shapefile")


cat("writing WLR prediction shapefile")
##WLR
Baseline.Factor <- factor(mydata$Baseline.Location)
WLR.final <- sapply_pb(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(WLR.predx[mydata$Baseline.Location == x]), y=c(WLR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
WLR.final2 <- SpatialLines(WLR.final)

WLR.tab <- data.frame(nYears=forecast.years,RateMeth="WLR",Source=data.path,Creator="R - AMBUR")

WLR.tab2 <-  WLR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(WLR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

WLR.final3 <- SpatialLinesDataFrame(WLR.final2, WLR.tab2)
#create shapefile and write it to the working directory
 projectionString <- proj4string(shapedata) # contains projection info
  proj4string(WLR.final3) <- projectionString

writeOGR(WLR.final3, ".", "WLR_forecast", driver="ESRI Shapefile")




#detach(mydata)

}

