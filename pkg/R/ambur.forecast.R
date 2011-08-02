ambur.forecast <-
function(years=50) {



#open ambur plotting file

winDialog("ok","Please select the 'results_stats.csv' file...")

data.path <- choose.files(default = "*.csv",multi = FALSE)

mydata <- read.csv(data.path, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(data.path)
setwd(dir.path)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_forecast", showWarnings=FALSE)
setwd("AMBUR_forecast")

dir.create(paste(time.stamp2," ","forecast",sep=""))
setwd(paste(time.stamp2," ","forecast",sep=""))


forecast.years <- years
offshore.cor <- ifelse(Baseline.Offshore == 1, -1,1)



EPR.predx <- ifelse(is.na(EPR) == FALSE,sin((Transect.Azimuth * pi/180)) * (EPR * forecast.years * offshore.cor) + Max.Date.Xcoord, Max.Date.Xcoord)
EPR.predy <- ifelse(is.na(EPR) == FALSE,cos((Transect.Azimuth * pi/180)) * (EPR * forecast.years * offshore.cor) + Max.Date.Ycoord, Max.Date.Ycoord)

LRR.predx <- ifelse(is.na(LRR.slope) == FALSE,sin((Transect.Azimuth * pi/180)) * (LRR.slope * forecast.years * offshore.cor) + Max.Date.Xcoord, Max.Date.Xcoord)
LRR.predy <- ifelse(is.na(LRR.slope) == FALSE,cos((Transect.Azimuth * pi/180)) * (LRR.slope * forecast.years * offshore.cor) + Max.Date.Ycoord, Max.Date.Ycoord)

WLR.predx <- ifelse(is.na(WLR.slope) == FALSE,sin((Transect.Azimuth * pi/180)) * (WLR.slope * forecast.years * offshore.cor) + Max.Date.Xcoord, Max.Date.Xcoord)
WLR.predy <- ifelse(is.na(WLR.slope) == FALSE,cos((Transect.Azimuth * pi/180)) * (WLR.slope * forecast.years * offshore.cor) + Max.Date.Ycoord, Max.Date.Ycoord)

#RLR.slope
#LMS.slope
#JK.avg
#JK.min
#JK.max

plot(EPR.predx,EPR.predy,type="l",asp="1",col="white",main="Shoreline Forecast",xlab="X",ylab="Y")
lines(Max.Date.Xcoord,Max.Date.Ycoord,col="gray")
lines(EPR.predx,EPR.predy,col="green")
lines(LRR.predx,LRR.predy,col="blue")
lines(WLR.predx,WLR.predy,col="red")


require(shapefiles)

id.field <- 1
dd <- data.frame(Id=c(id.field),X=EPR.predx,Y=EPR.predy)
ddTable <- data.frame(Id=id.field,nYears=forecast.years,RateMeth="EPR",Source=data.path,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("EPR_forecast",sep=""), arcgis=T)


id.field <- 1
dd <- data.frame(Id=c(id.field),X=LRR.predx,Y=LRR.predy)
ddTable <- data.frame(Id=id.field,nYears=forecast.years,RateMeth="LRR",Source=data.path,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("LRR_forecast",sep=""), arcgis=T)


id.field <- 1
dd <- data.frame(Id=c(id.field),X=WLR.predx,Y=WLR.predy)
ddTable <- data.frame(Id=id.field,nYears=forecast.years,RateMeth="WLR",Source=data.path,Creator="R - AMBUR")
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("WLR_forecast",sep=""), arcgis=T)



detach("package:shapefiles")
detach(mydata)

}

