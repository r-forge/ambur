ambur.forecast <-
function(years=50) {

require(tcltk)
 require(rgdal)

years=50    #for testing

#open ambur plotting file
tkmessageBox(message = "Please select the 'results_stats.csv' file...")


data.path <- tk_choose.files(default = "*.csv",multi = FALSE)

mydata <- read.csv(data.path, header=TRUE, sep=",")
attach(mydata)

dir.path <- dirname(data.path)
setwd(dir.path)

time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_forecast", showWarnings=FALSE)
setwd("AMBUR_forecast")

#dir.create(paste(time.stamp2," ","forecast",sep=""))
#setwd(paste(time.stamp2," ","forecast",sep=""))


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



##############################################################################

mydata$Baseline.Location <- ifelse(is.na(mydata$Baseline.Location) ==TRUE, 1, mydata$Baseline.Location)


##EPR
Baseline.Factor <- factor(mydata$Baseline.Location)
EPR.final <- sapply(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(EPR.predx[mydata$Baseline.Location == x]), y=c(EPR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
EPR.final2 <- SpatialLines(EPR.final)

EPR.tab <- data.frame(nYears=forecast.years,RateMeth="EPR",Source=data.path,Creator="R - AMBUR")

EPR.tab2 <-  EPR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(EPR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

EPR.final3 <- SpatialLinesDataFrame(EPR.final2, EPR.tab2)
#create shapefile and write it to the working directory
writeOGR(EPR.final3, ".", "EPR_forecast", driver="ESRI Shapefile")



##LRR
Baseline.Factor <- factor(mydata$Baseline.Location)
LRR.final <- sapply(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(LRR.predx[mydata$Baseline.Location == x]), y=c(LRR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
LRR.final2 <- SpatialLines(LRR.final)

LRR.tab <- data.frame(nYears=forecast.years,RateMeth="LRR",Source=data.path,Creator="R - AMBUR")

LRR.tab2 <-  LRR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(LRR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

LRR.final3 <- SpatialLinesDataFrame(LRR.final2, LRR.tab2)
#create shapefile and write it to the working directory
writeOGR(LRR.final3, ".", "LRR_forecast", driver="ESRI Shapefile")


##WLR
Baseline.Factor <- factor(mydata$Baseline.Location)
WLR.final <- sapply(levels(Baseline.Factor), function(x)
list(Lines(list(Line(list(x=c(WLR.predx[mydata$Baseline.Location == x]), y=c(WLR.predy[mydata$Baseline.Location == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
WLR.final2 <- SpatialLines(WLR.final)

WLR.tab <- data.frame(nYears=forecast.years,RateMeth="WLR",Source=data.path,Creator="R - AMBUR")

WLR.tab2 <-  WLR.tab[rep(1, length(unique(Baseline.Factor))),]
row.names(WLR.tab2) <- seq(0, length(unique(Baseline.Factor))-1,1)

WLR.final3 <- SpatialLinesDataFrame(WLR.final2, WLR.tab2)
#create shapefile and write it to the working directory
writeOGR(WLR.final3, ".", "WLR_forecast", driver="ESRI Shapefile")




detach(mydata)

}

