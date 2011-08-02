ambur.filtertran <-
function(winsize=5, indv=1) {

filterwinsize <- (winsize-1)/2
indv.base <- indv


require(shapefiles)
require(spatstat)

winDialog("ok","Please select the transects...")

path1 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path1)

path1a <- paste(substr(path1,1,max(del.ext)-4),".shp",sep="")

my.shapefile <- read.shp(path1a)

mydata <- convert.to.simple(my.shapefile)

path2 <- dirname(path1)
setwd(path2)

transects.dbf <- data.frame(read.dbf(paste(substr(path1,1,max(del.ext)-4),".dbf",sep="")))

colnames(transects.dbf) <- gsub("ID", "Id", colnames(transects.dbf))

transects.dbf$dbf.Id <- seq(1,length(transects.dbf$dbf.Id),by=1)


winDialog("ok","Please select the inner baseline...")
path3 <- choose.files(default = "*.shp",multi = FALSE)

del.ext <- nchar(path3)

path3a <- paste(substr(path3,1,max(del.ext)-4),".shp",sep="")

my.shapefile2 <- read.shp(path3a)

inner_base <- convert.to.simple(my.shapefile2)





time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)



dir.create(paste(time.stamp2," ","filtered",sep=""))
setwd(paste(time.stamp2," ","filtered",sep=""))

#set up master data table

trandata <- data.frame(transects.dbf)
colnames(trandata) <- gsub("dbf.", "", colnames(trandata))



#######################
n.adjrec <- 0

in.az <- trandata$Azimuth



#cast filter transects using recommended filter window size or user specified window for unclosed baselines
n.adj <- ifelse(filterwinsize > 0, filterwinsize, n.adjrec)

filter.az <- numeric(length(in.az))

for (i in 1:length(in.az)) {


if(i > n.adj && i < (length(in.az)-n.adj)) filter.window <- c(in.az[(i-1):(i-n.adj)],in.az[i],in.az[(i+1):(i+n.adj)]) else filter.window <- 0



filter.span <- sum(as.logical(filter.window[(filter.window > 90) && (filter.window < 270)]))*0   #find if there are transects radiating > 180 degrees

if(min(filter.window,na.rm = TRUE) < 90 && max(filter.window,na.rm = TRUE) > 270 && i > n.adj && i < (length(in.az)-n.adj) && filter.span ==0) filter.prep <- c(filter.window[filter.window >270]-360,filter.window[filter.window <= 270])  else  filter.prep <- filter.window



filter.prep2 <- ifelse(i > n.adj && i < (length(in.az)-n.adj),mean(filter.prep),in.az[i])



filter.az[i] <- ifelse(filter.prep2 < 0, filter.prep2 + 360, filter.prep2)

}

#########################################################################test to get closed baseline averages for windows at the end and beginning of the baseline

for (i in 1:length(in.az)) {

if(i <= n.adj) filter.window <- c(in.az[((length(in.az))-(n.adj-(i-1))+1):(length(in.az))],in.az[1:i],in.az[(i+1):(i+n.adj)]) else filter.window <- filter.az[i]

if(i >= (length(in.az)-n.adj)) filter.window <- c(in.az[(i-1):(i-n.adj)],in.az[i:(length(in.az))],in.az[1:(n.adj-(length(in.az)-i))]) else filter.window <- filter.window * 1


filter.span <- sum(as.logical(filter.window[(filter.window > 90) && (filter.window < 270)]))*0   #find if there are transects radiating > 180 degrees


if(min(filter.window,na.rm = TRUE) < 90 && max(filter.window,na.rm = TRUE) > 270 && i <= n.adj | i >= (length(in.az)-n.adj) && filter.span ==0) filter.prep <- c(filter.window[filter.window >270]-360,filter.window[filter.window <= 270])  else  filter.prep <- filter.az[i]

filter.prep2 <- mean(filter.prep)



filter.az[i] <- ifelse(filter.prep2 < 0, filter.prep2 + 360, filter.prep2)
}
#####################################################################################



filter.xy.az <- filter.az
filter.x <- sin((filter.az * pi/180)) * (trandata$TranDist) + trandata$StartX
filter.y <- cos((filter.az * pi/180)) * (trandata$TranDist) + trandata$StartY



##############experimental

az.cos <- cos(in.az * (pi/180))
az.sin <- sin(in.az * (pi/180))


n.adj <- ifelse(filterwinsize > 0, filterwinsize, n.adjrec)

filter.az2 <- numeric(length(in.az))
avg.cos <- numeric(length(in.az))
avg.sin <- numeric(length(in.az))
avg.val <- numeric(length(in.az))
avg.val.base <- numeric(length(in.az))
base.test <- numeric(length(in.az))
base.az  <- numeric(length(in.az))

for (i in 1:length(in.az)) {


if(i > n.adj && i < (length(az.cos)-n.adj)) filter.window1 <- c(az.cos[(i-1):(i-n.adj)],az.cos[i],az.cos[(i+1):(i+n.adj)]) else filter.window1 <- 0

if(i > n.adj && i < (length(az.sin)-n.adj)) filter.window2 <- c(az.sin[(i-1):(i-n.adj)],az.sin[i],az.sin[(i+1):(i+n.adj)]) else filter.window2 <- 0


if(i > n.adj && i < (length(az.cos)-n.adj)) filter.window3 <- c(trandata$BaseOrder[(i-1):(i-n.adj)],trandata$BaseOrder[i],trandata$BaseOrder[(i+1):(i+n.adj)]) else filter.window3 <- c(0,1)  #added to select for individual baselines

base.test[i] <- as.numeric(min(filter.window3*indv.base,na.rm=TRUE) == max(filter.window3*indv.base,na.rm=TRUE)) #added to select for individual baselines



avg.cos[i] <- ifelse(i > n.adj && i < (length(in.az)-n.adj),mean(filter.window1,na.rm=TRUE),az.cos[i])

avg.sin[i] <- ifelse(i > n.adj && i < (length(in.az)-n.adj) ,mean(filter.window2,na.rm=TRUE),az.sin[i])


avg.val[i] <- (atan2(avg.sin[i],avg.cos[i]))* (180/pi)
avg.val.base[i] <- (atan2(az.sin[i],az.cos[i]))* (180/pi) #added to select for individual baselines

base.az[i]  <- ifelse(base.test[i] > 0,avg.val[i],avg.val.base[i]) #added to select for individual baselines

filter.az2[i] <- ifelse(base.az[i] < 0, base.az[i] + 360, base.az[i])

#filter.az2[i] <- ifelse(avg.val[i] < 0, avg.val[i] + 360, avg.val[i])

}



filter2.x <- sin((filter.az2 * pi/180)) * (trandata$TranDist*2) + trandata$StartX
filter2.y <- cos((filter.az2 * pi/180)) * (trandata$TranDist*2) + trandata$StartY


########trim transects from 2nd filter
test.wx <- c(trandata$StartX,filter2.x)
test.wy <- c(trandata$StartY,filter2.y)

Cx <- inner_base$X
Cy <- inner_base$Y
Cx2 <- c(Cx[-1],Cx[length(Cx)])
Cy2 <- c(Cy[-1],Cy[length(Cy)])

Test.w <- owin()
Test.w <- owin(c(min(test.wx),max(test.wx)), c(min(test.wy),max(test.wy)))

TY.w <- owin()
TY.w <- owin(c(min(Cx),max(Cx)), c(min(Cy),max(Cy)))
TY <- psp(Cx,Cy,Cx2,Cy2,window=Test.w)

trim.x <- numeric(length(trandata$Transect))
trim.y  <- numeric(length(trandata$Transect))
trim.length  <- numeric(length(trandata$Transect))

for (i in 1:length(trandata$Transect)) {
 
b.x <- psp(trandata$StartX[i], trandata$StartY[i], filter2.x[i], filter2.y[i], window=Test.w)

 inner.trim <- crossing.psp(b.x,TY)

pts.dists <- ((inner.trim$x - trandata$StartX[i])^2 +  (inner.trim$y - trandata$StartY[i])^2)^(1/2)

 trim.x[i] <- ifelse(inner.trim$n == 0, filter2.x[i],inner.trim$x[which.min(pts.dists)] )
 trim.y[i] <- ifelse(inner.trim$n == 0, filter2.y[i],inner.trim$y[which.min(pts.dists)] )
 trim.length[i] <- ifelse(inner.trim$n == 0,trandata$TranDist,min(pts.dists))
 } 












########

plot(trandata$StartX,trandata$StartY,col="white",asp=1,xlab="X",ylab="Y")
segments(trandata$StartX,trandata$StartY,trim.x,trim.y,col="blue")

b <- max(trandata$BaseOrder,na.rm=TRUE)


filterdata <- trandata


#filterdata$EndX <- filter.x
#filterdata$EndY <- filter.y
#filterdata$Azimuth <- filter.az

#id.field <- filterdata$Id
#dd <- data.frame(Id=c(id.field,id.field),X=c(filterdata$StartX,filterdata$EndX),Y=c(filterdata$StartY,filterdata$EndY))
#ddTable <- data.frame(filterdata)
#ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
#write.shapefile(ddShapefile, paste("b",b,"transects_near_filter",sep=""), arcgis=T)




##write 2nd filter shapefile
#filterdata$EndX <- filter2.x
#filterdata$EndY <- filter2.y
#filterdata$Azimuth <- filter.az2

#id.field <- filterdata$Id
#dd <- data.frame(Id=c(id.field,id.field),X=c(filterdata$StartX,filterdata$EndX),Y=c(filterdata$StartY,filterdata$EndY))
#ddTable <- data.frame(filterdata)
#ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
#write.shapefile(ddShapefile, paste("b",b,"transects_near_filter2",sep=""), arcgis=T) 

##write 2nd filter trim shapefile
filterdata$EndX <- trim.x
filterdata$EndY <- trim.y
filterdata$TranDist <- trim.length



id.field <- filterdata$Id
dd <- data.frame(Id=c(id.field,id.field),X=c(filterdata$StartX,filterdata$EndX),Y=c(filterdata$StartY,filterdata$EndY))
ddTable <- data.frame(filterdata)
ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 3)
write.shapefile(ddShapefile, paste("b",b,"transects__filter",sep=""), arcgis=T) 

detach("package:shapefiles")
detach("package:spatstat")
}

