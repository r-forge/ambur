ambur.filtertran <-
function(winsize=5, indv=1) {

#require(tcltk)
#require(rgdal)
#require(rgeos)


#winsize=5
#indv=1


tkmessageBox(message = "Please select the transects shapefile...")
filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- if(interactive()) tk_choose.files(filter = filetype)
shapename <- gsub(".shp", "", basename(getdata))
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)

workingdir <- dirname(getdata)
setwd(workingdir)


tkmessageBox(message = "Please select the inner baseline shapefile...")
getdata2 <- tk_choose.files("","Choose file",multi = FALSE,filetype,1)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)


dir.create(paste(time.stamp2," ","filtered",sep=""))
setwd(paste(time.stamp2," ","filtered",sep=""))


pb <- tkProgressBar("AMBUR: progress bar", "Filtering...", 0, 100, 1)



#set up master data table
trandata <- attrtable

#repair transects and IDs to make certain they are in sequential order

trandata$Transect <- seq(1,length(trandata$Transect),1)
trandata$Id <- seq(1,length(trandata$Transect),1)

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





##establish moving average function
move.avg <- function(x,n=winsize){filter(x,rep(1/n,n), sides=2)}


###established azimuth filter function
filter.azimuths <- function(az.data){
cos.az <- cos(az.data * pi/180)
sin.az <- sin(az.data * pi/180)
avg.cosx <- move.avg(cos.az)
avg.siny <- move.avg(sin.az)
filter.az.rads <- atan2(avg.siny, avg.cosx)
filter.az.deg <- ifelse( (filter.az.rads *(180/pi)) < 0, (filter.az.rads *(180/pi)) + 360, (filter.az.rads *(180/pi)))
filter.az.final <- ifelse(is.na(filter.az.deg) == TRUE, trandata$Azimuth, filter.az.deg)
return(filter.az.final)}


setTkProgressBar(pb, 10 , "AMBUR: progress bar", "Step 1 of 10: calculating new transect azimuths")


### filter across all baselines
filter.az.all <- filter.azimuths(trandata$Azimuth)

### filter across individual baselines
Baseline.Factor <- factor(trandata$BaseOrder)
filter.az.indv <- sapply_pb(levels(Baseline.Factor), function(x) if (length(trandata$Azimuth[trandata$BaseOrder == x]) > winsize) filter.azimuths(trandata$Azimuth[trandata$BaseOrder == x]) else trandata$Azimuth[trandata$BaseOrder == x] ,simplify = TRUE)
filter.az.indv <- unlist(filter.az.indv, use.names = FALSE)

### get the appropriate filter azimuth based on user option
filter.az.sel <- if(indv == 1) filter.az.indv else filter.az.all


### get filter transects ending XY coordinates
filter2.x <- sin((filter.az.sel * pi/180)) * (trandata$TranDist*2) + trandata$StartX
filter2.y <- cos((filter.az.sel * pi/180)) * (trandata$TranDist*2) + trandata$StartY





setTkProgressBar(pb, 20 , "AMBUR: progress bar", "Step 2 of 10: building new transect lines")
#Pcnt.Complete <-  25
#info <- sprintf("%d%% Building new transect lines...", Pcnt.Complete)
#setTkProgressBar(pb, 25 , sprintf("AMBUR: Filter transects (%s)", info), info)



### build spatial lines for transects to intersect the baseline
Transect.Factor <- factor(trandata$Transect)
shape.prep <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(trandata$StartX[trandata$Transect == x], filter2.x[trandata$Transect == x]), y=c(trandata$StartY[trandata$Transect == x],filter2.y[trandata$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.prep2 <- SpatialLines(shape.prep)
shape.prep3 <- SpatialLinesDataFrame(shape.prep2, trandata)



#Pcnt.Complete <-  50
#info <- sprintf("%d%% Projecting to inner baseline...", Pcnt.Complete)
#setTkProgressBar(pb, 50 , sprintf("AMBUR: Filter transects (%s)", info), info)
setTkProgressBar(pb, 30 , "AMBUR: progress bar", "Step 3 of 10: projecting to inner baseline")

####intersect
##########################################


int <- gIntersects(shapedata2, shape.prep3, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

pb2 <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(seq(along=vec))), 50)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb2, i, sprintf("AMBUR: Capture baseline positions (%s)", info), info)

vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(shapedata2[i,], shape.prep3[int[,i],])  else 0 }

cond <- lapply(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]

out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))

close(pb2)

##############################

transID <- data.frame(nrn)[,2]
baseID <- data.frame(nrn)[,1]
INT_X <-  data.frame(coordinates(out))$x
INT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(1,length(INT_X),1)
inter.data <- data.frame(INT_X,INT_Y,transID,baseID,sortID)

tran.data <- data.frame(shapedata@data)
tran.data$Id <- as.numeric(row.names(tran.data))

tet <- merge(tran.data,inter.data , by.x = "Id", by.y = "transID", sort=FALSE, all.x=TRUE)
tet2 <- data.frame(tet[ order(tet[,"Transect"]) , ])


###added to correct for multiple interestions with the baseline  (8-20-2011) start:
Transect.Factor <- factor(tran.data$Transect)
tet2dist <- (((tet2[,"INT_X"]- tet2[,"StartX"])^2 +  (tet2[,"INT_Y"] - tet2[,"StartY"])^2)^(1/2))

setTkProgressBar(pb, 40 , "AMBUR: progress bar", "Step 4 of 10: calculating distance table")
tet2disttab <- data.frame(sapply_pb(levels(Transect.Factor), function(x) min(tet2dist[tet2$Transect == x],na.rm=FALSE) ,simplify = TRUE))

setTkProgressBar(pb, 50 , "AMBUR: progress bar", "Step 5 of 10: calculating secondary distance table")
tet2disttab2 <- data.frame(sapply_pb(levels(Transect.Factor), function(x) tet2dist[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))

setTkProgressBar(pb, 60 , "AMBUR: progress bar", "Step 6 of 10: adjusting x coordinates")
tet2intx <- data.frame(sapply_pb(levels(Transect.Factor), function(x) tet2$INT_X[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))

setTkProgressBar(pb, 70 , "AMBUR: progress bar", "Step 7 of 10: adjusting y coordinates")
tet2inty <- data.frame(sapply_pb(levels(Transect.Factor), function(x) tet2$INT_Y[tet2$Transect == x][tet2dist[tet2$Transect == x]== min(tet2dist[tet2$Transect == x],na.rm=FALSE)] ,simplify = TRUE))

tet3 <- data.frame(tran.data$Transect,tet2disttab,tet2disttab,tet2intx,tet2inty)
colnames(tet3) <- c("Transect","MinDist1","MinDist_check","INT_X","INT_Y")
##############end


setTkProgressBar(pb, 80 , "AMBUR: progress bar", "Step 8 of 10: building final data table")
### make new attribute table with filtered values
new_trandata <-  trandata
new_trandata[,"Azimuth"] <- ifelse(is.na(tet3$INT_X) == TRUE, as.numeric(tran.data$Azimuth), as.numeric(filter.az.sel))
new_trandata[,"EndX"]<- ifelse(is.na(tet3$INT_X) == TRUE, as.numeric(tran.data$EndX), as.numeric(tet3$INT_X))
new_trandata[,"EndY"] <- ifelse(is.na(tet3$INT_X) == TRUE, as.numeric(tran.data$EndY), as.numeric(tet3$INT_Y))
new_trandata[,"TranDist"] <- (((new_trandata[,"EndX"]- new_trandata[,"StartX"])^2 +  (new_trandata[,"EndY"] - new_trandata[,"StartY"])^2)^(1/2))





#Pcnt.Complete <-  75
#info <- sprintf("%d%% Finalizing filtered transects ...", Pcnt.Complete)
#setTkProgressBar(pb, 75 , sprintf("AMBUR: Filter transects (%s)", info), info)

setTkProgressBar(pb, 90 , "AMBUR: progress bar", "Step 9 of 10: building final shapefile")
### build spatial lines for final filtered transects shapefile
Transect.Factor <- factor(new_trandata$Transect)
shape.final <- sapply_pb(levels(Transect.Factor), function(x)
list(Lines(list(Line(list(x=c(new_trandata$StartX[new_trandata$Transect == x], new_trandata$EndX[new_trandata$Transect == x]), y=c(new_trandata$StartY[new_trandata$Transect == x],new_trandata$EndY[new_trandata$Transect == x])))), ID=(as.numeric(x)-1)))
,simplify = TRUE)
shape.final2 <- SpatialLines(shape.final)
shape.final3 <- SpatialLinesDataFrame(shape.final2, new_trandata)

   # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata2) # contains projection info
  
  proj4string(shape.final3) <- projectionString


#Pcnt.Complete <-  90
#info <- sprintf("%d%% Creating shapefile ...", Pcnt.Complete)
#setTkProgressBar(pb, 90 , sprintf("AMBUR: Filter transects (%s)", info), info)
setTkProgressBar(pb, 95 , "AMBUR: progress bar", "Step 10 of 10: writing final shapefile")
#create shapefile and write it to the working directory
writeOGR(shape.final3, ".", "filtered_transects", driver="ESRI Shapefile")
setTkProgressBar(pb, 100 , "AMBUR: progress bar", "finished!")



####plots
plot(c(trandata$StartX,trandata$EndX),c(trandata$StartY,trandata$EndY),col="white",asp=1,xlab="X",ylab="Y")
segments(trandata$StartX,trandata$StartY,trandata$EndX,trandata$EndY,col="gray")
segments(new_trandata$StartX,new_trandata$StartY,new_trandata$EndX,new_trandata$EndY,col="blue")

Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, 100 , sprintf("AMBUR: Filter transects (%s)", info), info)



##qc check on original azimuths
#orig.end.x <- sin((trandata$Azimuth * pi/180)) * (trandata$TranDist) + trandata$StartX
#orig.end.y <- cos((trandata$Azimuth * pi/180)) * (trandata$TranDist) + trandata$StartY
#round(orig.end.x - trandata$EndX,0)
#round(orig.end.y - trandata$EndY,0)

##for testing purposes:
#trandata$BaseOrder[300:max(length(trandata$BaseOrder))] <- 1
#trandata$BaseOrder
#filter.az.indv
#filter.az.all - filter.az.sel
#sum(filter.az.all - filter.az.indv)
#length(filter.az.all)
#length(filter.az.indv)
#indv = 2
#edit(trandata)


}