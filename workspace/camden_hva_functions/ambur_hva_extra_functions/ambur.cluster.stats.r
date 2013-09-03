ambur.adjacent.stats <-
function(d=100) {


# Establish the inputs
d <- d


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the Max Date Points shapefile...")
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
 pb <- tkProgressBar("AMBUR: progress bar", "Creating near transects...", 0, 100, 1)

Transect.Factor <- factor(shapedata$Transect)

setTkProgressBar(pb, 10 , "AMBUR: progress bar", "Calculating transect locations...")

#chc <- sapply_pb(levels(Transect.Factor), function(x) cutree(hclust(dist(data.frame(rownames=rownames(shapedata@data), x=coordinates(shapedata)[,1][shapedata$Transect == x], y=coordinates(shapedata)[,2][shapedata$Transect == x])), method="complete"),h=d)   ,simplify = FALSE)

#trans.indv <- sapply_pb(levels(Transect.Factor), function(x) data.frame("BaseOrder"=unique(shapedata$BaseOrder[shapedata$BaseOrder == x]),coordinates(spsample(shapedata[shapedata$BaseOrder == x,], round(sum(SpatialLinesLengths(shapedata[shapedata$BaseOrder == x,]))/sampledist,0), offset=0.000000, "regular") )) ,simplify = FALSE)




chc <- hclust(dist(data.frame(rownames=rownames(shapedata@data), x=coordinates(shapedata)[,1], y=coordinates(shapedata)[,2])), method="complete")
              
              
# Distance with a m threshold  
chc.d <- cutree(chc, h=d)               
 
 
shapedata@data <- data.frame(shapedata@data, Clust=chc.d) 


# Plot results
plot(shapedata, col=factor(shapedata@data$Clust), pch=19)
     box(col="black")
       title(main="Clustering")
        legend("topleft", legend=paste("Cluster", 1:4,sep=""),
                   col=palette()[1:4], pch=rep(19,4), bg="white")
                   
                   
                   
     
              
