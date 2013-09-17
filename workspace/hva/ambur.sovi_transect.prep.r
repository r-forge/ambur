ambur.sovi_transect.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the SoVI shapefile...")
filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata <- tk_choose.files(filter = filters,multi = FALSE)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata <- readOGR(getdata,layer=shapename)
attrtable <- data.frame(shapedata)


tkmessageBox(message = "Please select a post-AMBUR analysis transect shapefile generated from ambur.statshape function...")
getdata2 <- tk_choose.files(filter = filters,multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))

shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)

 ########
time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
setwd("AMBUR_HVA_data_prep")

#dir.create(paste(time.stamp2," ","shorepts",sep=""))
#setwd(paste(time.stamp2," ","shorepts",sep=""))

####create spatialLines only

#shapedataA <- as(shapedata, "SpatialLines")
#shapedataB <- as(shapedata2, "SpatialLines")

int <- gIntersects(shapedata2, shapedata, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

pb <- tkProgressBar("AMBUR: progress bar", "This might take a moment...", 0, max(length(seq(along=vec))), 50)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb, i, sprintf("AMBUR: Capture shoreline positions (%s)", info), info)

vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(shapedata2[i,], shapedata[int[,i],], byid=TRUE)  else 0 }






#cond <- lapply(vec, function(x) length(x) > 0)

cond <- lapply(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]

out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))

#testout <- data.frame(rn,SpatialLinesLengths(out))
#rownames(testout) <- rn
#outaa <- SpatialLinesDataFrame(out,testout)
#writeOGR(outaa, ".", "ambur_test", driver="ESRI Shapefile") 



Pcnt.Complete <-  50                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)



transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
#POINT_X <-  data.frame(coordinates(out))$x
#POINT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(0,length(rn)-1,1)

inter.data <- data.frame(transID,shoreID,sortID)
tran.data <- data.frame(shapedata2@data)
shore.data <- data.frame(shapedata@data)


tran.data$Id <- as.numeric(row.names(tran.data))

shore.data$mergeID <- as.numeric(row.names(shore.data))


tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)

tet3 <- tet2[ order(tet2[,"sortID"]) , ]

tet3$Id <- tet2[,"sortID"]

tet3$Distance <- SpatialLinesLengths(out)  #ambur.statshape generated transects puts underscores under Start_X and Start_Y to distinguish between pre analysis transects

#tet3$BASE_LOC <- ifelse(tet3$BASE_LOC == "NA", tet3$BaseOrder, tet3$BASE_LOC)

row.names(tet3) <- rn

outputdata <- SpatialLinesDataFrame(out,tet3)

 writeOGR(outputdata, ".", "ambur_sovi_lines_all", driver="ESRI Shapefile") 

#####################cull multiple intersections to get max length of transect line that is not water

mydata <- data.frame(outputdata)

mydata$Id2  <- seq(1,length(mydata[,"Transect"]),1)   #sets the Id field equal to that of the shape IDs in the spatial lines file

test <- data.frame(mydata[ ,"Transect"] , mydata[ ,"Distance"],mydata[ ,"Id2"])[ order(mydata[ ,"Transect"], mydata[ ,"Distance"]),]   #sort to get proper order

colnames(test) <- c("Transect", "Distance" ,"Id2")


testblah <- data.frame(test)

testblah2 <- aggregate(. ~ testblah$Transect,data = testblah[,2:3],FUN=function(x) c(mn =min(x), n=x[1] ) )  #get first occurence of the min distance ID

colnames(testblah2) <- c("TRANSECT", "DISTANCE" ,"ID2")

testblah3 <- testblah2[,3]

 valid.pts <-  as.numeric(testblah3[,1])

 shapedatacull <- outputdata[valid.pts,]
 shapetablecull <-  data.frame(shapedatacull)
 
 ######HVA index setup
hva_cat <- numeric(length(shapetablecull[,"Distance"] ))

hva_cat[shapetablecull[,"SoVI"] > 5] <- 5
hva_cat[shapetablecull[,"SoVI"] > 1 & shapetablecull[,"SoVI"] < 5] <- 4
hva_cat[shapetablecull[,"SoVI"] > -1 & shapetablecull[,"SoVI"] < 1] <- 3
hva_cat[shapetablecull[,"SoVI"] > -5 & shapetablecull[,"SoVI"] < -1] <- 2
hva_cat[shapetablecull[,"SoVI"] < -5] <- 1

shapetablecull$hva_cat <- hva_cat          


sovi_hva_table <- cbind(shapetablecull$Transect, hva_cat)

colnames(sovi_hva_table) <- c("Transect","hva_cat")

write.table(sovi_hva_table, file = "ambur_hva_sovi_data.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)

  plot(outputdata)
 plot(shapedatacull)

lns.output <- SpatialLinesDataFrame(shapedatacull,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(lns.output) <- projectionString

writeOGR(lns.output, ".", "ambur_hva_sovi_lns", driver="ESRI Shapefile")

 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)































}