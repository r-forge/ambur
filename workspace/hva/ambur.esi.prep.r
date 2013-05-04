ambur.esi.prep <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1


require(rgdal)
require(rgeos)
 require(tcltk)

tkmessageBox(message = "Please select the ESI shoreline shapefile...")
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


###########

####create spatialLines only

shapedataA <- as(shapedata, "SpatialLines")
shapedataB <- as(shapedata2, "SpatialLines")

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


Pcnt.Complete <-  50                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
POINT_X <-  data.frame(coordinates(out))$x
POINT_Y <-  data.frame(coordinates(out))$y


sortID <- seq(0,length(POINT_X)-1,1)

inter.data <- data.frame(POINT_X,POINT_Y,transID,shoreID,sortID)
tran.data <- data.frame(shapedata2@data)
shore.data <- data.frame(shapedata@data)


tran.data$Id <- as.numeric(row.names(tran.data))

shore.data$mergeID <- as.numeric(row.names(shore.data))


tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)

tet3 <- tet2[ order(tet2[,"sortID"]) , ]

tet3$Id <- tet2[,"sortID"]

tet3$Distance <- (((tet3$POINT_X - tet3$Start_X)^2 +  (tet3$Start_Y - tet3$POINT_Y)^2)^(1/2))   #ambur.statshape generated transects puts underscores under Start_X and Start_Y to distinguish between pre analysis transects

#tet3$BASE_LOC <- ifelse(tet3$BASE_LOC == "NA", tet3$BaseOrder, tet3$BASE_LOC)

 out2 <- SpatialPoints(cbind(tet3$POINT_X,tet3$POINT_Y)) #recalculate point coordinates in case they become out of sequence
outputdata <- SpatialPointsDataFrame(out2,tet3)



Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Extracting first intersection points (%s)", info), info)


 # Note that readOGR method reads the .prj file when it exists
   projectionString <- proj4string(shapedata) # contains projection info

  proj4string(outputdata) <- projectionString

# writeOGR(outputdata, ".", "ambur_hva_esi_all_pts", driver="ESRI Shapefile")    #write all intersection points to a shapefile





#####################cull multiple intersections to get first intersection from transect origin Start_X

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
ESI_values <- shapetablecull$MOSTSENSIT

ESI_values <- gsub("A", "", ESI_values)
ESI_values <- gsub("B", "", ESI_values)
ESI_values <- gsub("C", "", ESI_values)
ESI_values <- gsub("D", "", ESI_values)
ESI_values <- gsub("[+]", "", ESI_values)
ESI_values <- as.numeric(ESI_values)

ESIcat1 <- c(1,2)
ESIcat2 <- c(3,4)
ESIcat3 <- c(5,6)
ESIcat4 <- c(7,8,9)
ESIcat5 <- c(10)

hva_cat <- numeric(length(ESI_values))

hva_cat[ESI_values %in% ESIcat1] <- 1
hva_cat[ESI_values %in% ESIcat2] <- 2
hva_cat[ESI_values %in% ESIcat3] <- 3
hva_cat[ESI_values %in% ESIcat4] <- 4
hva_cat[ESI_values %in% ESIcat5] <- 5

shapetablecull$hva_cat <- hva_cat          

esi_hva_table <- cbind(shapetablecull$Transect, hva_cat)

colnames(esi_hva_table) <- c("Transect","hva_cat")

write.table(esi_hva_table, file = "ambur_hva_esi_data.csv", sep = ",", row.names = FALSE)


###########


 Pcnt.Complete <-  90
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)

  plot(outputdata)
 plot(shapedatacull)

 pts.output <- SpatialPointsDataFrame(shapedatacull,shapetablecull)

    projectionString <- proj4string(shapedata) # contains projection info

  proj4string(pts.output) <- projectionString

writeOGR(pts.output, ".", "ambur_hva_esi_pts", driver="ESRI Shapefile")

 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i * 1 , sprintf("AMBUR: Writing new points to shapefile (%s)", info), info)




}