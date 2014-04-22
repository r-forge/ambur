ambur.hva <-
function(loc="location",bufferdist=150,bufferdist2=10,removeanthro="no",file1="path",file2="path",file3="path",file4="path",dir1="path",scr=c(-1,-0.2,0.2,1),tscr=c(0.5,0.2),sscr=c(0.16,0.074),sovi=c(6.77,4.06,-1.37,-4.09)) {

# establish required packages

#require(rgdal)
#require(rgeos)
#require(tcltk)

#set variables
loc <- loc
bufferdist <- bufferdist
bufferdist2 <- bufferdist2
removeanthro <- removeanthro
getdata <-  file1
getdata2 <- file2
getdata3 <- file3
getdata4 <- file4
dir1 <- dir1
scr <- scr
tscr <- tscr
sscr <- sscr
sovi <- sovi


# set up monitoring functions
lapply_pb <- function(X, FUN, ...)
{
 env <- environment()
 pb_Total <- length(X)
 counter <- 0
 pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   

 # wrapper around FUN
 wrapper <- function(...){
   curVal <- get("counter", envir = env)
   assign("counter", curVal +1 ,envir=env)
   setTxtProgressBar(get("pb", envir=env), curVal +1)
   FUN(...)
 }
 res <- lapply(X, wrapper, ...)
 close(pb)
 res
}



# establish the inputs

filetype <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)

getdata <- ifelse(file1 == "path",tk_choose.files("","Select a post-AMBUR analysis envelope transect shapefile",multi = FALSE,filetype,1 ),file1)
shapename <- gsub(".shp", "", basename(getdata))
workingdir <- dirname(getdata)
setwd(workingdir)
shapedata1 <- readOGR(getdata,layer=shapename)
attrtable1 <- data.frame(shapedata1)

#set NA values in MaxClass to unclassified to assist with culling the anthropogenic data properly
attrtable1$MaxClass1 <- as.character(attrtable1$MaxClass1)
attrtable1$MaxClass1[is.na(attrtable1$MaxClass1)] <- "unclassified"

#cull the shapefile of anthropogenic transects depending on parameter option setting
if (removeanthro == "yes") {

shapedata <-  shapedata1[attrtable1[,"MaxClass1"] != "anthropogenic"   ,]
attrtable <- data.frame(shapedata)
 } else {
shapedata <-  shapedata1
attrtable <- data.frame(shapedata)
}

getdata2 <- ifelse(file2 == "path",tk_choose.files("","Select the FEMA Q3/FIRM shapefile",multi = FALSE,filetype,1),file2)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)


getdata3 <- ifelse(file3 == "path",tk_choose.files("","Please select the NOAA/NWS MOM-AGL or USACOE SLOSH shapefile",multi = FALSE,filetype,1),file3)
shapename3 <- gsub(".shp", "", basename(getdata3))
shapedata3 <- readOGR(getdata3,layer=shapename3)
attrtable3 <- data.frame(shapedata3)


getdata4 <- ifelse(file4 == "path",tk_choose.files("","Please select the NOAA SoVI shapefile",multi = FALSE,filetype,1),file4)
shapename4 <- gsub(".shp", "", basename(getdata4))
shapedata4 <- readOGR(getdata4,layer=shapename4)
attrtable4 <- data.frame(shapedata4)


# select the directory to store the results in and create necessary time-stamped folder to house results

#tkmessageBox(message = "Please select a directory to store AMBUR-HVA results in")
getdir <- ifelse(dir1 == "path",tk_choose.dir(default = "", caption = "Please select a directory to store AMBUR-HVA results in"),dir1)
workingdir <- getdir
setwd(workingdir)
#dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
#setwd("AMBUR_HVA_data_prep")
time.stamp1 <- as.character(Sys.time())
time.stamp2 <- gsub("[:]", "_", time.stamp1)
dir.create(paste(time.stamp2," ",loc,sep=""))
setwd(paste(time.stamp2," ",loc,sep=""))

##### write the parameter variables to a text file
  txt <- c("ambur.hva( ","loc=",paste('"',loc,'"',sep=""),", ","bufferdist=",bufferdist,", ","bufferdist2=",bufferdist2,", ","removeanthro=",paste('"',removeanthro,'"',sep=""),", ","file1=",paste('"',file1,'"',sep=""),", ","file2=",paste('"',file2,'"',sep=""),", ","file3=",paste('"',file3,'"',sep=""),", ","file4=",paste('"',file4,'"',sep=""),", ","dir1=",paste('"',dir1,'"',sep=""),", ","scr=","c(",paste(as.character(scr), sep="' '", collapse=", "),")", "," ,"tscr=","c(",paste(as.character(tscr), sep="' '", collapse=", "),")",",", "sscr=","c(",paste(as.character(sscr), sep="' '", collapse=", "),")", ",","sovi=","c(",paste(as.character(sovi), sep="' '", collapse=", "),")", ")")

writeLines(txt, "parameters.txt")
#####


###########################################################################################################################
# establish progress bar
pb <- tkProgressBar("AMBUR: progress bar", "Task 1 of 4: Shoreline Change Module...", 0, 100, 5)

# deconstruct shapefile from transects to points of the youngest shoreline and get coordinates
 out <- SpatialPoints(cbind(attrtable[,"Max_DateX"],attrtable[,"Max_DateY"]))



# update progress bar
Pcnt.Complete <-  10
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 1 of 4:Preparing data... (%s)", info), info)

# HVA parameter setup for shoreline change rate rank
EPR_values <- attrtable[,"EPR"]

scr_rank <- character(length(EPR_values ))

scr_rank[attrtable[,"EPR"] >= scr[4]] <- "very low"
scr_rank[attrtable[,"EPR"] > scr[3] & attrtable[,"EPR"] < scr[4]] <- "low"
scr_rank[attrtable[,"EPR"] >= scr[2] & attrtable[,"EPR"] <= scr[3]] <- "medium"
scr_rank[attrtable[,"EPR"] > scr[1] & attrtable[,"EPR"] < scr[2]] <- "high"
scr_rank[attrtable[,"EPR"] <= scr[1]] <- "very high"

attrtable$scr_score <- attrtable[,"EPR"]
attrtable$scr_rank <- scr_rank          


# HVA parameter setup for temporal shoreline change rate rank
EPR_SDEra_values <- attrtable[,"EPR_SDEra"]

tscrv_rank <- character(length(EPR_SDEra_values ))

tscrv_rank[attrtable[,"EPR_SDEra"] >= tscr[1]] <- "high"
tscrv_rank[attrtable[,"EPR_SDEra"] > tscr[2] & attrtable[,"EPR_SDEra"] < tscr[1]] <- "medium"
tscrv_rank[attrtable[,"EPR_SDEra"] <= tscr[2]] <- "low"

attrtable$tscrv_score <- attrtable[,"EPR_SDEra"]
attrtable$tscrv_rank <- tscrv_rank 
          

#################get spatial variability at a given buffer distance
outptsdf <- SpatialPointsDataFrame(out,attrtable)

out2 <- gBuffer(out, byid=TRUE, id=NULL, width=bufferdist, quadsegs=5, capStyle="ROUND",joinStyle="ROUND", mitreLimit=1.0)

out2a.df <- data.frame(attrtable$Transect)
colnames(out2a.df) <- "Transect" 
out2a <- SpatialPolygonsDataFrame(out2,out2a.df)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(out2a) <- projectionString
#writeOGR(out2a, ".", "ambur_hva_spatvar_polys", driver="ESRI Shapefile")  #write buffer polygons for error checking
out3 <- over(out2, outptsdf[,"EPR"], returnList = FALSE, fn = sd)

svar_values <- out3[,1]
sscrv_rank <- character(length(svar_values  ))

sscrv_rank[svar_values >= sscr[1]] <- "high"
sscrv_rank[svar_values > sscr[2] & svar_values < sscr[1]] <- "medium"
sscrv_rank[svar_values <= sscr[2]] <- "low"

attrtable$sscrv_score <- svar_values 
attrtable$sscrv_rank <- sscrv_rank


##################HVA Category

hva_cat <- numeric(length(svar_values  ))

hva_cat[scr_rank == "very high"] <- 5


hva_cat[scr_rank == "high" & tscrv_rank == "high" & sscrv_rank == "high"] <- 5

hva_cat[(scr_rank == "high") & (tscrv_rank == "high") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 4

hva_cat[(scr_rank == "high") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "high")] <- 4

hva_cat[(scr_rank == "high") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 4

hva_cat[scr_rank == "medium" & tscrv_rank == "high" & sscrv_rank == "high"] <- 4


hva_cat[(scr_rank == "medium") & (tscrv_rank == "high") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 3

hva_cat[(scr_rank == "medium") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "high")] <- 3

hva_cat[(scr_rank == "medium") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 3

hva_cat[scr_rank == "low" & tscrv_rank == "high" & sscrv_rank == "high"] <- 3


hva_cat[(scr_rank == "low") & (tscrv_rank == "high") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 2

hva_cat[(scr_rank == "low") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "high")] <- 2

hva_cat[(scr_rank == "low") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 2

hva_cat[scr_rank == "very low" & tscrv_rank == "high" & sscrv_rank == "high"] <- 2


hva_cat[(scr_rank == "very low") & (tscrv_rank == "high") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 1

hva_cat[(scr_rank == "very low") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "high")] <- 1

hva_cat[(scr_rank == "very low") & (tscrv_rank == "low" | tscrv_rank == "medium") & (sscrv_rank == "low" | sscrv_rank == "medium")] <- 1



attrtable$hva_cat <- hva_cat 

sc_hva_table <- data.frame(cbind(attrtable$Transect,attrtable[,"EPR"],attrtable[,"EPR_SDEra"],svar_values,scr_rank,tscrv_rank,sscrv_rank,hva_cat))

colnames(sc_hva_table) <- c("Transect","scr","tscrv","sscrv","scr_rank","tscrv_rank","sscrv_rank","hva_cat")
row.names(sc_hva_table) <- row.names(attrtable)  #to prepare table for a spatialdataframe

#write.table(sc_hva_table, file = "ambur_hva_shorechng_data.csv", sep = ",", row.names = FALSE)   #write CSV table for testing


 Pcnt.Complete <-  25
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 1 of 4:Writing final Shoreline Change HVA shapefile (%s)", info), info)

pts.output <- SpatialPointsDataFrame(out,attrtable) 
projectionString <- proj4string(shapedata) # contains projection info
proj4string(pts.output) <- projectionString
#writeOGR(pts.output, ".", "ambur_hva_shoreline_change_all", driver="ESRI Shapefile") #write a shapefile with all of the original data plus HVA data

transects.output <- SpatialLinesDataFrame(shapedata,sc_hva_table)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(transects.output) <- projectionString
writeOGR(transects.output, ".", "ambur_hva_shoreline_change", driver="ESRI Shapefile")

#write Google Earth KML file
#ambur_hva_shorechng_pts.kml <- spTransform(pts.output, CRS("+proj=longlat +datum=WGS84"))
#writeOGR(ambur_hva_shorechng_pts.kml["hva_cat"],"ambur_hva_shorechng_pts.kml","hva_cat","KML")
 
########################################################################################################################### 
 Pcnt.Complete <-  25
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 2 of 4: Starting Inundation Module... (%s)", info), info)

 
# set up FEMA Q3 Flood parameter ranking
q3_rank <- numeric(length(row.names(attrtable2) ))

q3_rank[attrtable2[,"FLD_ZONE"] == "0.2 PCT ANNUAL CHANCE FLOOD HAZARD"] <- 1
q3_rank[attrtable2[,"FLD_ZONE"] == "D"] <- 1
q3_rank[attrtable2[,"FLD_ZONE"] == "X"] <- 2
q3_rank[attrtable2[,"FLD_ZONE"] == "X_500"] <- 2
q3_rank[attrtable2[,"FLD_ZONE"] == "C"] <- 2
q3_rank[attrtable2[,"FLD_ZONE"] == "B"] <- 3
q3_rank[attrtable2[,"FLD_ZONE"] == "A"] <- 4
q3_rank[attrtable2[,"FLD_ZONE"] == "AE"] <- 4
q3_rank[attrtable2[,"FLD_ZONE"] == "AE (Floodway)"] <- 4
q3_rank[attrtable2[,"FLD_ZONE"] == "AH"] <- 4
q3_rank[attrtable2[,"FLD_ZONE"] == "VE"] <- 5
q3_rank[attrtable2[,"FLD_ZONE"] == "V"] <- 5
q3_rank[attrtable2[,"FLD_ZONE"] == "OPEN WATER"] <- 5

attrtable2$q3_rank <- q3_rank          


###########


 Pcnt.Complete <-  30
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 2 of 4: Assembling Q3 ranked polygons (%s)", info), info)

q3.plys.output <- SpatialPolygonsDataFrame(shapedata2,attrtable2)
projectionString <- proj4string(shapedata2) # contains projection info
proj4string(q3.plys.output) <- projectionString
#writeOGR(q3.plys.output, ".", "ambur_hva_q3_polys", driver="ESRI Shapefile")    #write polygon shapefile for error checking with all attribute data


# set up NOAA SLOSH MOM parameter ranking
slosh_rank <- numeric(length(row.names(attrtable3) ))

if ("c5_high" %in% colnames(attrtable3)) {

slosh_rank[attrtable3[,"c5_high"] > 0 & attrtable3[,"c5_high"] != 99.9] <- 1
slosh_rank[attrtable3[,"c4_high"] > 0 & attrtable3[,"c4_high"] != 99.9] <- 2
slosh_rank[attrtable3[,"c3_high"] > 0 & attrtable3[,"c3_high"] != 99.9] <- 3
slosh_rank[attrtable3[,"c2_high"] > 0 & attrtable3[,"c2_high"] != 99.9] <- 4
slosh_rank[attrtable3[,"c1_high"] > 0 & attrtable3[,"c1_high"] != 99.9] <- 5
 } else {
slosh_rank[attrtable3[,"slosh_cat"] == 1] <- 5
slosh_rank[attrtable3[,"slosh_cat"] == 2] <- 4
slosh_rank[attrtable3[,"slosh_cat"] == 3] <- 3
slosh_rank[attrtable3[,"slosh_cat"] == 4] <- 2
slosh_rank[attrtable3[,"slosh_cat"] == 5] <- 1
slosh_rank[attrtable3[,"slosh_cat"] == 0] <- 1
}
attrtable3$slosh_rank <- slosh_rank          


###########
Pcnt.Complete <-  35
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 2 of 4: Assembling ranked SLOSH polygons (%s)", info), info)

slosh.plys.output <- SpatialPolygonsDataFrame(shapedata3,attrtable3)
projectionString <- proj4string(shapedata3) # contains projection info
proj4string(slosh.plys.output) <- projectionString
#writeOGR(slosh.plys.output, ".", "ambur_hva_slosh_polys", driver="ESRI Shapefile")    #write polygon shapefile for error checking with all attribute data

## begin analysis

shapedataA <- as(q3.plys.output, "SpatialPolygons")
shapedataB <- as(slosh.plys.output, "SpatialPolygons")

int <- gIntersects(q3.plys.output, slosh.plys.output, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

Pcnt.Complete <-  40
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 2 of 4:Intersections (%s)", info), info)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i-1)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb, Pcnt.Complete, sprintf("AMBUR: Task 2 of 4:Intersecting Q3 & SLOSH (%s)", info), info)

vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(q3.plys.output[i,], slosh.plys.output[int[,i],], byid=TRUE)  else 0 }
cond <- lapply_pb(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]
out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))


Pcnt.Complete <-  45                                   
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 2 of 4: Creating output shapefile (%s)", info), info)

transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
sortID <- seq(1,length(transID),1)
inter.data <- data.frame(transID,shoreID,sortID)
tran.data <- data.frame(q3.plys.output@data)
shore.data <- data.frame(slosh.plys.output@data)
tran.data$Id <- as.numeric(row.names(tran.data))
shore.data$mergeID <- as.numeric(row.names(shore.data))
tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)
tet3 <- tet2[ order(tet2[,"sortID"]) , ]
tet3$Id <- tet2[,"sortID"]
row.names(tet3) <- seq(1,length(tet3$Id),1)

#change the shape IDs to match the tet3 IDs
out2 <- spChFIDs(out, as.character(as.character(row.names(tet3))))
hva.calc <- sqrt((tet3$q3_rank * tet3$slosh_rank)/2)
tet3$hva_score <- hva.calc

####scale hva.calc values back to 1 to 5
n.parameters <- 2
cat1 <- sqrt((1^n.parameters)/n.parameters)
cat2 <- sqrt((2^n.parameters)/n.parameters)
cat3 <- sqrt((3^n.parameters)/n.parameters)
cat4 <- sqrt((4^n.parameters)/n.parameters)
cat5 <- sqrt((5^n.parameters)/n.parameters)
break.pts <- c(cat1,cat2,cat3,cat4,cat5)
break.pts2 <- ((break.pts[-1] - break.pts[-5])/2) + break.pts[-5]
break.pts3 <- c(break.pts[1],break.pts2,break.pts[5])


hva_scaled <- numeric(length(hva.calc ))

hva_scaled[hva.calc <= break.pts3[2]] <- 1
hva_scaled[hva.calc > break.pts3[2] & hva.calc <= break.pts3[3]] <- 2
hva_scaled[hva.calc > break.pts3[3] & hva.calc <= break.pts3[4]] <- 3
hva_scaled[hva.calc > break.pts3[4] & hva.calc <= break.pts3[5]] <- 4
hva_scaled[hva.calc > break.pts3[5]] <- 5

tet3$hva_cat <- hva_scaled

inun_hva_table <- data.frame(cbind(tet3$q3_rank,tet3$slosh_rank,tet3$hva_cat))

colnames(inun_hva_table) <- c("q3_rank","slosh_rank","hva_cat")
row.names(inun_hva_table) <- row.names(tet3)  #to prepare table for a spatialdataframe




Pcnt.Complete <-  50
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Task 2 of 4:Writing final Inundation HVA shapefile (%s)", info), info)

outputdata2 <- SpatialPolygonsDataFrame(out2,tet3)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(outputdata2) <- projectionString
#writeOGR(outputdata2, ".", "ambur_hva_inundation_all", driver="ESRI Shapefile") #write all data for error checking

outputdata2a <- SpatialPolygonsDataFrame(out2,inun_hva_table)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(outputdata2a) <- projectionString
writeOGR(outputdata2a, ".", "ambur_hva_inundation", driver="ESRI Shapefile") #write just rankings and HVA data


########################################################################################################################### 
 Pcnt.Complete <-  55
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 3 of 4: SLOSH Inundation & SoVI Module... (%s)", info), info)

# HVA parameter setup for SoVI rank
SoVi_values <- attrtable4[,"SoVI"]

sovi_rank <- character(length(row.names(attrtable4) ))

sovi_rank[attrtable4[,"SoVI"] >= sovi[1]] <- 5
sovi_rank[attrtable4[,"SoVI"] > sovi[2] & attrtable4[,"SoVI"] < sovi[1]] <- 4
sovi_rank[attrtable4[,"SoVI"] >= sovi[3] & attrtable4[,"SoVI"] <= sovi[2]] <- 3
sovi_rank[attrtable4[,"SoVI"] > sovi[4] & attrtable4[,"SoVI"] < sovi[3]] <-2
sovi_rank[attrtable4[,"SoVI"] <= sovi[4]] <- 1

attrtable4$sovi_score <- attrtable4[,"SoVI"]
attrtable4$sovi_rank <- sovi_rank 

sovi.plys.output <- SpatialPolygonsDataFrame(shapedata4,attrtable4)
projectionString <- proj4string(shapedata4) # contains projection info
proj4string(sovi.plys.output) <- projectionString
writeOGR(sovi.plys.output, ".", "ambur_hva_sovi_scaled", driver="ESRI Shapefile")    #write polygon shapefile for error checking with all attribute data

## begin analysis

shapedataA <- as(sovi.plys.output, "SpatialPolygons")
shapedataB <- as(outputdata2a, "SpatialPolygons")

int <- gIntersects(sovi.plys.output, outputdata2a, byid=TRUE)
vec <- vector(mode="list", length=dim(int)[2])

Pcnt.Complete <-  60
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 3 of 4: Intersecting Innundation & SoVI data (%s)", info), info)

for (i in seq(along=vec)) {
Pcnt.Complete <-  round(((i-1)/ length(seq(along=vec))) * 100, 0) 
Pcnt.Complete2 <- paste(Pcnt.Complete," ","%",sep="") 
info <- sprintf("%1.0f percent done", Pcnt.Complete)   
setTkProgressBar(pb, Pcnt.Complete, sprintf("AMBUR: Task 3 of 4: Intersecting Innundation & SoVI data (%s)", info), info)

vec[[i]] <- if (sum(int[,i]) != 0) gIntersection(sovi.plys.output[i,], outputdata2a[int[,i],], byid=TRUE)  else 0 }
cond <- lapply_pb(vec, function(x) class(x) != "numeric")
vec2 <- vec[unlist(cond)]
out <- do.call("rbind", vec2)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))


Pcnt.Complete <-  60                                  
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Task 3 of 4: Creating new polygons for shapefile (%s)", info), info)

transID <- data.frame(nrn)[,1]
shoreID <- data.frame(nrn)[,2]
sortID <- seq(1,length(transID),1)
inter.data <- data.frame(transID,shoreID,sortID)
tran.data <- data.frame(sovi.plys.output@data)
shore.data <- data.frame(outputdata2a@data)
tran.data$Id <- as.numeric(row.names(tran.data))
shore.data$mergeID <- as.numeric(row.names(shore.data))
tet <- merge(inter.data,tran.data , by.x = "transID", by.y = "Id", sort=FALSE)
tet2 <- merge(tet,shore.data, by.x = "shoreID", by.y = "mergeID", sort=FALSE)
tet3 <- tet2[ order(tet2[,"sortID"]) , ]
tet3$Id <- tet2[,"sortID"]
row.names(tet3) <- seq(1,length(tet3$Id),1)

#change the shape IDs to match the tet3 IDs
out2 <- spChFIDs(out, as.character(as.character(row.names(tet3))))
hva.calc <- sqrt((as.numeric(tet3$sovi_rank) * as.numeric(tet3$hva_cat))/2)
tet3$hva_score <- hva.calc

####scale hva.calc values back to 1 to 5
n.parameters <- 2
cat1 <- sqrt((1^n.parameters)/n.parameters)
cat2 <- sqrt((2^n.parameters)/n.parameters)
cat3 <- sqrt((3^n.parameters)/n.parameters)
cat4 <- sqrt((4^n.parameters)/n.parameters)
cat5 <- sqrt((5^n.parameters)/n.parameters)
break.pts <- c(cat1,cat2,cat3,cat4,cat5)
break.pts2 <- ((break.pts[-1] - break.pts[-5])/2) + break.pts[-5]
break.pts3 <- c(break.pts[1],break.pts2,break.pts[5])


hva_scaled <- numeric(length(hva.calc ))

hva_scaled[hva.calc <= break.pts3[2]] <- 1
hva_scaled[hva.calc > break.pts3[2] & hva.calc <= break.pts3[3]] <- 2
hva_scaled[hva.calc > break.pts3[3] & hva.calc <= break.pts3[4]] <- 3
hva_scaled[hva.calc > break.pts3[4] & hva.calc <= break.pts3[5]] <- 4
hva_scaled[hva.calc > break.pts3[5]] <- 5

tet3$hva_cat_c1 <- hva_scaled

inun_sovi_hva_table <- data.frame(cbind(tet3$q3_rank,tet3$slosh_rank,tet3$hva_cat,tet3$sovi_rank,tet3$hva_cat_c1))

colnames(inun_sovi_hva_table) <- c("q3_rank","slosh_rank","inun_cat","sovi_rank","hva_cat")
row.names(inun_sovi_hva_table) <- row.names(tet3)  #to prepare table for a spatialdataframe

Pcnt.Complete <-  75
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Task 3 of 4: Writing final Inundation/SoVI HVA shapefile (%s)", info), info)

outputdata3 <- SpatialPolygonsDataFrame(out2,tet3)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(outputdata3) <- projectionString
#writeOGR(outputdata3, ".", "ambur_hva_inundation_sovi_all", driver="ESRI Shapefile") #write all data for error checking

outputdata3a <- SpatialPolygonsDataFrame(out2,inun_sovi_hva_table)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(outputdata3a) <- projectionString
writeOGR(outputdata3a, ".", "ambur_hva_inundation_sovi", driver="ESRI Shapefile") #write just rankings and HVA data

########################################################################################################################### 
 Pcnt.Complete <-  80
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 4 of 4: Shoreline Change/Inundation/SoVI Module... (%s)", info), info)


 # buffer transects based on their widths or transect spacing


out2 <- gBuffer(transects.output, byid=TRUE, id=NULL, width=bufferdist2, quadsegs=5, capStyle="ROUND",joinStyle="MITRE", mitreLimit=1.0)

out2a.df <- data.frame(sc_hva_table)
 out2a <- SpatialPolygonsDataFrame(out2,sc_hva_table)
projectionString <- proj4string(shapedata) # contains projection info
proj4string(out2a) <- projectionString
#writeOGR(out2a, ".", "ambur_hva_transect_buffers", driver="ESRI Shapefile")  #write buffer polygons for error checking


#reg4 <- gUnaryUnion(out2a, hva_cat)
#reg4.tab <- data.frame(hva_cat=row.names(reg4))
#row.names(reg4.tab) <- row.names(reg4)
#out2b <- SpatialPolygonsDataFrame(reg4,reg4.tab)
#writeOGR(out2b, ".", "ambur_hva_transect_buffers2", driver="ESRI Shapefile")


## begin analysis for composite HVA
Pcnt.Complete <-  85
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR: Task 4 of 4: Intersecting Shoreline Change & Innundation & SoVI data (%s)", info), info)


shapedataA <- as(out2a, "SpatialPolygons")
shapedataB <- as(outputdata3a, "SpatialPolygons")

shapedataAt <- out2a
shapedataBt <- outputdata3a

#merge polygons: 
 
shapeAmiss <- colnames(out2a@data)[!colnames(out2a@data) %in% colnames(outputdata3a@data)]
shapeBmiss <- colnames(outputdata3a@data)[!colnames(outputdata3a@data) %in% colnames(out2a@data)]

shapedataAt@data[shapeBmiss] <- NA
shapedataBt@data[shapeAmiss] <- NA

#assign unique row.names to 2nd shapefile in order to set up merging
shapedataBt2 <- spChFIDs(shapedataBt,as.character((max(length(shapedataAt@data[,1]))+1):(max(length(shapedataBt@data[,1]))+max(length(shapedataAt@data[,1]))))) 


Pcnt.Complete <-  90                                  
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.50 , sprintf("AMBUR: Task 4 of 4: Creating new polygons for output shapefile (%s)", info), info)

#merge data: 
shapedataCC <- rbind(shapedataBt2,shapedataAt)


Pcnt.Complete <-  95
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, i *0.75 , sprintf("AMBUR: Task 4 of 4: Writing final Inundation/SoVI HVA shapefile (%s)", info), info)

 

outputdata4a <- shapedataCC
projectionString <- proj4string(shapedata) # contains projection info
proj4string(outputdata4a) <- projectionString
writeOGR(outputdata4a, ".", "ambur_hva_sc_inundation_sovi", driver="ESRI Shapefile") #write just rankings and HVA data


########################################################################################################################### 
 Pcnt.Complete <-  100
info <- sprintf("%d%% done", Pcnt.Complete)
setTkProgressBar(pb, Pcnt.Complete , sprintf("AMBUR-HVA: All task are complete! (%s)", info), info)


############################################  RASTER EXPERIMENTS   #############################################
##commented lines to be activated only have one "#" symbol
##establish the dimensions of the grid/raster

#Cx <- coordinates(outputdata4a)[,1]
#Cy <- coordinates(outputdata4a)[,2]
#cellsize <- 50
#Cx2 <- length(seq(min(Cx),max(Cx),by=cellsize))
#Cy2 <- length(seq(min(Cy),max(Cy),by=cellsize))

#r <- raster()
#extent(r) <- extent(outputdata4a)
#ncol(r) <- Cx2
#nrow(r) <- Cy2
#test111 <-  rasterize(outputdata4a,r, as.numeric(as.character(outputdata4a@data$hva_cat)), fun='last', background=NA, mask=FALSE, update=TRUE, updateValue='all', filename="", getCover=FALSE, silent=FALSE, progress='text' )

##replace NA values with 5 for open water
##test111[is.na(test111)] <- 5  ## Keep this line removed
#writeRaster(test111, filename="test.img", format="HFA", overwrite=TRUE, progress='text' )

##############################  build interpolated rastera
#xy <- data.frame(xyFromCell(test111, 1:ncell(test111)))
#vals <- getValues(test111)
#r1 <- na.omit(cbind(xy,vals))

####################### inverse distance weighted (IDW)
## require(gstats)  ### keep this line removed and add GSTATS to package requirements
#r <- test111
#mg <- gstat(id = "vals", formula = r1$vals~1, locations = ~x+y, data=r1, nmax=7, set=list(idp = .5))
#z <- interpolate(r, mg, progress='text')
##z <- mask(z, r) ##keep this line removed
#writeRaster(z, filename="test_IDW.img", format="HFA", overwrite=TRUE, progress='text' )

################ kriging
#r2 <- r1
#coordinates(r2) <- ~x+y
#projection(r2) <- projection(shapedata)
#v <- variogram(r2$vals~1, r2, progress='text'  )
#m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
#gOK <- gstat(id="vals", formula=r2$vals~1, data=r2, model=m, nmax=20)
#OK <- interpolate(r, gOK, progress='text' )
##OK <- mask(OK, r)  ## keep this line removed
#writeRaster(OK, filename="test_Kriging.img", format="HFA", overwrite=TRUE, progress='text' )

############### test sampling polygons with points ##try  as.numeric(as.character(outputdata4a@data$hva_cat))
#blah1a <- spsample(outputdata4a, cellsize=20, type='regular')
#b1out <- over(blah1a,outputdata4a[,"hva_cat"], returnList = FALSE, fn = max)
#hvacat_values <- b1out[,1]
#blah1at <- data.frame(coordinates(blah1a),hvacat_values)
#colnames(blah1at) <- c("x_coord","y_coord","hva_cat")
#blah1out <- SpatialPointsDataFrame(blah1a,blah1at)
#projectionString <- proj4string(shapedata) # contains projection info
#proj4string(blah1out) <- projectionString
#writeOGR(blah1out, ".", "ambur_hva_sample_pts", driver="ESRI Shapefile")  #write buffer polygons for error checking

}
