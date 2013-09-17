
ambur.inundationHVA_noSoVI <-
function(userinput1=1) {


# Establish the inputs
nothing <- userinput1



require(rgdal)
require(rgeos)
require(tcltk)


tkmessageBox(message = "Please select the AMBUR_HVA_data_prep directory...")
getdata <- tk_choose.dir(default = "", caption = "Select directory")
workingdir <- getdata
setwd(workingdir)

tkmessageBox(message = "Please select a post-AMBUR analysis transect shapefile generated from ambur.statshape function...")
 filters <- matrix(c("Shapefile", ".shp"), 1, 2, byrow = TRUE)
getdata2 <- tk_choose.files(filter = filters,multi = FALSE)
shapename2 <- gsub(".shp", "", basename(getdata2))
shapedata2 <- readOGR(getdata2,layer=shapename2)
attrtable2 <- data.frame(shapedata2)

#write.table(attrtable2, file = "ambur_hva__data.csv", sep = ",", row.names = FALSE)






time.stamp1 <- as.character(Sys.time())

time.stamp2 <- gsub("[:]", "_", time.stamp1)

#dir.create("AMBUR_HVA_data_prep", showWarnings=FALSE)
#setwd("AMBUR_HVA_data_prep")

#######################


multmerge = function(mypath){
filenames=list.files(path=mypath, full.names=TRUE, pattern = "*.csv" )
datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
Reduce(function(x,y) {merge(x,y,by.x = "Transect", by.y = "Transect", all= TRUE, sort=TRUE)}, datalist)
 }

 #list.files(path=getdata, full.names=TRUE, pattern = "*.csv" )
# lf <- list.files(path=getdata, full.names=TRUE)
# lapply(lf, function(x){read.csv(file=x,header=T)})
 
 merge.tab <- multmerge(getdata)

merged.files <- list.files(path=getdata, full.names=FALSE, pattern = "*.csv" )[-1]

rename1 <- gsub("ambur_hva_", "", merged.files)
rename2 <- gsub("data.csv", "cat", rename1)

merge.tab2 <- as.data.frame(merge.tab)    

colnames(merge.tab2)[(length(colnames(merge.tab2)) - (length(rename2)-1)):(length(colnames(merge.tab2)))] <- rename2

hva_columns <- merge.tab2[,-1]

merge.tab2[,-1][merge.tab2[,-1]==0] <- NA

hva.product <- apply(merge.tab2[,-1],1,prod,na.rm=TRUE)

test.variables <-  ifelse(is.na(merge.tab2[,-1]) == TRUE, 0, 1)      ###get number of variables present
n.variables <- apply(test.variables,1,sum)

hva.calc <- sqrt(hva.product/n.variables)

####scale hva.calc values back to 1 to 5
n.parameters <- length(merged.files)
cat1 <- sqrt((1^n.parameters)/n.parameters)
cat2 <- sqrt((2^n.parameters)/n.parameters)
cat3 <- sqrt((3^n.parameters)/n.parameters)
cat4 <- sqrt((4^n.parameters)/n.parameters)
cat5 <- sqrt((5^n.parameters)/n.parameters)

break.pts <- c(cat1,cat2,cat3,cat4,cat5)

break.pts2 <- ((break.pts[-1] - break.pts[-5])/2) + break.pts[-5]

break.pts3 <- c(break.pts[1],break.pts2,break.pts[5])

#plot(break.pts3,asp=1)
#points(break.pts,col="blue")


hva_scaled <- numeric(length(hva.calc ))

hva_scaled[hva.calc <= break.pts3[2]] <- 1
hva_scaled[hva.calc > break.pts3[2] & hva.calc <= break.pts3[3]] <- 2
hva_scaled[hva.calc > break.pts3[3] & hva.calc <= break.pts3[4]] <- 3
hva_scaled[hva.calc > break.pts3[4] & hva.calc <= break.pts3[5]] <- 4
hva_scaled[hva.calc > break.pts3[5]] <- 5





merge.tab2$hva_calc <- hva.calc
merge.tab2$hva_scaled <- hva_scaled

final.tab <- merge(attrtable2,merge.tab2,by.x = "Transect", by.y = "Transect", all.x= TRUE, sort=TRUE)

row.names(final.tab) <- as.numeric(row.names(attrtable2))

lns.output <- SpatialLinesDataFrame(shapedata2,final.tab)

    projectionString <- proj4string(shapedata2) # contains projection info

  proj4string(lns.output) <- projectionString

writeOGR(lns.output, ".", "ambur_inundation_hva_transects", driver="ESRI Shapefile")

pts.prep <- SpatialPoints(cbind(final.tab$Max_DateX,final.tab$Max_DateY))
pts.output <- SpatialPointsDataFrame(pts.prep,final.tab)
proj4string(pts.output) <- projectionString
writeOGR(pts.output, ".", "ambur.inundationHVA_noSoVI_pts", driver="ESRI Shapefile")


}